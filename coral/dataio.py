"""
This python module contains functions for reading Synthetic Aperture
Radar images

Supports:
  - GAMMA flat binary MLI images (+ .par file)
  - Plain GeoTIFF (+ GAMMA-style .mli.par file)
  - SNAP BEAM-DIMAP products (.dim + .data folder), e.g. the output of
    radiometric calibration on a Sentinel-1 IW SLC scene
"""
import numpy as np
import rasterio
from rasterio.windows import Window
import xml.etree.ElementTree as ET
import os


def readfile(file, sub_im, cr, band=None):
    """Function to read an image file.

    Parameters
    ----------
    file : str
        Path to the image. Extension selects the reader:
        '.dim' -> SNAP BEAM-DIMAP product, '.tif' -> GeoTIFF,
        anything else -> GAMMA flat binary MLI.
    sub_im : int
        Half-width (pixels) of the subset window around the target.
    cr : tuple(int, int)
        (column, row) pixel coordinates of the target.
    band : str, optional
        Only used for BEAM-DIMAP products. Name (or substring) of the
        band to read, e.g. 'Sigma0_VV'. If omitted, the first band
        matching a calibrated-intensity naming pattern is used.
    """

    root, ext = os.path.splitext(file)

    if ext == '.dim':
        print('Reading SNAP BEAM-DIMAP image:', file)
        data, rho_r, rho_a, theta = readdimap(file, sub_im, cr, band)
        return data, rho_r, rho_a, theta

    if ext == '.tif':
        print('Reading tiff image:', file)
        par = readpar(root + '.mli.par')
        data = readtiff(file, sub_im, cr)

    else: # must be GAMMA flat binary float format
        print('Reading flat binary image', file)
        par = readpar(root + ext + '.par')
        data = readmli(file, par, sub_im, cr)

    # extract relevant metadata
    rho_r = float(par['range_pixel_spacing'].split()[0])
    rho_a = float(par['azimuth_pixel_spacing'].split()[0])
    theta = float(par['incidence_angle'].split()[0])

    return data, rho_r, rho_a, theta


def readpar(file):
    """Function to read a GAMMA 'par' file into a dictionary"""
    par={}
    with open(file) as f:
        for line in f:
            if "Gamma" or " " in line:
                break # ignore header line
        for line in f:
            line=line.rstrip() # remove blank lines and whitespace
            if line and not "title" in line:
                (key, val) = line.split(":")
                par[str(key)] = val
    return par


def readmli(datafile, par, sub_im, cr):
    """Function to read a GAMMA mli file and provide a subsetted image"""
    ct = int(par['range_samples']) * int(par['azimuth_lines'])

    dt = np.dtype('>f4') # GAMMA files are big endian 32 bit float

    d = np.fromfile(datafile, dtype=dt, count=ct)

    d = d.reshape(int(par['azimuth_lines']), int(par['range_samples']))
    #print("Number of elements and size of the array is",d.size, d.shape)
    #d[d==0]= np.nan # convert zeros to nan
    return d[cr[1]-sub_im:cr[1]+sub_im,cr[0]-sub_im:cr[0]+sub_im]


def readtiff(datafile, sub_im, cr):
    """Function to read a tiff and provide a subsetted image"""

    with rasterio.open(datafile) as src:
        d = src.read(1, window=Window(cr[0]-sub_im, cr[1]-sub_im, sub_im*2, sub_im*2))

    #print("Number of elements and size of the array is",d.size, d.shape)
    #d[d==0]= np.nan # convert zeros to nan
    return d


# ---------------------------------------------------------------------
# SNAP BEAM-DIMAP support
# ---------------------------------------------------------------------

def readdimap(dim_file, sub_im, cr, band=None):
    """Read a SNAP BEAM-DIMAP product (.dim + .data folder), typically
    the output of Radiometric Calibration on a Sentinel-1 IW SLC scene,
    and provide a subsetted intensity image plus the metadata needed by
    corner_reflector.py.

    The pixel data (.img files inside <product>.data/) are plain ENVI
    raster files, which GDAL/rasterio can open directly - so the actual
    subsetting re-uses the same rasterio.Window approach as readtiff().
    The incidence angle is not a single scalar in SNAP products; it is
    stored as a coarser-resolution tie-point grid
    ('incidenceAngleFromEllipsoid' or similar), which is bilinearly
    interpolated to the exact target pixel location here.

    Parameters
    ----------
    dim_file : str
        Path to the .dim XML file. Its companion .data folder must sit
        alongside it, e.g. 'S1_calibrated.dim' + 'S1_calibrated.data/'.
    sub_im : int
        Half-width of the subset window (pixels) around the target.
    cr : tuple(int, int)
        (column, row) pixel coordinates of the target in the full image.
    band : str, optional
        Name (or substring) of the band to read, e.g. 'Sigma0_VV'. If
        not given, the first band matching a calibrated-intensity
        pattern (Sigma0, Gamma0, Beta0, Intensity, Amplitude) is used.

    Returns
    -------
    data : ndarray
        Subsetted intensity image (sub_im*2 x sub_im*2).
    rho_r, rho_a : float
        Range and azimuth pixel spacing (m), read from Abstracted_Metadata.
    theta : float
        Incidence angle (deg) at the target location, bilinearly
        interpolated from the incidence-angle tie-point grid.
    """
    root, _ = os.path.splitext(dim_file)
    data_dir = root + '.data'

    xroot = ET.parse(dim_file).getroot()

    rho_r, rho_a = _read_dimap_spacing(xroot, dim_file)

    band_name, band_path = _find_dimap_band(xroot, data_dir, dim_file, band)
    print('  using band:', band_name)

    with rasterio.open(band_path) as src:
        d = src.read(1, window=Window(cr[0]-sub_im, cr[1]-sub_im,
                                       sub_im*2, sub_im*2))

    theta = _read_dimap_incidence_angle(xroot, data_dir, cr)

    return d, rho_r, rho_a, theta


def _read_dimap_incidence_angle(xroot, data_dir, cr, name_hint='incid'):
    """Get the incidence angle (deg) at pixel location cr = (col, row).

    Terrain-corrected (TC) BEAM-DIMAP products store the incidence angle
    as a full-resolution raster band (e.g. 'projectedLocalIncidenceAngle')
    rather than a coarse tie-point grid. This tries that route first,
    reading the single pixel value directly, and falls back to the
    tie-point-grid interpolation (used by non-TC products) if no such
    band is found.
    """
    band_path = None
    for df in xroot.iter('Data_File'):
        path_el = df.find('DATA_FILE_PATH')
        if path_el is None:
            continue
        href = path_el.get('href')
        bname = os.path.splitext(os.path.basename(href))[0]
        if name_hint.lower() in bname.lower():
            band_path = os.path.join(data_dir, os.path.basename(href))
            break

    if band_path is not None:
        # apply the same .hdr -> actual data file fix used in
        # _find_dimap_band, in case this band's href also points at .hdr
        if band_path.lower().endswith('.hdr'):
            no_ext = band_path[:-4]
            if os.path.exists(no_ext):
                band_path = no_ext
            elif os.path.exists(no_ext + '.img'):
                band_path = no_ext + '.img'

        with rasterio.open(band_path) as src:
            val = src.read(1, window=Window(cr[0], cr[1], 1, 1))
        return float(val[0, 0])

    # fall back to tie-point grid (non-TC products)
    return _read_dimap_tie_point(xroot, data_dir, cr, name_hint=name_hint)


def _read_dimap_spacing(xroot, dim_file):
    """Extract range/azimuth pixel spacing (m) from a BEAM-DIMAP
    Abstracted_Metadata block."""
    rho_r = rho_a = None
    for attr in xroot.iter('MDATTR'):
        name = (attr.get('name') or '').lower()
        if name == 'range_spacing':
            rho_r = float(attr.text)
        elif name == 'azimuth_spacing':
            rho_a = float(attr.text)

    if rho_r is None or rho_a is None:
        raise ValueError(
            f"Could not find range_spacing/azimuth_spacing metadata in "
            f"{dim_file}. Run inspect_dimap('{dim_file}') to list the "
            "available MDATTR names and adjust _read_dimap_spacing() if "
            "your SNAP version uses different attribute names."
        )
    return rho_r, rho_a


def _find_dimap_band(xroot, data_dir, dim_file, band=None):
    """Locate the .img file for the requested (or best-guess) intensity
    band inside a BEAM-DIMAP .data folder."""
    band_files = {}
    for df in xroot.iter('Data_File'):
        path_el = df.find('DATA_FILE_PATH')
        if path_el is None:
            continue
        href = path_el.get('href')
        bname = os.path.splitext(os.path.basename(href))[0]
        band_files[bname] = os.path.join(data_dir, os.path.basename(href))

    if band is not None:
        matches = [b for b in band_files if band.lower() in b.lower()]
    else:
        patterns = ('sigma0', 'gamma0', 'beta0', 'intensity', 'amplitude')
        matches = [b for b in band_files
                   if any(p in b.lower() for p in patterns)]

    if not matches:
        raise ValueError(
            f"No matching band found in {dim_file}. Available bands: "
            f"{list(band_files.keys())}. Pass e.g. readfile(file, sub_im, "
            "cr, band='Sigma0_VV')."
        )

    band_name = matches[0]
    band_path = band_files[band_name]

    # perbaikan: kalau href mengarah ke .hdr, cari file data aslinya
    if band_path.lower().endswith('.hdr'):
        no_ext = band_path[:-4]           # buang ".hdr"
        if os.path.exists(no_ext):
            band_path = no_ext            # file data tanpa ekstensi
        elif os.path.exists(no_ext + '.img'):
            band_path = no_ext + '.img'   # file data dengan .img

    return band_name, band_path


def _read_dimap_tie_point(xroot, data_dir, cr, name_hint='incid'):
    """Read a BEAM-DIMAP tie-point grid (e.g. incidenceAngleFromEllipsoid)
    and bilinearly interpolate its value at pixel location cr = (col, row)
    in full-resolution image coordinates."""

    info = None
    for tpg in xroot.iter('Tie_Point_Grid_Info'):
        name_el = tpg.find('TIE_POINT_GRID_NAME')
        if name_el is not None and name_hint.lower() in name_el.text.lower():
            info = tpg
            break

    if info is None:
        raise ValueError(
            f"Could not find a tie-point grid matching '{name_hint}' in "
            "the .dim file. Run inspect_dimap(dim_file) to list the "
            "available tie-point grid names and pass the correct one via "
            "name_hint."
        )

    grid_name = info.find('TIE_POINT_GRID_NAME').text
    offset_x = float(info.find('TIE_POINT_OFFSET_X').text)
    offset_y = float(info.find('TIE_POINT_OFFSET_Y').text)
    subsamp_x = float(info.find('TIE_POINT_SUBSAMPLING_X').text)
    subsamp_y = float(info.find('TIE_POINT_SUBSAMPLING_Y').text)

    # locate the tie-point grid's raster file
    tp_path = None
    for tpf in xroot.iter('Tie_Point_Grid_File'):
        path_el = tpf.find('TIE_POINT_GRID_FILE_PATH')
        if path_el is None:
            continue
        href = path_el.get('href')
        if os.path.splitext(os.path.basename(href))[0] == grid_name:
            tp_path = os.path.join(data_dir, os.path.basename(href))
            break

    if tp_path is None:
        # tie-point grids are conventionally stored in a
        # 'tie_point_grids' sub-folder even when not explicitly listed
        candidate = os.path.join(data_dir, 'tie_point_grids', grid_name + '.img')
        if os.path.exists(candidate):
            tp_path = candidate
        else:
            raise ValueError(
                f"Could not locate the raster file for tie-point grid "
                f"'{grid_name}'."
            )

    with rasterio.open(tp_path) as src:
        grid = src.read(1).astype(float)

    # map full-resolution pixel coords -> tie-point grid coords
    tp_x = (cr[0] - offset_x) / subsamp_x
    tp_y = (cr[1] - offset_y) / subsamp_y

    x0 = int(np.floor(tp_x))
    y0 = int(np.floor(tp_y))
    x1 = min(x0 + 1, grid.shape[1] - 1)
    y1 = min(y0 + 1, grid.shape[0] - 1)
    x0 = max(x0, 0)
    y0 = max(y0, 0)

    fx = tp_x - x0
    fy = tp_y - y0

    v00, v01 = grid[y0, x0], grid[y0, x1]
    v10, v11 = grid[y1, x0], grid[y1, x1]

    theta = (v00*(1-fx)*(1-fy) + v01*fx*(1-fy) +
             v10*(1-fx)*fy    + v11*fx*fy)

    return float(theta)


def inspect_dimap(dim_file):
    """Debug helper: print the bands, tie-point grids, and
    Abstracted_Metadata attributes found in a BEAM-DIMAP .dim file, so
    you can confirm/adjust the names used by readdimap() for your
    specific SNAP version/workflow."""
    xroot = ET.parse(dim_file).getroot()

    print('Bands (Data_File):')
    for df in xroot.iter('Data_File'):
        path_el = df.find('DATA_FILE_PATH')
        if path_el is not None:
            print('  -', path_el.get('href'))

    print('\nTie-point grids:')
    for tpg in xroot.iter('Tie_Point_Grid_Info'):
        name_el = tpg.find('TIE_POINT_GRID_NAME')
        if name_el is not None:
            print('  -', name_el.text)

    print('\nAbstracted_Metadata attributes containing "spac" or "incid":')
    for attr in xroot.iter('MDATTR'):
        name = (attr.get('name') or '')
        if 'spac' in name.lower() or 'incid' in name.lower():
            print(f'  - {name}: {attr.text}')
