# CoRAL-Sentinel1-TCR

A modified extension of **CoRAL** (Corner Reflector Analysis Library), originally developed by Geoscience Australia, adapted for radiometric analysis (RCS/SCR) of a Trihedral Corner Reflector (TCR) using Sentinel-1 SAR imagery.

> **Note on authorship.** This is **not** an original implementation. It is a modification and extension of the original CoRAL library (Copyright 2018 Geoscience Australia, Apache License 2.0), developed as part of an undergraduate/graduate thesis project. All original CoRAL source code, theoretical formulation, and licensing terms are retained and attributed to the original authors (Garthwaite, 2017; Geoscience Australia). See [`LICENSE`](./LICENSE) for the full Apache-2.0 license text.

## Original Reference

This work builds on the corner reflector analysis methodology and formulation described in:

> Garthwaite, M. C. (2017). *On the Design of Radar Corner Reflectors for Deformation Monitoring in Multi-Frequency InSAR*. Remote Sensing, 9(7), 648. https://doi.org/10.3390/rs9070648

and the CoRAL software library:

> Geoscience Australia. *CoRAL: Corner Reflector Analysis Library*. https://github.com/GeoscienceAustralia/CoRAL

## Summary of Modifications

The original CoRAL library computes the Radar Cross Section (RCS) and Signal-to-Clutter Ratio (SCR) of a corner reflector by defining the clutter region as a simple square annulus (a large clutter window minus a smaller, co-centred target window). This repository introduces the following modifications, implemented **without altering** the original CoRAL formulation (Garthwaite, 2017, Eq. 2, 7, 8) or the original `corner_reflector.py` / `dataio.py` source files:

### 1. SNAP BEAM-DIMAP (`.dim`) input support

The original CoRAL library only supports GAMMA flat-binary and plain GeoTIFF inputs. This modification (`coral/dataio.py`) adds a reader for SNAP BEAM-DIMAP products (`.dim` + `.data` folder) — the standard output format of ESA SNAP after radiometric calibration and terrain correction — including:

- automatic band resolution for a requested polarisation (e.g. `Sigma0_VV`), with fallback handling for both full-swath and single-subswath (`split`) Sentinel-1 IW products
- incidence-angle extraction from either a full-resolution raster band (Terrain-Corrected products) or a tie-point grid with bilinear interpolation (non-Terrain-Corrected products)
- range/azimuth pixel spacing extraction from BEAM-DIMAP `Abstracted_Metadata`, ensuring the RCS illuminated-area term (Garthwaite, 2017, Eq. 2) is computed from the true, scene-specific pixel spacing rather than an assumed constant

### 2. Cross-exclusion clutter window definition

`coral/cross_exclusion.py` (new module) reimplements the clutter-region definition following the cross-exclusion window described in Garthwaite et al. (2015, *The Design of Radar Corner Reflectors for the Australian Geophysical Observing System*, Geoscience Australia Record 2015/03, Fig. 4.11), in which a plus-shaped ("cross") region spanning the full width and height of the clutter window is excluded prior to averaging the remaining four corner quadrants as clutter. This is a stricter definition than the original CoRAL square-annulus clutter window, intended to exclude residual mainlobe/sidelobe contamination from the clutter estimate. All downstream energy, SCR, and RCS formulas (`calc_total_energy`, `calc_scr`, `calc_rcs`) are reused unmodified from the original `corner_reflector.py`.

### 3. Global clutter median normalisation (optional)

An optional processing variant (see `RCS_Corner_Reflector_Analysis_CrossExclusion_ClutterNormalized.ipynb`) additionally computes a single global median across cross-excluded clutter pixels (linear sigma-nought domain), pooled across acquisitions in a time series, and subtracts this median from the cross-excluded clutter pixels of each individual scene prior to the RCS/SCR computation (floored at zero to avoid non-physical negative linear values). This is implemented entirely from the notebook via a `readdimap` monkey-patch, without modifying `coral/dataio.py`. Target-window pixels are never affected by this step.

## Repository Structure

# CoRAL - Corner Reflector Analysis Library

[![Travis](https://img.shields.io/travis/GeoscienceAustralia/CoRAL/master.svg?label=Travis%20CI)](https://travis-ci.org/GeoscienceAustralia/CoRAL)
[![Coverage Status](https://coveralls.io/repos/github/GeoscienceAustralia/CoRAL/badge.svg?branch=master)](https://coveralls.io/github/GeoscienceAustralia/CoRAL?branch=master)

CoRAL is a library of tools for analysing the response of corner reflectors (or other targets) in Synthetic Aperture Radar (SAR) images. The library started as a set of Python scripts developed during corner reflector prototyping experiments at Geoscience Australia. This work, and the implemented algorithms, are documented in [Garthwaite 2017](https://doi.org/10.3390/rs9070648).

### Dependencies

This package requires the following:

```
Python==3.6
NumPy
Matplotlib
Rasterio
Jupyterlab (to use the Jupyter Notebook)
```

### Running

An example script `Calc_CR_response.py` is provided, and can be executed as follows: `python3 Calc_CR_response.py`. This example will calculate the Radar Cross Section (RCS) and Signal to Clutter Ratio (SCR) for a single corner reflector. Results can be viewed in the generated image files (copies are lodged in the `data/` directory).

### Example data

A stack of 7 example SAR images is contained in the `data/` directory. The example SAR images are a set of multi-looked Sigma-0 intensity images acquired by Sentinel-1 in Interferometric Wide Swath mode. These images have been processed from Single Look Complex data using the GAMMA Remote Sensing software with 4 Looks in range and 1 Look in azimuth. Each tiff file (`*.tif`) contains 200*200 floating point values. Each data file has an accompanying text file (`*.mli.par`) containing metadata generated by the GAMMA Remote Sensing software. The example images, in descending-pass radar geometry, cover a small area cropped from the full images containing a test site "SERF" at Samford, Queensland, Australia. A corner reflector target was deployed at SERF for teaching purposes in August 2018.

### Jupyter Notebook

A Jupyter notebook is included in the package that demonstrates an example CoRAL analysis using the included SAR data.

### Testing

To run unit tests from within the CoRAL repository directory:

```
python -m unittest discover tests/
```

### Authors

* **Matt Garthwaite** - [mcgarth](https://github.com/mcgarth)

## License

Copyright 2018 Geoscience Australia

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
