# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[Unreleased]](https://github.com/kovihome/pmutil) - transient search
### Added
- Using precalibrated FITS (SeeStar FITS images)
- New algorithm for aligning and stacking
- Save graphs with --saveGraph option
- Manage flat library
- Archivate raw images and result files

### Modified
- FITS headers contain instrument and camera name
- Use observer namecode from FITS header
- Store mag limit in catalog files
- Store standard coefficients in catalog files
- Preserve records even mag is below limit
- Make sequence and combined photometry together 
- Update reference catalog
- Reference catalog contains field stars
- Make saved command file executable

## [[1.1] - 2021-01-22](https://github.com/kovihome/pmutil/releases/tag/pmutil-v1.1.0) - Standardization
### Added
- Calculate magnitude limit with linear fit
- Select best reference image for aligning
- Create histogram equalized jpg
- Create monochrome jpg
- .pmlib folder for configuration and common files
- Save calibrated image in jpg
- More magnitude calculation methods
- Standardization
- Manage hot/dead pixels

## [[1.0] - 2020-01-26](https://github.com/kovihome/pmutil/releases/tag/pmutil-v1.0.0) - Original release
### Added
- Image calibration with FITSH
- Photometry with Astromerty.net and SExtractor
- Using stored flats
- Reference catalog for variables and comparison stars
- Observation report in AAVSO format
