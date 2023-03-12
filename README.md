## pmutil

**pmutil** is several bash and python scripts, which helps perform full fotometry pipeline. It is based on mainly the FITSH, Astrometry.net and SExtractor pakcages.

# Preconditions

These scripts run on Linux, and they have some dependecies:

- **python3** 3.6+ for .py scripts
- python packages: **astroalign**, **astropy**, **astroquery**, **matplotlib**, **numpy**, **Pillow**, **photutils**, **scipy**, **xmltodict**
- **FITSH** 0.9.3+ (https://fitsh.net)
- **RawTran** 0.3.8+ (http://integral.physics.muni.cz/rawtran/)
- **ExifTool** 10.80+ (https://sourceforge.net/projects/exiftool/)
- **Astrometry.net** 0.77+ (http://astrometry.net)
- **SExtractor** 2.19.5+ (http://www.astromatic.net/software/sextractor)
- **wcstools** 3.9.5+ (http://tdc-www.harvard.edu/software/wcstools/)

# Download and install pmutil

Download the pmutil source from https://github.com/kovihome/pmutil/archive/master.zip to folder <pmutil_install>

Run the <pmutil_install>/src/main/configure script to install the package.

Some configuration files of third party programs placed into the ~/bin folder, a few configuration values need to be changed.

astrometry.cfg: lines containing 'add_path'
sex.cgf: 'PARAMETERS_NAME' value

# Usage

TBD.
