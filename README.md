## pmutil

**pmutil** performs full fotometry pipeline. It is based on mainly the FITSH, Astrometry.net and SExtractor packages and run on Linux based machines.

# Preconditions

Before installing pmutil, you need the packages to be installed:

- **python** 3.8+ for .py scripts
- python packages: **astroalign**, **astropy**, **astroquery**, **matplotlib**, **numpy**, **pillow**, **photutils**, **rawpy**, **scipy**, **xmltodict**
- **FITSH** 0.9.3+ (https://fitsh.net)
- **Astrometry.net** 0.81+ (http://astrometry.net)
- **SExtractor** 2.19.5+ (http://www.astromatic.net/software/sextractor)
- **ExifTool** 10.82+ (https://sourceforge.net/projects/exiftool/)

---
# Download and install pmutil

Download the pmutil source from https://github.com/kovihome/pmutil/archive/master.zip to folder `<pmutil_install>`.

Run the `<pmutil_install>/src/main/configure` script to install the package.

Some configuration files of third party programs placed into the `~/.pmlib` folder, a few configuration values need to be changed.

- astrometry.cfg: lines containing 'add_path'
- sex.cgf: 'PARAMETERS_NAME' value

---
# Usage

For more detailed description see the [pmutil wiki pages](https://github.com/kovihome/pmutil/wiki).

### Calibration
Prepare and process raw data:
```bash
  ppl-calibration <options> image-folder
```

### Photometry
Perform photometric measurements:
```bash
  ppl-photometry <options> image-folder
```

### Transient Search
Identify transient objects:
```bash
  ppl-transient <options> image-folder
```

### Helper Scripts
- `ppl-refcat`: Create a reference catalog.
- `ppl-colorize`: Create color image from calibrated images.
- `ppl-clean`: Clean up work folders.

---
## Contributing
Contributions are welcome! Please follow these steps:
1. Fork the repository.
2. Create a feature branch.
3. Submit a pull request with your changes.

---
## **Support This Project**

If you find this project useful, please consider supporting it:

<!--
- [GitHub Sponsors](https://github.com/sponsors/yourusername)
- [Patreon](https://www.patreon.com/yourproject)
- [Open Collective](https://opencollective.com/yourproject)
- [Liberapay]()
- [Buy me a coffee](https://buymeacoffee.com)
-->
- **_[PayPal](https://paypal.me/kovihome?country.x=HU&locale.x=hu_HU)_**

---
## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---
## Support
If you encounter issues, please open an issue on the [GitHub Issues](https://github.com/your-repository/pmutil/issues) page or contact the maintainer.

---
## Acknowledgements
Special thanks to the developers of [FITSH](https://fitsh.net/), [Astrometry.net](http://astrometry.net/), and [SExtractor](https://www.astromatic.net/software/sextractor) for providing essential tools for this pipeline.

---
## **Contact**

If you have any questions, feel free to contact us at:
- **Email**: [kovihome86@gmail.com](mailto:kovihome86@gmail.com)
- **GitHub Issues**: [Issues Page](https://github.com/kovihome/ReqSmith/issues)
