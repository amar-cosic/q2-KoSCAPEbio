# Changelog

All notable changes to this project will be documented in this file.

## [1.0.2] – 2026-06-11

### Changed
- Updated main README with proper QIIME 2 installation workflow.
- Removed `requirements.txt` for the QIIME plugin to avoid dependency conflicts.
- Improved database curation README with clearer dependency instructions and standalone usage.

### Fixed
- Addressed Biopython dependency vulnerability reported by Dependabot.
- Improved overall documentation clarity and consistency.


## [1.0.1] – 2025-05-08
### Fixed
- Handled missing values in relative abundance output by filling NaN with 0 after normalization.
- Corrected Figures in both README files.

## [1.0.0] – 2025-02-26
### Changed
- Updated example instruction in README (`--i-` → `--p-`).

## [1.0.0] – 2025-02-25
### Added
- Zenodo DOI link in README.
- Initial release of `koscapebio`, a QIIME 2 plugin for computing species-level relative abundances from QIIME.
