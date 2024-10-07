<--
Copyright (C) 2020-2024 CERN and UCLouvain.
Licensed under the GNU Lesser General Public License (version 3 or later).
Created by: A. Valassi (Oct 2023) for the MG5aMC CUDACPP plugin.
Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.
-->

# Changelog

All notable changes to this project will be documented in this file.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com).

--------------------------------------------------------------------------------

## [Unreleased] - 2024-10-06

### Added

- Infrastructure issues
  - AV ([#945]) Added cudacpp_bldall runcard to produce multi-backend gridpacks.
  - AV ([#700]) Added cudacpp_helinl and cudacpp_hrdcod runcards to support HELINL=1 and HRDCOD=1 builds.
  - AV ([#957]) In internal "tlau" tests, instrumented python code in gridpacks to provide timing profiles for event generation.
  - AV Enhanced internal "tlau" tests and added a few test logs and various scripts to analyse them.
  - AV ([#972]) Added new RDTSC-based timers (faster than the legacy chrono-based timers) and adapted all internal APIs.

### Changed

- Updated cudacpp version to 1.00.01.

- Infrastructure issues
  - AV ([#700]) Renamed the floating_type runcard as cudacpp_fptype for consistency with other cudacpp runcards.

### Fixed

- Platform-specific issues
  - AV ([#1011]) Added workaround for Floating Point Exceptions in vxxxxx in the HIP backend.

- Infrastructure issues
  - AV ([#1013]) Fixed release scripts to create 'v1.00.01' tags from a '(1,0,1)' python tuple.
  - AV ([#1015]) Removed add_input_for_banner from output.py (plugin_run_card is not needed in cudacpp).
  - AV ([#995]) In cudacpp_config.mk moved default FPTYPE from 'd' to 'm' (already the default cudacpp_backend in run_card.dat).

--------------------------------------------------------------------------------

## [1.00.00] - 2024-10-03

### Added

- (OM+AV+SR+SH+ZW+JT+DM) First release of the MG5aMC CUDACPP plugin.
  - Validated and released for MG5aMC version 3.6.0.
  - Hosted in the https://github.com/madgraph5/madgraph4gpu original repo.
  - Repo uses the original directory structure (plugin is epochX/cudacpp/CODEGEN/PLUGIN/CUDACPP_SA_OUTPUT).

### Known issues

- This section lists some of the main new issues identified in release v1.00.00.

- General issues
  - ([#959]) Cross-section instabilities when changing vector size between 32 and 16384.
  - ([#993]) LHE file mismatch (fortran vs cudacpp) in the experimental multi-backend gridpacks.

- Platform-specific issues
  - ([#1011]) Floating Point Exceptions in vxxxxx in the HIP backend.

- Physics-process-specific issues
  - ([#944]) Cross-section mismatch (fortran vs cudacpp) in Drell-Yan plus 4 jets.
  - ([#942]) Floating Point Exceptions in Drell-Yan plus 0 to 2 jets (workaround: `CUDACPP_RUNTIME_DISABLEFPE=1`).
  - ([#846]) ME mismatch (HRDCOD=1 vs HRDCOD=1) in EWdim6 models.
  - ([#601]) Builds fail with very complex final states (e.g. gg to ttgggg).

--------------------------------------------------------------------------------

[1.00.00]: https://github.com/madgraph5/madgraph4gpu/releases/tag/cudacpp_for3.6.0_v1.00.00
[Unreleased]: https://github.com/madgraph5/madgraph4gpu/releases/compare/cudacpp_for3.6.0_v1.00.00...HEAD

[#601]: https://github.com/madgraph5/madgraph4gpu/issues/601
[#700]: https://github.com/madgraph5/madgraph4gpu/issues/700
[#846]: https://github.com/madgraph5/madgraph4gpu/issues/846
[#942]: https://github.com/madgraph5/madgraph4gpu/issues/942
[#944]: https://github.com/madgraph5/madgraph4gpu/issues/944
[#945]: https://github.com/madgraph5/madgraph4gpu/issues/945
[#957]: https://github.com/madgraph5/madgraph4gpu/issues/957
[#959]: https://github.com/madgraph5/madgraph4gpu/issues/959
[#972]: https://github.com/madgraph5/madgraph4gpu/issues/972
[#993]: https://github.com/madgraph5/madgraph4gpu/issues/993
[#995]: https://github.com/madgraph5/madgraph4gpu/issues/995
[#1011]: https://github.com/madgraph5/madgraph4gpu/issues/1011
[#1013]: https://github.com/madgraph5/madgraph4gpu/issues/1013
[#1015]: https://github.com/madgraph5/madgraph4gpu/issues/1015
