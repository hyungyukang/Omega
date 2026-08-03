#!/usr/bin/env python3

"""
Validate the configuration files used to build Omega's runtime config

Every entry in the packaged configuration files is validated, rather than only
the entries needed by a single mesh. This ensures the configuration files are
valid and complete, when called from ./case.setup.

The following files are validated:

  components/omega/cime_config/omega_buildnml/data/input_files.yaml
    Checks that every mesh assigns an input file to each of the required
    IOStreams (HorzMeshIn, InitialVertCoord, InitialState), and that no
    IOStream is assigned more than one file. IOStream names are checked
    against the IOStreams Omega defines in configs/Default.yml.

  components/omega/cime_config/omega_buildnml/data/config_overrides.yaml
    Checks that the coupled overrides are present, and that every mesh with
    overrides is a mesh defined in input_files.yaml. The coupled and mesh
    specific options are checked against configs/Default.yml, so that an
    override cannot silently add a new option, rather than setting an existing
    one. IOStreams overrides are not checked, as they are allowed to define
    new streams.

  components/omega/cime_config/omega_buildnml/validate.py
    Checked that KNOWN_STREAMS has not drifted from the IOStreams defined in
    configs/Default.yml.

Exits non-zero, and reports all the problems found, if any file is invalid.
"""

import argparse

from omega_buildnml import (
    DEFAULT_CONFIG_PATH,
    read_config_overrides,
    read_default_config,
    read_input_files_config,
)


def main() -> None:
    """
    Validate all of Omega's packaged configuration files.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.parse_args()

    # validate all the mesh entries in input_files.yaml
    read_input_files_config()

    # validate the coupled and all the mesh entries in config_overrides.yaml
    read_config_overrides()

    # check KNOWN_STREAMS has not drifted from the IOStreams in Default.yml
    read_default_config(DEFAULT_CONFIG_PATH, check_streams=True)

    print("PASS: Omega configuration files are valid")


if __name__ == "__main__":
    main()
