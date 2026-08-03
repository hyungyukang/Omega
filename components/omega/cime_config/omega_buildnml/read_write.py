from importlib import resources
from pathlib import Path
from typing import IO

import yaml

from ._types import PathLike, YamlMapping


def read_default_config(path: PathLike) -> YamlMapping:
    """
    Read the default configuration file from the package resources.

    Parameters:
    -----------
    path : PathLike
        Path to default config file (i.e. components/omega/config/Defaults.yml)

    Returns:
    --------
    dict[str, Any]
        The read configuration as a dictionary.
    """
    path = Path(path)

    if not path.is_file():
        err_msg = f"{path} does not exist or is not a file"
        raise FileNotFoundError(err_msg)

    with path.open("r", encoding="utf-8") as f:
        config = _read_yaml_file(f)

    if "Omega" not in config:
        err_msg = f"{path} does not contain a top-level 'Omega' section"
        raise ValueError(err_msg)

    # TODO: check that the config is valid (e.g., required keys are present)
    return config["Omega"]


def read_input_files_config() -> YamlMapping:
    """
    Read the input_files.yaml configuration file from the package resources.

    Returns:
    --------
    dict[str, Any]
        The read configuration as a dictionary.
    """
    return _read_packaged_yaml_file("input_files.yaml")


def read_config_overrides() -> YamlMapping:
    """
    Read config_overrides.yaml configuration file from the package resources.

    Returns:
    --------
    dict[str, Any]
        The read configuration as a dictionary.
    """
    return _read_packaged_yaml_file("config_overrides.yaml")


def write_yaml_mapping(
    mapping: YamlMapping, file_path: PathLike
) -> None:
    """
    Wite a mapping to a YAML file.

    Parameters:
    -----------
    mapping : dict[str, Any]
        The mapping to write to the YAML file.
    file_path : PathLike
        The path to the output YAML file.

    Returns:
    --------
    None
    """
    with Path(file_path).open("w", encoding='utf-8') as f:
        yaml.safe_dump(mapping, f, sort_keys=False, default_flow_style=False)


def write_input_data_list(
    streams_files: dict[str, str], casebuild: PathLike
) -> None:
    """
    Build omega.input_data_list

    Enables automatic retrival of missing input files

    Parameters
    ----------
    streams_files : dict[str, str]
        Dict mapping stream names to their corresponding input file paths.
    casebuild : Path
        Path to the case build directory (i.e. CASEBUILD)

    Returns
    -------
    None
    """

    unique_files = list(dict.fromkeys(streams_files.values()))

    input_data_list = [
        f"omega_input_{index} = {filename}"
        for index, filename in enumerate(unique_files, start=1)
    ]

    path = Path(casebuild) / "omega.input_data_list"
    with path.open("w", encoding='utf-8') as f:
        for input_file in input_data_list:
            f.write(f"{input_file}\n")


def _read_yaml_file(f: IO[str]) -> YamlMapping:
    """
    Read a YAML mapping from a file.

    Parameters
    ----------
    f : TextIO
        File-like object to read the YAML mapping from.
    Returns:
    -------
    dict[str, Any]
        Read YAML mapping.
    """

    # TODO: Reject duplicate key mappings
    return yaml.safe_load(f)


def _read_packaged_yaml_file(file_name: str) -> YamlMapping:
    """
    Read a YAML mapping packaged within config_builder/data

    Parameters
    ----------
    file_path : str
        Name of the packaged YAML file to read.

    Returns:
    --------
    dict[str, Any]
        Read YAML mapping.
    """
    resource = resources.files("omega_buildnml.data").joinpath(file_name)

    if not resource.is_file():
        err_msg = f"Packaged configuration file '{file_name}' not found."
        raise FileNotFoundError(err_msg)

    with resource.open("r", encoding="utf-8") as f:
        return _read_yaml_file(f)
