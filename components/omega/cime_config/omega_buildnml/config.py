from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

from ._types import PathLike, YamlMapping


def build_omega_config(
    defaults: YamlMapping,
    coupled_overrides: YamlMapping,
    mesh_overrides: YamlMapping,
    runtime_overrides: YamlMapping,
) -> YamlMapping:
    """
    Build the Omega configuration dictionary.

    Parameters:
    -----------
    defaults : dict[str, Any]
        Default configuration values, loaded from config/Defaults.yaml.
    coupled_overrides : dict[str, Any]
        Coupled model overrides, loaded from cime_config/config_overrides.yaml
    mesh_overrides: dict[str, Any]
        Mesh specific overrides, loaded from cime_config/mesh_overrides.yaml
    runtime_overrides : dict[str, Any]
        Runtime specific overrides, based on CIME case configuration.

    Returns:
    --------
    dict[str, Any]:
        Final Omega configuration dictionary.
    """
    if "Omega" in defaults:
        defaults = defaults["Omega"]

    config = deepcopy(defaults)

    config = _deep_merge(config, coupled_overrides)
    config = _deep_merge(config, mesh_overrides)
    config = _deep_merge(config, runtime_overrides)

    return {"Omega": config}


def resolve_streams_files(
    input_files: YamlMapping, mesh_name: str, din_loc_root: PathLike
) -> dict[str, str]:
    """
    Resolve input filenames for all Omega IOStreams

    Uses mesh name and input_files.yaml to determine the correct input files
    for each IOStream.

    Parameters:
    -----------
    input_files : dict[str, str]
        Parsed content of cime_config/input_files.yaml.
    mesh_name : str
        CIME ocean grid name
    din_loc_root : Path
        Path to root of E3SM inpute data directory.

    Returns:
    --------
    dict[str, str]
        Mapping of Omega IOStream names to resolved input filenames.
    """

    # NOTE: This function has **very** pedantic error checking as a way of
    # ensuring the entries in input_files.yaml are correct and complete

    err_suffix = (
        "Please check your setting in `components/omega/cime_config/"
        "input_files.yaml`"
    )
    meshes: YamlMapping = input_files.get("meshes", {})

    if mesh_name not in meshes:
        err_msg = (
            f"Unsupported OCN_GRID for Omega: {mesh_name}. \n" + err_suffix
        )
        raise ValueError(err_msg)

    mesh_definition: YamlMapping = meshes[mesh_name]

    inputs = mesh_definition.get("inputs")
    if not inputs:
        err_msg = (
            f"No input files defined for: {mesh_name}. \n" + err_suffix
        )
        raise ValueError(err_msg)

    mesh_dir = Path(din_loc_root) / "ocn" / "omega" / mesh_name
    streams_files = {}

    for index, input_group in enumerate(inputs):

        err_msg = (
            "Missing {key} in input group {index} for mesh: {mesh_name}. \n" +
            err_suffix
        )
        if 'file' not in input_group:
            _err_msg = err_msg.format(
                key='file', index=index, mesh_name=mesh_name
            )
            raise ValueError(err_msg)
        if 'streams' not in input_group:
            _err_msg = err_msg.format(
                key='streams', index=index, mesh_name=mesh_name
            )
            raise ValueError(err_msg)

        file_name = input_group['file']
        streams = input_group['streams']

        if not isinstance(file_name, str) or not file_name:
            _err_msg = err_msg.format(
                key='file', index=index, mesh_name=mesh_name
            )
            raise ValueError(_err_msg)

        if not isinstance(streams, list) or not streams:
            _err_msg = err_msg.format(
                key='streams', index=index, mesh_name=mesh_name
            )
            raise ValueError(_err_msg)

        resolved_file_path = mesh_dir / file_name

        for stream in streams:
            if stream in streams_files:
                err_msg = (
                    f"Stream '{stream}' is assigned more than once for mesh: "
                    f"{mesh_name}. \n" + err_suffix
                )
                raise ValueError(err_msg)

            streams_files[stream] = str(resolved_file_path)

    return streams_files


def build_runtime_overrides(
    calendar: str,
    continue_run: bool,
    case_name: str,
    streams_files: dict[str, str],
) -> YamlMapping:
    """
    Build the runtime override dictionary from the CIME case configuration.

    Parameters:
    -----------
    calendar : str
        CIME calendar type (i.e., "NO_LEAP", "GREGORIAN").
    continue_run : bool
        Whether to continue a previous run.
    case_name : str
        CIME case name.
    streams_files : dict[str, str]
        Resolved input filenames keyed by Omega IOStream name.

    Returns:
    --------
    dict[str, Any]
        Runtime overrides dictionary.
    """
    required_streams = {"HorzMeshIn", "InitialVertCoord", "InitialState"}

    missing_streams = required_streams - set(streams_files)
    if missing_streams:
        raise ValueError(
            f"Missing required input streams: {', '.join(missing_streams)}"
        )

    io_overrides = {
        stream_name: {"Filename": filename}
        for stream_name, filename in streams_files.items()
    }

    if continue_run:
        io_overrides["InitialState"]["FreqUnits"] = "Never"
        io_overrides["RestartRead"] = {"FreqUnits": "OnStartup"}
    else:
        io_overrides["InitialState"]["FreqUnits"] = "OnStartup"
        io_overrides["RestartRead"] = {"FreqUnits": "Never"}

    io_overrides["RestartWrite"] = {
        "Filename": f"{case_name}.omega.r.$Y-$M-$D_$h.$m.$s"
    }

    return {
        "TimeIntegration": {"CalendarType": _to_omega_calendar(calendar)},
        "IOStreams": io_overrides,
    }


def _deep_merge(
    base: Mapping[str, Any],
    override: Mapping[str, Any],
) -> YamlMapping:
    """
    Recursively merge two mappings without modifying either input.

    Nested mappings are merged. All other values, including lists, are replaced
    by the override value.

    Parameters:
    -----------
    base (dict[str, Any])
        The base mapping
    override (dict[str, Any]):
        The mapping with the override values

    Returns:
        dict[str, Any]: Merged mapping
    """
    merged: YamlMapping = deepcopy(dict(base))

    for key, override_value in override.items():
        base_value = merged.get(key)

        if (
            isinstance(base_value, Mapping) and
            isinstance(override_value, Mapping)
        ):
            merged[key] = _deep_merge(base_value, override_value)
        else:
            merged[key] = deepcopy(override_value)

    return merged


def _to_omega_calendar(calendar: str) -> str:
    """
    Convert a CIME calendar string to corresponding Omega calendar string

    Args:
        calendar (str): CIME calendar string.

    Returns:
        str: Omega calendar string.
    """
    if calendar == "NO_LEAP":
        return "No Leap"
    elif calendar == "GREGORIAN":
        return "Gregorian"
    else:
        raise ValueError(f"Unsupported calendar type: {calendar}")
