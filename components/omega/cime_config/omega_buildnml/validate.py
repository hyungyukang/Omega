from typing import Optional

from ._types import YamlMapping

DATA_PATH = "components/omega/cime_config/omega_buildnml/data"

INPUT_FILES_PATH = f"{DATA_PATH}/input_files.yaml"

OVERRIDES_PATH = f"{DATA_PATH}/config_overrides.yaml"

#: IOStreams defined in ``components/omega/configs/Default.yml``
KNOWN_STREAMS = frozenset(
    {
        "HorzMeshIn",
        "InitialVertCoord",
        "InitialState",
        "Forcing",
        "RestartRead",
        "RestartWrite",
        "History",
        "Highfreq",
    }
)

#: IOStreams that every mesh must provide an input file for
REQUIRED_STREAMS = frozenset(
    {"HorzMeshIn", "InitialVertCoord", "InitialState"}
)

#: Keys allowed in a single mesh entry
MESH_KEYS = frozenset({"inputs"})

#: Keys allowed in a single entry of a mesh's ``inputs`` list
INPUT_GROUP_KEYS = frozenset({"file", "streams"})

#: Keys allowed at the top level of ``config_overrides.yaml``
OVERRIDES_KEYS = frozenset({"coupled", "meshes"})


def validate_input_files_config(
    input_files: YamlMapping, mesh_name: Optional[str] = None
) -> YamlMapping:
    """
    Validate the contents of the ``input_files.yaml`` configuration.

    All problems found are collected and reported together, rather than
    raising on the first one encountered.

    Parameters:
    -----------
    input_files : dict[str, Any]
        Parsed content of ``cime_config/omega_buildnml/data/input_files.yaml``
    mesh_name : str, optional
        The name of the mesh to validate. If not provided, all mesh entries
        will be validated.

    Returns:
    --------
    dict[str, Any]
        The validated configuration.

    Raises:
    -------
    ValueError
        If any required keys are missing or if any values are invalid.
    """
    if not isinstance(input_files, dict) or not input_files:
        err_msg = f"`{INPUT_FILES_PATH}` is empty or is not a mapping."
        raise ValueError(err_msg)

    unknown_keys = sorted(set(input_files) - {"meshes"})
    if unknown_keys:
        _raise(
            [f"Unknown top-level key(s): {', '.join(unknown_keys)}."],
            INPUT_FILES_PATH,
        )

    meshes: YamlMapping = input_files.get("meshes", {})
    if not isinstance(meshes, dict) or not meshes:
        _raise(
            ["`meshes` is missing, empty, or is not a mapping."],
            INPUT_FILES_PATH,
        )

    if mesh_name is None:
        errors = []
        for name in meshes:
            errors.extend(_validate_input_files_entry(input_files, name))
        _raise(errors, INPUT_FILES_PATH)
        return input_files

    if mesh_name not in meshes:
        err_msg = (
            f"Unsupported OCN_GRID for Omega: {mesh_name}. \n"
            f"Could not find entry in `{INPUT_FILES_PATH}`"
        )
        raise ValueError(err_msg)

    _raise(_validate_input_files_entry(input_files, mesh_name),
           INPUT_FILES_PATH)

    return input_files


def validate_config_overrides(
    config_overrides: YamlMapping,
    input_files: YamlMapping,
    mesh_name: Optional[str] = None,
) -> YamlMapping:
    """
    Validate the contents of the ``config_overrides.yaml`` configuration.

    All problems found are collected and reported together, rather than
    raising on the first one encountered.

    Mesh specific overrides are optional, so a mesh without an entry is not
    an error.

    Parameters:
    -----------
    config_overrides : dict[str, Any]
        Parsed content of
        ``cime_config/omega_buildnml/data/config_overrides.yaml``
    input_files : dict[str, Any]
        Parsed content of ``cime_config/omega_buildnml/data/input_files.yaml``,
        which defines the meshes Omega supports.
    mesh_name : str, optional
        The name of the mesh to validate. If not provided, all mesh entries
        will be validated.

    Returns:
    --------
    dict[str, Any]
        The validated configuration.

    Raises:
    -------
    ValueError
        If any required keys are missing or if any values are invalid.
    """
    if not isinstance(config_overrides, dict) or not config_overrides:
        err_msg = f"`{OVERRIDES_PATH}` is empty or is not a mapping."
        raise ValueError(err_msg)

    errors: list[str] = []

    unknown_keys = sorted(set(config_overrides) - OVERRIDES_KEYS)
    if unknown_keys:
        errors.append(f"Unknown top-level key(s): {', '.join(unknown_keys)}.")

    coupled = config_overrides.get("coupled")
    if not isinstance(coupled, dict) or not coupled:
        errors.append("`coupled` is missing, empty, or is not a mapping.")

    meshes: YamlMapping = config_overrides.get("meshes", {})
    if not isinstance(meshes, dict):
        errors.append("`meshes` is not a mapping.")
        _raise(errors, OVERRIDES_PATH)

    if mesh_name is None:
        for name in meshes:
            errors.extend(
                _validate_config_overrides_entry(
                    config_overrides, input_files, name
                )
            )
    elif mesh_name in meshes:
        errors.extend(
            _validate_config_overrides_entry(
                config_overrides, input_files, mesh_name
            )
        )

    _raise(errors, OVERRIDES_PATH)

    return config_overrides


def _raise(errors: list[str], config_path: str) -> None:
    """
    Raise a single ``ValueError`` describing all accumulated errors.

    Does nothing when ``errors`` is empty.

    Parameters:
    -----------
    errors : list[str]
        Error messages collected during validation.
    config_path : str
        Path of the configuration file the errors were found in.

    Raises:
    -------
    ValueError
        If ``errors`` is non-empty.
    """
    if not errors:
        return

    details = "\n".join(f"  - {error}" for error in errors)
    err_msg = (
        f"Invalid Omega configuration:\n{details}\n"
        f"Please check your setting in `{config_path}`"
    )
    raise ValueError(err_msg)


def _validate_input_files_entry(
    input_files: YamlMapping, mesh_name: str
) -> list[str]:
    """
    Validate that the specified mesh has a valid configuration in input_files.

    Parameters:
    -----------
    input_files : dict[str, Any]
        Parsed content of ``cime_config/omega_buildnml/data/input_files.yaml``
    mesh_name : str
        The name of the mesh to validate.

    Returns:
    --------
    list[str]
        Error messages describing any problems found. Empty when the entry is
        valid.
    """
    meshes: YamlMapping = input_files["meshes"]
    mesh: YamlMapping = meshes[mesh_name]

    if not isinstance(mesh, dict) or not mesh:
        return [f"Mesh: {mesh_name} is empty or is not a mapping."]

    errors: list[str] = []
    streams_files = {}

    unknown_keys = sorted(set(mesh) - MESH_KEYS)
    if unknown_keys:
        errors.append(
            f"Unknown key(s) {', '.join(unknown_keys)} for "
            f"mesh: {mesh_name}."
        )

    inputs = mesh.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        errors.append(
            f"Missing inputs, or inputs is not a non-empty list, for "
            f"mesh: {mesh_name}."
        )
        return errors

    missing_msg = "Missing {key} in input group {index} for mesh: {mesh_name}."

    for index, input_group in enumerate(inputs):

        if not isinstance(input_group, dict):
            errors.append(
                f"Input group {index} is not a mapping for "
                f"mesh: {mesh_name}."
            )
            continue

        unknown_keys = sorted(set(input_group) - INPUT_GROUP_KEYS)
        if unknown_keys:
            errors.append(
                f"Unknown key(s) {', '.join(unknown_keys)} in input group "
                f"{index} for mesh: {mesh_name}."
            )

        file_name = input_group.get('file')
        streams = input_group.get('streams')

        if not isinstance(file_name, str) or not file_name:
            errors.append(
                missing_msg.format(
                    key='file', index=index, mesh_name=mesh_name
                )
            )

        if not isinstance(streams, list) or not streams:
            errors.append(
                missing_msg.format(
                    key='streams', index=index, mesh_name=mesh_name
                )
            )
            continue

        for stream in streams:
            if not isinstance(stream, str) or not stream:
                errors.append(
                    f"Stream names must be non-empty strings, got "
                    f"'{stream}' in input group {index} for "
                    f"mesh: {mesh_name}."
                )
                continue

            if stream not in KNOWN_STREAMS:
                errors.append(
                    f"Unknown IOStream '{stream}' in input group {index} for "
                    f"mesh: {mesh_name}. Valid IOStreams are: "
                    f"{', '.join(sorted(KNOWN_STREAMS))}."
                )
                continue

            if stream in streams_files:
                errors.append(
                    f"Stream '{stream}' is assigned more than once for "
                    f"mesh: {mesh_name}."
                )
                continue

            # just store filename; the full path will be resolved later
            # must store something to test for duplicates
            streams_files[stream] = str(file_name)

    missing_streams = REQUIRED_STREAMS - set(streams_files)
    if missing_streams:
        errors.append(
            f"Missing required IOStream(s) "
            f"{', '.join(sorted(missing_streams))} for mesh: {mesh_name}."
        )

    return errors


def _validate_config_overrides_entry(
    config_overrides: YamlMapping, input_files: YamlMapping, mesh_name: str
) -> list[str]:
    """
    Validate the overrides of a single mesh in config_overrides.

    Parameters:
    -----------
    config_overrides : dict[str, Any]
        Parsed content of
        ``cime_config/omega_buildnml/data/config_overrides.yaml``
    input_files : dict[str, Any]
        Parsed content of ``cime_config/omega_buildnml/data/input_files.yaml``,
        which defines the meshes Omega supports.
    mesh_name : str
        The name of the mesh to validate.

    Returns:
    --------
    list[str]
        Error messages describing any problems found. Empty when the entry is
        valid.
    """
    meshes: YamlMapping = config_overrides["meshes"]
    supported_meshes: YamlMapping = input_files.get("meshes", {})
    overrides: YamlMapping = meshes[mesh_name]

    errors: list[str] = []

    if not isinstance(overrides, dict) or not overrides:
        errors.append(f"Mesh: {mesh_name} is empty or is not a mapping.")

    if mesh_name not in supported_meshes:
        errors.append(
            f"Unsupported mesh: {mesh_name}. Could not find entry in "
            f"`{INPUT_FILES_PATH}`."
        )

    return errors
