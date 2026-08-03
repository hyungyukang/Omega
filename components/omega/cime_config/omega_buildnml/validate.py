from typing import Optional

from ._types import YamlMapping

CONFIG_PATH = (
    "components/omega/cime_config/omega_buildnml/data/input_files.yaml"
)

ERR_SUFFIX = f"Please check your setting in `{CONFIG_PATH}`"

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
        "GlobalStats",
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
        err_msg = f"`{CONFIG_PATH}` is empty or is not a mapping."
        raise ValueError(err_msg)

    unknown_keys = sorted(set(input_files) - {"meshes"})
    if unknown_keys:
        _raise([f"Unknown top-level key(s): {', '.join(unknown_keys)}."])

    meshes: YamlMapping = input_files.get("meshes", {})
    if not isinstance(meshes, dict) or not meshes:
        _raise(["`meshes` is missing, empty, or is not a mapping."])

    if mesh_name is None:
        errors = []
        for name in meshes:
            errors.extend(_validate_input_files_entry(input_files, name))
        _raise(errors)
        return input_files

    if mesh_name not in meshes:
        err_msg = (
            f"Unsupported OCN_GRID for Omega: {mesh_name}. \n"
            f"Could not find entry in `{CONFIG_PATH}`"
        )
        raise ValueError(err_msg)

    _raise(_validate_input_files_entry(input_files, mesh_name))

    return input_files


def _raise(errors: list[str]) -> None:
    """
    Raise a single ``ValueError`` describing all accumulated errors.

    Does nothing when ``errors`` is empty.

    Parameters:
    -----------
    errors : list[str]
        Error messages collected during validation.

    Raises:
    -------
    ValueError
        If ``errors`` is non-empty.
    """
    if not errors:
        return

    details = "\n".join(f"  - {error}" for error in errors)
    err_msg = (
        f"Invalid Omega input file configuration:\n{details}\n{ERR_SUFFIX}"
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
