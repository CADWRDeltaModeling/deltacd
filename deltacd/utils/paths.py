"""Path resolution helpers shared across DETAW and DCD."""

from __future__ import annotations

from pathlib import Path


def convert_to_absolute_path(filename: str | Path, dir_input_base: str | Path) -> str:
    """Resolve a path against ``dir_input_base`` when it is relative.

    Parameters
    ----------
    filename : str or pathlib.Path
        Input or output file path from a YAML configuration.
    dir_input_base : str or pathlib.Path
        Directory containing the main YAML configuration file.

    Returns
    -------
    str
        Resolved path string.
    """
    path = Path(filename)
    if not path.is_absolute():
        path = Path(dir_input_base) / path
    return str(path)
