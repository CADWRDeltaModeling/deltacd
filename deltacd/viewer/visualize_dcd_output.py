"""Interactive Panel/HoloViews viewer for DeltaCD output NetCDF files.

The app displays two tabs:
  - Area Output  – time-series plots for DETAW subarea variables.
  - Node Output  – time-series plots for DSM2 node variables.

All input paths are read from a YAML config file (default: config.yaml next to
this script). Run with:
    panel serve visualize_dcd_output.py --show
    panel serve visualize_dcd_output.py --show --args --config path/to/config.yaml
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import json

import yaml
import holoviews as hv
import numpy as np
import pandas as pd
import panel as pn
import xarray as xr


hv.extension("bokeh")
pn.extension(sizing_mode="stretch_width")

# To run the app, use the command with a custom config file:
# panel serve visualize_dcd_output.py --show --args --config path/to/config.yaml


# Shared layout constants used by both the area and node tabs.
PLOT_HEIGHT = 480      # pixel height of HoloViews plots
CONTROLS_WIDTH = 440   # pixel width of the left-hand controls column

# These globals are populated at module load time once the config has been
# parsed (see the bootstrap block at the bottom of the file).  They are
# declared here so that the function definitions below can reference them.
FACTOR_PATH_BY_VARIABLE: dict[str, Path] = {}
SUBAREAS_GEOJSON: dict | None = None
NODES_GEOJSON: dict | None = None


def _utm_to_latlon(easting: float, northing: float, zone: int = 10, northern: bool = True) -> tuple[float, float]:
    """Convert UTM Zone N coordinates to WGS84 (lat, lon) in decimal degrees."""
    # WGS84 ellipsoid parameters
    a = 6378137.0            # semi-major axis (m)
    f = 1 / 298.257223563    # flattening
    b = a * (1 - f)          # semi-minor axis
    e2 = 1 - (b / a) ** 2   # first eccentricity squared
    ep2 = (a / b) ** 2 - 1  # second eccentricity squared
    k0 = 0.9996              # central scale factor

    # Remove false easting; adjust for southern hemisphere if needed
    x = easting - 500000.0
    y = northing if northern else northing - 10000000.0
    lon0 = math.radians((zone - 1) * 6 - 180 + 3)  # central meridian of zone

    n = (a - b) / (a + b)
    M = y / k0
    mu = M / (a * (1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256))

    # Iteratively solve for the footpoint latitude (phi1) using a series expansion
    phi1 = mu
    for _ in range(5):
        phi1 = (
            mu
            + (3 * n / 2 - 27 * n ** 3 / 32) * math.sin(2 * phi1)
            + (21 * n ** 2 / 16 - 55 * n ** 4 / 32) * math.sin(4 * phi1)
            + (151 * n ** 3 / 96) * math.sin(6 * phi1)
            + (1097 * n ** 4 / 512) * math.sin(8 * phi1)
        )

    # Intermediate quantities at the footpoint latitude
    sin_p = math.sin(phi1)
    cos_p = math.cos(phi1)
    tan_p = math.tan(phi1)
    N1 = a / math.sqrt(1 - e2 * sin_p ** 2)          # radius of curvature in prime vertical
    T1 = tan_p ** 2
    C1 = ep2 * cos_p ** 2
    R1 = a * (1 - e2) / (1 - e2 * sin_p ** 2) ** 1.5  # radius of curvature in meridian
    D = x / (N1 * k0)

    lat = phi1 - (N1 * tan_p / R1) * (
        D ** 2 / 2
        - (5 + 3 * T1 + 10 * C1 - 4 * C1 ** 2 - 9 * ep2) * D ** 4 / 24
        + (61 + 90 * T1 + 298 * C1 + 45 * T1 ** 2 - 252 * ep2 - 3 * C1 ** 2) * D ** 6 / 720
    )
    lon = lon0 + (
        D
        - (1 + 2 * T1 + C1) * D ** 3 / 6
        + (5 - 2 * C1 + 28 * T1 - 3 * C1 ** 2 + 8 * ep2 + 24 * T1 ** 2) * D ** 5 / 120
    ) / cos_p

    return math.degrees(lat), math.degrees(lon)


def _reproject_coords(coords: list, zone: int = 10, northern: bool = True) -> list:
    """Recursively reproject a GeoJSON coordinate array from UTM to WGS84."""
    if not coords:
        return coords
    if isinstance(coords[0], (int, float)):
        lat, lon = _utm_to_latlon(coords[0], coords[1], zone, northern)
        return [lon, lat]
    return [_reproject_coords(c, zone, northern) for c in coords]


def _reproject_geojson(geojson: dict, zone: int = 10, northern: bool = True) -> dict:
    """Return a new GeoJSON FeatureCollection with coordinates in WGS84."""
    import copy
    result = copy.deepcopy(geojson)
    result.pop("crs", None)  # remove projected CRS declaration
    for feature in result.get("features", []):
        geom = feature.get("geometry")
        if geom and "coordinates" in geom:
            geom["coordinates"] = _reproject_coords(geom["coordinates"], zone, northern)
    return result


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Interactive viewer for DeltaCD output.")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).resolve().parent / "config.yaml",
        help="Path to the YAML configuration file.",
    )
    args, _ = parser.parse_known_args()
    return args


def _load_config(config_path: Path) -> dict:
    """Load and return the YAML configuration as a plain dict."""
    with config_path.open(encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _resolve_path(base: Path, value: str) -> Path:
    """Resolve a path string relative to *base* unless it is already absolute."""
    p = Path(value)
    return p if p.is_absolute() else (base / p).resolve()


# ── Shared utilities ──────────────────────────────────────────────────────────


def to_string_list(values: np.ndarray) -> list[str]:
    """Convert a numpy array of coordinate labels to plain Python strings.

    NetCDF4 string variables are sometimes read back as byte strings; this
    handles that case transparently alongside regular str/int/float labels.
    """
    result: list[str] = []
    for value in np.asarray(values).tolist():
        if isinstance(value, bytes):
            result.append(value.decode("utf-8"))
        else:
            result.append(str(value))
    return result


def to_timestamp_index(ds: xr.Dataset) -> pd.DatetimeIndex:
    """Return the dataset's time coordinate as a pandas DatetimeIndex.

    Three-level fallback handles the variety of time encodings produced by
    different versions of xarray/cftime:
      1. Use the xarray index directly (covers standard numpy datetime64).
      2. Call to_datetimeindex() for CFTime indexes.
      3. Parse raw values via pd.to_datetime or manual strftime as a last resort.
    """
    try:
        index = ds.indexes["time"]
        if hasattr(index, "to_datetimeindex"):  # CFTimeIndex (xarray < 2024)
            return pd.DatetimeIndex(index.to_datetimeindex())
        return pd.DatetimeIndex(index)
    except Exception:
        values = ds["time"].values

    try:
        return pd.to_datetime(values)
    except Exception:
        # Last resort: objects that have strftime (e.g. cftime.datetime)
        converted = []
        for value in values:
            if hasattr(value, "strftime"):
                converted.append(pd.Timestamp(value.strftime("%Y-%m-%d")))
            else:
                converted.append(pd.Timestamp(str(value)))
        return pd.DatetimeIndex(converted)


def sort_numeric_strings(values: list[str]) -> list[str]:
    """Sort a list of strings, placing purely numeric strings first in numeric order."""
    return sorted(
        values,
        key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value),
    )


def aggregate_frame(frame: pd.DataFrame, aggregation: str) -> pd.DataFrame:
    if aggregation == "Daily":
        return frame
    if aggregation == "Monthly mean":
        return frame.resample("MS").mean()
    if aggregation == "Calendar-year mean":
        return frame.resample("YS").mean()
    if aggregation == "Water-year mean":
        return frame.resample("YS-OCT").mean()
    raise ValueError(f"Unsupported aggregation: {aggregation}")


# ── Map tab ───────────────────────────────────────────────────────────────────


def make_map(
    selected_areas: list[str],
    selected_nodes: list[str],
    subareas_geojson: dict,
    nodes_geojson: dict,
) -> pn.pane.HTML | pn.pane.Markdown:
    try:
        import folium
    except ImportError:
        return pn.pane.Markdown(
            "Install `folium` to enable the map: `pip install folium`"
        )

    # Centre the map on the Sacramento–San Joaquin Delta.
    m = folium.Map(location=[38.05, -121.55], zoom_start=10, tiles="cartodb positron")

    # Convert to sets for O(1) membership checks inside the style callback.
    selected_area_set = {str(a) for a in (selected_areas or [])}
    selected_node_set = {str(n) for n in (selected_nodes or [])}

    def subarea_style(feature: dict) -> dict:
        new_sub = str(feature["properties"]["NEW_SUB"])
        if new_sub in selected_area_set:
            return {"fillColor": "#0b6e4f", "color": "#0b6e4f", "weight": 2, "fillOpacity": 0.6}
        return {"fillColor": "#aad9c3", "color": "#555", "weight": 1, "fillOpacity": 0.2}

    folium.GeoJson(
        subareas_geojson,
        style_function=subarea_style,
        tooltip=folium.GeoJsonTooltip(
            fields=["NEW_SUB", "SUB_NAME"],
            aliases=["Subarea ID:", "Name:"],
        ),
        name="Subareas",
    ).add_to(m)

    for feature in nodes_geojson["features"]:
        node_id = str(feature["properties"]["id"])
        lon, lat = feature["geometry"]["coordinates"]
        is_selected = node_id in selected_node_set
        folium.CircleMarker(
            location=[lat, lon],
            radius=7 if is_selected else 4,
            color="#cc0000" if is_selected else "#555",
            fill=True,
            fill_color="#ee3333" if is_selected else "#999",
            fill_opacity=0.9 if is_selected else 0.5,
            tooltip=f"Node: {node_id}",
        ).add_to(m)

    folium.LayerControl().add_to(m)
    html_str = m._repr_html_()
    # folium wraps the map in a percentage-height div; override it with a fixed pixel height
    html_str = html_str.replace(
        "height:0;padding-bottom:60%;",
        f"height:{PLOT_HEIGHT}px;",
    )
    return pn.pane.HTML(html_str, height=PLOT_HEIGHT, sizing_mode="stretch_width")


# ── Area tab ──────────────────────────────────────────────────────────────────


def area_open_dataset(path: Path) -> xr.Dataset:
    """Open the area NetCDF and ensure area coordinate labels are plain strings."""
    ds = xr.open_dataset(path, decode_times=True)
    # Normalise any of the recognised area coordinate names to plain Python
    # strings so that label-based selection works consistently.
    for coord_name in ("area", "area_id", "subarea", "subarea_id"):
        if coord_name in ds.variables:
            ds = ds.assign_coords({coord_name: to_string_list(ds[coord_name].values)})
    return ds


def area_detect_dimension(ds: xr.Dataset) -> str:
    """Return the name of the spatial dimension in the area dataset.

    Checks a list of known names first; falls back to inspecting data
    variables for any non-time dimension when none of the known names match.
    """
    for candidate in ("area", "area_id", "subarea", "subarea_id"):
        if candidate in ds.dims:
            return candidate
    # Generic fallback: take the first non-time dimension found.
    for _name, data_var in ds.data_vars.items():
        dims = tuple(dim for dim in data_var.dims if dim != "time")
        if dims:
            return dims[0]
    raise ValueError("Could not detect an area dimension in the dataset.")


def area_dataset_summary(
    ds: xr.Dataset, time_index: pd.DatetimeIndex, area_dim: str
) -> dict[str, object]:
    variables = [
        name for name, data_var in ds.data_vars.items() if data_var.dims == ("time", area_dim)
    ]
    return {
        "variables": variables,
        "area_dim": area_dim,
        "area_count": ds.sizes.get(area_dim, 0),
        "time_count": ds.sizes.get("time", 0),
        "start": time_index.min(),
        "end": time_index.max(),
    }


def area_to_dataframe(
    ds: xr.Dataset,
    variable: str,
    area_dim: str,
    areas: list[str],
    time_index: pd.DatetimeIndex,
) -> pd.DataFrame:
    """Extract *variable* for the requested *areas* as a (time × area) DataFrame."""
    # Transpose to guarantee time is always the row axis regardless of storage order.
    data_array = ds[variable].sel({area_dim: areas}).transpose("time", area_dim)
    frame = data_array.to_pandas()
    # A single-area selection returns a Series; promote it to a one-column DataFrame.
    if isinstance(frame, pd.Series):
        frame = frame.to_frame(name=areas[0])
    frame.index = time_index
    frame.columns = to_string_list(np.asarray(frame.columns))
    return frame


def area_make_plot(
    ds: xr.Dataset,
    time_index: pd.DatetimeIndex,
    area_dim: str,
    variable: str,
    areas: list[str],
    date_start: object,
    date_end: object,
    aggregation: str,
    rolling_window: int,
) -> hv.Overlay | hv.Curve:
    if not areas:
        return hv.Curve([]).opts(
            height=PLOT_HEIGHT, responsive=True, title="Select at least one subarea"
        )
    frame = area_to_dataframe(ds, variable, area_dim, areas, time_index)
    # Clip to the user-selected date range; fall back to full extent if unset.
    start = pd.Timestamp(date_start) if date_start is not None else time_index.min()
    end = pd.Timestamp(date_end) if date_end is not None else time_index.max()
    filtered = aggregate_frame(frame.loc[start:end], aggregation)
    if rolling_window > 1:  # apply optional smoothing after aggregation
        filtered = filtered.rolling(window=rolling_window, min_periods=1).mean()
    if filtered.empty:
        return hv.Curve([]).opts(
            height=PLOT_HEIGHT, responsive=True, title="No data in the selected date range"
        )
    # Build one Curve per selected area and overlay them.
    curves = [
        hv.Curve((filtered.index, filtered[column]), kdims="time", vdims=variable, label=column)
        for column in filtered.columns
    ]
    return hv.Overlay(curves).opts(
        height=PLOT_HEIGHT,
        legend_position="right",
        responsive=True,
        show_grid=True,
        tools=["hover"],
        title=f"{variable.replace('_', ' ').title()} by subarea",
        xlabel="Time",
        ylabel=variable,
    )


def area_make_stats_table(
    ds: xr.Dataset,
    time_index: pd.DatetimeIndex,
    area_dim: str,
    variable: str,
    areas: list[str],
    date_start: object,
    date_end: object,
    aggregation: str,
) -> pn.viewable.Viewable:
    if not areas:
        return pn.pane.Markdown("Select one or more subareas to see summary statistics.")
    frame = area_to_dataframe(ds, variable, area_dim, areas, time_index)
    start = pd.Timestamp(date_start) if date_start is not None else time_index.min()
    end = pd.Timestamp(date_end) if date_end is not None else time_index.max()
    filtered = aggregate_frame(frame.loc[start:end], aggregation)
    if filtered.empty:
        return pn.pane.Markdown("No values are available for the current selection.")
    summary = pd.DataFrame(
        {
            "mean": filtered.mean(),
            "min": filtered.min(),
            "max": filtered.max(),
            "latest": filtered.iloc[-1],
        }
    ).round(3)
    summary.index.name = area_dim
    return pn.pane.DataFrame(summary, sizing_mode="stretch_width", height=280)


def build_area_tab(path: Path) -> pn.Row:
    if not path.exists():
        return pn.Row(
            pn.pane.Markdown(f"### Dataset not found\n\nExpected file at `{path}`."),
            sizing_mode="stretch_width",
        )

    ds = area_open_dataset(path)
    time_index = to_timestamp_index(ds)
    area_dim = area_detect_dimension(ds)
    summary = area_dataset_summary(ds, time_index, area_dim)

    area_options = to_string_list(ds[area_dim].values)
    if area_options and all(option.isdigit() for option in area_options):
        area_options = sort_numeric_strings(area_options)

    available_variables = summary["variables"]
    if not available_variables:
        return pn.Row(
            pn.pane.Markdown(
                f"### No variables found\n\nNo time-by-{area_dim} variables in `{path.name}`."
            ),
        )

    variable = pn.widgets.Select(
        name="Variable", options=available_variables, value=available_variables[0]
    )
    areas = pn.widgets.MultiChoice(
        name="Subareas",
        options=area_options,
        value=area_options[: min(4, len(area_options))],
        delete_button=True,
        placeholder="Choose one or more subareas",
    )
    date_start = pn.widgets.DatePicker(
        name="Start date",
        value=summary["start"].date(),
    )
    date_end = pn.widgets.DatePicker(
        name="End date",
        value=summary["end"].date(),
    )
    aggregation = pn.widgets.Select(
        name="Aggregation",
        options=["Daily", "Monthly mean", "Calendar-year mean", "Water-year mean"],
        value="Monthly mean",
    )
    rolling_window = pn.widgets.IntSlider(name="Rolling window", start=1, end=90, step=1, value=1)

    notes = pn.pane.Markdown(
        "\n".join(
            [
                f"- File: `{path.name}`",
                f"- Subarea dimension: `{area_dim}`",
                f"- Subareas: `{summary['area_count']}`",
                f"- Time steps: `{summary['time_count']}`",
                f"- Variables: `{', '.join(available_variables)}`",
                f"- Period: `{summary['start'].date()}` to `{summary['end'].date()}`",
            ]
        ),
        sizing_mode="stretch_width",
    )
    if SUBAREAS_GEOJSON is not None and NODES_GEOJSON is not None:
        def _area_map(sel_areas: list[str]) -> object:
            return make_map(sel_areas, [], SUBAREAS_GEOJSON, NODES_GEOJSON)
        map_pane: pn.viewable.Viewable = pn.panel(pn.bind(_area_map, areas))
    else:
        map_pane = pn.pane.Markdown("GeoJSON files not found.")

    info_tabs = pn.Tabs(
        ("Map", map_pane),
        ("Dataset Overview", notes),
        active=0,
        tabs_location="above",
        sizing_mode="stretch_width",
    )

    date_row = pn.Row(date_start, date_end, sizing_mode="stretch_width")

    plot = pn.bind(
        area_make_plot, ds, time_index, area_dim, variable, areas, date_start, date_end, aggregation, rolling_window
    )
    stats = pn.bind(area_make_stats_table, ds, time_index, area_dim, variable, areas, date_start, date_end, aggregation)

    controls = pn.Column(
        info_tabs, variable, areas, date_row, aggregation, rolling_window, width=CONTROLS_WIDTH
    )
    main = pn.Column(
        pn.panel(plot),
        pn.pane.Markdown("### Summary statistics"),
        pn.panel(stats),
        sizing_mode="stretch_width",
    )
    return pn.Row(controls, main, sizing_mode="stretch_width")


# ── Node tab ──────────────────────────────────────────────────────────────────


def node_open_dataset(path: Path) -> xr.Dataset:
    """Open the node NetCDF and ensure the 'node' coordinate is plain strings."""
    ds = xr.open_dataset(path, decode_times=True)
    # Normalise node labels to plain strings whether 'node' is currently a
    # coordinate or just a variable (both cases require the same treatment).
    if "node" not in ds.coords and "node" in ds.variables:
        ds = ds.assign_coords(node=("node", to_string_list(ds["node"].values)))
    elif "node" in ds.coords:
        ds = ds.assign_coords(node=("node", to_string_list(ds["node"].values)))
    return ds


def node_dataset_summary(ds: xr.Dataset, time_index: pd.DatetimeIndex) -> dict[str, object]:
    # The crop variable is a scalar or single-element array; extract its value
    # for display in the dataset overview panel.
    crop_value = "Unknown"
    if "crop" in ds.variables:
        crop_data = np.asarray(ds["crop"].values)
        if crop_data.size:
            crop_value = str(crop_data.reshape(-1)[0])
    variables = [
        name for name, data_var in ds.data_vars.items() if data_var.dims == ("time", "node")
    ]
    return {
        "crop": crop_value,
        "variables": variables,
        "node_count": ds.sizes.get("node", 0),
        "time_count": ds.sizes.get("time", 0),
        "start": time_index.min(),
        "end": time_index.max(),
    }


def node_to_dataframe(
    ds: xr.Dataset, variable: str, nodes: list[str], time_index: pd.DatetimeIndex
) -> pd.DataFrame:
    data_array = ds[variable].sel(node=nodes).transpose("time", "node")
    frame = data_array.to_pandas()
    if isinstance(frame, pd.Series):
        frame = frame.to_frame(name=nodes[0])
    frame.index = time_index
    frame.columns = to_string_list(np.asarray(frame.columns))
    return frame


def node_make_plot(
    ds: xr.Dataset,
    time_index: pd.DatetimeIndex,
    variable: str,
    nodes: list[str],
    date_start: object,
    date_end: object,
    aggregation: str,
    rolling_window: int,
) -> hv.Overlay | hv.Curve:
    if not nodes:
        return hv.Curve([]).opts(
            height=PLOT_HEIGHT, responsive=True, title="Select at least one node"
        )
    frame = node_to_dataframe(ds, variable, nodes, time_index)
    start = pd.Timestamp(date_start) if date_start is not None else time_index.min()
    end = pd.Timestamp(date_end) if date_end is not None else time_index.max()
    filtered = aggregate_frame(frame.loc[start:end], aggregation)
    if rolling_window > 1:
        filtered = filtered.rolling(window=rolling_window, min_periods=1).mean()
    if filtered.empty:
        return hv.Curve([]).opts(
            height=PLOT_HEIGHT, responsive=True, title="No data in the selected date range"
        )
    curves = [
        hv.Curve((filtered.index, filtered[column]), kdims="time", vdims=variable, label=column)
        for column in filtered.columns
    ]
    return hv.Overlay(curves).opts(
        height=PLOT_HEIGHT,
        legend_position="right",
        responsive=True,
        show_grid=True,
        tools=["hover"],
        title=f"{variable.replace('_', ' ').title()} by node",
        xlabel="Time",
        ylabel=f"{variable} (cfs)",
    )


def node_make_stats_table(
    ds: xr.Dataset,
    time_index: pd.DatetimeIndex,
    variable: str,
    nodes: list[str],
    date_start: object,
    date_end: object,
    aggregation: str,
) -> pn.viewable.Viewable:
    if not nodes:
        return pn.pane.Markdown("Select one or more nodes to see summary statistics.")
    frame = node_to_dataframe(ds, variable, nodes, time_index)
    start = pd.Timestamp(date_start) if date_start is not None else time_index.min()
    end = pd.Timestamp(date_end) if date_end is not None else time_index.max()
    filtered = aggregate_frame(frame.loc[start:end], aggregation)
    if filtered.empty:
        return pn.pane.Markdown("No values are available for the current selection.")
    summary = pd.DataFrame(
        {
            "mean": filtered.mean(),
            "min": filtered.min(),
            "max": filtered.max(),
            "latest": filtered.iloc[-1],
        }
    ).round(3)
    summary.index.name = "node"
    return pn.pane.DataFrame(summary, sizing_mode="stretch_width", height=280)


def load_factor_table(path: Path) -> pd.DataFrame:
    """Load a diversion/drainage factor CSV and return unique (area_id, node) pairs.

    The factor tables map each DETAW subarea to one or more DSM2 nodes so that
    selecting a subarea in the node tab automatically highlights the relevant
    nodes.  Only the two lookup columns are kept.
    """
    if not path.exists():
        return pd.DataFrame(columns=["area_id", "node"])
    frame = pd.read_csv(path, usecols=["area_id", "node"])
    frame = frame.dropna(subset=["area_id", "node"]).copy()
    frame["area_id"] = frame["area_id"].astype(str)
    frame["node"] = frame["node"].astype(str)
    return frame.drop_duplicates().sort_values(["area_id", "node"])


def build_node_tab(path: Path) -> pn.Row:
    if not path.exists():
        return pn.Row(
            pn.pane.Markdown(f"### Dataset not found\n\nExpected file at `{path}`."),
            sizing_mode="stretch_width",
        )

    ds = node_open_dataset(path)
    time_index = to_timestamp_index(ds)
    summary = node_dataset_summary(ds, time_index)
    node_options = to_string_list(ds["node"].values)
    node_option_set = set(node_options)
    available_variables = summary["variables"]
    # Pre-load all factor tables once so callbacks don't re-read files on each interaction.
    factor_tables = {
        variable_name: load_factor_table(factor_path)
        for variable_name, factor_path in FACTOR_PATH_BY_VARIABLE.items()
    }

    variable = pn.widgets.Select(
        name="Variable", options=available_variables, value=available_variables[0]
    )
    subareas = pn.widgets.MultiChoice(
        name="Subareas",
        options=[],
        value=[],
        delete_button=True,
        placeholder="Choose one or more subareas",
    )
    nodes = pn.widgets.MultiChoice(
        name="Nodes",
        options=node_options,
        value=node_options[: min(4, len(node_options))],
        delete_button=True,
        placeholder="Choose one or more nodes",
    )

    def resolve_factor_table(variable_name: str) -> pd.DataFrame:
        """Return the factor table for *variable_name*, filtered to nodes present in the dataset."""
        frame = factor_tables.get(variable_name, pd.DataFrame(columns=["area_id", "node"]))
        if frame.empty:
            return frame
        # Drop any rows whose node doesn't exist in the loaded dataset.
        return frame[frame["node"].isin(node_option_set)]

    def update_subarea_and_node_options(*_: object) -> None:
        """Rebuild subarea and node widget options when the selected variable changes.

        Different variables use different factor tables (e.g. diversion vs
        drainage), so the set of mappable subareas and nodes must be refreshed.
        Previously selected values are retained where they remain valid.
        """
        frame = resolve_factor_table(variable.value)
        previous_nodes = list(nodes.value)
        previous_subareas = list(subareas.value)

        if frame.empty:
            subareas.options = []
            subareas.value = []
            nodes.options = node_options
            nodes.value = [node for node in previous_nodes if node in node_option_set]
            return

        available_subareas = sort_numeric_strings(frame["area_id"].drop_duplicates().tolist())
        available_nodes = sort_numeric_strings(frame["node"].drop_duplicates().tolist())

        subareas.options = available_subareas
        retained_subareas = [area for area in previous_subareas if area in available_subareas]
        subareas.value = retained_subareas

        nodes.options = available_nodes
        if retained_subareas:
            selected_nodes = sort_numeric_strings(
                frame.loc[frame["area_id"].isin(retained_subareas), "node"]
                .drop_duplicates()
                .tolist()
            )
            nodes.value = selected_nodes
        else:
            nodes.value = [node for node in previous_nodes if node in set(available_nodes)]

    def update_nodes_from_subareas(event: object) -> None:
        """When the user picks subareas, auto-select the corresponding nodes."""
        del event
        selected_subareas = list(subareas.value)
        if not selected_subareas:
            return
        frame = resolve_factor_table(variable.value)
        if frame.empty:
            return
        selected_nodes = sort_numeric_strings(
            frame.loc[frame["area_id"].isin(selected_subareas), "node"].drop_duplicates().tolist()
        )
        nodes.value = selected_nodes

    variable.param.watch(update_subarea_and_node_options, "value")
    subareas.param.watch(update_nodes_from_subareas, "value")
    update_subarea_and_node_options()

    date_start = pn.widgets.DatePicker(
        name="Start date",
        value=summary["start"].date(),
    )
    date_end = pn.widgets.DatePicker(
        name="End date",
        value=summary["end"].date(),
    )
    aggregation = pn.widgets.Select(
        name="Aggregation",
        options=["Daily", "Monthly mean", "Calendar-year mean", "Water-year mean"],
        value="Monthly mean",
    )
    rolling_window = pn.widgets.IntSlider(name="Rolling window", start=1, end=90, step=1, value=1)

    notes = pn.pane.Markdown(
        "\n".join(
            [
                f"- File: `{path.name}`",
                f"- Crop: `{summary['crop']}`",
                f"- Nodes: `{summary['node_count']}`",
                f"- Time steps: `{summary['time_count']}`",
                f"- Variables: `{', '.join(available_variables)}`",
                "- Subarea source: diversion/seepage/depletion_from_waterbody use diversion factors, drainage uses drainage factors",
                f"- Period: `{summary['start'].date()}` to `{summary['end'].date()}`",
            ]
        ),
        sizing_mode="stretch_width",
    )
    if SUBAREAS_GEOJSON is not None and NODES_GEOJSON is not None:
        def _node_map(sel_areas: list[str], sel_nodes: list[str]) -> object:
            return make_map(sel_areas, sel_nodes, SUBAREAS_GEOJSON, NODES_GEOJSON)
        map_pane: pn.viewable.Viewable = pn.panel(pn.bind(_node_map, subareas, nodes))
    else:
        map_pane = pn.pane.Markdown("GeoJSON files not found.")

    info_tabs = pn.Tabs(
        ("Map", map_pane),
        ("Dataset Overview", notes),
        active=0,
        tabs_location="above",
        sizing_mode="stretch_width",
    )

    selector_row = pn.Row(
        pn.Column(subareas, sizing_mode="stretch_width"),
        pn.Column(nodes, sizing_mode="stretch_width"),
        sizing_mode="stretch_width",
    )
    date_row = pn.Row(date_start, date_end, sizing_mode="stretch_width")

    plot = pn.bind(
        node_make_plot, ds, time_index, variable, nodes, date_start, date_end, aggregation, rolling_window
    )
    stats = pn.bind(node_make_stats_table, ds, time_index, variable, nodes, date_start, date_end, aggregation)

    controls = pn.Column(
        info_tabs, variable, selector_row, date_row, aggregation, rolling_window, width=CONTROLS_WIDTH
    )
    main = pn.Column(
        pn.panel(plot),
        pn.pane.Markdown("### Summary statistics"),
        pn.panel(stats),
        sizing_mode="stretch_width",
    )
    return pn.Row(controls, main, sizing_mode="stretch_width")


# ── Combined app ──────────────────────────────────────────────────────────────


def build_app(area_path: Path, node_path: Path) -> pn.template.FastListTemplate:
    area_content = build_area_tab(area_path)
    node_content = build_node_tab(node_path)

    # Inject custom CSS to style the top-level area/node tabs to match the
    # DWR green colour scheme.
    tab_styles = [
        """
        .bk-tab {
            font-size: 14px;
            font-weight: 600;
            padding: 10px 28px;
            background-color: #e8f5ef;
            color: #0b6e4f;
            border: 2px solid #0b6e4f;
            border-bottom: none;
            border-radius: 6px 6px 0 0;
            margin-right: 6px;
            cursor: pointer;
            transition: background-color 0.15s, color 0.15s;
        }
        .bk-tab:hover {
            background-color: #c2e8d5;
        }
        .bk-tab.bk-active {
            background-color: #0b6e4f;
            color: #ffffff;
        }
        .bk-tabs-header {
            margin-bottom: 0;
        }
        """
    ]

    tabs = pn.Tabs(
        ("Area Output", area_content),
        ("Node Output", node_content),
        sizing_mode="stretch_width",
        tabs_location="above",
        active=0,
        stylesheets=tab_styles,
    )

    return pn.template.FastListTemplate(
        title="DeltaCD Output Viewer",
        main=[tabs],
        accent_base_color="#0b6e4f",
        header_background="#0b6e4f",
    )


# ── Bootstrap ─────────────────────────────────────────────────────────────────
# This block runs at module load time (both via `panel serve` and direct
# execution).  It:
#   1. Parses the --config argument.
#   2. Loads the YAML config and resolves all paths relative to the config file.
#   3. Populates the module-level globals used by the tab builders.
#   4. Constructs the Panel app and registers it as servable.

ARGS = parse_args()
_config = _load_config(ARGS.config)
_base = ARGS.config.resolve().parent  # all relative paths in the config are relative to this dir

# Map each variable type to the CSV that holds its area→node factor mapping.
# Seepage and depletion share the diversion factor file.
FACTOR_PATH_BY_VARIABLE = {
    "diversion": _resolve_path(_base, _config["diversion_factors_path"]),
    "drainage": _resolve_path(_base, _config["drainage_factors_path"]),
    "seepage": _resolve_path(_base, _config["diversion_factors_path"]),
    "depletion_from_waterbody": _resolve_path(_base, _config["diversion_factors_path"]),
}

# Load and reproject GeoJSON layers from UTM Zone 10N → WGS84 for folium.
# If either file is missing the map tab will show a fallback message.
_subareas_geojson_path = _resolve_path(_base, _config["subareas_geojson_path"])
_nodes_geojson_path = _resolve_path(_base, _config["nodes_geojson_path"])
SUBAREAS_GEOJSON = (
    _reproject_geojson(json.loads(_subareas_geojson_path.read_text(encoding="utf-8")))
    if _subareas_geojson_path.exists()
    else None
)
NODES_GEOJSON = (
    _reproject_geojson(json.loads(_nodes_geojson_path.read_text(encoding="utf-8")))
    if _nodes_geojson_path.exists()
    else None
)

APP = build_app(
    _resolve_path(_base, _config["area_path"]),
    _resolve_path(_base, _config["node_path"]),
)
APP.servable()


if __name__ == "__main__":
    pn.serve(APP, show=True, title="DeltaCD Output Viewer")
    