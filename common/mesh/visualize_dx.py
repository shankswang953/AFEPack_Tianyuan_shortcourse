#!/usr/bin/env python3
"""Render a two-dimensional AFEPack/OpenDX triangular mesh.

Purpose:
    Create a PNG, SVG, or PDF preview from an AFEPack-generated ``.dx`` file.
    Files containing only positions/connections are drawn as meshes. Files
    containing a nodal or element scalar array are drawn as colored fields.

Usage:
    python3 visualize_dx.py INPUT.dx [-o OUTPUT.png]
    python3 visualize_dx.py INPUT.dx --mesh-only --show-nodes

Arguments:
    INPUT.dx              OpenDX file to read.
    -o, --output PATH     Figure path. Defaults to INPUT with a ``.png`` suffix.
    --mesh-only           Ignore scalar data and draw only the triangulation.
    --show-nodes          Draw a marker at every mesh node.
    --no-edges            Hide triangle edges in scalar-field plots.
    --cmap NAME           Matplotlib color map for scalar data.
    --title TEXT          Custom figure title.
    --dpi INTEGER         Raster output resolution.

Output:
    The parent directory of OUTPUT is created automatically. The output format
    is selected from the filename suffix, for example ``.png`` or ``.svg``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
import os
from pathlib import Path
import shlex
import tempfile
from typing import Sequence


# Use a writable, machine-independent cache location before importing
# Matplotlib. This is useful on clusters and in read-only teaching directories.
_CACHE_ROOT = Path(tempfile.gettempdir()) / "afepack-dx-plot-cache"
_CACHE_ROOT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


@dataclass(frozen=True)
class DXArray:
    """One inline OpenDX array."""

    object_id: str
    value_type: str
    rank: int
    shape: tuple[int, ...]
    item_count: int
    rows: list[list[float | int]]


@dataclass(frozen=True)
class DXMesh:
    """Triangular geometry and optional scalar data extracted from OpenDX."""

    points: list[list[float]]
    triangles: list[list[int]]
    scalar_values: list[float] | None
    scalar_location: str | None


def _parse_array_header(line: str, source: Path) -> tuple[str, str, int, tuple[int, ...], int] | None:
    """Return the metadata from an inline OpenDX array header."""

    try:
        tokens = shlex.split(line, comments=True)
    except ValueError as error:
        raise ValueError(f"{source}: invalid quoted text in OpenDX header") from error

    if len(tokens) < 8 or tokens[0] != "object" or "class" not in tokens:
        return None

    class_index = tokens.index("class")
    if class_index + 1 >= len(tokens) or tokens[class_index + 1] != "array":
        return None
    if len(tokens) < 2:
        return None

    try:
        value_type = tokens[tokens.index("type") + 1]
        rank = int(tokens[tokens.index("rank") + 1])
        try:
            item_index = tokens.index("item")
        except ValueError:
            item_index = tokens.index("items")
        item_count = int(tokens[item_index + 1])
        data_index = tokens.index("data")
    except (ValueError, IndexError) as error:
        raise ValueError(f"{source}: unsupported OpenDX array header: {line.strip()}") from error

    if data_index + 1 >= len(tokens) or tokens[data_index + 1] != "follows":
        raise ValueError(f"{source}: only inline 'data follows' arrays are supported")

    if rank == 0:
        shape: tuple[int, ...] = ()
    else:
        if "shape" not in tokens:
            raise ValueError(f"{source}: rank-{rank} array has no shape")
        shape_index = tokens.index("shape")
        try:
            shape = tuple(int(value) for value in tokens[shape_index + 1 : item_index])
        except ValueError as error:
            raise ValueError(f"{source}: invalid OpenDX array shape") from error
        if len(shape) != rank:
            raise ValueError(
                f"{source}: rank {rank} does not match shape {shape!r}"
            )

    return tokens[1], value_type, rank, shape, item_count


def _parse_number(token: str, integer: bool, source: Path) -> float | int:
    try:
        return int(token) if integer else float(token)
    except ValueError as error:
        kind = "integer" if integer else "number"
        raise ValueError(f"{source}: expected an OpenDX {kind}, found {token!r}") from error


def _read_arrays_and_components(source: Path) -> tuple[dict[str, DXArray], dict[str, str]]:
    """Read inline arrays and field-component references from an OpenDX file."""

    lines = source.read_text(encoding="utf-8").splitlines()
    arrays: dict[str, DXArray] = {}
    components: dict[str, str] = {}
    line_index = 0

    while line_index < len(lines):
        line = lines[line_index]
        stripped = line.strip()

        if stripped.startswith("component "):
            tokens = shlex.split(stripped, comments=True)
            if len(tokens) >= 4 and tokens[0] == "component" and tokens[2] == "value":
                components[tokens[1]] = tokens[3]
            line_index += 1
            continue

        metadata = _parse_array_header(line, source)
        if metadata is None:
            line_index += 1
            continue

        object_id, value_type, rank, shape, item_count = metadata
        width = math.prod(shape) if shape else 1
        expected_values = item_count * width
        integer = value_type.lower() in {
            "byte",
            "short",
            "int",
            "integer",
            "long",
            "ubyte",
            "ushort",
            "uint",
            "ulong",
        }
        values: list[float | int] = []
        line_index += 1

        while line_index < len(lines) and len(values) < expected_values:
            data_text = lines[line_index].split("#", 1)[0].strip()
            if data_text:
                data_tokens = data_text.split()
                remaining = expected_values - len(values)
                if len(data_tokens) > remaining:
                    raise ValueError(
                        f"{source}: object {object_id} contains more data than declared"
                    )
                values.extend(
                    _parse_number(token, integer, source) for token in data_tokens
                )
            line_index += 1

        if len(values) != expected_values:
            raise ValueError(
                f"{source}: object {object_id} declares {expected_values} values "
                f"but contains {len(values)}"
            )

        rows = [values[offset : offset + width] for offset in range(0, len(values), width)]
        arrays[object_id] = DXArray(
            object_id=object_id,
            value_type=value_type,
            rank=rank,
            shape=shape,
            item_count=item_count,
            rows=rows,
        )

    return arrays, components


def _get_component(
    arrays: dict[str, DXArray],
    components: dict[str, str],
    name: str,
    fallback_id: str,
    source: Path,
) -> DXArray:
    object_id = components.get(name, fallback_id)
    try:
        return arrays[object_id]
    except KeyError as error:
        raise ValueError(
            f"{source}: field component {name!r} refers to missing object {object_id!r}"
        ) from error


def read_dx_mesh(source: Path, ignore_scalar_data: bool = False) -> DXMesh:
    """Read a 2-D triangular AFEPack/OpenDX mesh and optional scalar array."""

    arrays, components = _read_arrays_and_components(source)
    positions = _get_component(arrays, components, "positions", "1", source)
    connections = _get_component(arrays, components, "connections", "2", source)

    if positions.rank != 1 or positions.shape != (2,):
        raise ValueError(
            f"{source}: positions must have rank 1 and shape 2, "
            f"found rank {positions.rank} and shape {positions.shape!r}"
        )
    if connections.rank != 1 or connections.shape != (3,):
        raise ValueError(
            f"{source}: connections must be triangular (rank 1, shape 3)"
        )

    points = [[float(value) for value in row] for row in positions.rows]
    triangles = [[int(value) for value in row] for row in connections.rows]
    if triangles:
        smallest_index = min(min(triangle) for triangle in triangles)
        largest_index = max(max(triangle) for triangle in triangles)
        if smallest_index < 0 or largest_index >= len(points):
            raise ValueError(
                f"{source}: triangle index range {smallest_index}..{largest_index} "
                f"is invalid for {len(points)} positions"
            )

    scalar_values: list[float] | None = None
    scalar_location: str | None = None
    data_id = components.get("data")
    if data_id is not None and not ignore_scalar_data:
        try:
            data = arrays[data_id]
        except KeyError as error:
            raise ValueError(
                f"{source}: data component refers to missing object {data_id!r}"
            ) from error
        if data.rank != 0:
            raise ValueError(
                f"{source}: only scalar OpenDX data can be plotted; "
                f"object {data_id} has rank {data.rank}"
            )
        scalar_values = [float(row[0]) for row in data.rows]
        if len(scalar_values) == len(points):
            scalar_location = "nodal"
        elif len(scalar_values) == len(triangles):
            scalar_location = "element"
        else:
            raise ValueError(
                f"{source}: scalar data has {len(scalar_values)} items; expected "
                f"{len(points)} nodal or {len(triangles)} element values"
            )

    return DXMesh(
        points=points,
        triangles=triangles,
        scalar_values=scalar_values,
        scalar_location=scalar_location,
    )


def render_dx(
    mesh: DXMesh,
    source: Path,
    output: Path,
    *,
    title: str | None,
    cmap: str,
    dpi: int,
    show_nodes: bool,
    show_edges: bool,
) -> None:
    """Render a parsed mesh to the requested figure file."""

    x_coordinates = [point[0] for point in mesh.points]
    y_coordinates = [point[1] for point in mesh.points]
    figure, axis = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)

    if mesh.scalar_values is None:
        axis.triplot(
            x_coordinates,
            y_coordinates,
            mesh.triangles,
            color="#334155",
            linewidth=0.55,
        )
        plot_kind = "mesh"
    else:
        if mesh.scalar_location == "nodal":
            field = axis.tripcolor(
                x_coordinates,
                y_coordinates,
                mesh.triangles,
                mesh.scalar_values,
                shading="gouraud",
                cmap=cmap,
            )
        else:
            field = axis.tripcolor(
                x_coordinates,
                y_coordinates,
                mesh.triangles,
                facecolors=mesh.scalar_values,
                shading="flat",
                cmap=cmap,
            )
        colorbar = figure.colorbar(field, ax=axis, shrink=0.86, pad=0.025)
        colorbar.set_label("value")
        if show_edges:
            axis.triplot(
                x_coordinates,
                y_coordinates,
                mesh.triangles,
                color="#0f172a",
                linewidth=0.25,
                alpha=0.45,
            )
        plot_kind = f"{mesh.scalar_location} scalar field"

    if show_nodes:
        axis.scatter(
            x_coordinates,
            y_coordinates,
            s=7,
            color="#dc2626",
            edgecolors="none",
            zorder=3,
        )

    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_title(title if title is not None else f"{source.name}: {plot_kind}")

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=dpi)
    plt.close(figure)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Render a 2-D triangular AFEPack/OpenDX file as a mesh or scalar field."
        )
    )
    parser.add_argument("input", type=Path, help="OpenDX input file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="figure path; default: INPUT with a .png suffix",
    )
    parser.add_argument(
        "--mesh-only",
        action="store_true",
        help="ignore scalar data and draw only the triangulation",
    )
    parser.add_argument(
        "--show-nodes",
        action="store_true",
        help="draw a marker at every mesh node",
    )
    parser.add_argument(
        "--no-edges",
        action="store_true",
        help="hide triangle edges in scalar-field plots",
    )
    parser.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib color map for scalar data (default: viridis)",
    )
    parser.add_argument("--title", help="custom figure title")
    parser.add_argument(
        "--dpi",
        type=int,
        default=180,
        help="raster output resolution (default: 180)",
    )
    return parser


def main(arguments: Sequence[str] | None = None) -> int:
    parser = build_argument_parser()
    options = parser.parse_args(arguments)

    if not options.input.is_file():
        parser.error(f"input file does not exist: {options.input}")
    if options.dpi <= 0:
        parser.error("--dpi must be positive")

    output = options.output if options.output is not None else options.input.with_suffix(".png")
    if not output.suffix:
        output = output.with_suffix(".png")

    try:
        mesh = read_dx_mesh(options.input, ignore_scalar_data=options.mesh_only)
        render_dx(
            mesh,
            options.input,
            output,
            title=options.title,
            cmap=options.cmap,
            dpi=options.dpi,
            show_nodes=options.show_nodes,
            show_edges=not options.no_edges,
        )
    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")

    if mesh.scalar_values is None:
        data_description = "mesh only"
    else:
        data_description = f"{mesh.scalar_location} scalar field"
    print(
        f"Wrote {output} ({len(mesh.points)} nodes, "
        f"{len(mesh.triangles)} triangles, {data_description})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
