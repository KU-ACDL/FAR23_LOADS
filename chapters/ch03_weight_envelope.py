from __future__ import annotations

import json
from itertools import combinations
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = PROJECT_ROOT / "inputs" / "weight_db.json"
GEOMETRY_OUTPUT_PATH = PROJECT_ROOT / "outputs" / "cases" / "ch5_geometry_output.json"
OUTPUT_PATH = PROJECT_ROOT / "outputs" / "cases" / "ch3_weight_envelope_output.json"
REPORT_PATH = PROJECT_ROOT / "outputs" / "reports" / "ch3_weight_envelope_report.txt"


# ============================================================
# Chapter 3 - WTENV.BAS
# Input : inputs/weight_db.json
# Output: outputs/cases/ch3_weight_envelope_output.json
# ============================================================


@dataclass
class WeightItem:
    NS: str
    W: float
    X: float
    Y: float = 0.0
    Z: float = 0.0
    C: int = 2
    INDEX: int = 0
    IIX: float = 0.0
    IIY: float = 0.0
    IIZ: float = 0.0
    IIXZ: float = 0.0


@dataclass
class EnvelopePoint:
    ADDED: str
    XBAR: float
    ZBAR: float
    WEIGHT: float


DEFAULT_PROGRAM_DATA = {
    "M": 100,
    "ID": "6 PLACE AIRPLANE",
    "F": "WTAFTCG.INP",
    "Q": True,
}


def read_json(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8-sig") as f:
        return json.load(f)


def write_json(data: dict[str, Any], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def write_text(text: str, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write(text)


def project_relative(path: str | Path) -> str:
    path = Path(path)
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def project_path(path: str | Path) -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    return PROJECT_ROOT / path


def validate_component(item: dict[str, Any], index: int) -> None:
    for key in ["W", "X"]:
        if key not in item:
            raise ValueError(f"Missing component key '{key}' in row {index}.")
    if "NS" not in item and "N" not in item:
        raise ValueError(f"Missing component name key 'NS' or 'N' in row {index}.")


def make_weight_item(item: dict[str, Any], index: int) -> WeightItem:
    validate_component(item, index)
    return WeightItem(
        NS=str(item.get("NS", item.get("N"))),
        W=float(item["W"]),
        X=float(item["X"]),
        Y=float(item.get("Y", 0.0)),
        Z=float(item.get("Z", 0.0)),
        C=int(item.get("C", 2)),
        INDEX=index,
        IIX=float(item.get("IIX", item.get("IXX", 0.0))),
        IIY=float(item.get("IIY", item.get("IYY", 0.0))),
        IIZ=float(item.get("IIZ", item.get("IZZ", 0.0))),
        IIXZ=float(item.get("IIXZ", item.get("IXZ", 0.0))),
    )


def load_case(
    path: str | Path,
    case_name: str | None = None,
) -> tuple[str, list[WeightItem], dict[str, Any]]:
    data = read_json(path)
    components = data.get("components")
    selected_case_name = str(data.get("name", "weight_db"))

    if components is None:
        cases = data.get("cases")
        if not isinstance(cases, dict) or not cases:
            raise ValueError("Input DB must contain a non-empty 'components' list.")

        selected_case_name = case_name or next(iter(cases))
        if selected_case_name not in cases:
            raise ValueError(f"Case '{selected_case_name}' not found in input DB.")
        components = cases[selected_case_name]

    if not isinstance(components, list) or not components:
        raise ValueError("Input DB must contain a non-empty component list.")

    ch3_config = data.get("ch3", {})
    if ch3_config and not isinstance(ch3_config, dict):
        raise ValueError("Input DB 'ch3' section must be an object.")

    return (
        selected_case_name,
        [make_weight_item(item, index) for index, item in enumerate(components, start=1)],
        ch3_config,
    )


def cg_point(items: list[WeightItem], added: str) -> EnvelopePoint:
    total_weight = sum(item.W for item in items)
    if total_weight <= 0:
        raise ValueError("Total weight must be positive.")

    return EnvelopePoint(
        ADDED=added,
        XBAR=sum(item.W * item.X for item in items) / total_weight,
        ZBAR=sum(item.W * item.Z for item in items) / total_weight,
        WEIGHT=total_weight,
    )


def compute_edge(
    base_items: list[WeightItem],
    discretionary_items: list[WeightItem],
    reverse: bool,
) -> list[dict[str, Any]]:
    active_items = list(base_items)
    active_loads: list[str] = []
    first_point = point_to_dict(cg_point(active_items, "MINIMUM WEIGHT"))
    first_point["LOADS"] = []
    points = [first_point]

    if reverse:
        ordered_items = sorted(discretionary_items, key=lambda component: (-component.X, component.INDEX))
    else:
        ordered_items = sorted(discretionary_items, key=lambda component: (component.X, -component.INDEX))

    for item in ordered_items:
        active_items.append(item)
        active_loads.append(item.NS)
        point = point_to_dict(cg_point(active_items, item.NS))
        point["LOADS"] = list(active_loads)
        points.append(point)

    return points


def compute_basic_points(
    empty_items: list[WeightItem],
    minimum_items: list[WeightItem],
) -> list[EnvelopePoint]:
    basic_points = [cg_point(empty_items, "EMPTY WEIGHT")]
    active_items = list(empty_items)
    minimum_additions = [item for item in minimum_items if item not in empty_items]

    for item in minimum_additions:
        active_items.append(item)
        if item.NS.upper().startswith("PILOT"):
            basic_points.append(cg_point(active_items, item.NS))

    minimum_point = cg_point(minimum_items, "MINIMUM WEIGHT")
    last_point = basic_points[-1]
    if (
        abs(last_point.XBAR - minimum_point.XBAR) > 1.0e-9
        or abs(last_point.WEIGHT - minimum_point.WEIGHT) > 1.0e-9
    ):
        basic_points.append(minimum_point)

    return basic_points


def compute_combination_points(
    base_items: list[WeightItem],
    discretionary_items: list[WeightItem],
) -> list[dict[str, Any]]:
    points: list[dict[str, Any]] = []
    for count in range(len(discretionary_items) + 1):
        for selected_items in combinations(discretionary_items, count):
            selected = list(selected_items)
            active_items = base_items + selected
            label = selected[-1].NS if selected else "MINIMUM WEIGHT"
            point = point_to_dict(cg_point(active_items, label))
            point["LOADS"] = [item.NS for item in selected]
            points.append(point)
    return points


def assign_load_point_ids(
    forward_edge: list[dict[str, Any]],
    aft_edge: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    load_points: list[dict[str, Any]] = []
    point_ids_by_loads: dict[tuple[str, ...], int] = {}

    for edge_name, edge in [("forward_edge", forward_edge), ("aft_edge", aft_edge)]:
        for point in edge:
            loads = tuple(sorted(str(load) for load in point.get("LOADS", [])))
            if loads not in point_ids_by_loads:
                point_id = len(load_points) + 1
                point_ids_by_loads[loads] = point_id
                load_points.append(
                    {
                        "POINT_ID": point_id,
                        "ADDED": point["ADDED"],
                        "XBAR": point["XBAR"],
                        "ZBAR": point["ZBAR"],
                        "WEIGHT": point["WEIGHT"],
                        "S1": point["S1"],
                        "LOADS": list(point.get("LOADS", [])),
                        "EDGES": [edge_name],
                    }
                )
            else:
                point_id = point_ids_by_loads[loads]
                load_point = load_points[point_id - 1]
                if edge_name not in load_point["EDGES"]:
                    load_point["EDGES"].append(edge_name)
            point["POINT_ID"] = point_id

    return load_points


def point_to_dict(point: EnvelopePoint) -> dict[str, float | str]:
    return {
        "ADDED": point.ADDED,
        "XBAR": point.XBAR,
        "ZBAR": point.ZBAR,
        "WEIGHT": point.WEIGHT,
        "S1": point.WEIGHT,
    }


def load_wing_geometry(path: str | Path) -> dict[str, float]:
    geometry_path = project_path(path)
    if not geometry_path.exists():
        raise FileNotFoundError(
            f"Chapter 5 geometry output not found: {project_relative(geometry_path)}. "
            "Run chapters/ch05_geometry.py before chapters/ch03_weight_envelope.py."
        )

    geometry_output = read_json(geometry_path)
    for surface in geometry_output.get("surfaces", []):
        surface_name = str(surface.get("N$", "")).upper()
        surface_output = surface.get("Surface Output", {})
        if "WING" in surface_name and "MAC" in surface_output and "XMACLE" in surface_output:
            return {
                "MAC": float(surface_output["MAC"]),
                "XMACLE": float(surface_output["XMACLE"]),
                "source": project_relative(geometry_path),
            }

    raise ValueError(
        f"Could not find wing MAC/XMACLE in Chapter 5 output: {project_relative(geometry_path)}"
    )


def compute_structural_limits(ch3_config: dict[str, Any]) -> dict[str, Any] | None:
    structural = ch3_config.get("structural_limits")
    if not structural:
        return None

    required = [
        "FWDGRL",
        "AFTGRL",
        "FWDREDL",
        "GROSS_WEIGHT",
        "REDUCED_WEIGHT",
    ]
    for key in required:
        if key not in structural:
            raise ValueError(f"Missing ch3.structural_limits key: {key}")

    geometry = load_wing_geometry(ch3_config.get("geometry_output_path", GEOMETRY_OUTPUT_PATH))
    mac = geometry["MAC"]
    xmacle = geometry["XMACLE"]
    gross_weight = float(structural["GROSS_WEIGHT"])
    reduced_weight = float(structural["REDUCED_WEIGHT"])

    fwd_gross = xmacle + float(structural["FWDGRL"]) / 100.0 * mac
    aft_gross = xmacle + float(structural["AFTGRL"]) / 100.0 * mac
    fwd_reduced = xmacle + float(structural["FWDREDL"]) / 100.0 * mac

    points = {
        "FWDGROSS": {
            "label": "Fwd Gross",
            "XBAR": fwd_gross,
            "WEIGHT": gross_weight,
            "S1": gross_weight,
        },
        "AFTGROSS": {
            "label": "Aft Gross",
            "XBAR": aft_gross,
            "WEIGHT": gross_weight,
            "S1": gross_weight,
        },
        "FWDRED": {
            "label": "Fwd Reduced",
            "XBAR": fwd_reduced,
            "WEIGHT": reduced_weight,
            "S1": reduced_weight,
        },
    }

    return {
        "input": {
            "MAC": mac,
            "XMACLE": xmacle,
            "geometry_source": geometry["source"],
            "FWDGRL": float(structural["FWDGRL"]),
            "AFTGRL": float(structural["AFTGRL"]),
            "FWDREDL": float(structural["FWDREDL"]),
            "GROSS_WEIGHT": gross_weight,
            "REDUCED_WEIGHT": reduced_weight,
        },
        "points": points,
        "plot_lines": [
            {
                "name": "forward_reduced_to_gross",
                "points": [
                    {"x": fwd_reduced, "y": reduced_weight},
                    {"x": fwd_gross, "y": gross_weight},
                ],
            },
            {
                "name": "gross_weight_limit",
                "points": [
                    {"x": fwd_gross, "y": gross_weight},
                    {"x": aft_gross, "y": gross_weight},
                ],
            },
        ],
        "plot_polygon": [
            {"x": fwd_reduced, "y": reduced_weight, "label": "Fwd Reduced"},
            {"x": fwd_gross, "y": gross_weight, "label": "Fwd Gross"},
            {"x": aft_gross, "y": gross_weight, "label": "Aft Gross"},
        ],
    }


def compute_weight_envelope(
    items: list[WeightItem],
    exclude_names: set[str] | None = None,
    ch3_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    ch3_config = ch3_config or {}
    exclude_names = exclude_names or set()
    empty_items = [item for item in items if item.C == 0 and item.W > 0]
    minimum_items = [item for item in items if item.C <= 1 and item.W > 0]
    discretionary_items = [
        item
        for item in items
        if item.C >= 2 and item.W > 0 and item.NS.upper() not in exclude_names
    ]
    excluded_items = [
        item
        for item in items
        if item.C >= 2 and item.W > 0 and item.NS.upper() in exclude_names
    ]

    if not empty_items:
        raise ValueError("At least one empty-weight item with C=0 is required.")
    if not minimum_items:
        raise ValueError("At least one minimum-weight item with C<=1 is required.")

    forward_edge = compute_edge(
        base_items=minimum_items,
        discretionary_items=discretionary_items,
        reverse=False,
    )
    aft_edge = compute_edge(
        base_items=minimum_items,
        discretionary_items=discretionary_items,
        reverse=True,
    )

    load_points = assign_load_point_ids(forward_edge, aft_edge)
    polygon_points = forward_edge + list(reversed(aft_edge[:-1]))
    structural_limits = compute_structural_limits(ch3_config)
    basic_points = compute_basic_points(empty_items, minimum_items)
    combination_points = compute_combination_points(minimum_items, discretionary_items)

    return {
        "empty_weight": point_to_dict(cg_point(empty_items, "EMPTY WEIGHT")),
        "minimum_weight": point_to_dict(cg_point(minimum_items, "MINIMUM WEIGHT")),
        "basic_points": [point_to_dict(point) for point in basic_points],
        "load_points": load_points,
        "combination_points": combination_points,
        "load_groups": {
            "empty_weight_items": [item.NS for item in empty_items],
            "minimum_weight_items": [item.NS for item in minimum_items],
            "discretionary_load_items": [item.NS for item in discretionary_items],
            "excluded_items": [item.NS for item in excluded_items],
        },
        "forward_edge": forward_edge,
        "aft_edge": aft_edge,
        "plot_data": {
            "title": ch3_config.get(
                "title",
                "USEFUL LOAD ENVELOPE AND STRUCTURAL LIMITS",
            ),
            "x_axis": ch3_config.get("x_axis", "Fuselage Station"),
            "y_axis": ch3_config.get("y_axis", "Weight"),
            "x_value": "XBAR",
            "y_value": "WEIGHT",
            "basic_points": [
                {"x": point.XBAR, "y": point.WEIGHT, "label": point.ADDED}
                for point in basic_points
            ],
            "combination_points": [
                {
                    "x": point["XBAR"],
                    "y": point["WEIGHT"],
                    "label": point["ADDED"],
                    "loads": point["LOADS"],
                }
                for point in combination_points
            ],
            "forward_edge": [
                {
                    "x": point["XBAR"],
                    "y": point["WEIGHT"],
                    "label": point["ADDED"],
                    "point_id": point["POINT_ID"],
                    "loads": point["LOADS"],
                }
                for point in forward_edge
            ],
            "aft_edge": [
                {
                    "x": point["XBAR"],
                    "y": point["WEIGHT"],
                    "label": point["ADDED"],
                    "point_id": point["POINT_ID"],
                    "loads": point["LOADS"],
                }
                for point in aft_edge
            ],
            "envelope_polygon": [
                {
                    "x": point["XBAR"],
                    "y": point["WEIGHT"],
                    "label": point["ADDED"],
                    "point_id": point["POINT_ID"],
                    "loads": point["LOADS"],
                }
                for point in polygon_points
            ],
            "structural_limits": structural_limits,
        },
        "structural_limits": structural_limits,
    }


def build_report(database_name: str, result: dict[str, Any]) -> str:
    lines = [
        "Weight CG Envelope",
        "",
        f"SOURCE DB: {database_name}",
        "",
        "FORWARD EDGE",
        "NO  ADDED                       XBAR       ZBAR     WEIGHT",
    ]

    for point in result["envelope"]["forward_edge"]:
        lines.append(
            f"{point['POINT_ID']:>2}  "
            f"{point['ADDED']:<26}"
            f"{point['XBAR']:>10.3f}"
            f"{point['ZBAR']:>11.3f}"
            f"{point['WEIGHT']:>11.3f}"
        )

    lines.extend(
        [
            "",
            "AFT EDGE",
            "NO  ADDED                       XBAR       ZBAR     WEIGHT",
        ]
    )

    for point in result["envelope"]["aft_edge"]:
        lines.append(
            f"{point['POINT_ID']:>2}  "
            f"{point['ADDED']:<26}"
            f"{point['XBAR']:>10.3f}"
            f"{point['ZBAR']:>11.3f}"
            f"{point['WEIGHT']:>11.3f}"
        )

    structural_limits = result["envelope"].get("structural_limits")
    if structural_limits:
        lines.extend(
            [
                "",
                "STRUCTURAL LIMIT POINTS",
                "POINT                         XBAR     WEIGHT",
            ]
        )
        for key in ["FWDRED", "FWDGROSS", "AFTGROSS"]:
            point = structural_limits["points"][key]
            lines.append(
                f"{point['label']:<26}"
                f"{point['XBAR']:>10.3f}"
                f"{point['WEIGHT']:>11.3f}"
            )

    return "\n".join(lines) + "\n"


def run_chapter3(
    input_path: str | Path = INPUT_PATH,
    output_path: str | Path = OUTPUT_PATH,
    report_path: str | Path = REPORT_PATH,
    case_name: str | None = None,
) -> dict[str, Any]:
    selected_database_name, items, ch3_config = load_case(input_path, case_name=case_name)
    program_data = DEFAULT_PROGRAM_DATA | ch3_config.get("program_data", {})
    plot_input = {
        "title": ch3_config.get("title", "USEFUL LOAD ENVELOPE AND STRUCTURAL LIMITS"),
        "x_axis": ch3_config.get("x_axis", "Fuselage Station"),
        "y_axis": ch3_config.get("y_axis", "Weight"),
        "geometry_output_path": project_relative(
            project_path(ch3_config.get("geometry_output_path", GEOMETRY_OUTPUT_PATH))
        ),
        "structural_limits": ch3_config.get("structural_limits"),
    }
    exclude_names = {name.upper() for name in ch3_config.get("exclude_names", [])}
    envelope = compute_weight_envelope(
        items,
        exclude_names=exclude_names,
        ch3_config=ch3_config,
    )

    output = {
        "chapter": 3,
        "program": "WTENV.BAS",
        "input": {
            "path": project_relative(input_path),
            "database": selected_database_name,
            "program_data": program_data,
            "plot_input": plot_input,
            "exclude_names": sorted(exclude_names),
        },
        "units": {
            "XBAR": "in",
            "ZBAR": "in",
            "WEIGHT": "lb",
        },
        "envelope": envelope,
    }

    write_json(output, output_path)
    write_text(build_report(selected_database_name, output), report_path)
    return output


if __name__ == "__main__":
    result = run_chapter3()
    print(json.dumps(result, indent=2, ensure_ascii=False))
