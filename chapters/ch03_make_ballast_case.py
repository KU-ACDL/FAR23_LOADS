from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from ch03_weight_envelope import (
    PROJECT_ROOT,
    WeightItem,
    cg_point,
    compute_structural_limits,
    load_case,
    project_path,
    project_relative,
)


INPUT_PATH = PROJECT_ROOT / "inputs" / "weight_db.json"
CONFIG_PATH = PROJECT_ROOT / "inputs" / "ballast_case_config.json"
ENVELOPE_OUTPUT_PATH = PROJECT_ROOT / "outputs" / "cases" / "ch3_weight_envelope_output.json"
DEFAULT_CASE_NAME = "aft_gross"


def read_json(path: str | Path) -> dict[str, Any]:
    path = project_path(path)
    with path.open("r", encoding="utf-8-sig") as f:
        return json.load(f)


def write_json(data: dict[str, Any], path: str | Path) -> None:
    path = project_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def write_text(text: str, path: str | Path) -> None:
    path = project_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write(text)


def component_to_dict(component: WeightItem, name_override: str | None = None) -> dict[str, Any]:
    data = asdict(component)
    if name_override:
        data["NS"] = name_override
    if abs(float(data.get("IIXZ", 0.0))) < 1.0e-12:
        data.pop("IIXZ", None)
    return data


def get_ballast_case_config(config: dict[str, Any], case_name: str) -> dict[str, Any]:
    cases = config.get("cases", {})
    if case_name not in cases:
        raise ValueError(f"Missing ballast case configuration: {case_name}")
    return cases[case_name]


def find_components_by_name(items: list[WeightItem], names: list[str]) -> list[WeightItem]:
    by_name = {item.NS.upper(): item for item in items}
    selected: list[WeightItem] = []
    missing: list[str] = []

    for name in names:
        item = by_name.get(name.upper())
        if item is None:
            missing.append(name)
        else:
            selected.append(item)

    if missing:
        raise ValueError(f"Selected load(s) not found in DB: {missing}")
    return selected


def load_selected_point(
    point_id: int,
    envelope_output_path: str | Path,
) -> dict[str, Any]:
    envelope_output = read_json(envelope_output_path)
    for point in envelope_output.get("envelope", {}).get("load_points", []):
        if int(point["POINT_ID"]) == point_id:
            return point
    raise ValueError(f"Selected point_id {point_id} not found in {project_relative(envelope_output_path)}")


def make_ballast_component(
    name: str,
    weight: float,
    x: float,
    z: float,
) -> WeightItem:
    return WeightItem(
        NS=name,
        W=weight,
        X=x,
        Y=0.0,
        Z=z,
        C=2,
        IIX=0.0,
        IIY=0.0,
        IIZ=0.0,
        IIXZ=0.0,
    )


def build_report(summary: dict[str, Any]) -> str:
    target = summary["target"]
    selected = summary["selected_point"]
    ballast = summary["ballast"]
    final = summary["final_point"]
    outputs = summary["outputs"]

    lines = [
        "Ballast Case Generation",
        "",
        f"CASE: {summary['case']}",
        f"SOURCE INPUT: {summary['source_input']}",
        f"SOURCE CONFIG: {summary['source_config']}",
        f"SOURCE ENVELOPE: {summary['source_envelope_output']}",
        f"SOURCE DB: {summary['source_db']}",
        "",
        "TARGET",
        f"NAME                 {target['name']}",
        f"WEIGHT               {target['WEIGHT']:.3f}",
        f"XBAR                 {target['XBAR']:.3f}",
        f"GEOMETRY SOURCE      {target['geometry_source']}",
        "",
        "SELECTED LOAD POINT",
        f"POINT ID             {selected['POINT_ID']}",
        f"ADDED                {selected['ADDED']}",
        f"WEIGHT               {selected['WEIGHT']:.3f}",
        f"XBAR                 {selected['XBAR']:.3f}",
        f"ZBAR                 {selected['ZBAR']:.3f}",
        "LOADS                " + ", ".join(selected["loads"]),
        "",
        "BALLAST",
        f"NAME                 {ballast['NS']}",
        f"WEIGHT               {ballast['W']:.3f}",
        f"X                    {ballast['X']:.3f}",
        f"Y                    {ballast['Y']:.3f}",
        f"Z                    {ballast['Z']:.3f}",
        "",
        "FINAL POINT",
        f"WEIGHT               {final['WEIGHT']:.3f}",
        f"XBAR                 {final['XBAR']:.3f}",
        f"ZBAR                 {final['ZBAR']:.3f}",
        "",
        "OUTPUT",
        f"CH4 INPUT DB         {outputs['ch4_input']}",
    ]
    return "\n".join(lines) + "\n"


def make_ballast_case(
    input_path: str | Path = INPUT_PATH,
    config_path: str | Path = CONFIG_PATH,
    case_name: str = DEFAULT_CASE_NAME,
) -> dict[str, Any]:
    selected_database_name, items, ch3_config = load_case(input_path, case_name=None)
    config = read_json(config_path)
    ballast_config = get_ballast_case_config(config, case_name)
    envelope_output_path = ballast_config.get("envelope_output_path", ENVELOPE_OUTPUT_PATH)

    target_name = str(ballast_config.get("target", "AFTGROSS")).upper()
    structural_limits = compute_structural_limits(ch3_config)
    if not structural_limits:
        raise ValueError("Structural limits are required to make a ballast case.")
    if target_name not in structural_limits["points"]:
        raise ValueError(f"Target structural limit '{target_name}' not found.")

    target = structural_limits["points"][target_name]
    target_weight = float(target["WEIGHT"])
    target_x = float(target["XBAR"])

    base_items = [item for item in items if item.C <= 1 and item.W > 0]
    selected_point = None
    if "selected_point_id" in ballast_config:
        selected_point = load_selected_point(
            int(ballast_config["selected_point_id"]),
            envelope_output_path,
        )
        selected_load_names = list(selected_point.get("LOADS", []))
    else:
        selected_load_names = list(ballast_config.get("selected_loads", []))
    selected_loads = find_components_by_name(items, selected_load_names)
    selected_items = base_items + selected_loads

    selected_cg = cg_point(selected_items, "SELECTED LOAD")
    ballast_weight = target_weight - selected_cg.WEIGHT
    if ballast_weight <= 0:
        raise ValueError(
            f"Selected load is already at or above target weight: "
            f"{selected_cg.WEIGHT:.3f} >= {target_weight:.3f}"
        )

    ballast_x = (
        target_weight * target_x - selected_cg.WEIGHT * selected_cg.XBAR
    ) / ballast_weight
    ballast_z = float(ballast_config.get("ballast_z", selected_cg.ZBAR))
    ballast = make_ballast_component(
        name=str(ballast_config.get("ballast_name", "BALLAST")),
        weight=ballast_weight,
        x=ballast_x,
        z=ballast_z,
    )

    final_items = selected_items + [ballast]
    final_point = cg_point(final_items, target_name)
    name_overrides = {
        str(key).upper(): str(value)
        for key, value in ballast_config.get("output_name_overrides", {}).items()
    }
    output_components = [
        component_to_dict(item, name_overrides.get(item.NS.upper()))
        for item in final_items
    ]

    output_input_path = ballast_config.get(
        "output_input_path",
        f"inputs/ch4_{case_name}_input.json",
    )
    output_report_path = ballast_config.get(
        "output_report_path",
        f"outputs/reports/ch3_{case_name}_ballast_case_report.txt",
    )

    ch4_input = {
        "cases": {
            case_name: output_components,
        }
    }
    summary = {
        "chapter": 3,
        "program": "WTENV.BAS",
        "source_input": project_relative(input_path),
        "source_config": project_relative(config_path),
        "source_envelope_output": project_relative(envelope_output_path),
        "source_db": selected_database_name,
        "case": case_name,
        "target": {
            "name": target_name,
            "WEIGHT": target_weight,
            "XBAR": target_x,
            "geometry_source": structural_limits["input"]["geometry_source"],
        },
        "selected_point": {
            "POINT_ID": selected_point.get("POINT_ID") if selected_point else None,
            "ADDED": selected_point.get("ADDED") if selected_point else "SELECTED LOAD",
            "WEIGHT": selected_cg.WEIGHT,
            "XBAR": selected_cg.XBAR,
            "ZBAR": selected_cg.ZBAR,
            "loads": selected_load_names,
        },
        "ballast": {
            "NS": ballast.NS,
            "W": ballast.W,
            "X": ballast.X,
            "Y": ballast.Y,
            "Z": ballast.Z,
        },
        "final_point": {
            "WEIGHT": final_point.WEIGHT,
            "XBAR": final_point.XBAR,
            "ZBAR": final_point.ZBAR,
        },
        "outputs": {
            "ch4_input": project_relative(project_path(output_input_path)),
            "report": project_relative(project_path(output_report_path)),
        },
    }

    write_json(ch4_input, output_input_path)
    write_text(build_report(summary), output_report_path)
    return summary


if __name__ == "__main__":
    result = make_ballast_case()
    print(json.dumps(result, indent=2, ensure_ascii=False))
