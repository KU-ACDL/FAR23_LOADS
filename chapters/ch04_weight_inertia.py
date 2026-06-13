import json
from dataclasses import dataclass, asdict
from math import atan2, cos, sin, degrees
from pathlib import Path
from typing import List, Dict, Any, Union


BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_PATH = BASE_DIR / "inputs" / "ch4_input.json"
OUTPUT_PATH = BASE_DIR / "outputs" / "ch4_output.json"


@dataclass
class Component:
    name: str
    weight: float  # [lb]
    x: float       # [in]
    y: float       # [in]
    z: float       # [in]
    type: int = 0
    ixx_cg: float = 0.0  # [lb-in^2]
    iyy_cg: float = 0.0  # [lb-in^2]
    izz_cg: float = 0.0  # [lb-in^2]
    ixz_cg: float = 0.0  # [lb-in^2]


@dataclass
class CG:
    x: float  # [in]
    y: float  # [in]
    z: float  # [in]


@dataclass
class AircraftAxesInertia:
    ixx: float  # [lb-in^2]
    iyy: float  # [lb-in^2]
    izz: float  # [lb-in^2]
    ixz: float  # [lb-in^2]


@dataclass
class PrincipalAxesInertia:
    principal_angle_deg: float
    ixp: float  # [lb-in^2]
    iyp: float  # [lb-in^2]
    izp: float  # [lb-in^2]


@dataclass
class Chapter4Result:
    chapter: int
    case_name: str
    total_weight: float
    cg: CG
    inertia_aircraft_axes: AircraftAxesInertia
    inertia_principal_axes: PrincipalAxesInertia


def compute_cg(components: List[Component]) -> tuple[float, CG]:
    total_weight = sum(c.weight for c in components)

    if total_weight <= 0:
        raise ValueError("Total weight must be positive.")

    x_cg = sum(c.weight * c.x for c in components) / total_weight
    y_cg = sum(c.weight * c.y for c in components) / total_weight
    z_cg = sum(c.weight * c.z for c in components) / total_weight

    return total_weight, CG(x=x_cg, y=y_cg, z=z_cg)


def compute_aircraft_axes_inertia(
    components: List[Component],
    cg: CG
) -> AircraftAxesInertia:
    """
    Aircraft axes 기준 관성 계산.
    """
    ixx = 0.0
    iyy = 0.0
    izz = 0.0
    ixz = 0.0

    for c in components:
        dx = c.x - cg.x
        dy = c.y - cg.y
        dz = c.z - cg.z

        ixx += c.weight * (dy**2 + dz**2) + c.ixx_cg
        iyy += c.weight * (dx**2 + dz**2) + c.iyy_cg
        izz += c.weight * (dx**2 + dy**2) + c.izz_cg

        # 중요: BASIC 출력과 맞추는 IXZ 부호
        ixz += c.weight * dx * dz + c.ixz_cg

    return AircraftAxesInertia(
        ixx=ixx,
        iyy=iyy,
        izz=izz,
        ixz=ixz
    )


def compute_principal_axes_inertia(
    inertia: AircraftAxesInertia
) -> PrincipalAxesInertia:
    """
    Principal axes 기준 관성 계산.
    """
    ixx = inertia.ixx
    iyy = inertia.iyy
    izz = inertia.izz
    ixz = inertia.ixz

    # tan(2 theta) = 2 Ixz / (Izz - Ixx)
    angle_rad = 0.5 * atan2(2.0 * ixz, izz - ixx)

    c = cos(angle_rad)
    s = sin(angle_rad)

    ixp = ixx * c**2 + izz * s**2 - ixz * sin(2.0 * angle_rad)
    izp = ixx * s**2 + izz * c**2 + ixz * sin(2.0 * angle_rad)
    iyp = iyy

    return PrincipalAxesInertia(
        principal_angle_deg=degrees(angle_rad),
        ixp=ixp,
        iyp=iyp,
        izp=izp
    )


def compute_case(case_name: str, components: List[Component]) -> Chapter4Result:
    total_weight, cg = compute_cg(components)
    aircraft_inertia = compute_aircraft_axes_inertia(components, cg)
    principal_inertia = compute_principal_axes_inertia(aircraft_inertia)

    return Chapter4Result(
        chapter=4,
        case_name=case_name,
        total_weight=total_weight,
        cg=cg,
        inertia_aircraft_axes=aircraft_inertia,
        inertia_principal_axes=principal_inertia
    )


def compute_all_cases(cases: Dict[str, List[Component]]) -> Dict[str, Any]:
    results = {}

    for case_name, components in cases.items():
        result = compute_case(case_name, components)
        results[case_name] = asdict(result)

    return {
        "chapter": 4,
        "cases": results
    }


def load_cases(path: Union[str, Path]) -> Dict[str, List[Component]]:
    with open(path, "r", encoding="utf-8-sig") as f:
        data = json.load(f)

    return {
        case_name: [Component(**component) for component in components]
        for case_name, components in data["cases"].items()
    }


def run_chapter4(
    input_path: Union[str, Path] = INPUT_PATH,
    output_path: Union[str, Path] = OUTPUT_PATH
) -> Dict[str, Any]:
    output = compute_all_cases(load_cases(input_path))

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    return output


if __name__ == "__main__":
    result = run_chapter4()
    print(json.dumps(result, indent=2, ensure_ascii=False))
