import json
from dataclasses import asdict, dataclass
from math import atan2, cos, degrees, sin
from pathlib import Path
from typing import Any, Dict, List, Union


BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_PATH = BASE_DIR / "inputs" / "ch4_input.json"
OUTPUT_PATH = BASE_DIR / "outputs" / "ch4_output.json"
COMPONENT_OUTPUT_PATH = BASE_DIR / "outputs" / "ch4_components_output.json"


@dataclass
class Component:
    NS: str
    W: float  # [lb]
    X: float  # [in]
    Y: float  # [in]
    Z: float  # [in]
    C: int = 0
    IIX: float = 0.0  # [lb-in^2]
    IIY: float = 0.0  # [lb-in^2]
    IIZ: float = 0.0  # [lb-in^2]
    IIXZ: float = 0.0  # [lb-in^2]


@dataclass
class CG:
    XBAR: float  # [in]
    YBAR: float  # [in]
    ZBAR: float  # [in]


@dataclass
class AirplaneAxesInertia:
    IX: float  # [lb-in^2]
    IY: float  # [lb-in^2]
    IZ: float  # [lb-in^2]
    IXZ: float  # [lb-in^2]


@dataclass
class ComponentContribution:
    NS: str
    W: float
    X: float
    Y: float
    Z: float
    C: int
    DX: float
    DY: float
    DZ: float
    IX_ITEM: float
    IY_ITEM: float
    IZ_ITEM: float
    IXZ_ITEM: float
    IIX: float
    IIY: float
    IIZ: float
    IIXZ: float
    IX_TOTAL: float
    IY_TOTAL: float
    IZ_TOTAL: float
    IXZ_TOTAL: float


@dataclass
class PrincipalAxesInertia:
    THETA_DEG: float
    PXI: float  # [lb-in^2]
    PYI: float  # [lb-in^2]
    PZI: float  # [lb-in^2]


@dataclass
class Chapter4Result:
    CHAPTER: int
    CASE_NAME: str
    S1: float
    XBAR: float
    YBAR: float
    ZBAR: float
    IX: float
    IY: float
    IZ: float
    IXZ: float
    THETA_DEG: float
    PXI: float
    PYI: float
    PZI: float
    COMPONENTS: List[ComponentContribution]


def compute_cg(components: List[Component]) -> tuple[float, CG]:
    S1 = sum(component.W for component in components)

    if S1 <= 0:
        raise ValueError("Total weight must be positive.")

    XBAR = sum(component.W * component.X for component in components) / S1
    YBAR = sum(component.W * component.Y for component in components) / S1
    ZBAR = sum(component.W * component.Z for component in components) / S1

    return S1, CG(XBAR=XBAR, YBAR=YBAR, ZBAR=ZBAR)


def compute_component_contributions(
    components: List[Component],
    cg: CG
) -> List[ComponentContribution]:
    contributions = []

    for component in components:
        DX = component.X - cg.XBAR
        DY = component.Y - cg.YBAR
        DZ = component.Z - cg.ZBAR

        IX_ITEM = component.W * (DY**2 + DZ**2)
        IY_ITEM = component.W * (DX**2 + DZ**2)
        IZ_ITEM = component.W * (DX**2 + DY**2)
        IXZ_ITEM = component.W * DX * DZ

        contributions.append(
            ComponentContribution(
                NS=component.NS,
                W=component.W,
                X=component.X,
                Y=component.Y,
                Z=component.Z,
                C=component.C,
                DX=DX,
                DY=DY,
                DZ=DZ,
                IX_ITEM=IX_ITEM,
                IY_ITEM=IY_ITEM,
                IZ_ITEM=IZ_ITEM,
                IXZ_ITEM=IXZ_ITEM,
                IIX=component.IIX,
                IIY=component.IIY,
                IIZ=component.IIZ,
                IIXZ=component.IIXZ,
                IX_TOTAL=IX_ITEM + component.IIX,
                IY_TOTAL=IY_ITEM + component.IIY,
                IZ_TOTAL=IZ_ITEM + component.IIZ,
                IXZ_TOTAL=IXZ_ITEM + component.IIXZ
            )
        )

    return contributions


def sum_airplane_axes_inertia(
    contributions: List[ComponentContribution]
) -> AirplaneAxesInertia:
    IX = 0.0
    IY = 0.0
    IZ = 0.0
    IXZ = 0.0

    for contribution in contributions:
        IX += contribution.IX_TOTAL
        IY += contribution.IY_TOTAL
        IZ += contribution.IZ_TOTAL
        IXZ += contribution.IXZ_TOTAL

    return AirplaneAxesInertia(
        IX=IX,
        IY=IY,
        IZ=IZ,
        IXZ=IXZ
    )


def compute_principal_axes_inertia(
    inertia: AirplaneAxesInertia
) -> PrincipalAxesInertia:
    IX = inertia.IX
    IY = inertia.IY
    IZ = inertia.IZ
    IXZ = inertia.IXZ

    # BASIC: THETA = ATN(2 * IXZ / (IZ - IX)) / 2.
    THETA = 0.5 * atan2(2.0 * IXZ, IZ - IX)

    C = cos(THETA)
    S = sin(THETA)

    PXI = IX * C**2 + IZ * S**2 - IXZ * sin(2.0 * THETA)
    PZI = IX * S**2 + IZ * C**2 + IXZ * sin(2.0 * THETA)
    PYI = IY

    return PrincipalAxesInertia(
        THETA_DEG=degrees(THETA),
        PXI=PXI,
        PYI=PYI,
        PZI=PZI
    )


def compute_case(case_name: str, components: List[Component]) -> Chapter4Result:
    S1, cg = compute_cg(components)
    component_contributions = compute_component_contributions(components, cg)
    airplane_axes_inertia = sum_airplane_axes_inertia(component_contributions)
    principal_axes_inertia = compute_principal_axes_inertia(airplane_axes_inertia)

    return Chapter4Result(
        CHAPTER=4,
        CASE_NAME=case_name,
        S1=S1,
        XBAR=cg.XBAR,
        YBAR=cg.YBAR,
        ZBAR=cg.ZBAR,
        IX=airplane_axes_inertia.IX,
        IY=airplane_axes_inertia.IY,
        IZ=airplane_axes_inertia.IZ,
        IXZ=airplane_axes_inertia.IXZ,
        THETA_DEG=principal_axes_inertia.THETA_DEG,
        PXI=principal_axes_inertia.PXI,
        PYI=principal_axes_inertia.PYI,
        PZI=principal_axes_inertia.PZI,
        COMPONENTS=component_contributions
    )


def select_components_by_C(components: List[Component], max_C: int) -> List[Component]:
    return [component for component in components if component.C <= max_C]


def compute_all_cases(cases: Dict[str, List[Component]]) -> tuple[Dict[str, Any], Dict[str, Any]]:
    results = {}
    component_results = {}
    weight_cases = {
        "C0_EMPTY_WEIGHT": 0,
        "C1_MINIMUM_WEIGHT": 1,
        "C2_GROSS_WEIGHT": 2
    }

    for case_name, components in cases.items():
        invalid_C = sorted({component.C for component in components if component.C not in weight_cases.values()})
        if invalid_C:
            raise ValueError(f"Unsupported C value(s) in {case_name}: {invalid_C}")

        results[case_name.upper()] = {}
        component_results[case_name.upper()] = {}
        for weight_case_name, max_C in weight_cases.items():
            selected_components = select_components_by_C(components, max_C)
            result = compute_case(weight_case_name, selected_components)
            result_data = asdict(result)
            component_results[case_name.upper()][weight_case_name] = {
                "CHAPTER": result_data["CHAPTER"],
                "CASE_NAME": result_data["CASE_NAME"],
                "COMPONENTS": result_data.pop("COMPONENTS")
            }
            results[case_name.upper()][weight_case_name] = result_data

    output = {
        "CHAPTER": 4,
        "CASES": results
    }
    component_output = {
        "CHAPTER": 4,
        "CASES": component_results
    }

    return output, component_output


def load_cases(path: Union[str, Path]) -> Dict[str, List[Component]]:
    with open(path, "r", encoding="utf-8-sig") as f:
        data = json.load(f)

    return {
        case_name: [Component(**component) for component in components]
        for case_name, components in data["cases"].items()
    }


def run_chapter4(
    input_path: Union[str, Path] = INPUT_PATH,
    output_path: Union[str, Path] = OUTPUT_PATH,
    component_output_path: Union[str, Path] = COMPONENT_OUTPUT_PATH
) -> Dict[str, Any]:
    output, component_output = compute_all_cases(load_cases(input_path))

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    Path(component_output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(component_output_path, "w", encoding="utf-8") as f:
        json.dump(component_output, f, indent=2, ensure_ascii=False)

    return output


if __name__ == "__main__":
    result = run_chapter4()
    print(json.dumps(result, indent=2, ensure_ascii=False))
