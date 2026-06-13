import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = PROJECT_ROOT / "inputs" / "ch6_input.json"
OUTPUT_PATH = PROJECT_ROOT / "outputs" / "ch6_output.json"


@dataclass
class Chapter6Input:
    CAT: str
    W: float
    S: float
    VH: float
    VSTALL: float
    V5FSTALL: float
    HMACH: float

    VC: float = 1.0
    VD: float = 1.0
    VA: float = 1.0
    VF: float = 1.0
    N1: float = 1.0
    N2: float = 1.0


def positive_load_factor_min(CAT: str, W: float) -> float:
    CAT = CAT.upper()

    if CAT == "N":
        return min(2.1 + 24000.0 / (W + 10000.0), 3.8)

    if CAT == "U":
        return 4.4

    if CAT == "A":
        return 6.0

    raise ValueError("CAT must be one of 'N', 'U', or 'A'")


def negative_load_factor_min(CAT: str, N1: float) -> float:
    CAT = CAT.upper()

    if CAT in ["N", "U"]:
        return -0.4 * N1

    if CAT == "A":
        return -0.5 * N1

    raise ValueError("CAT must be one of 'N', 'U', or 'A'")


def vc_factor(CAT: str, WOS: float) -> float:
    CAT = CAT.upper()

    if CAT in ["N", "U"]:
        K1 = 33.0
        K1_100 = 28.6
    elif CAT == "A":
        K1 = 36.0
        K1_100 = 28.6
    else:
        raise ValueError("CAT must be one of 'N', 'U', or 'A'")

    if WOS <= 20.0:
        return K1

    if WOS >= 100.0:
        return K1_100

    return K1 - (K1 - K1_100) * (WOS - 20.0) / 80.0


def vd_factor(CAT: str, WOS: float) -> float:
    CAT = CAT.upper()

    if CAT == "N":
        K2 = 1.40
        K2_100 = 1.35
    elif CAT == "U":
        K2 = 1.50
        K2_100 = 1.35
    elif CAT == "A":
        K2 = 1.55
        K2_100 = 1.35
    else:
        raise ValueError("CAT must be one of 'N', 'U', or 'A'")

    if WOS <= 20.0:
        return K2

    if WOS >= 100.0:
        return K2_100

    return K2 - (K2 - K2_100) * (WOS - 20.0) / 80.0


def adjust_up_to_minimum(chosen: float, minimum: float) -> float:
    if chosen <= 1.0:
        return minimum

    return max(chosen, minimum)


def adjust_negative_load_factor(N2: float, NEGMAN: float) -> float:
    # STRSPEED.BAS lines 430-440: only less-conservative values are adjusted.
    if N2 > NEGMAN:
        return NEGMAN

    return N2


def mach_limits(VC: float, VD: float, HMACH: float) -> Dict[str, float]:
    T = 59.0 - 0.003566 * HMACH
    A = 29.02 * math.sqrt(T + 459.4)

    if HMACH > 35332.0:
        A = 575.0

    SIGMA = (1.0 - 0.000006879 * HMACH) ** 4.258

    if HMACH > 35332.0:
        SIGMA = (
            0.00072725
            * math.exp(-0.00004778 * (HMACH - 35332.0))
        ) / 0.002378

    VCTRUE = VC / math.sqrt(SIGMA)
    VDTRUE = VD / math.sqrt(SIGMA)
    MC = VCTRUE / A
    MD = VDTRUE / A

    return {
        "T": T,
        "A": A,
        "SIGMA": SIGMA,
        "VCTRUE": VCTRUE,
        "VDTRUE": VDTRUE,
        "MC": MC,
        "MD": MD,
    }


def analyze_chapter6(inp: Chapter6Input) -> Dict[str, Any]:
    CAT = inp.CAT.upper()
    W = inp.W
    S = inp.S
    VH = inp.VH
    VSTALL = inp.VSTALL
    V5FSTALL = inp.V5FSTALL
    HMACH = inp.HMACH

    WOS = W / S

    NMAN = positive_load_factor_min(CAT, W)
    N1 = adjust_up_to_minimum(inp.N1, NMAN)

    NEGMAN = negative_load_factor_min(CAT, N1)
    N2 = adjust_negative_load_factor(inp.N2, NEGMAN)

    K1 = vc_factor(CAT, WOS)
    VCMIN_RAW = K1 * math.sqrt(WOS)
    VCMIN = VCMIN_RAW
    if VCMIN > 0.9 * VH:
        VCMIN = 0.9 * VH

    VC = adjust_up_to_minimum(inp.VC, VCMIN)

    K2 = vd_factor(CAT, WOS)
    VDMIN_BY_FACTOR = K2 * VCMIN
    VDMIN_BY_MARGIN = 1.25 * VC
    VDMIN = VDMIN_BY_FACTOR
    if VDMIN < VDMIN_BY_MARGIN:
        VDMIN = VDMIN_BY_MARGIN

    VD = adjust_up_to_minimum(inp.VD, VDMIN)

    VAMIN = VSTALL * math.sqrt(N1)
    VA = adjust_up_to_minimum(inp.VA, VAMIN)

    VFMIN = max(1.4 * VSTALL, 1.8 * V5FSTALL)
    VF = adjust_up_to_minimum(inp.VF, VFMIN)

    mach = mach_limits(VC, VD, HMACH)

    adjusted_items = []
    if N1 != inp.N1:
        adjusted_items.append("N1")
    if N2 != inp.N2:
        adjusted_items.append("N2")
    if VC != inp.VC:
        adjusted_items.append("VC")
    if VD != inp.VD:
        adjusted_items.append("VD")
    if VA != inp.VA:
        adjusted_items.append("VA")
    if VF != inp.VF:
        adjusted_items.append("VF")

    return {
        "chapter": 6,
        "output": {
            "VC": VC,
            "VD": VD,
            "VA": VA,
            "VF": VF,
            "N1": N1,
            "N2": N2,
            "MC": mach["MC"],
            "MD": mach["MD"],
            "INPUT": {
                "CAT": CAT,
                "W": W,
                "S": S,
                "VH": VH,
                "VSTALL": VSTALL,
                "V5FSTALL": V5FSTALL,
                "HMACH": HMACH,
            },
            "MINIMUMS": {
                "WOS": WOS,
                "K1": K1,
                "VCMIN": VCMIN,
                "VCMIN_RAW": VCMIN_RAW,
                "K2": K2,
                "VDMIN": VDMIN,
                "VDMIN_BY_FACTOR": VDMIN_BY_FACTOR,
                "VDMIN_BY_MARGIN": VDMIN_BY_MARGIN,
                "VAMIN": VAMIN,
                "VFMIN": VFMIN,
                "NMAN": NMAN,
                "NEGMAN": NEGMAN,
            },
            "MACH": mach,
            "adjusted_items": adjusted_items,
        },
    }


def analyze_chapter6_from_json(input_json: Dict[str, Any]) -> Dict[str, Any]:
    raw = input_json["input"]
    return analyze_chapter6(Chapter6Input(**raw))


def run_chapter6(
    input_path: str | Path = INPUT_PATH,
    output_path: str | Path = OUTPUT_PATH,
) -> Dict[str, Any]:
    with open(input_path, "r", encoding="utf-8-sig") as f:
        input_data = json.load(f)

    result = analyze_chapter6_from_json(input_data)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    return result


if __name__ == "__main__":
    result = run_chapter6()
    print(json.dumps(result, indent=2, ensure_ascii=False))
