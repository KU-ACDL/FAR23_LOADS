import json
import math
from dataclasses import dataclass
from typing import Dict, Any


INPUT_PATH = "ch6_input.json"
OUTPUT_PATH = "ch6_output.json"


@dataclass
class Chapter6Input:
    category: str                  # "N", "U", "A"
    weight: float
    wing_area: float
    vh: float
    vs_clean: float
    vs_flap: float
    mach_shoulder_altitude: float

    chosen_vc: float = 1.0
    chosen_vd: float = 1.0
    chosen_va: float = 1.0
    chosen_vf: float = 1.0
    chosen_n_pos: float = 1.0
    chosen_n_neg: float = 1.0


def speed_of_sound(altitude: float) -> float:
    """
    ISA 기준 음속 계산.
    결과 단위: knots EAS 근사값
    """
    gamma = 1.4
    r = 287.05
    t0 = 288.15
    lapse = 0.0065

    altitude_m = altitude * 0.3048
    temperature_k = t0 - lapse * altitude_m

    a_ms = math.sqrt(gamma * r * temperature_k)
    a_kt = a_ms * 1.943844

    return a_kt


def positive_load_factor_min(category: str, weight: float) -> float:
    category = category.upper()

    if category == "N":
        return min(2.1 + 24000.0 / (weight + 10000.0), 3.8)

    if category == "U":
        return 4.4

    if category == "A":
        return 6.0

    raise ValueError("category는 'N', 'U', 'A' 중 하나여야 함.")


def negative_load_factor_min(category: str, n_pos: float) -> float:
    category = category.upper()

    if category in ["N", "U"]:
        return -0.4 * n_pos

    if category == "A":
        return -0.5 * n_pos

    raise ValueError("category는 'N', 'U', 'A' 중 하나여야 함.")


def vc_factor(category: str, wing_loading: float) -> float:
    """
    VC_min = K * sqrt(W/S)
    """
    category = category.upper()

    if category in ["N", "U"]:
        k_low = 33.0
        k_high = 28.6
    elif category == "A":
        k_low = 36.0
        k_high = 28.6
    else:
        raise ValueError("category는 'N', 'U', 'A' 중 하나여야 함.")

    if wing_loading <= 20.0:
        return k_low

    if wing_loading >= 100.0:
        return k_high

    return k_low - (k_low - k_high) * (wing_loading - 20.0) / 80.0


def vd_factor(category: str, wing_loading: float) -> float:
    """
    VD_min = K * VC_min
    """
    category = category.upper()

    if category == "N":
        k_low = 1.40
        k_high = 1.35
    elif category == "U":
        k_low = 1.50
        k_high = 1.35
    elif category == "A":
        k_low = 1.55
        k_high = 1.35
    else:
        raise ValueError("category는 'N', 'U', 'A' 중 하나여야 함.")

    if wing_loading <= 20.0:
        return k_low

    if wing_loading >= 100.0:
        return k_high

    return k_low - (k_low - k_high) * (wing_loading - 20.0) / 80.0


def select_value(chosen: float, minimum: float) -> float:
    """
    BASIC 출력 기준:
    chosen 값이 1이면 사용자가 지정하지 않은 것으로 보고 minimum 사용.
    chosen 값이 minimum보다 작으면 minimum으로 보정.
    """
    if chosen <= 1.0:
        return minimum

    return max(chosen, minimum)

def mach_limits(vc, vd, shoulder_altitude):
    h = shoulder_altitude

    t = 59.0 - 0.003566 * h
    a = 29.02 * math.sqrt(t + 459.4)

    if h > 35332.0:
        a = 575.0

    sigma = (1.0 - 0.000006879 * h) ** 4.258

    if h > 35332.0:
        sigma = (
            0.00072725
            * math.exp(-0.00004778 * (h - 35332.0))
        ) / 0.002378

    vc_true = vc / math.sqrt(sigma)
    vd_true = vd / math.sqrt(sigma)

    mc = vc_true / a
    md = vd_true / a

    return mc, md


def analyze_chapter6(inp: Chapter6Input) -> Dict[str, Any]:
    category = inp.category.upper()

    wing_loading = inp.weight / inp.wing_area

    # 1. Load factor
    n_pos_min = positive_load_factor_min(category, inp.weight)
    n_pos = select_value(inp.chosen_n_pos, n_pos_min)

    n_neg_min = negative_load_factor_min(category, n_pos)

    if inp.chosen_n_neg <= 1.0:
        n_neg = n_neg_min
    else:
        n_neg = min(inp.chosen_n_neg, n_neg_min)

    # 2. VC
    vc_min_raw = vc_factor(category, wing_loading) * math.sqrt(wing_loading)
    vc_min = min(vc_min_raw, 0.9 * inp.vh)
    vc = select_value(inp.chosen_vc, vc_min)

    # 3. VD
    vd_min_by_factor = vd_factor(category, wing_loading) * vc_min
    vd_min_by_margin = 1.25 * vc
    vd_min = max(vd_min_by_factor, vd_min_by_margin)
    vd = select_value(inp.chosen_vd, vd_min)

    # 4. VA
    va_min = inp.vs_clean * math.sqrt(n_pos)
    va = select_value(inp.chosen_va, va_min)
    va = min(va, vc)

    # 5. VF
    vf_min = max(
        1.4 * inp.vs_clean,
        1.8 * inp.vs_flap
    )
    vf = select_value(inp.chosen_vf, vf_min)

    # 6. Mach limit
    mc, md = mach_limits(vc, vd, inp.mach_shoulder_altitude)

    adjusted_items = []

    if n_pos != inp.chosen_n_pos:
        adjusted_items.append("n_positive")

    if n_neg != inp.chosen_n_neg:
        adjusted_items.append("n_negative")

    if vc != inp.chosen_vc:
        adjusted_items.append("vc")

    if vd != inp.chosen_vd:
        adjusted_items.append("vd")

    if va != inp.chosen_va:
        adjusted_items.append("va")

    if vf != inp.chosen_vf:
        adjusted_items.append("vf")

    return {
        "chapter": 6,
        "output": {
            "category": category,
            "wing_loading": wing_loading,
            "load_factors": {
                "n_positive": n_pos,
                "n_negative": n_neg
            },
            "speeds": {
                "vc": vc,
                "vd": vd,
                "va": va,
                "vf": vf
            },
            "mach_limits": {
                "mc": mc,
                "md": md,
                "shoulder_altitude": inp.mach_shoulder_altitude
            },
            "adjusted_items": adjusted_items
        }
    }

def analyze_chapter6_from_json(input_json: Dict[str, Any]) -> Dict[str, Any]:
    return analyze_chapter6(Chapter6Input(**input_json["input"]))


def run_chapter6(
    input_path: str = INPUT_PATH,
    output_path: str = OUTPUT_PATH
) -> Dict[str, Any]:
    with open(input_path, "r", encoding="utf-8-sig") as f:
        input_data = json.load(f)

    result = analyze_chapter6_from_json(input_data)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    return result


if __name__ == "__main__":
    result = run_chapter6()
    print(json.dumps(result, indent=2, ensure_ascii=False))
