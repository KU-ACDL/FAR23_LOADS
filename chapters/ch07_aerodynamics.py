import json
import math
from copy import deepcopy
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = PROJECT_ROOT / "inputs" / "ch7_input.json"
CH5_GEOMETRY_PATH = PROJECT_ROOT / "outputs" / "cases" / "ch5_geometry_output.json"
CH8_OUTPUT_PATH = PROJECT_ROOT / "outputs" / "cases" / "ch8_output.json"
OUTPUT_PATH = PROJECT_ROOT / "outputs" / "cases" / "ch7_aerodynamics_output.json"
REPORT_PATH = PROJECT_ROOT / "outputs" / "reports" / "ch7_aerodynamics_report.txt"


PI_BASIC = 3.1416
DEG_PER_RAD = 57.3
# Appendix A uses sea-level Reynolds number at 70 mph. The restored BASIC
# listing omits this inverse kinematic-viscosity factor in line 2155.
REYNOLDS_PER_FOOT_SECOND = 6380.0


class Chapter7InputError(ValueError):
    pass


def select_airloads_method(inp: Dict[str, Any]) -> Dict[str, Any]:
    mach = _as_float(_get(inp, "MACH", 0.0, required=False), "MACH")
    sweep = _as_float(_get(inp, "SBA", 0.0, required=False), "SBA")

    if mach < 0 or mach >= 1:
        raise Chapter7InputError("MACH must be between 0 and 1")

    use_airload4 = abs(sweep) > 15.0 or mach >= 0.4
    gc_input = _get(inp, "GC", None, required=False)
    if gc_input is None:
        gc = 1.0 / math.sqrt(1.0 - mach ** 2) if mach >= 0.4 else 1.0
    else:
        gc = _as_float(gc_input, "GC")
    if gc < 1.0:
        raise Chapter7InputError("GC must be at least 1.0")

    return {
        "METHOD": "AIRLOAD4.BAS" if use_airload4 else "AIRLOADS.BAS",
        "MACH": mach,
        "SBA": sweep,
        "GC": gc,
        "USE_AIRLOAD4": use_airload4,
    }


def _get_group(data: Dict[str, Any], *names: str, required: bool = True) -> Dict[str, Any]:
    for name in names:
        if name in data and data[name] is not None:
            return data[name]
    if required:
        raise Chapter7InputError(f"Missing input group: {names[0]}")
    return {}


def _get(group: Dict[str, Any], key: str, default: Any = None, required: bool = True) -> Any:
    if key in group:
        return group[key]
    if required:
        raise Chapter7InputError(f"Missing input variable: {key}")
    return default


def _as_float_list(values: List[Any], key: str) -> List[float]:
    if values is None:
        raise Chapter7InputError(f"Missing list: {key}")
    return [float(v) for v in values]


def _as_float(value: Any, key: str) -> float:
    try:
        return float(value)
    except Exception as exc:
        raise Chapter7InputError(f"Invalid float for {key}: {value}") from exc


def _as_int(value: Any, key: str) -> int:
    try:
        return int(value)
    except Exception as exc:
        raise Chapter7InputError(f"Invalid integer for {key}: {value}") from exc


def _left1(value: Any) -> str:
    return str(value).strip().upper()[:1]


def _check_len(values: List[Any], expected: int, key: str) -> None:
    if len(values) != expected:
        raise Chapter7InputError(f"{key} length must be {expected}, got {len(values)}")


def _interp_linear(x: float, stations: List[float], values: List[float]) -> float:
    """
    AIRLOADS.BAS의 station별 선형 보간을 Python으로 옮긴 함수.
    입력 station 범위를 약간 벗어나면 첫 구간 또는 마지막 구간으로 외삽한다.
    """
    if len(stations) != len(values):
        raise Chapter7InputError("stations and values length mismatch")
    if len(stations) < 2:
        raise Chapter7InputError("at least two stations are required")

    if x <= stations[0]:
        i = 0
    elif x >= stations[-1]:
        i = len(stations) - 2
    else:
        i = 0
        for k in range(len(stations) - 1):
            if stations[k] <= x <= stations[k + 1]:
                i = k
                break

    x0 = stations[i]
    x1 = stations[i + 1]
    y0 = values[i]
    y1 = values[i + 1]

    if x1 == x0:
        raise Chapter7InputError("duplicate station causes division by zero")

    return (x - x0) * (y1 - y0) / (x1 - x0) + y0


def calculate_tau(inp: Dict[str, Any]) -> Dict[str, float]:
    """TAU.BAS correction for non-elliptical wing planform."""
    TAPR = _as_float(_get(inp, "TAPR"), "TAPR")
    TIPR = _as_float(_get(inp, "TIPR"), "TIPR")

    if TAPR < 0:
        raise Chapter7InputError("TAPR must be non-negative")
    if TIPR < 0 or TIPR > 1:
        raise Chapter7InputError("TIPR must be between 0 and 1")

    TAU0 = (
        0.206209
        - 1.26146 * TAPR
        + 3.05385 * TAPR ** 2
        - 2.8027 * TAPR ** 3
        + 0.976801 * TAPR ** 4
    )
    TAU1 = (
        0.112203
        - 0.577843 * TAPR
        + 1.08306 * TAPR ** 2
        - 0.596856 * TAPR ** 3
        + 0.194241 * TAPR ** 4
    )
    TAU2 = (
        0.0302789
        + 0.0294027 * TAPR
        - 0.470926 * TAPR ** 2
        + 0.880983 * TAPR ** 3
        - 0.394766 * TAPR ** 4
    )
    TAU3 = 0.0

    if TIPR <= 0.1:
        TAU = TAU0 + TIPR * (TAU1 - TAU0) / 0.1
    elif TIPR <= 0.2:
        TAU = TAU1 + (TIPR - 0.1) * (TAU2 - TAU1) / 0.1
    else:
        TAU = TAU2 + (TIPR - 0.2) * (TAU3 - TAU2) / 0.8

    return {
        "TAPR": TAPR,
        "TIPR": TIPR,
        "TAU0": TAU0,
        "TAU1": TAU1,
        "TAU2": TAU2,
        "TAU3": TAU3,
        "TAU": TAU,
    }


def _polyfit_basic_style(XDATA: List[float], YDATA: List[float], NP: int = 4) -> List[float]:
    """
    AIRLOADS.BAS 4090~4490의 정규방정식 기반 polynomial fit.
    반환값은 [constant, 1st, 2nd, ...] 순서.
    """
    HH = len(XDATA)
    if HH < NP + 1:
        raise Chapter7InputError(f"Polynomial fit requires at least {NP + 1} points")

    SX99 = [0.0 for _ in range(NP * 2 + 1)]
    SYX = [0.0 for _ in range(NP + 1)]

    for j in range(NP * 2 + 1):
        SX99[j] = sum(x ** j for x in XDATA)

    for j in range(NP + 1):
        SYX[j] = sum(y * (x ** j) for x, y in zip(XDATA, YDATA))

    n = NP + 1
    A2 = [[0.0 for _ in range(n + 1)] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            A2[i][j] = SX99[i + j]
        A2[i][n] = SYX[i]

    # Gaussian elimination with partial pivoting for numerical safety.
    for k in range(n - 1):
        pivot_row = max(range(k, n), key=lambda r: abs(A2[r][k]))
        if abs(A2[pivot_row][k]) < 1e-15:
            raise Chapter7InputError("Singular matrix in polynomial fit")
        if pivot_row != k:
            A2[k], A2[pivot_row] = A2[pivot_row], A2[k]

        pivot = A2[k][k]
        for j in range(k, n + 1):
            A2[k][j] /= pivot

        for i in range(k + 1, n):
            factor = A2[i][k]
            for j in range(k, n + 1):
                A2[i][j] -= factor * A2[k][j]

    if abs(A2[n - 1][n - 1]) < 1e-15:
        raise Chapter7InputError("Singular matrix in polynomial fit")

    coef = [0.0 for _ in range(n)]
    coef[n - 1] = A2[n - 1][n] / A2[n - 1][n - 1]

    for k in range(n - 2, -1, -1):
        value = A2[k][n]
        for j in range(k + 1, n):
            value -= A2[k][j] * coef[j]
        coef[k] = value

    return coef


@dataclass
class Chapter7State:
    # Wing geometry calculation data
    H: int
    B: float
    DY: float
    YE: List[float]
    XF: List[float]
    XA: List[float]
    CE: List[float]
    DA: List[float]
    CX25: List[float]
    C50X: List[float]
    A: float
    AREA_TOTAL: float
    SC2: float
    SAYE: float
    SBARXC: float
    MAC: float
    YBAR: float
    XMACLE: float
    ZAR: float

    # Filled by later calculations
    MO: Optional[List[float]] = None
    SMOCDY: Optional[float] = None
    BARMO: Optional[float] = None
    S: Optional[float] = None
    S6: Optional[float] = None
    CCLA1: Optional[List[float]] = None
    CLA1: Optional[List[float]] = None
    CHECKCL1: Optional[float] = None

    REFANG: Optional[List[float]] = None
    WSUM: Optional[float] = None
    XSUM: Optional[float] = None
    AWO: Optional[float] = None
    AO: Optional[List[float]] = None
    FCLB: Optional[List[float]] = None
    VCCLB: Optional[List[float]] = None
    L: Optional[float] = None
    YLIN: Optional[float] = None
    YLOB: Optional[float] = None
    AVECLB: Optional[float] = None
    THETA_FAIR: Optional[List[Optional[float]]] = None
    UCLB: Optional[List[float]] = None
    BCHECK: Optional[float] = None

    AA: Optional[List[float]] = None
    R3N: Optional[List[float]] = None
    C3LMAX: Optional[List[float]] = None
    CM: Optional[List[float]] = None
    WING_STALL_CL: Optional[float] = None
    STALLCL: Optional[float] = None

    MM: Optional[float] = None
    ALPHA: Optional[float] = None
    A_local: Optional[List[float]] = None
    A1: Optional[List[float]] = None
    KCL: Optional[List[float]] = None
    CID: Optional[List[float]] = None
    CPD: Optional[List[float]] = None
    ZCOEFM: Optional[List[float]] = None
    CD: Optional[List[float]] = None
    G1: Optional[float] = None
    G2: Optional[float] = None
    G3: Optional[float] = None
    G7: Optional[float] = None
    G4CLW: Optional[float] = None
    G5CDW: Optional[float] = None
    G6CMW: Optional[float] = None
    ANRW2WL: Optional[float] = None
    AIRLOAD_METHOD: str = "AIRLOADS.BAS"
    GC: float = 1.0
    SBA: float = 0.0
    USE_AIRLOAD4: bool = False


def load_wing_geometry(inp: Dict[str, Any]) -> Chapter7State:
    element_station = _get_group(
        _get_group(inp, "Intermediate"),
        "Element Station",
    )
    surface_sum = _get_group(
        _get_group(inp, "Intermediate"),
        "Surface Sum",
    )
    surface_output = _get_group(inp, "Surface Output")
    element_output = _get(inp, "Element Output")
    if not isinstance(element_output, list):
        raise Chapter7InputError("Element Output must be a JSON array")

    H = _as_int(_get(inp, "H"), "H")
    if H < 2 or H > 100:
        raise Chapter7InputError("H must be between 2 and 100")
    if len(element_output) != H:
        raise Chapter7InputError(
            f"Element Output length must be {H}, got {len(element_output)}"
        )

    DY = _as_float(_get(element_station, "DY"), "DY")
    A = _as_float(_get(surface_output, "A"), "A")
    MAC = _as_float(_get(surface_output, "MAC"), "MAC")
    YBAR = _as_float(_get(surface_output, "YBAR"), "YBAR")
    XMACLE = _as_float(_get(surface_output, "XMACLE"), "XMACLE")
    ZAR = _as_float(_get(surface_output, "ZAR"), "ZAR")
    SC2 = _as_float(_get(surface_sum, "SC2"), "SC2")
    SAYE = _as_float(_get(surface_sum, "SAYE"), "SAYE")
    SBARXC = _as_float(_get(surface_sum, "SBARXC"), "SBARXC")

    XF = [_as_float(_get(row, "XF(I)"), "XF(I)") for row in element_output]
    XA = [_as_float(_get(row, "XA(I)"), "XA(I)") for row in element_output]
    YE = [_as_float(_get(row, "YE(I)"), "YE(I)") for row in element_output]
    CE = [_as_float(_get(row, "CE(I)"), "CE(I)") for row in element_output]
    DA = [
        _as_float(_get(row, "CE(I)*DY"), "CE(I)*DY")
        for row in element_output
    ]

    AREA_TOTAL = 2.0 * A
    B = math.sqrt(ZAR * AREA_TOTAL)
    CX25 = [xf + 0.25 * ce for xf, ce in zip(XF, CE)]
    C50X = [xf + 0.5 * ce for xf, ce in zip(XF, CE)]

    return Chapter7State(
        H=H,
        B=B,
        DY=DY,
        YE=YE,
        XF=XF,
        XA=XA,
        CE=CE,
        DA=DA,
        CX25=CX25,
        C50X=C50X,
        A=A,
        AREA_TOTAL=AREA_TOTAL,
        SC2=SC2,
        SAYE=SAYE,
        SBARXC=SBARXC,
        MAC=MAC,
        YBAR=YBAR,
        XMACLE=XMACLE,
        ZAR=ZAR,
    )


def calculate_additive_lift_distribution(state: Chapter7State, inp: Dict[str, Any]) -> Dict[str, Any]:
    NO = _as_int(_get(inp, "NO"), "NO")
    WS = _as_float_list(_get(inp, "WS"), "WS")
    M_input = _as_float_list(_get(inp, "M"), "M")
    GC = _as_float(_get(inp, "GC", 1.0, required=False), "GC")
    use_airload4 = bool(_get(inp, "USE_AIRLOAD4", False, required=False))
    M_corrected = [value * GC for value in M_input]

    _check_len(WS, NO, "WS")
    _check_len(M_corrected, NO, "M")

    MO = [_interp_linear(y, WS, M_corrected) for y in state.YE]
    SMOCDY = sum(MO[j] * state.DA[j] for j in range(state.H))
    BARMO = SMOCDY / state.A
    S = 2.0 * state.A

    CCLA1 = []
    CLA1 = []
    S6 = 0.0

    for k in range(state.H):
        ellipse_term = (4.0 * S / (PI_BASIC * state.B)) * math.sqrt(max(0.0, 1.0 - (2.0 * state.YE[k] / state.B) ** 2))
        CCLA1_K = 0.5 * (MO[k] * state.CE[k] / BARMO + ellipse_term)
        CLA1_K = CCLA1_K / state.CE[k]
        CCLA1.append(CCLA1_K)
        CLA1.append(CLA1_K)
        S6 += CLA1_K * state.DA[k]

    CHECKCL1 = S6 / state.A
    if use_airload4 and abs(CHECKCL1 - 1.0) >= 0.0001:
        CCLA1 = [value / CHECKCL1 for value in CCLA1]
        CLA1 = [value / CHECKCL1 for value in CLA1]
        S6 = sum(CLA1[j] * state.DA[j] for j in range(state.H))
        CHECKCL1 = S6 / state.A

    state.MO = MO
    state.SMOCDY = SMOCDY
    state.BARMO = BARMO
    state.S = S
    state.S6 = S6
    state.CCLA1 = CCLA1
    state.CLA1 = CLA1
    state.CHECKCL1 = CHECKCL1

    return {
        "GC": GC,
        "M_CORRECTED(I)": M_corrected,
        "MO(I)": MO,
        "SMOCDY": SMOCDY,
        "BARMO": BARMO,
        "S": S,
        "S6": S6,
        "CCLA1(I)": CCLA1,
        "CLA1(I)": CLA1,
        "CHECKCL1": CHECKCL1,
    }


def calculate_basic_lift_distribution(state: Chapter7State, inp: Dict[str, Any]) -> Dict[str, Any]:
    if state.MO is None:
        raise Chapter7InputError("Additive lift distribution must be calculated first")

    NN = _as_int(_get(inp, "NN"), "NN")
    ZWS = _as_float_list(_get(inp, "ZWS"), "ZWS")
    RANG = _as_float_list(_get(inp, "RANG"), "RANG")
    YDIS = _as_float(_get(inp, "YDIS"), "YDIS")

    _check_len(ZWS, NN, "ZWS")
    _check_len(RANG, NN, "RANG")

    REFANG = [_interp_linear(y, ZWS, RANG) for y in state.YE]
    WSUM = 0.0
    XSUM = 0.0

    for l in range(state.H):
        WCEARMO = REFANG[l] * state.MO[l] * state.DA[l]
        XCMO = state.MO[l] * state.DA[l]
        WSUM += WCEARMO
        XSUM += XCMO

    AWO = WSUM / XSUM

    AO = []
    FCLB = []
    VCCLB = []
    for i in range(state.H):
        AO_I = REFANG[i] - AWO
        FCLB_I = AO_I * state.MO[i] / 2.0
        VCCLB_I = state.CE[i] * FCLB_I
        AO.append(AO_I)
        FCLB.append(FCLB_I)
        VCCLB.append(VCCLB_I)

    L = None
    YLIN = None
    YLOB = None
    AVECLB = None
    THETA_FAIR: List[Optional[float]] = [None for _ in range(state.H)]
    UCLB = [0.0 for _ in range(state.H)]

    if YDIS == 0:
        UCLB = FCLB.copy()
    else:
        if (state.B / 2.0 - YDIS) > YDIS:
            L = (state.B / 2.0 - YDIS) / 2.0
        else:
            L = YDIS / 2.0

        if L == 0:
            raise Chapter7InputError("YDIS makes discontinuity fairing length zero")

        YLIN = YDIS - L
        YLOB = YDIS + L

        for i in range(state.H - 1):
            if state.YE[i] <= YDIS <= state.YE[i + 1]:
                AVECLB = (VCCLB[i] + VCCLB[i + 1]) / 2.0
                break

        if AVECLB is None:
            raise Chapter7InputError("YDIS is not located between adjacent YE elements")

        for i in range(state.H):
            if YLIN <= state.YE[i] <= YLOB:
                THETA = PI_BASIC * (state.YE[i] - YLIN) / (2.0 * L)
                VCCLB[i] = (VCCLB[i] - AVECLB) * abs(math.cos(THETA)) + AVECLB
                THETA_FAIR[i] = THETA
            UCLB[i] = VCCLB[i] / state.CE[i]

    S7 = sum(UCLB[i] * state.DA[i] for i in range(state.H))
    BCHECK = S7 / state.A

    state.REFANG = REFANG
    state.WSUM = WSUM
    state.XSUM = XSUM
    state.AWO = AWO
    state.AO = AO
    state.FCLB = FCLB
    state.VCCLB = VCCLB
    state.L = L
    state.YLIN = YLIN
    state.YLOB = YLOB
    state.AVECLB = AVECLB
    state.THETA_FAIR = THETA_FAIR
    state.UCLB = UCLB
    state.BCHECK = BCHECK

    return {
        "REFANG(I)": REFANG,
        "WSUM": WSUM,
        "XSUM": XSUM,
        "AWO": AWO,
        "AO(I)": AO,
        "FCLB(I)": FCLB,
        "VCCLB(I)": VCCLB,
        "L": L,
        "YLIN": YLIN,
        "YLOB": YLOB,
        "AVECLB": AVECLB,
        "THETA(I)": THETA_FAIR,
        "UCLB(I)": UCLB,
        "BCHECK": BCHECK,
    }


def calculate_stall_cl(state: Chapter7State, inp: Dict[str, Any]) -> Dict[str, Any]:
    if state.UCLB is None or state.CLA1 is None:
        raise Chapter7InputError("Basic and additive lift distributions must be calculated first")

    MR = _as_int(_get(inp, "MR"), "MR")
    WWS = _as_float_list(_get(inp, "WWS"), "WWS")
    C1LMAX = _as_float_list(_get(inp, "C1LMAX"), "C1LMAX")
    R1N = _as_float_list(_get(inp, "R1N"), "R1N")
    C2LMAX = _as_float_list(_get(inp, "C2LMAX"), "C2LMAX")
    R2N = _as_float_list(_get(inp, "R2N"), "R2N")
    CHORD = _as_float_list(_get(inp, "CHORD"), "CHORD")

    for key, values in {
        "WWS": WWS,
        "C1LMAX": C1LMAX,
        "R1N": R1N,
        "C2LMAX": C2LMAX,
        "R2N": R2N,
        "CHORD": CHORD,
    }.items():
        _check_len(values, MR, key)

    AA = []
    R3N = []
    C3LMAX = []

    for i in range(MR):
        AA_I = (math.log(C1LMAX[i]) - math.log(C2LMAX[i])) / (math.log(R1N[i]) - math.log(R2N[i]))
        R3N_I = 70.0 * 1.467 * CHORD[i] / 12.0 * REYNOLDS_PER_FOOT_SECOND
        C3LMAX_I = C2LMAX[i] * (R3N_I / R2N[i]) ** AA_I
        AA.append(AA_I)
        R3N.append(R3N_I)
        C3LMAX.append(C3LMAX_I)

    CM = [_interp_linear(y, WWS, C3LMAX) for y in state.YE]

    QCL = 0.2
    STALLCL = 0.0
    KCL = [0.0 for _ in range(state.H)]

    # BASIC 코드처럼 0.01씩 증가시키며 첫 stall 조건을 찾는다.
    guard = 0
    while True:
        guard += 1
        if guard > 10000:
            raise Chapter7InputError("STALL CL iteration did not converge")

        STALLCL = 0.0
        for j in range(state.H):
            KCL[j] = state.UCLB[j] + QCL * state.CLA1[j]
            if KCL[j] >= CM[j]:
                STALLCL = QCL

        if STALLCL > 0.2:
            break
        QCL += 0.01

    state.AA = AA
    state.R3N = R3N
    state.C3LMAX = C3LMAX
    state.CM = CM
    state.WING_STALL_CL = QCL
    state.STALLCL = STALLCL
    state.KCL = KCL.copy()

    return {
        "AA(I)": AA,
        "R3N(I)": R3N,
        "C3LMAX(I)": C3LMAX,
        "CM(I)": CM,
        "WING_STALL_CL": QCL,
        "QCL": QCL,
        "KCL(I)": KCL.copy(),
        "STALLCL": STALLCL,
    }


def _calculate_kcl_cid_cd_for_qcl(
    state: Chapter7State,
    QCL: float,
    GC: Optional[float] = None,
    SBA: Optional[float] = None,
    use_airload4: Optional[bool] = None,
) -> Tuple[float, List[float], List[float], List[float], List[float], List[float]]:
    if state.BARMO is None or state.UCLB is None or state.CLA1 is None or state.MO is None or state.REFANG is None or state.AWO is None or state.CPD is None:
        raise Chapter7InputError("Required distributions are missing")
    if state.MM is None:
        raise Chapter7InputError("MM is missing")

    GC = state.GC if GC is None else GC
    SBA = state.SBA if SBA is None else SBA
    use_airload4 = state.USE_AIRLOAD4 if use_airload4 is None else use_airload4

    ALPHA = QCL / (state.MM / DEG_PER_RAD)
    KCL = []
    A_local = []
    A1 = []
    CID = []
    CD = []

    for j in range(state.H):
        KCL_J = state.UCLB[j] + QCL * state.CLA1[j]
        A_J = ALPHA - state.AWO + state.REFANG[j]
        A1_J = A_J - KCL_J / state.MO[j]
        CID_J = KCL_J * A1_J / DEG_PER_RAD
        if use_airload4:
            CID_J *= GC
        CD_J = state.CPD[j] + CID_J

        KCL.append(KCL_J)
        A_local.append(A_J)
        A1.append(A1_J)
        CID.append(CID_J)
        CD.append(CD_J)

    if use_airload4 and abs(SBA) > 15.0:
        col19 = []
        for j in range(state.H):
            col16 = state.CE[j] * KCL[j] / (state.MAC * QCL)
            col18 = (1.0 - 2.0 * state.YE[j] / state.B) ** 2 * (
                1.0 - math.cos(SBA / DEG_PER_RAD)
            )
            col19.append(col16 - col18)

        cl_col19 = sum(
            col19[j] * state.DA[j] for j in range(state.H)
        ) / state.A
        if abs(cl_col19) < 1e-15:
            raise Chapter7InputError("AIRLOAD4 sweep correction normalization is zero")

        for j in range(state.H):
            col20 = col19[j] / cl_col19
            KCL[j] = col20 * state.MAC / state.CE[j] * QCL
            A1[j] = A_local[j] - KCL[j] / state.MO[j]
            CID[j] = KCL[j] * A1[j] / DEG_PER_RAD
            CD[j] = state.CPD[j] + CID[j]

    return ALPHA, KCL, A_local, A1, CID, CD


def calculate_wing_aero_distribution(state: Chapter7State, inp: Dict[str, Any]) -> Dict[str, Any]:
    if state.BARMO is None or state.UCLB is None or state.CLA1 is None or state.MO is None or state.REFANG is None or state.AWO is None:
        raise Chapter7InputError("Previous lift distribution calculations must be completed")

    TAU = _as_float(_get(inp, "TAU"), "TAU")
    GC = _as_float(_get(inp, "GC", 1.0, required=False), "GC")
    SBA = _as_float(_get(inp, "SBA", 0.0, required=False), "SBA")
    use_airload4 = bool(_get(inp, "USE_AIRLOAD4", False, required=False))
    QCL = _as_float(_get(inp, "QCL"), "QCL")
    MZ = _as_int(_get(inp, "MZ"), "MZ")
    WIS = _as_float_list(_get(inp, "WIS"), "WIS")
    COD = _as_float_list(_get(inp, "COD"), "COD")
    MC = _as_int(_get(inp, "MC"), "MC")
    STA = _as_float_list(_get(inp, "STA"), "STA")
    MMC = _as_float_list(_get(inp, "MMC"), "MMC")

    _check_len(WIS, MZ, "WIS")
    _check_len(COD, MZ, "COD")
    _check_len(STA, MC, "STA")
    _check_len(MMC, MC, "MMC")

    if TAU < 0 or TAU > 1:
        raise Chapter7InputError("TAU must be between 0 and 1")

    # Appendix A results place the Peery TAU correction in the finite-wing
    # denominator. This reproduces M8 = 0.080358 for the example airplane.
    MM = (state.BARMO * DEG_PER_RAD) / (
        1.0 + (state.BARMO * DEG_PER_RAD) * (1.0 + TAU) / (PI_BASIC * state.ZAR)
    )
    state.MM = MM

    if use_airload4:
        COD = [value * GC for value in COD]
        MMC = [value * GC for value in MMC]

    CPD = [_interp_linear(y, WIS, COD) for y in state.YE]
    ZCOEFM = [_interp_linear(y, STA, MMC) for y in state.YE]
    state.CPD = CPD
    state.ZCOEFM = ZCOEFM

    ALPHA, KCL, A_local, A1, CID, CD = _calculate_kcl_cid_cd_for_qcl(
        state, QCL, GC, SBA, use_airload4
    )

    G1 = sum(KCL[j] * state.DA[j] for j in range(state.H))
    G2 = sum(CD[j] * state.DA[j] for j in range(state.H))
    G3 = sum(ZCOEFM[j] * state.CE[j] * state.DA[j] for j in range(state.H))
    G7 = state.SC2

    G4CLW = G1 / state.A
    G5CDW = G2 / state.A
    G6CMW = G3 / G7
    ANRW2WL = ALPHA - state.AWO

    state.ALPHA = ALPHA
    state.KCL = KCL
    state.A_local = A_local
    state.A1 = A1
    state.CID = CID
    state.CD = CD
    state.G1 = G1
    state.G2 = G2
    state.G3 = G3
    state.G7 = G7
    state.G4CLW = G4CLW
    state.G5CDW = G5CDW
    state.G6CMW = G6CMW
    state.ANRW2WL = ANRW2WL
    state.AIRLOAD_METHOD = "AIRLOAD4.BAS" if use_airload4 else "AIRLOADS.BAS"
    state.GC = GC
    state.SBA = SBA
    state.USE_AIRLOAD4 = use_airload4

    return {
        "METHOD": "AIRLOAD4.BAS" if use_airload4 else "AIRLOADS.BAS",
        "GC": GC,
        "SBA": SBA,
        "TAU": TAU,
        "MM": MM,
        "ALPHA": ALPHA,
        "A(I)": A_local,
        "A1(I)": A1,
        "KCL(I)": KCL,
        "CID(I)": CID,
        "CPD(I)": CPD,
        "ZCOEFM(I)": ZCOEFM,
        "CD(I)": CD,
        "G1": G1,
        "G2": G2,
        "G3": G3,
        "G7": G7,
        "G4CLW": G4CLW,
        "G5CDW": G5CDW,
        "G6CMW": G6CMW,
        "ANRW2WL": ANRW2WL,
    }


def calculate_landing_gear(state: Chapter7State, inp: Dict[str, Any], S_ft2: float) -> Dict[str, Any]:
    ANOSEGEAR = _as_float(_get(inp, "ANOSEGEAR"), "ANOSEGEAR")
    AMAINGEAR = _as_float(_get(inp, "AMAINGEAR"), "AMAINGEAR")
    XNG = _as_float(_get(inp, "XNG"), "XNG")
    ZNG = _as_float(_get(inp, "ZNG"), "ZNG")
    XMG = _as_float(_get(inp, "XMG"), "XMG")
    ZMG = _as_float(_get(inp, "ZMG"), "ZMG")
    CDNOSEGEAR = _as_float(_get(inp, "CDNOSEGEAR"), "CDNOSEGEAR")
    CDMAINGEAR = _as_float(_get(inp, "CDMAINGEAR"), "CDMAINGEAR")
    ZCG = _as_float(_get(inp, "ZCG"), "ZCG")

    CDNG = CDNOSEGEAR * ANOSEGEAR / S_ft2
    CDMG = CDMAINGEAR * AMAINGEAR / S_ft2
    CDLG = CDNG + CDMG

    # XNG, XMG는 원 BASIC에서도 입력받지만 moment 계산에는 사용하지 않는다.
    _ = XNG
    _ = XMG

    CMNG = -CDNOSEGEAR * ANOSEGEAR * (ZCG / 12.0 - ZNG / 12.0) / (state.MAC / 12.0 * S_ft2)
    CMMG = -CDMAINGEAR * AMAINGEAR * (ZCG / 12.0 - ZMG / 12.0) / (state.MAC / 12.0 * S_ft2)
    CMLG = CMNG + CMMG

    return {
        "CDNG": CDNG,
        "CDMG": CDMG,
        "CDLG": CDLG,
        "CMNG": CMNG,
        "CMMG": CMMG,
        "CMLG": CMLG,
    }


def calculate_airplane_less_tail(
    state: Chapter7State,
    inp: Dict[str, Any],
    landing_gear_inp: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    if state.MM is None or state.CPD is None or state.ZCOEFM is None or state.AWO is None:
        raise Chapter7InputError("Wing aero distribution must be calculated first")

    WF = _as_float(_get(inp, "WF"), "WF")
    LF = _as_float(_get(inp, "LF"), "LF")
    PC = _as_float(_get(inp, "PC"), "PC")
    FA = _as_float(_get(inp, "FA"), "FA")
    AF = _as_float(_get(inp, "AF"), "AF")
    THETA = _as_float(_get(inp, "THETA"), "THETA")
    LGEXT = _left1(_get(inp, "LGEXT", _get(inp, "LGEXT$", "N", required=False), required=False))
    TA = _as_float(_get(inp, "TA"), "TA")
    NEGST = _as_float(_get(inp, "NEGST"), "NEGST")
    ST = _as_float(_get(inp, "ST"), "ST")
    DELTA = _as_float(_get(inp, "DELTA"), "DELTA")

    if DELTA <= 0:
        raise Chapter7InputError("DELTA must be positive")

    KF = 3.338167e-03 + 4.358083e-05 * PC + 5.75194e-07 * PC ** 2 + 1.850009e-07 * PC ** 3
    S_ft2 = 2.0 * state.A / 144.0
    SLOPCM = FA * KF * WF ** 2 * LF / (S_ft2 * state.MAC * state.MM / 12.0 / DEG_PER_RAD)

    AD = 0.105 * AF
    CFD = AD / S_ft2
    AT = 0.012 * 2.0 * TA
    CTD = AT / S_ft2

    landing = {
        "CDNG": 0.0,
        "CDMG": 0.0,
        "CDLG": 0.0,
        "CMNG": 0.0,
        "CMMG": 0.0,
        "CMLG": 0.0,
    }

    if LGEXT == "Y":
        if not landing_gear_inp:
            raise Chapter7InputError("Landing Gear Option Data is required when LGEXT$ is 'Y'")
        landing = calculate_landing_gear(state, landing_gear_inp, S_ft2)

    CDLG = landing["CDLG"]
    CMLG = landing["CMLG"]

    curve_rows = []
    XDATA = []
    YDATA = []

    qcl_values = []
    q = NEGST
    # BASIC FOR loop에 맞춰 end 포함 처리.
    while q <= ST + DELTA * 1e-9:
        qcl_values.append(q)
        q += DELTA

    last_G4CLW = None
    last_G5CDW = None
    last_G6CMW = None
    last_U8CMF0 = None
    last_G8CMF = None
    last_G9CDF = None
    last_ANRW2WL = None

    for QCL in qcl_values:
        ALPHA, KCL, A_local, A1, CID, CD = _calculate_kcl_cid_cd_for_qcl(state, QCL)
        ANRW2WL = ALPHA - state.AWO

        G1 = sum(KCL[j] * state.DA[j] for j in range(state.H))
        G2 = sum(CD[j] * state.DA[j] for j in range(state.H))
        G3 = sum(
            state.ZCOEFM[j] * state.CE[j] * state.DA[j]
            for j in range(state.H)
        )
        G7 = state.SC2

        G4CLW = G1 / state.A
        G5CDW = G2 / state.A
        G6CMW = G3 / G7
        U8CMF0 = THETA * SLOPCM * (state.MM / DEG_PER_RAD)
        G8CMF = SLOPCM * QCL
        G9CDF = CFD + CTD + CDLG

        CDAPT = G9CDF + G5CDW
        CMFUS = G8CMF + U8CMF0
        WFCM = G6CMW + G8CMF + U8CMF0
        WFLGCM = G6CMW + G8CMF + U8CMF0 + CMLG

        curve_rows.append({
            "QCL": QCL,
            "ANRW2WL": ANRW2WL,
            "CLW": G4CLW,
            "CDW": G5CDW,
            "CDFT": G9CDF,
            "CDFTLG": G9CDF,
            "CDAPT": CDAPT,
            "CMW": G6CMW,
            "CMFUS": CMFUS,
            "WFCM": WFCM,
            "WFLGCM": WFLGCM,
        })

        XDATA.append(G4CLW)
        YDATA.append(CDAPT)

        last_G4CLW = G4CLW
        last_G5CDW = G5CDW
        last_G6CMW = G6CMW
        last_U8CMF0 = U8CMF0
        last_G8CMF = G8CMF
        last_G9CDF = G9CDF
        last_ANRW2WL = ANRW2WL

    CD_COEF = _polyfit_basic_style(XDATA, YDATA, NP=4)

    M8 = state.MM / DEG_PER_RAD
    M9 = M8 * state.AWO

    # Original BASIC uses current G6CMW after loop. It is constant for this model.
    G6CMW_const = last_G6CMW if last_G6CMW is not None else 0.0
    U8CMF0_const = THETA * SLOPCM * (state.MM / DEG_PER_RAD)

    M10 = G6CMW_const + U8CMF0_const + CMLG + SLOPCM * M9
    M11 = SLOPCM * M8

    return {
        "S": S_ft2,
        "KF": KF,
        "SLOPCM": SLOPCM,
        "AD": AD,
        "CFD": CFD,
        "AT": AT,
        "CTD": CTD,
        "CDLG": CDLG,
        "CMLG": CMLG,
        "G8CMF": last_G8CMF,
        "U8CMF0": last_U8CMF0,
        "G9CDF": last_G9CDF,
        "Y(I)": YDATA,
        "X(I)": XDATA,
        "X(1)": CD_COEF[0],
        "X(2)": CD_COEF[1],
        "X(3)": CD_COEF[2],
        "X(4)": CD_COEF[3],
        "X(5)": CD_COEF[4],
        "M8": M8,
        "M9": M9,
        "M10": M10,
        "M11": M11,
        "curve_rows": curve_rows,
        "CL": {
            "C0": M9,
            "C1": M8,
            "C2": 0.0,
            "C3": 0.0,
            "C4": 0.0,
        },
        "CD": {
            "D0": CD_COEF[0],
            "D1": CD_COEF[1],
            "D2": CD_COEF[2],
            "D3": CD_COEF[3],
            "D4": CD_COEF[4],
        },
        "CM": {
            "M0": M10,
            "M1": M11,
            "M2": 0.0,
            "M3": 0.0,
            "M4": 0.0,
        },
        "CDNG": landing["CDNG"],
        "CDMG": landing["CDMG"],
        "CMNG": landing["CMNG"],
        "CMMG": landing["CMMG"],
        "last_CLW": last_G4CLW,
        "last_CDW": last_G5CDW,
        "last_ANRW2WL": last_ANRW2WL,
    }


def calculate_airloads(state: Chapter7State, inp: Dict[str, Any]) -> Dict[str, Any]:
    if state.KCL is None or state.CD is None or state.ZCOEFM is None or state.ANRW2WL is None:
        raise Chapter7InputError("Wing aero distribution must be calculated before airloads")

    V = _as_float(_get(inp, "V"), "V")
    AIRL = str(_get(inp, "AIRL$", "", required=False))
    WL = _as_float(_get(inp, "WL"), "WL")
    SL = _as_float(_get(inp, "SL"), "SL")

    Q = V ** 2 / 295.0
    AN = state.ANRW2WL

    L_load = []
    D_load = []
    ML = []
    LZ = []
    DX = []

    for j in range(state.H):
        L_J = state.KCL[j] * state.DA[j] * Q / 144.0
        D_J = state.CD[j] * state.DA[j] * Q / 144.0
        ML_J = state.ZCOEFM[j] * state.CE[j] * state.DA[j] * Q / 144.0
        LZ_J = L_J * math.cos(AN / DEG_PER_RAD) + D_J * math.sin(AN / DEG_PER_RAD)
        DX_J = D_J * math.cos(AN / DEG_PER_RAD) - L_J * math.sin(AN / DEG_PER_RAD)

        L_load.append(L_J)
        D_load.append(D_J)
        ML.append(ML_J)
        LZ.append(LZ_J)
        DX.append(DX_J)

    SZ = [0.0 for _ in range(state.H)]
    MXX = [0.0 for _ in range(state.H)]
    TYY = [0.0 for _ in range(state.H)]

    SZ[state.H - 1] = LZ[state.H - 1]
    SZ_acc = SZ[state.H - 1]
    MXX_acc = 0.0
    TYY_acc = 0.0

    for i in range(state.H - 2, -1, -1):
        SZ_acc += LZ[i]
        SZ[i] = SZ_acc
        MXX_acc += SZ[i + 1] * state.DY
        MXX[i] = MXX_acc
        TYY_acc -= SZ[i + 1] * (state.CX25[i + 1] - state.CX25[i])
        TYY[i] = TYY_acc

    Z = [WL + math.tan(SL / DEG_PER_RAD) * y for y in state.YE]

    SX = [0.0 for _ in range(state.H)]
    MZZ = [0.0 for _ in range(state.H)]
    TVYY = [0.0 for _ in range(state.H)]

    SX[state.H - 1] = DX[state.H - 1]
    SX_acc = SX[state.H - 1]
    MZZ_acc = 0.0
    TVYY_acc = 0.0

    for i in range(state.H - 2, -1, -1):
        SX_acc += DX[i]
        SX[i] = SX_acc
        MZZ_acc += SX[i + 1] * (state.YE[i + 1] - state.YE[i])
        MZZ[i] = MZZ_acc
        TVYY_acc += SX[i + 1] * (Z[i + 1] - Z[i])
        TVYY[i] = TVYY_acc

    TRQ = [0.0 for _ in range(state.H)]
    TRQ_acc = 0.0
    for i in range(state.H - 1, -1, -1):
        TRQ_acc += ML[i]
        TRQ[i] = TRQ_acc

    TMYY = [TYY[i] + TVYY[i] + TRQ[i] for i in range(state.H)]

    return {
        "AIRL$": AIRL,
        "Q": Q,
        "AN": AN,
        "L(I)": L_load,
        "D(I)": D_load,
        "ML(I)": ML,
        "LZ(I)": LZ,
        "DX(I)": DX,
        "Z(I)": Z,
        "SX(I)": SX,
        "SZ(I)": SZ,
        "MXX(I)": MXX,
        "TYY(I)": TYY,
        "TVYY(I)": TVYY,
        "TRQ(I)": TRQ,
        "TMYY(I)": TMYY,
        "MZZ(I)": MZZ,
    }


def state_wing_geometry_output(state: Chapter7State) -> Dict[str, Any]:
    return {
        "B": state.B,
        "DY": state.DY,
        "YE(I)": state.YE,
        "XF(I)": state.XF,
        "XA(I)": state.XA,
        "CE(I)": state.CE,
        "CE(I)*DY": state.DA,
        "CX25(I)": state.CX25,
        "C50X(I)": state.C50X,
        "A": state.A,
        "AREA_TOTAL": state.AREA_TOTAL,
        "SC2": state.SC2,
        "SAYE": state.SAYE,
        "SBARXC": state.SBARXC,
        "MAC": state.MAC,
        "YBAR": state.YBAR,
        "XMACLE": state.XMACLE,
        "ZAR": state.ZAR,
        "H": state.H,
    }


def _analyze_single_configuration(input_json: Dict[str, Any]) -> Dict[str, Any]:
    wing_geometry_inp = _get_group(input_json, "Wing geometry data", "wing_geometry")
    method_inp = _get_group(input_json, "AIRLOAD Method Data", "airload_method", required=False)
    additive_inp = _get_group(input_json, "Additive Lift Distribution Data", "additive_lift_distribution")
    basic_inp = _get_group(input_json, "Basic Lift Distribution Data", "basic_lift_distribution")
    stall_inp = _get_group(input_json, "Stall CL Data", "stall_cl", required=False)
    wing_aero_inp = _get_group(input_json, "Wing Aero Coefficient Distribution Data", "wing_aero_coefficient_distribution")
    tau_direct_inp = _get_group(input_json, "TAU Input Data", "tau_input", required=False)
    tau_inp = _get_group(input_json, "TAU Calculation Data", "tau_calculation", required=False)
    airplane_less_tail_inp = _get_group(input_json, "Airplane Less Tail Data", "airplane_less_tail", required=False)
    landing_gear_inp = _get_group(input_json, "Landing Gear Option Data", "landing_gear_option", required=False)
    airloads_inp = _get_group(input_json, "Airloads Option Data", "airloads_option", required=False)

    state = load_wing_geometry(wing_geometry_inp)

    intermediate: Dict[str, Any] = {
        "Wing Geometry Calculation Data": state_wing_geometry_output(state),
    }
    output: Dict[str, Any] = {
        "Wing Geometry Calculation Data": {
            "YE(I)": state.YE,
            "CE(I)": state.CE,
            "CE(I)*DY": state.DA,
            "AREA_TOTAL": state.AREA_TOTAL,
            "MAC": state.MAC,
            "YBAR": state.YBAR,
            "XMACLE": state.XMACLE,
            "ZAR": state.ZAR,
            "B": state.B,
            "H": state.H,
            "DY": state.DY,
        }
    }

    method_result = select_airloads_method(method_inp)
    intermediate["AIRLOAD Method Selection Data"] = method_result
    output["AIRLOAD Method Data"] = {
        "METHOD": method_result["METHOD"],
        "USE_AIRLOAD4": method_result["USE_AIRLOAD4"],
        "MACH": method_result["MACH"],
        "SBA": method_result["SBA"],
        "GC": method_result["GC"],
    }
    additive_inp = _deep_merge(
        additive_inp,
        {
            "GC": method_result["GC"],
            "USE_AIRLOAD4": method_result["USE_AIRLOAD4"],
        },
    )

    intermediate["Additive Lift Distribution Calculation Data"] = calculate_additive_lift_distribution(state, additive_inp)
    output["Additive Lift Distribution Data"] = {
        "CHECKCL1": state.CHECKCL1,
        "CCLA1(I)": state.CCLA1,
        "CLA1(I)": state.CLA1,
    }

    intermediate["Basic Lift Distribution Calculation Data"] = calculate_basic_lift_distribution(state, basic_inp)
    output["Basic Lift Distribution Data"] = {
        "AWO": state.AWO,
        "REFANG(I)": state.REFANG,
        "AO(I)": state.AO,
        "VCCLB(I)": state.VCCLB,
        "FCLB(I)": state.FCLB,
        "UCLB(I)": state.UCLB,
        "BCHECK": state.BCHECK,
    }

    if stall_inp:
        intermediate["Stall CL Calculation Data"] = calculate_stall_cl(state, stall_inp)
        output["Stall CL Data"] = {
            "R3N(I)": state.R3N,
            "C3LMAX(I)": state.C3LMAX,
            "CM(I)": state.CM,
            "STALLCL": state.WING_STALL_CL,
            "KCL(I)": state.KCL,
        }

    if method_result["USE_AIRLOAD4"]:
        if not tau_inp:
            raise Chapter7InputError(
                "tau_calculation with TAPR and TIPR is required for AIRLOAD4.BAS"
            )
        tau_result = calculate_tau(tau_inp)
        intermediate["TAU Calculation Data"] = tau_result
        wing_aero_inp = _deep_merge(wing_aero_inp, {"TAU": tau_result["TAU"]})
        output["AIRLOAD Method Data"].update({
            "TAPR": tau_result["TAPR"],
            "TIPR": tau_result["TIPR"],
            "TAU": tau_result["TAU"],
        })
    else:
        tau_value = _as_float(_get(tau_direct_inp, "TAU"), "TAU")
        intermediate["TAU Input Data"] = {"TAU": tau_value}
        wing_aero_inp = _deep_merge(wing_aero_inp, {"TAU": tau_value})
        output["AIRLOAD Method Data"]["TAU"] = tau_value

    wing_aero_inp = _deep_merge(
        wing_aero_inp,
        {
            "GC": method_result["GC"],
            "SBA": method_result["SBA"],
            "USE_AIRLOAD4": method_result["USE_AIRLOAD4"],
        },
    )
    intermediate["Wing Aero Coefficient Distribution Calculation Data"] = calculate_wing_aero_distribution(state, wing_aero_inp)
    output["Wing Aero Coefficient Distribution Data"] = {
        "KCL(I)": state.KCL,
        "CID(I)": state.CID,
        "CPD(I)": state.CPD,
        "CD(I)": state.CD,
        "ZCOEFM(I)": state.ZCOEFM,
        "ALPHA": state.ALPHA,
        "ANRW2WL": state.ANRW2WL,
        "CLW": state.G4CLW,
        "CDW": state.G5CDW,
        "CMW": state.G6CMW,
    }

    if airplane_less_tail_inp:
        airplane_less_tail = calculate_airplane_less_tail(state, airplane_less_tail_inp, landing_gear_inp)
        intermediate["Airplane Less Tail Calculation Data"] = airplane_less_tail
        lgext = _get(
            airplane_less_tail_inp,
            "LGEXT",
            _get(airplane_less_tail_inp, "LGEXT$", "N", required=False),
            required=False,
        )
        if landing_gear_inp and _left1(lgext) == "Y":
            intermediate["Landing Gear Calculation Data"] = {
                "CDNG": airplane_less_tail["CDNG"],
                "CDMG": airplane_less_tail["CDMG"],
                "CDLG": airplane_less_tail["CDLG"],
                "CMNG": airplane_less_tail["CMNG"],
                "CMMG": airplane_less_tail["CMMG"],
                "CMLG": airplane_less_tail["CMLG"],
            }

        curve_last = airplane_less_tail["curve_rows"][-1] if airplane_less_tail["curve_rows"] else {}
        output["Airplane Less Tail Data"] = {
            "CDFT": curve_last.get("CDFT"),
            "CDFTLG": curve_last.get("CDFTLG"),
            "CDAPT": curve_last.get("CDAPT"),
            "CMFUS": curve_last.get("CMFUS"),
            "WFCM": curve_last.get("WFCM"),
            "WFLGCM": curve_last.get("WFLGCM"),
            "STALLCL": state.WING_STALL_CL,
            "CL": airplane_less_tail["CL"],
            "CD": airplane_less_tail["CD"],
            "CM": airplane_less_tail["CM"],
            "CDNG": airplane_less_tail["CDNG"],
            "CDMG": airplane_less_tail["CDMG"],
            "CDLG": airplane_less_tail["CDLG"],
            "CMNG": airplane_less_tail["CMNG"],
            "CMMG": airplane_less_tail["CMMG"],
            "CMLG": airplane_less_tail["CMLG"],
            "curve_rows": airplane_less_tail["curve_rows"],
        }

    if airloads_inp:
        airloads = calculate_airloads(state, airloads_inp)
        intermediate["Airloads Calculation Data"] = airloads
        output["Airloads Option Data"] = {
            "CX25(I)": state.CX25,
            "YE(I)": state.YE,
            "Z(I)": airloads["Z(I)"],
            "DX(I)": airloads["DX(I)"],
            "LZ(I)": airloads["LZ(I)"],
            "ML(I)": airloads["ML(I)"],
            "SX(I)": airloads["SX(I)"],
            "SZ(I)": airloads["SZ(I)"],
            "MXX(I)": airloads["MXX(I)"],
            "TMYY(I)": airloads["TMYY(I)"],
            "MZZ(I)": airloads["MZZ(I)"],
        }

    return {
        "chapter": 7,
        "intermediate": intermediate,
        "output": output,
    }


def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    merged = deepcopy(base)
    for key, value in override.items():
        if key == "extends":
            continue
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def _get_ch5_wing_geometry(ch5_output: Dict[str, Any]) -> Dict[str, Any]:
    surfaces = ch5_output.get("surfaces")
    if not isinstance(surfaces, list):
        raise Chapter7InputError(
            "ch5_geometry_output.json must contain a 'surfaces' array"
        )

    for surface in surfaces:
        if not isinstance(surface, dict):
            continue
        surface_name = str(surface.get("N$", "")).upper()
        if "WING SURFACE" in surface_name:
            return deepcopy(surface)

    raise Chapter7InputError(
        "WING SURFACE was not found in ch5_geometry_output.json"
    )


def _merge_ch5_geometry(
    user_input: Dict[str, Any],
    ch5_output: Dict[str, Any],
) -> Dict[str, Any]:
    merged = deepcopy(user_input)
    common = merged.setdefault("common", {})
    if not isinstance(common, dict):
        raise Chapter7InputError("'common' must be a JSON object")
    common["wing_geometry"] = _get_ch5_wing_geometry(ch5_output)
    return merged


def _resolve_configuration(
    name: str,
    configurations: Dict[str, Any],
    common: Dict[str, Any],
    resolving: Optional[List[str]] = None,
) -> Dict[str, Any]:
    if name not in configurations:
        raise Chapter7InputError(f"Unknown configuration: {name}")

    resolving = resolving or []
    if name in resolving:
        chain = " -> ".join(resolving + [name])
        raise Chapter7InputError(f"Circular configuration inheritance: {chain}")

    raw = configurations[name]
    if not isinstance(raw, dict):
        raise Chapter7InputError(f"Configuration '{name}' must be a JSON object")

    parent_name = raw.get("extends")
    if parent_name is None:
        resolved = deepcopy(common)
    else:
        resolved = _resolve_configuration(
            str(parent_name),
            configurations,
            common,
            resolving + [name],
        )

    return _deep_merge(resolved, raw)


def _select_ch8_cases(
    ch8_output: Dict[str, Any],
    requested_cases: Any,
) -> List[Dict[str, Any]]:
    results = ch8_output.get("results")
    if not isinstance(results, list):
        raise Chapter7InputError("ch8_output.json must contain a 'results' array")

    cases_by_number: Dict[int, Dict[str, Any]] = {}
    for row in results:
        if not isinstance(row, dict):
            continue
        case_number = _as_int(_get(row, "CASE"), "CASE")
        if case_number in cases_by_number:
            raise Chapter7InputError(
                f"Duplicate CASE {case_number} in ch8_output.json"
            )
        cases_by_number[case_number] = row

    if isinstance(requested_cases, str):
        if requested_cases.lower() != "all":
            raise Chapter7InputError(
                "airloads_settings.cases must be 'all' or a list of CASE numbers"
            )
        return [cases_by_number[number] for number in sorted(cases_by_number)]

    if not isinstance(requested_cases, list) or not requested_cases:
        raise Chapter7InputError(
            "airloads_settings.cases must be a non-empty list or 'all'"
        )

    selected = []
    seen = set()
    for value in requested_cases:
        case_number = _as_int(value, "airloads_settings.cases")
        if case_number in seen:
            raise Chapter7InputError(
                f"Duplicate CASE {case_number} in airloads_settings.cases"
            )
        if case_number not in cases_by_number:
            raise Chapter7InputError(
                f"CASE {case_number} was not found in ch8_output.json"
            )
        seen.add(case_number)
        selected.append(cases_by_number[case_number])
    return selected


def _calculate_ch8_airloads_cases(
    input_json: Dict[str, Any],
    configurations: Dict[str, Dict[str, Any]],
    common: Dict[str, Any],
    ch8_output: Optional[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    settings = input_json.get("airloads_settings")
    if settings is None:
        return []
    if not isinstance(settings, dict):
        raise Chapter7InputError("'airloads_settings' must be a JSON object")
    if ch8_output is None:
        raise Chapter7InputError(
            "ch8_output.json is required when airloads_settings is present"
        )

    WL = _as_float(_get(settings, "WL"), "airloads_settings.WL")
    SL = _as_float(_get(settings, "SL"), "airloads_settings.SL")
    selected_cases = _select_ch8_cases(
        ch8_output,
        _get(settings, "cases"),
    )

    output_cases = []
    for case_row in selected_cases:
        case_number = _as_int(_get(case_row, "CASE"), "CASE")
        condition = str(_get(case_row, "CONDITION"))
        configuration_name = str(_get(case_row, "CONFIG")).upper()
        QCL = _as_float(_get(case_row, "CL"), "CL")
        V = _as_float(_get(case_row, "V"), "V")

        resolved_input = _resolve_configuration(
            configuration_name,
            configurations,
            common,
        )
        wing_aero_input = _get_group(
            resolved_input,
            "Wing Aero Coefficient Distribution Data",
            "wing_aero_coefficient_distribution",
        )
        resolved_input["wing_aero_coefficient_distribution"] = _deep_merge(
            wing_aero_input,
            {"QCL": QCL},
        )

        airload_description = f"CASE {case_number} {condition}"
        resolved_input["airloads_option"] = {
            "V": V,
            "AIRL$": airload_description,
            "WL": WL,
            "SL": SL,
        }

        calculation = _analyze_single_configuration(resolved_input)
        airloads = calculation["intermediate"].get("Airloads Calculation Data")
        if not isinstance(airloads, dict):
            raise Chapter7InputError(
                f"Airloads calculation did not produce data for CASE {case_number}"
            )

        output_cases.append({
            "CASE": case_number,
            "CONDITION": condition,
            "CONFIG": configuration_name,
            "CG": case_row.get("CG"),
            "ALT": case_row.get("ALT"),
            "WT": case_row.get("WT"),
            "XCG": case_row.get("XCG"),
            "ZCG": case_row.get("ZCG"),
            "QCL": QCL,
            "V": V,
            "AIRL$": airload_description,
            "WL": WL,
            "SL": SL,
            "airloads": airloads,
        })

    return output_cases


def analyze_chapter7_from_json(
    input_json: Dict[str, Any],
    ch8_output: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    configurations_input = input_json.get("configurations")
    if not isinstance(configurations_input, list) or not configurations_input:
        raise Chapter7InputError("'configurations' must be a non-empty JSON array")

    common = input_json.get("common", {})
    if not isinstance(common, dict):
        raise Chapter7InputError("'common' must be a JSON object")

    configurations: Dict[str, Dict[str, Any]] = {}
    configuration_order: List[str] = []
    for item in configurations_input:
        if not isinstance(item, dict):
            raise Chapter7InputError("Each configuration must be a JSON object")
        name = str(_get(item, "CNF")).upper()
        if name in configurations:
            raise Chapter7InputError(f"Duplicate configuration: {name}")
        configurations[name] = item
        configuration_order.append(name)

    requested = input_json.get("run_configurations", configuration_order)
    if isinstance(requested, str):
        requested_names = configuration_order if requested.lower() == "all" else [requested.upper()]
    elif isinstance(requested, list):
        requested_names = [str(name).upper() for name in requested]
    else:
        raise Chapter7InputError("'run_configurations' must be 'all', a name, or a list of names")

    if not requested_names:
        raise Chapter7InputError("'run_configurations' cannot be empty")

    output_configurations = []
    geometry_output = None
    for name in requested_names:
        resolved_input = _resolve_configuration(name, configurations, common)
        result = _analyze_single_configuration(resolved_input)
        calculation_output = result["output"]
        aero = calculation_output.get("Airplane Less Tail Data", {})
        method = calculation_output.get("AIRLOAD Method Data", {})
        wing_geometry = calculation_output.get("Wing Geometry Calculation Data", {})
        if geometry_output is None:
            geometry_output = {
                "MAC": wing_geometry.get("MAC"),
                "XW": (
                    wing_geometry.get("XMACLE") + 0.25 * wing_geometry.get("MAC")
                    if wing_geometry.get("XMACLE") is not None and wing_geometry.get("MAC") is not None
                    else None
                ),
                "S": (
                    wing_geometry.get("AREA_TOTAL") / 144.0
                    if wing_geometry.get("AREA_TOTAL") is not None
                    else None
                ),
            }
        input_configuration = configurations[name]
        output_configurations.append({
            "CNF": name,
            "STALLCL": aero.get("STALLCL"),
            "NE": input_configuration.get("NE"),
            "CL": aero.get("CL"),
            "CD": aero.get("CD"),
            "CM": aero.get("CM"),
            "weight_cases": deepcopy(input_configuration.get("weight_cases", [])),
        })

    geometry_input = input_json.get("geometry", {})
    geometry = {
        "MAC": (geometry_output or {}).get("MAC"),
        "XTC": geometry_input.get("XTC"),
        "XTF": geometry_input.get("XTF"),
        "XW": (geometry_output or {}).get("XW"),
        "ZW": geometry_input.get("ZW"),
        "S": (geometry_output or {}).get("S"),
    }
    analysis_input = input_json.get("analysis", {})
    analysis = {
        "MN": _as_float(_get_group(common, "airload_method").get("MACH", 0.0), "MACH"),
        "altitudes": deepcopy(analysis_input.get("altitudes", [])),
    }
    airloads_cases = _calculate_ch8_airloads_cases(
        input_json,
        configurations,
        common,
        ch8_output,
    )

    return {
        "geometry": geometry,
        "category": deepcopy(input_json.get("category", {})),
        "design_speeds": deepcopy(input_json.get("design_speeds", {})),
        "limit_loads": deepcopy(input_json.get("limit_loads", {})),
        "analysis": analysis,
        "configurations": output_configurations,
        "airloads_cases": airloads_cases,
    }


def _format_ch7_report(input_json: Dict[str, Any]) -> str:
    common = input_json["common"]
    configurations = {
        str(item["CNF"]).upper(): item for item in input_json["configurations"]
    }
    requested = [
        str(name).upper()
        for name in input_json.get("run_configurations", configurations.keys())
    ]
    lines: List[str] = []

    for name in requested:
        flap_deg = configurations[name].get("FLAP_DEG")
        result = _analyze_single_configuration(
            _deep_merge(common, configurations[name])
        )
        data = result["intermediate"]
        geometry = data["Wing Geometry Calculation Data"]
        additive = data["Additive Lift Distribution Calculation Data"]
        basic = data["Basic Lift Distribution Calculation Data"]
        stall = data.get("Stall CL Calculation Data", {})
        tau = data.get("TAU Calculation Data", {})
        wing = data["Wing Aero Coefficient Distribution Calculation Data"]
        airplane = data.get("Airplane Less Tail Calculation Data", {})

        lines.extend([
            (
                f"WING AERODYNAMIC COEFFICIENTS, {name} "
                f"(FLAP {flap_deg:g} DEG)"
                if flap_deg is not None
                else f"WING AERODYNAMIC COEFFICIENTS, {name}"
            ),
            "=" * 72,
            "",
            "WING GEOMETRY CALCULATIONS",
            f"AREA TOTAL = {geometry['AREA_TOTAL']:.4f}",
            f"MAC = {geometry['MAC']:.6f}  YBAR = {geometry['YBAR']:.6f}  "
            f"XMACLE = {geometry['XMACLE']:.6f}",
            f"ASPECT RATIO = {geometry['ZAR']:.6f}  SPAN = {geometry['B']:.3f}  "
            f"ELEMENTS = {geometry['H']}  DY = {geometry['DY']:.5f}",
            "",
            "ELEM          YE          CE",
        ])
        for j, (ye, ce) in enumerate(zip(geometry["YE(I)"], geometry["CE(I)"]), 1):
            lines.append(f"{j:4d} {ye:11.5f} {ce:11.5f}")

        lines.extend([
            "",
            "ADDITIVE LIFT DISTRIBUTION",
            f"CHECKCL1 = {additive['CHECKCL1']:.8f}",
            "ELEM          YE      CC(LA1)       C(LA1)",
        ])
        for j, (ye, cc, cl) in enumerate(
            zip(geometry["YE(I)"], additive["CCLA1(I)"], additive["CLA1(I)"]), 1
        ):
            lines.append(f"{j:4d} {ye:11.5f} {cc:12.5f} {cl:12.7f}")

        lines.extend([
            "",
            "BASIC LIFT DISTRIBUTION",
            f"AWO = {basic['AWO']:.7f}",
            "ELEM     REF ANGLE           AO        CC1B         C1B",
        ])
        for j, values in enumerate(zip(
            basic["REFANG(I)"], basic["AO(I)"], basic["VCCLB(I)"], basic["UCLB(I)"]
        ), 1):
            refang, ao, cc1b, c1b = values
            lines.append(
                f"{j:4d} {refang:13.5f} {ao:12.5f} {cc1b:12.5f} {c1b:11.5f}"
            )
        lines.append(f"BCHECK = {basic['BCHECK']:.8g}")

        if stall:
            lines.extend([
                "",
                "STALL CALCULATIONS",
                f"WING STALL CL = {stall['STALLCL']:.7f}",
                "ELEM          YE   ELEMENT CLMAX    CL AT STALL",
            ])
            for j, values in enumerate(zip(
                geometry["YE(I)"], stall["CM(I)"], stall["KCL(I)"]
            ), 1):
                ye, clmax, kcl = values
                lines.append(f"{j:4d} {ye:11.5f} {clmax:15.6f} {kcl:14.6f}")

        lines.extend([
            "",
            "WING AERO COEFFICIENT DISTRIBUTIONS",
            f"METHOD = {wing['METHOD']}",
            f"TAU = {wing['TAU']:.8f}"
            + (
                f"  TAPR = {tau['TAPR']:.8f}  TIPR = {tau['TIPR']:.8f}"
                if tau else ""
            ),
            f"GC = {wing['GC']:.8f}  SBA = {wing['SBA']:.5f} DEG",
            f"SPANWISE DISTRIBUTIONS FOR CL = {wing['G4CLW']:.6f}",
            "ELEM          CL         CDI         CPD          CD          CM",
        ])
        for j, values in enumerate(zip(
            wing["KCL(I)"], wing["CID(I)"], wing["CPD(I)"],
            wing["CD(I)"], wing["ZCOEFM(I)"]
        ), 1):
            cl, cdi, cpd, cd, cm = values
            lines.append(
                f"{j:4d} {cl:11.6f} {cdi:11.6f} {cpd:11.6f} "
                f"{cd:11.6f} {cm:11.6f}"
            )
        lines.extend([
            f"ALPHA = {wing['ALPHA']:.7f}",
            f"ANRW2WL = {wing['ANRW2WL']:.7f}",
            f"CL(WING) = {wing['G4CLW']:.7f}",
            f"CD(WING) = {wing['G5CDW']:.7f}",
            f"CM(WING) = {wing['G6CMW']:.7f}",
        ])

        if airplane:
            cl = airplane["CL"]
            cd = airplane["CD"]
            cm = airplane["CM"]
            lines.extend([
                "",
                "AIRPLANE LESS TAIL AERO COEFFICIENTS",
                "ANGLE(RW TO WL)  CL(WING)  CD(WING)  CD(F+T)  CD(A/P)  "
                "CM(WING)  CM(FUS)  CM(W+F)",
            ])
            for row in airplane["curve_rows"]:
                lines.append(
                    f"{row['ANRW2WL']:15.5f} {row['CLW']:9.5f} "
                    f"{row['CDW']:9.5f} {row['CDFT']:8.5f} "
                    f"{row['CDAPT']:8.5f} {row['CMW']:9.5f} "
                    f"{row['CMFUS']:8.5f} {row['WFCM']:9.5f}"
                )
            lines.extend([
                "",
                "EQUATIONS FOR AERO COEFFICIENTS FOR AIRPLANE LESS TAIL",
                f"CL={cl['C0']:+.6f}{cl['C1']:+.6f}*ANGLE(RW TO WL)",
                f"CD={cd['D0']:+.6f}{cd['D1']:+.6f}*CL"
                f"{cd['D2']:+.6f}*CL^2{cd['D3']:+.6f}*CL^3"
                f"{cd['D4']:+.6f}*CL^4",
                f"CM={cm['M0']:+.6f}{cm['M1']:+.6f}*ANGLE(RW TO WL)",
            ])

        lines.extend(["", "\f", ""])

    return "\n".join(lines)


def run_chapter7(
    input_path: str = INPUT_PATH,
    ch5_geometry_path: str = CH5_GEOMETRY_PATH,
    ch8_output_path: str = CH8_OUTPUT_PATH,
    output_path: str = OUTPUT_PATH,
    report_path: str = REPORT_PATH,
) -> Dict[str, Any]:
    with open(input_path, "r", encoding="utf-8-sig") as f:
        user_input = json.load(f)

    with open(ch5_geometry_path, "r", encoding="utf-8-sig") as f:
        ch5_output = json.load(f)

    with open(ch8_output_path, "r", encoding="utf-8-sig") as f:
        ch8_output = json.load(f)

    input_data = _merge_ch5_geometry(user_input, ch5_output)

    result = analyze_chapter7_from_json(input_data, ch8_output)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    with open(report_path, "w", encoding="utf-8") as f:
        f.write(_format_ch7_report(input_data))

    return result


if __name__ == "__main__":
    run_chapter7()
