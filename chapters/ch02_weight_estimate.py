from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = PROJECT_ROOT / "inputs" / "initial_aircraft.json"
OUTPUT_PATH = PROJECT_ROOT / "outputs" / "cases" / "ch2_weight_estimate_output.json"
REPORT_PATH = PROJECT_ROOT / "outputs" / "reports" / "ch2_weight_estimate_report.txt"


# ============================================================
# Chapter 2 - WTESTIMA.BAS
# Input : inputs/initial_aircraft.json
# Output: outputs/cases/ch2_weight_estimate_output.json
# ============================================================


ENGTYPE_DESCRIPTION = {
    "RF": "RECIPROCAL 4 CYCLE ENGINE",
    "RT": "RECIPROCAL 2 CYCLE ENGINE",
    "TC": "TURBOCHARGED RECIPROCAL ENGINE",
    "LC": "LIQUID COOLED RECIPROCAL ENGINE",
    "TP": "TURBOPROP ENGINE",
}


def read_json(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    with path.open("r", encoding="utf-8-sig") as f:
        return json.load(f)


def write_json(data: dict[str, Any], path: str | Path) -> None:
    path = Path(path)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def write_text(text: str, path: str | Path) -> None:
    path = Path(path)
    with path.open("w", encoding="utf-8") as f:
        f.write(text)


def round0(x: float) -> int:
    """
    BASIC 출력값처럼 정수 표시용.
    Python round()는 bankers rounding이라 여기서는 일반 반올림으로 처리.
    """
    return int(math.floor(x))


def validate_input(data: dict[str, Any]) -> None:
    required = ["NOENGS", "HP", "ENGTYPE$", "SEATS", "BAG", "P$", "HOURS"]

    for key in required:
        if key not in data:
            raise ValueError(f"Missing input key: {key}")

    if int(data["NOENGS"]) < 1 or int(data["NOENGS"]) > 6:
        raise ValueError("NOENGS must be 1~6")

    if float(data["HP"]) < 1 or float(data["HP"]) > 3000:
        raise ValueError("HP must be 1~3000")

    if int(data["SEATS"]) < 1 or int(data["SEATS"]) > 12:
        raise ValueError("SEATS must be 1~12")

    if float(data["HOURS"]) <= 0 or float(data["HOURS"]) > 10:
        raise ValueError("HOURS must be > 0 and <= 10")

    if float(data["BAG"]) < 0 or float(data["BAG"]) > 6250:
        raise ValueError("BAG must be 0~6250")

    engtype = str(data["ENGTYPE$"]).upper()
    if engtype not in ENGTYPE_DESCRIPTION:
        raise ValueError(
            f"ENGTYPE$ must be one of {list(ENGTYPE_DESCRIPTION.keys())}, got {engtype}"
        )

    p = str(data["P$"]).upper()
    if p not in ["P", "U", "Y", "N"]:
        raise ValueError("P$ must be P/U or Y/N")


def normalize_P(P_dollar: str) -> str:
    """Normalize legacy Y/N input to the original WTESTIMA.BAS P/U symbols."""
    p = P_dollar.upper()
    if p == "Y":
        return "P"
    if p == "N":
        return "U"
    return p


def compute_K(NOENGS: int, SEATS: int, ENGTYPE: str, P_dollar: str) -> float:
    """
    Wempty/Wto 비율 K 계산.

    기본:
        K = 0.62

    조건 보정:
        SEATS = 1       -> K - 0.04
        P$ = Y          -> K + 0.02
        NOENGS >= 2     -> K + 0.01
        ENGTYPE$ = LC   -> K + 0.01
        ENGTYPE$ = TC   -> K + 0.01
        ENGTYPE$ = TP   -> K - 0.05
        ENGTYPE$ = RT   -> K - 0.01
    """

    K = 0.62

    if SEATS == 1:
        K -= 0.04

    if P_dollar == "P":
        K += 0.02

    if NOENGS >= 2:
        K += 0.01

    if ENGTYPE == "LC":
        K += 0.01

    if ENGTYPE == "TC":
        K += 0.01

    if ENGTYPE == "TP":
        K -= 0.05

    if ENGTYPE == "RT":
        K -= 0.01

    return K


def compute_WTFUEL(HP: float, HOURS: float, ENGTYPE: str) -> float:
    """
    연료중량 계산.

    RF, TC, LC:
        WTFUEL = 0.75 * HP * 0.5 * HOURS

    RT:
        WTFUEL = 0.75 * HP * 0.7 * HOURS

    TP:
        WTFUEL = HP * 0.55 * HOURS
    """

    if ENGTYPE in ["RF", "TC", "LC"]:
        return 0.75 * HP * 0.5 * HOURS

    if ENGTYPE == "RT":
        return 0.75 * HP * 0.7 * HOURS

    if ENGTYPE == "TP":
        return HP * 0.55 * HOURS

    raise ValueError(f"Unknown ENGTYPE$: {ENGTYPE}")


def compute_WTENGINSTALLED(HP: float, NOENGS: int, ENGTYPE: str) -> float:
    """
    장착 엔진 중량 계산.

    문서 예제 RF, HP=265, NOENGS=1에서 490 lb가 나와야 하므로
    reciprocating 계열은 아래 근사식을 사용.

        WTENGINSTALLED = NOENGS * 2.575 * (HP / NOENGS)^0.950

    TP는 별도 계수가 필요할 수 있음.
    일단 동일 식으로 계산하되 ENGTYPE별로 분리 가능하게 함수로 둠.
    """

    hp_per_engine = HP / NOENGS

    if ENGTYPE == "RF":
        return NOENGS * (
            105.8439 + 1.448059 * hp_per_engine + 6.31254e-6 * hp_per_engine**2
        )

    if ENGTYPE == "RT":
        return NOENGS * (
            26.08666 + 0.650924 * hp_per_engine - 4.431183e-3 * hp_per_engine**2
        )

    if ENGTYPE == "TC":
        return NOENGS * (
            155.6418 + 1.4689 * hp_per_engine + 3.37101e-4 * hp_per_engine**2
        )

    if ENGTYPE == "LC":
        return NOENGS * (
            387.2534 + 1.02973 * hp_per_engine - 4.09947e-4 * hp_per_engine**2
        )

    if ENGTYPE == "TP":
        return 0.48 * HP * 1.3

    raise ValueError(f"Unknown ENGTYPE$: {ENGTYPE}")


def compute_ch2(aircraft_data: dict[str, Any]) -> dict[str, Any]:
    validate_input(aircraft_data)

    # ========================================================
    # INPUT
    # ========================================================
    NOENGS = int(aircraft_data["NOENGS"])
    HP = float(aircraft_data["HP"])
    ENGTYPE_dollar = str(aircraft_data["ENGTYPE$"]).upper()
    SEATS = int(aircraft_data["SEATS"])
    BAG = float(aircraft_data["BAG"])
    P_dollar = normalize_P(str(aircraft_data["P$"]))
    HOURS = float(aircraft_data["HOURS"])
    AIRPLANE_dollar = str(
        aircraft_data.get("AIRPLANE$", aircraft_data.get("AIRPLANE_NAME", "AIRPLANE"))
    )

    # ========================================================
    # INTERMEDIATE - Useful Load
    # ========================================================
    WTSEATS = SEATS * 170.0
    WTFUEL = compute_WTFUEL(HP=HP, HOURS=HOURS, ENGTYPE=ENGTYPE_dollar)
    USEFUL = WTFUEL + WTSEATS + BAG

    # ========================================================
    # INTERMEDIATE - Weight Estimation
    # ========================================================
    K = compute_K(
        NOENGS=NOENGS,
        SEATS=SEATS,
        ENGTYPE=ENGTYPE_dollar,
        P_dollar=P_dollar,
    )

    WTO = USEFUL / (1.0 - K)

    # ========================================================
    # INTERMEDIATE - Structure Weight
    # ========================================================
    WTFUS = 0.0982 * WTO
    WTWING = 0.1036 * WTO
    WTTAIL = 0.0234 * WTO
    WTNAC = 0.0146 * WTO
    WTLANDGEAR = 0.0571 * WTO
    WTCONTROLS = 0.0150 * WTO

    WTSTRUCT = (
        WTFUS
        + WTWING
        + WTTAIL
        + WTNAC
        + WTLANDGEAR
        + WTCONTROLS
    )

    # ========================================================
    # INTERMEDIATE - Powerplant Weight
    # ========================================================
    WTENGINSTALLED = compute_WTENGINSTALLED(
        HP=HP,
        NOENGS=NOENGS,
        ENGTYPE=ENGTYPE_dollar,
    )

    WTFUELSYS = 0.1068 * WTENGINSTALLED

    if NOENGS >= 2:
        WTEXHAUST = 0.251 * WTENGINSTALLED
    else:
        WTEXHAUST = 0.147 * WTENGINSTALLED

    WTENGOTHER = 0.1757 * WTENGINSTALLED

    WTPROP = NOENGS * 0.2515 * ((HP / NOENGS) ** 1.04)

    WTPPGROUP = (
        WTENGINSTALLED
        + WTFUELSYS
        + WTEXHAUST
        + WTENGOTHER
    )

    # ========================================================
    # INTERMEDIATE - System Weight
    # ========================================================
    if NOENGS >= 2:
        WTNAVEQUIP = 0.0118 * WTO
        WTPNEUMATIC = 0.0 * WTO
        WTELECTRICAL = 0.0269 * WTO
        WTELECTRONIC = 0.0024 * WTO
        WTFURNEQUIP = 0.0458 * WTO
        WTENVIRANTIICE = 0.0118 * WTO
        WTMISC = 0.0079 * WTO
        WTTOTALSYS = 0.119 * WTO

        misc_output_name = "WTMISC"
        misc_output_value = WTMISC

    else:
        WTNAVEQUIP = 0.0044 * WTO
        WTPNEUMATIC = 0.00099 * WTO
        WTELECTRICAL = 0.0241 * WTO
        WTELECTRONIC = 0.0 * WTO
        WTFURNEQUIP = 0.0441 * WTO
        WTENVIRANTIICE = 0.0031 * WTO
        WTMISC = 0.0022 * WTO
        WTTOTALSYS = 0.0774 * WTO

        misc_output_name = "WTMISC"
        misc_output_value = WTMISC

    # ========================================================
    # INTERMEDIATE - Component Weight Sum / Check
    # ========================================================
    SUMWTS = WTSTRUCT + WTPPGROUP + WTTOTALSYS

    # Options & miscellaneous
    OPTMISC = WTO - USEFUL - SUMWTS
    OPTMISC_ITERATIONS = 1

    while OPTMISC < 0:
        OPTMISC_ITERATIONS += 1
        if OPTMISC_ITERATIONS > 1000:
            raise RuntimeError("OPTMISC did not converge after 1000 iterations.")

        WTO *= 1.01

        WTFUS = 0.0982 * WTO
        WTWING = 0.1036 * WTO
        WTTAIL = 0.0234 * WTO
        WTNAC = 0.0146 * WTO
        WTLANDGEAR = 0.0571 * WTO
        WTCONTROLS = 0.0150 * WTO
        WTSTRUCT = (
            WTFUS
            + WTWING
            + WTTAIL
            + WTNAC
            + WTLANDGEAR
            + WTCONTROLS
        )

        WTENGINSTALLED = compute_WTENGINSTALLED(
            HP=HP,
            NOENGS=NOENGS,
            ENGTYPE=ENGTYPE_dollar,
        )
        WTFUELSYS = 0.1068 * WTENGINSTALLED
        if NOENGS >= 2:
            WTEXHAUST = 0.251 * WTENGINSTALLED
        else:
            WTEXHAUST = 0.147 * WTENGINSTALLED
        WTENGOTHER = 0.1757 * WTENGINSTALLED
        WTPROP = NOENGS * 0.2515 * ((HP / NOENGS) ** 1.04)
        WTPPGROUP = (
            WTENGINSTALLED
            + WTFUELSYS
            + WTEXHAUST
            + WTENGOTHER
        )

        if NOENGS >= 2:
            WTNAVEQUIP = 0.0118 * WTO
            WTPNEUMATIC = 0.0 * WTO
            WTELECTRICAL = 0.0269 * WTO
            WTELECTRONIC = 0.0024 * WTO
            WTFURNEQUIP = 0.0458 * WTO
            WTENVIRANTIICE = 0.0118 * WTO
            WTMISC = 0.0079 * WTO
            WTTOTALSYS = 0.119 * WTO
            misc_output_name = "WTMISC"
            misc_output_value = WTMISC
        else:
            WTNAVEQUIP = 0.0044 * WTO
            WTPNEUMATIC = 0.00099 * WTO
            WTELECTRICAL = 0.0241 * WTO
            WTELECTRONIC = 0.0 * WTO
            WTFURNEQUIP = 0.0441 * WTO
            WTENVIRANTIICE = 0.0031 * WTO
            WTMISC = 0.0022 * WTO
            WTTOTALSYS = 0.0774 * WTO
            misc_output_name = "WTMISC"
            misc_output_value = WTMISC

        SUMWTS = WTSTRUCT + WTPPGROUP + WTTOTALSYS
        OPTMISC = WTO - USEFUL - SUMWTS

    # EMPTY, ETR은 원 코드 변수라기보다는 출력식/구현용 값
    EMPTY = WTO - USEFUL
    ETR = EMPTY / WTO
    useful_load_items = {
        "OPTMISC": OPTMISC,
        "PILOT": 170.0,
        **{f"PASSENGER NO. {i}": 170.0 for i in range(2, SEATS + 1)},
        "BAGGAGE": BAG,
        "FUEL": WTFUEL,
    }
    useful_load_items_rounded = {
        "OPTMISC": round0(OPTMISC),
        "PILOT": 170,
        **{f"PASSENGER NO. {i}": 170 for i in range(2, SEATS + 1)},
        "BAGGAGE": round0(BAG),
        "FUEL": round0(WTFUEL),
    }

    # ========================================================
    # OUTPUT
    # ========================================================
    output = {
        "chapter": 2,
        "program": "WTESTIMA.BAS",
        "units": {
            "weight": "lb",
            "power": "hp",
            "time": "hr",
        },
        "INPUT": {
            "Aircraft Data": {
                "NOENGS": NOENGS,
                "HP": HP,
                "ENGTYPE$": ENGTYPE_dollar,
                "ENGTYPE_DESCRIPTION": ENGTYPE_DESCRIPTION[ENGTYPE_dollar],
                "SEATS": SEATS,
                "BAG": BAG,
                "P$": P_dollar,
                "P_DESCRIPTION": "PRESSURIZED" if P_dollar == "P" else "UNPRESSURIZED",
                "HOURS": HOURS,
                "AIRPLANE$": AIRPLANE_dollar,
            }
        },
        "INTERMEDIATE": {
            "Useful Load": {
                "WTSEATS": WTSEATS,
                "WTFUEL": WTFUEL,
                "USEFUL": USEFUL,
            },
            "Weight Estimation": {
                "K": K,
                "WTO": WTO,
                "EMPTY": EMPTY,
                "OPTMISC_ITERATIONS": OPTMISC_ITERATIONS,
            },
            "Structure Weight": {
                "WTFUS": WTFUS,
                "WTWING": WTWING,
                "WTTAIL": WTTAIL,
                "WTNAC": WTNAC,
                "WTLANDGEAR": WTLANDGEAR,
                "WTCONTROLS": WTCONTROLS,
                "WTSTRUCT": WTSTRUCT,
            },
            "Powerplant Weight": {
                "WTENGINSTALLED": WTENGINSTALLED,
                "WTPROP": WTPROP,
                "WTFUELSYS": WTFUELSYS,
                "WTEXHAUST": WTEXHAUST,
                "WTENGOTHER": WTENGOTHER,
                "WTPPGROUP": WTPPGROUP,
            },
            "System Weight": {
                "WTNAVEQUIP": WTNAVEQUIP,
                "WTPNEUMATIC": WTPNEUMATIC,
                "WTELECTRICAL": WTELECTRICAL,
                "WTELECTRONIC": WTELECTRONIC,
                "WTFURNEQUIP": WTFURNEQUIP,
                "WTENVIRANTIICE": WTENVIRANTIICE,
                misc_output_name: misc_output_value,
                "WTTOTALSYS": WTTOTALSYS,
            },
            "Component Weight Sum": {
                "WTSTRUCT": WTSTRUCT,
                "WTPPGROUP": WTPPGROUP,
                "WTTOTALSYS": WTTOTALSYS,
                "SUMWTS": SUMWTS,
                "WTENGINSTALLED": WTENGINSTALLED,
            },
            "Check": {
                "OPTMISC": OPTMISC,
                "OPTMISC_ITERATIONS": OPTMISC_ITERATIONS,
            },
        },
        "OUTPUT": {
            "Weight Summary": {
                "WTO": WTO,
                "USEFUL": USEFUL,
                "EMPTY": EMPTY,
                "ETR": ETR,
            },
            "Structure Weight": {
                "WTWING": WTWING,
                "WTFUS": WTFUS,
                "WTTAIL": WTTAIL,
                "WTNAC": WTNAC,
                "WTLANDGEAR": WTLANDGEAR,
                "WTCONTROLS": WTCONTROLS,
                "WTSTRUCT": WTSTRUCT,
            },
            "Powerplant Weight": {
                "WTENGINSTALLED": WTENGINSTALLED,
                "WTPROP": WTPROP,
                "WTFUELSYS": WTFUELSYS,
                "WTEXHAUST": WTEXHAUST,
                "WTENGOTHER": WTENGOTHER,
                "WTPPGROUP": WTPPGROUP,
            },
            "System Weight": {
                "WTNAVEQUIP": WTNAVEQUIP,
                "WTPNEUMATIC": WTPNEUMATIC,
                "WTELECTRICAL": WTELECTRICAL,
                "WTELECTRONIC": WTELECTRONIC,
                "WTFURNEQUIP": WTFURNEQUIP,
                "WTENVIRANTIICE": WTENVIRANTIICE,
                misc_output_name: misc_output_value,
                "WTTOTALSYS": WTTOTALSYS,
            },
            "Optional / Useful Load Items": {
                **useful_load_items,
            },
        },
        "OUTPUT_ROUNDED": {
            "Weight Summary": {
                "WTO": round0(WTO),
                "USEFUL": round0(USEFUL),
                "EMPTY": round0(EMPTY),
                "ETR": math.floor(100 * ETR) / 100,
            },
            "Structure Weight": {
                "WTWING": round0(WTWING),
                "WTFUS": round0(WTFUS),
                "WTTAIL": round0(WTTAIL),
                "WTNAC": round0(WTNAC),
                "WTLANDGEAR": round0(WTLANDGEAR),
                "WTCONTROLS": round0(WTCONTROLS),
                "WTSTRUCT": round0(WTSTRUCT),
            },
            "Powerplant Weight": {
                "WTENGINSTALLED": round0(WTENGINSTALLED),
                "WTPROP": round0(WTPROP),
                "WTFUELSYS": round0(WTFUELSYS),
                "WTEXHAUST": round0(WTEXHAUST),
                "WTENGOTHER": round0(WTENGOTHER),
                "WTPPGROUP": round0(WTPPGROUP),
            },
            "System Weight": {
                "WTNAVEQUIP": round0(WTNAVEQUIP),
                "WTPNEUMATIC": round0(WTPNEUMATIC),
                "WTELECTRICAL": round0(WTELECTRICAL),
                "WTELECTRONIC": round0(WTELECTRONIC),
                "WTFURNEQUIP": round0(WTFURNEQUIP),
                "WTENVIRANTIICE": round0(WTENVIRANTIICE),
                misc_output_name: round0(misc_output_value),
                "WTTOTALSYS": round0(WTTOTALSYS),
            },
            "Optional / Useful Load Items": {
                **useful_load_items_rounded,
            },
        },
    }

    return output


def report_line(label: str, value: Any | None = None, value_col: int = 28) -> str:
    if value is None:
        return label
    return f"{label:<{value_col}}{value}"


def build_ch2_report(result: dict[str, Any]) -> str:
    aircraft = result["INPUT"]["Aircraft Data"]
    rounded = result["OUTPUT_ROUNDED"]
    summary = rounded["Weight Summary"]
    structure = rounded["Structure Weight"]
    powerplant = rounded["Powerplant Weight"]
    systems = rounded["System Weight"]
    useful_items = rounded["Optional / Useful Load Items"]
    airplane_name = aircraft.get("AIRPLANE$", aircraft.get("AIRPLANE_NAME", "AIRPLANE"))

    lines: list[str] = [
        "Estimated Weight",
        "",
        f"ESTIMATED WEIGHT DATA FOR {airplane_name}",
        "",
        "INPUT",
        "",
        report_line("MAX CONTINUOUS HP", round0(aircraft["HP"])),
        report_line("NUMBER OF ENGINES", aircraft["NOENGS"]),
        report_line("NUMBER OF SEATS", aircraft["SEATS"]),
        report_line("HOURS AT CRUISE POWER", round0(aircraft["HOURS"])),
        report_line("MAX BAGGAGE WEIGHT", round0(aircraft["BAG"])),
        aircraft["P_DESCRIPTION"],
        aircraft["ENGTYPE_DESCRIPTION"],
        "",
        "OUTPUT",
        "",
        report_line("MAX TAKE OFF WT", summary["WTO"]),
        report_line("USEFUL LOAD", summary["USEFUL"]),
        report_line("EMPTY WEIGHT", summary["EMPTY"]),
        report_line("W(EMPTY)/W(TO)", summary["ETR"]),
        "",
        report_line("WING", structure["WTWING"]),
        report_line("FUSELAGE", structure["WTFUS"]),
        report_line("TAIL", structure["WTTAIL"]),
        report_line("NACELLE", structure["WTNAC"]),
        report_line("LANDING GEAR", structure["WTLANDGEAR"]),
        report_line("CONTROLS", structure["WTCONTROLS"]),
        report_line("TOTAL STRUCTURE", structure["WTSTRUCT"], value_col=42),
        "",
        report_line("ENGINE INSTALLED", powerplant["WTENGINSTALLED"]),
        report_line("(INCLUDES PROPELLER(S))", f"( {powerplant['WTPROP']} )"),
        report_line("FUEL SYSTEM", powerplant["WTFUELSYS"]),
        report_line("EXHAUST", powerplant["WTEXHAUST"]),
        report_line("OTHER ENGINE DETAILS", powerplant["WTENGOTHER"]),
        report_line("TOTAL POWERPLANT", powerplant["WTPPGROUP"], value_col=42),
        "",
        report_line("INSTRUMENTS & NAV EQUIP", systems["WTNAVEQUIP"]),
        report_line("PNEUMATICS", systems["WTPNEUMATIC"]),
        report_line("ELECTRICAL", systems["WTELECTRICAL"]),
        report_line("ELECTRONICS", systems["WTELECTRONIC"]),
        report_line("FURNISHINGS & EQUIPMENT", systems["WTFURNEQUIP"]),
        report_line("ENVIRONMENTAL & ANTI-ICE", systems["WTENVIRANTIICE"]),
    ]

    lines.extend(
        [
            report_line("MISC OTHER SYSTEM WT", systems["WTMISC"]),
            report_line("TOTAL SYSTEMS WEIGHT", systems["WTTOTALSYS"], value_col=42),
            "",
            report_line("OPTIONS & MISCELLANEOUS", useful_items["OPTMISC"]),
            "",
            report_line("EMPTY WEIGHT", summary["EMPTY"], value_col=52),
            "",
            report_line("PILOT", useful_items["PILOT"]),
        ]
    )

    for passenger_no in range(2, int(aircraft["SEATS"]) + 1):
        label = f"PASSENGER NO. {passenger_no}"
        lines.append(report_line(label, useful_items[label]))

    lines.extend(
        [
            report_line("BAGGAGE", useful_items["BAGGAGE"]),
            report_line("FUEL", useful_items["FUEL"]),
            "",
            "",
            report_line("USEFUL LOAD", summary["USEFUL"], value_col=52),
        ]
    )

    return "\n".join(lines) + "\n"


def run_ch2(
    input_json_path: str | Path = INPUT_PATH,
    output_json_path: str | Path = OUTPUT_PATH,
    report_txt_path: str | Path = REPORT_PATH,
) -> None:
    db = read_json(input_json_path)

    if "Aircraft Data" in db:
        aircraft_data = db["Aircraft Data"]
    else:
        aircraft_data = db

    result = compute_ch2(aircraft_data)
    Path(output_json_path).parent.mkdir(parents=True, exist_ok=True)
    Path(report_txt_path).parent.mkdir(parents=True, exist_ok=True)
    write_json(result, output_json_path)
    write_text(build_ch2_report(result), report_txt_path)

    print(f"Saved output JSON: {output_json_path}")
    print(f"Saved report TXT: {report_txt_path}")


if __name__ == "__main__":
    run_ch2()
