from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = PROJECT_ROOT / "inputs" / "wing_db.json"
OUTPUT_PATH = PROJECT_ROOT / "outputs" / "cases" / "ch5_geometry_output.json"
REPORT_PATH = PROJECT_ROOT / "outputs" / "reports" / "ch5_geometry_report.txt"


# ============================================================
# Chapter 5 - WINGGEOM.BAS
# Input : inputs/wing_db.json
# Output: outputs/cases/ch5_geometry_output.json
# ============================================================


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


def validate_surface(surface: dict[str, Any]) -> None:
    """Validate Chapter 5 surface input lengths."""

    name = surface["N$"]

    M = int(surface["M"])
    N = int(surface["N"])
    H = int(surface["H"])

    XLE = surface["XLE(I)"]
    YLE = surface["YLE(I)"]
    XTE = surface["XTE(I)"]
    YTE = surface["YTE(I)"]

    if len(XLE) != M:
        raise ValueError(f"{name}: M={M}, but len(XLE(I))={len(XLE)}")

    if len(YLE) != M:
        raise ValueError(f"{name}: M={M}, but len(YLE(I))={len(YLE)}")

    if len(XTE) != N:
        raise ValueError(f"{name}: N={N}, but len(XTE(I))={len(XTE)}")

    if len(YTE) != N:
        raise ValueError(f"{name}: N={N}, but len(YTE(I))={len(YTE)}")

    if H <= 0:
        raise ValueError(f"{name}: H must be positive.")


def interp_x_from_y(
    target_y: float,
    X_points: list[float],
    Y_points: list[float],
    surface_name: str,
    edge_name: str,
) -> float:
    """
    Linear interpolation along one WINGGEOM edge.

    Edge coordinates are given as:
        XLE(I), YLE(I)
        XTE(I), YTE(I)

    Returns the X coordinate at the requested YE(I) station.
    """

    for k in range(len(Y_points) - 1):
        y1 = float(Y_points[k])
        y2 = float(Y_points[k + 1])
        x1 = float(X_points[k])
        x2 = float(X_points[k + 1])

        # Use the segment containing target_y.
        if (y1 <= target_y <= y2) or (y2 <= target_y <= y1):
            if math.isclose(y2, y1):
                return x1

            ratio = (target_y - y1) / (y2 - y1)
            return x1 + ratio * (x2 - x1)

    raise ValueError(
        f"{surface_name}: YE(I)={target_y} is outside {edge_name} range. "
        f"Y range = {min(Y_points)} ~ {max(Y_points)}"
    )


def compute_surface(surface: dict[str, Any]) -> dict[str, Any]:
    """
    Chapter 5 WINGGEOM calculation.

    Variable naming follows the original WINGGEOM.BAS listing:
        INPUT:
            N$, Q2$, M, N, H, Q4$, C3$
            XLE(I), YLE(I), XTE(I), YTE(I)

        INTERMEDIATE:
            DY, YE(I), XF(I), XA(I), CE(I), DA, A(I), A
            C2DY, SC2, ZDAYE, SAYE, DBARXC, SBARXC
            XBAR, YBAR, MAC, XMACLE, ZAR

        OUTPUT:
            Coordinate Print:
                XLE(I), YLE(I), XTE(I), YTE(I)

            Surface Output:
                A, MAC, YBAR, XMACLE, ZAR

            Element Output:
                XF(I), XA(I), YE(I), CE(I), CE(I)*DY
    """

    validate_surface(surface)

    # -----------------------------
    # INPUT
    # -----------------------------
    N_dollar = surface["N$"]
    Q2_dollar = str(surface["Q2$"]).upper()
    M = int(surface["M"])
    N = int(surface["N"])
    H = int(surface["H"])
    Q4_dollar = str(surface["Q4$"]).upper()
    C3_dollar = str(surface["C3$"]).upper()

    XLE = [float(v) for v in surface["XLE(I)"]]
    YLE = [float(v) for v in surface["YLE(I)"]]
    XTE = [float(v) for v in surface["XTE(I)"]]
    YTE = [float(v) for v in surface["YTE(I)"]]

    # -----------------------------
    # Initialize intermediate sums.
    # -----------------------------
    DY = (YLE[-1] - YLE[0]) / H

    A = 0.0
    SC2 = 0.0
    SAYE = 0.0
    SBARXC = 0.0

    element_rows: list[dict[str, float]] = []

    # -----------------------------
    # Element calculation.
    # -----------------------------
    for idx in range(H):
        # YE(I): element midpoint span station.
        YE_i = YLE[0] + DY * (idx + 0.5)

        # XF(I), XA(I): leading/trailing edge X at YE(I).
        XF_i = interp_x_from_y(
            target_y=YE_i,
            X_points=XLE,
            Y_points=YLE,
            surface_name=N_dollar,
            edge_name="leading edge",
        )

        XA_i = interp_x_from_y(
            target_y=YE_i,
            X_points=XTE,
            Y_points=YTE,
            surface_name=N_dollar,
            edge_name="trailing edge",
        )

        # CE(I), DA, A(I)
        CE_i = XA_i - XF_i
        DA = CE_i * DY
        A_i = DA

        # Surface Sum
        C2DY = CE_i**2 * DY
        ZDAYE = DA * YE_i
        DBARXC = DA * (XF_i + XA_i) / 2.0

        A += DA
        SC2 += C2DY
        SAYE += ZDAYE
        SBARXC += DBARXC

        element_rows.append(
            {
                "XF(I)": XF_i,
                "XA(I)": XA_i,
                "YE(I)": YE_i,
                "CE(I)": CE_i,
                "DA": DA,
                "A(I)": A_i,
                "C2DY": C2DY,
                "ZDAYE": ZDAYE,
                "DBARXC": DBARXC,
                "CE(I)*DY": CE_i * DY,
            }
        )

    if math.isclose(A, 0.0):
        raise ValueError(f"{N_dollar}: A is zero. Check geometry input.")

    # -----------------------------
    # Final Geometry
    # -----------------------------
    XBAR = SBARXC / A
    YBAR = SAYE / A
    MAC = SC2 / A
    XMACLE = XBAR - MAC / 2.0

    # -----------------------------
    # Aspect Ratio
    # -----------------------------
    # Q2$ = Y: surface is symmetric about Y=0.
    #   total span = 2 * max Y
    #   total area = 2 * A
    #   ZAR = (2Y)^2 / (2A)
    #
    # Q2$ = N: non-symmetric or one-sided surface.
    #   span = max(Y) - min(Y)
    #   area = A
    #   ZAR = span^2 / A
    # -----------------------------
    if Q2_dollar == "Y":
        span_total = 2.0 * max(max(YLE), max(YTE))
        area_total = 2.0 * A
        ZAR = span_total**2 / area_total
    else:
        span_total = max(max(YLE), max(YTE)) - min(min(YLE), min(YTE))
        area_total = A
        ZAR = span_total**2 / area_total

    # -----------------------------
    # OUTPUT JSON
    # -----------------------------
    result = {
        "N$": N_dollar,
        "Q2$": Q2_dollar,
        "M": M,
        "N": N,
        "H": H,
        "Q4$": Q4_dollar,
        "C3$": C3_dollar,
        "Coordinate Print": {
            "XLE(I)": XLE,
            "YLE(I)": YLE,
            "XTE(I)": XTE,
            "YTE(I)": YTE,
        },
        "Intermediate": {
            "Element Station": {
                "DY": DY,
            },
            "Surface Sum": {
                "A": A,
                "SC2": SC2,
                "SAYE": SAYE,
                "SBARXC": SBARXC,
            },
            "Final Geometry": {
                "XBAR": XBAR,
                "YBAR": YBAR,
                "MAC": MAC,
                "XMACLE": XMACLE,
                "ZAR": ZAR,
            },
        },
        "Surface Output": {
            "A": A,
            "MAC": MAC,
            "YBAR": YBAR,
            "XMACLE": XMACLE,
            "ZAR": ZAR,
        },
        "Element Output": [
            {
                "XF(I)": row["XF(I)"],
                "XA(I)": row["XA(I)"],
                "YE(I)": row["YE(I)"],
                "CE(I)": row["CE(I)"],
                "CE(I)*DY": row["CE(I)*DY"],
            }
            for row in element_rows
        ],
    }

    return result


def format_number(value: float, decimals: int = 3, trim: bool = False) -> str:
    rounded = round(float(value), decimals)
    if decimals == 0:
        return str(int(round(rounded)))
    text = f"{rounded:.{decimals}f}"
    if trim:
        text = text.rstrip("0").rstrip(".")
    return text


def surface_report_title(surface_name: str) -> str:
    name = surface_name.strip()
    upper_name = name.upper()

    if upper_name.startswith("WING "):
        base = "WING"
    elif " FOR " in upper_name:
        base = name[: upper_name.index(" FOR ")]
    else:
        base = name

    words = [word.capitalize() for word in base.replace("$", "").split()]
    title = " ".join(words)
    if "Geometry" not in title:
        title = f"{title} Geometry"
    return title


def format_coordinate_rows(x_values: list[float], y_values: list[float]) -> list[str]:
    rows = []
    for idx, (x_value, y_value) in enumerate(zip(x_values, y_values), start=1):
        rows.append(
            f"{idx:>4}"
            f" {format_number(x_value, 5, trim=True):>14}"
            f" {format_number(y_value, 5, trim=True):>14}"
        )
    return rows


def format_surface_report(surface: dict[str, Any]) -> str:
    coordinates = surface["Coordinate Print"]
    surface_output = surface["Surface Output"]

    lines = [
        surface_report_title(surface["N$"]),
        " " * 28 + "AERODYNAMIC SURFACE GEOMETRIC PROPERTIES",
        "",
        surface["N$"],
        "",
        "SYM ABOUT CL" if surface["Q2$"] == "Y" else "NOT SYM ABOUT CL",
        "",
        "COORDINATES OF LEADING EDGE",
        "POINT NO        FUS STA (XLE)   WING STA (YLE)",
        *format_coordinate_rows(coordinates["XLE(I)"], coordinates["YLE(I)"]),
        "",
        "COORDINATES OF TRAILING EDGE",
        "POINT NO        FUS STA (XTE)   WING STA (YTE)",
        *format_coordinate_rows(coordinates["XTE(I)"], coordinates["YTE(I)"]),
        "",
        "AREA/SIDE   MAC        YLE(MAC)    XLE(MAC)    ASPECT RATIO",
        (
            f"{format_number(surface_output['A'], 0):>8}"
            f"{format_number(surface_output['MAC']):>10}"
            f"{format_number(surface_output['YBAR']):>12}"
            f"{format_number(surface_output['XMACLE']):>12}"
            f"{format_number(surface_output['ZAR']):>14}"
        ),
    ]

    if surface.get("Q4$") == "Y":
        lines.extend(
            [
                "",
                "ELEMENT DATA",
                "ELEM        XLE        XTE          Y          C       AREA",
            ]
        )
        for idx, row in enumerate(surface["Element Output"], start=1):
            lines.append(
                f"{idx:>4}"
                f"{format_number(row['XF(I)']):>11}"
                f"{format_number(row['XA(I)']):>11}"
                f"{format_number(row['YE(I)']):>11}"
                f"{format_number(row['CE(I)']):>11}"
                f"{format_number(row['CE(I)*DY']):>11}"
            )

    return "\n".join(lines)


def build_geometry_report(output: dict[str, Any]) -> str:
    sections = [format_surface_report(surface) for surface in output["surfaces"]]
    return "\n\n\n".join(sections) + "\n"


def run_ch5(
    input_json_path: str | Path = INPUT_PATH,
    output_json_path: str | Path = OUTPUT_PATH,
    report_txt_path: str | Path = REPORT_PATH,
) -> None:
    db = read_json(input_json_path)

    surfaces = db.get("surfaces", [])
    if not surfaces:
        raise ValueError("No surfaces found in input DB.")

    output = {
        "chapter": db.get("chapter", 5),
        "program": db.get("program", "WINGGEOM.BAS"),
        "units": {
            "A": "in^2",
            "MAC": "in",
            "YBAR": "in",
            "XMACLE": "in",
            "ZAR": "-",
            "DY": "in",
            "XF(I)": "in",
            "XA(I)": "in",
            "YE(I)": "in",
            "CE(I)": "in",
            "CE(I)*DY": "in^2",
        },
        "surfaces": [],
    }

    for surface in surfaces:
        output["surfaces"].append(compute_surface(surface))

    Path(output_json_path).parent.mkdir(parents=True, exist_ok=True)
    Path(report_txt_path).parent.mkdir(parents=True, exist_ok=True)
    write_json(output, output_json_path)
    write_text(build_geometry_report(output), report_txt_path)

    print(f"Saved output JSON: {output_json_path}")
    print(f"Saved report TXT: {report_txt_path}")


if __name__ == "__main__":
    run_ch5()
