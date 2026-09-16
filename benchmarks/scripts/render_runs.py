"""The gradient and optimization records as a table, and as a pdf.

    python render_runs.py ../data/gradient/2026-09-16_m4max_caffeine.json

Writes <name>.md and <name>.pdf beside the data, so the numbers and the table made
from them travel together and a table is never hand-typed. The suite is read from
the file rather than named on the command line.

Speedups are computed here and not stored: they are read against the `full` row of
the same molecule, basis and functional in the same file, so a re-measured row can
never leave a stale ratio behind it.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from render import provenance_line as _provenance_line


def provenance_line(run):
    """The shared line, less a thread count that was never recorded.

    The first gradient and optimization suites passed None for the threads, so
    the shared line renders "1 rank of None threads". The runs are measured and
    will not be repeated to record it, so it is dropped from the line rather
    than printed as a word.
    """
    return _provenance_line(run).replace(" of None threads", "")

NAMES = {"full": "four-centre", "ri_jk_simd": "RI-JK simd"}
MODES = {"in_memory": "in memory", "direct": "direct"}

SPECS = {
    # suite -> (header, per-row cells, explicit column widths)
    "gradient": (
        ["functional", "basis", "nao", "naux", "method", "gradient (s)",
         "speedup", "SCF (s)", "vs four-centre"],
        [7, 9, 4, 5, 15, 8, 7, 7, 10],
    ),
    "optimization": (
        ["functional", "basis", "nao", "naux", "method", "total (s)", "speedup",
         "steps", "s/step", "energy (a.u.)"],
        [7, 9, 4, 5, 15, 8, 7, 5, 7, 13],
    ),
    "tda": (
        ["functional", "basis", "nao", "naux", "method", "TDA (s)", "speedup",
         "iter", "s/iter", "SCF (s)", "max dE (a.u.)"],
        [7, 9, 4, 5, 15, 8, 7, 5, 7, 7, 10],
    ),
}


def method_name(row):
    name = NAMES.get(row["method"], row["method"])
    mode = MODES.get(row.get("ri_mode"))
    return f"{name}, {mode}" if mode else name


def _label(doc):
    """What the run solved for, which the suite alone does not say.

    The excited state suite runs two solvers into records of one shape, so the
    column which carries the time is named from the rows and not from the suite.
    A record written before the solver was recorded is the Tamm-Dancoff one.
    """
    if doc["suite"] != "tda":
        return doc["suite"].upper()
    return doc["rows"][0].get("solver", "tda").upper()


def _headed(doc, header):
    """The header with the time column named for the solver."""
    return [h.replace("TDA (s)", f"{_label(doc)} (s)") for h in header]


def _cost(row, suite):
    if suite == "gradient":
        return row["grad_wall"]
    if suite == "tda":
        return row["tda_wall"]
    return row["wall"]


def cells_of(doc):
    """One list of cells per row, grouped so a basis reads as a block."""
    suite = doc["suite"]
    rows = doc["rows"]
    out = []

    for functional in dict.fromkeys(r["functional"] for r in rows):
        for basis in dict.fromkeys(r["basis"] for r in rows
                                   if r["functional"] == functional):
            picked = [r for r in rows if r["functional"] == functional
                      and r["basis"] == basis]
            ref = next((r for r in picked if r["method"] == "full"), None)

            # NOTE: from whichever row of the group carries it. The heading cells
            # are written on the first row, which is the four-centre one, and that
            # row resolves no identity and has no auxiliary basis.
            naux = next((r["naux"] for r in picked if r.get("naux")), None)

            for n, r in enumerate(picked):
                speed = (f'{_cost(ref, suite) / _cost(r, suite):.2f}'
                         if ref else "--")
                head = [functional, basis, str(r["nao"]),
                        str(naux or "--")] if n == 0 else ["", "", "", ""]

                if suite == "gradient":
                    if r["method"] == "full" or ref is None:
                        agree = "--"
                    else:
                        d = (np.array(r["_gradient"]) - np.array(ref["_gradient"]))
                        agree = f'{np.abs(d).max():.1e}'
                    out.append(head + [method_name(r), f'{r["grad_wall"]:.2f}',
                                       speed, f'{r["scf_wall"]:.2f}', agree])
                elif suite == "tda":
                    if r["method"] == "full" or ref is None:
                        agree = "--"
                    else:
                        d = (np.array(r["excitation_energies"]) -
                             np.array(ref["excitation_energies"]))
                        agree = f'{np.abs(d).max():.1e}'
                    # NOTE: the time per iteration beside the total. Two runs which
                    # converge in different numbers of iterations are not compared
                    # on speed by their totals alone.
                    out.append(head + [method_name(r), f'{r["tda_wall"]:.2f}',
                                       speed, str(r["iterations"]),
                                       f'{r["tda_wall"] / max(r["iterations"], 1):.2f}',
                                       f'{r["scf_wall"]:.2f}', agree])
                else:
                    out.append(head + [method_name(r), f'{r["wall"]:.1f}', speed,
                                       str(r["steps"]),
                                       f'{r["wall_per_step"]:.2f}',
                                       f'{r["energy"]:.8f}'])
    return out


def scaling(doc):
    """The exponent p of cost proportional to nao**p, per functional and method.

    Fitted over every basis of the file, which is the only thing in it that varies
    the size. The auxiliary basis does not vary with it -- one fitting set serves
    every orbital basis here -- so the exponent of the resolution of the identity
    is in the orbital dimension alone and is not the scaling of the method with
    the problem.
    """
    suite = doc["suite"]
    series = {}

    for functional in dict.fromkeys(r["functional"] for r in doc["rows"]):
        for method in dict.fromkeys(r["method"] for r in doc["rows"]):
            picked = sorted((r for r in doc["rows"]
                             if r["functional"] == functional
                             and r["method"] == method),
                            key=lambda r: r["nao"])
            if len(picked) < 2:
                continue
            nao = np.array([r["nao"] for r in picked], dtype=float)
            cost = np.array([_cost(r, suite) for r in picked], dtype=float)
            slope = np.polyfit(np.log(nao), np.log(cost), 1)[0]
            series[(functional, method)] = (nao, cost, slope)

    return series


def _select(doc, basis):
    """The run with only the rows of one basis, so a table can be made of it."""
    if basis is None:
        return doc
    rows = [r for r in doc["rows"] if r["basis"] == basis]
    if not rows:
        raise SystemExit(f"no rows of basis {basis} in this run")
    return {**doc, "rows": rows}


def render(path, basis=None):
    doc = _select(json.loads(Path(path).read_text()), basis)
    header = _headed(doc, SPECS[doc["suite"]][0])
    lines = [f'## {_label(doc).lower()}: {doc["rows"][0]["molecule"]}', "",
             provenance_line(doc["run"]), "",
             "| " + " | ".join(header) + " |",
             "| " + " | ".join("---" for _ in header) + " |"]
    lines += ["| " + " | ".join(c for c in row) + " |" for row in cells_of(doc)]

    # NOTE: an exponent needs two sizes to fit. A run of one basis is a
    # measurement and not a scaling, and says so by leaving the section out
    # rather than by printing an empty table.
    if doc["suite"] == "gradient" and scaling(doc):
        lines += ["", "### Scaling, cost proportional to nao to the p", "",
                  "| functional | method | p |", "| --- | --- | ---: |"]
        for (f, m), (_, _, p) in scaling(doc).items():
            lines.append(f'| {f} | {NAMES.get(m, m)} | {p:.2f} |')

    return "\n".join(lines)


def render_pdf(path, out, basis=None):
    import matplotlib
    matplotlib.use("pdf")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    doc = _select(json.loads(Path(path).read_text()), basis)
    header, widths = SPECS[doc["suite"]]
    header = _headed(doc, header)
    cells = cells_of(doc)
    line = provenance_line(doc["run"]).replace("`", "")
    title = f'{_label(doc)}: {doc["rows"][0]["molecule"]}'

    with PdfPages(out) as pdf:
        height = 1.6 + 0.22 * len(cells)
        fig, ax = plt.subplots(figsize=(12.5, height))
        ax.axis("off")
        ax.text(0.0, 1.0, title, transform=ax.transAxes, fontsize=13,
                fontweight="bold")
        ax.text(0.0, 1.0 - 0.9 / height, line, transform=ax.transAxes,
                fontsize=7.5, color="#475569")

        # a band per basis, so the methods of one basis read as a group
        shades, band = [], 0
        for row in cells:
            if row[0] or row[1]:
                band += 1
            shades.append(["#ffffff" if band % 2 else "#f1f5f9"] * len(header))

        # NOTE: explicit widths. Left to itself the table divides the row evenly
        # and the columns carrying words are clipped mid-name.
        norm = [w / sum(widths) for w in widths]
        table = ax.table(cellText=cells, colLabels=header, cellColours=shades,
                         colWidths=norm, cellLoc="right", loc="upper center",
                         bbox=[0, 0, 1, 1.0 - 1.35 / height])
        table.auto_set_font_size(False)
        table.set_fontsize(7.0)

        speed_col = header.index("speedup")
        for (r, c), cell in table.get_celld().items():
            cell.set_edgecolor("#cbd5e1")
            cell.set_linewidth(0.4)
            if r == 0:
                cell.set_facecolor("#1e293b")
                cell.set_text_props(color="white", fontweight="bold")
            elif c == speed_col and cells[r - 1][c] not in ("--", "1.00"):
                cell.set_text_props(fontweight="bold")
            if c in (0, 1, 4):
                cell.set_text_props(ha="left")
                cell._text.set_x(0.03)

        pdf.savefig(fig, bbox_inches="tight", pad_inches=0.3)
        plt.close(fig)

        if doc["suite"] == "gradient" and scaling(doc):
            _scaling_page(pdf, doc, title, line)

    return out


def _scaling_page(pdf, doc, title, line):
    import matplotlib.pyplot as plt

    series = scaling(doc)

    fig, (ax, tx) = plt.subplots(1, 2, figsize=(12.5, 5.2),
                                 gridspec_kw={"width_ratios": [3, 2]})
    fig.suptitle(f'{title} -- scaling with the orbital basis', fontsize=13,
                 fontweight="bold", x=0.09, ha="left")
    fig.text(0.09, 0.90, line, fontsize=7.5, color="#475569", ha="left")

    styles = {("HF", "full"): ("#b91c1c", "o-"),
              ("HF", "ri_jk_simd"): ("#1d4ed8", "o-"),
              ("B3LYP", "full"): ("#f97316", "s--"),
              ("B3LYP", "ri_jk_simd"): ("#0891b2", "s--")}

    for (f, m), (nao, cost, p) in series.items():
        colour, style = styles.get((f, m), ("#475569", "o-"))
        ax.loglog(nao, cost, style, color=colour, markersize=4.5, linewidth=1.4,
                  label=f'{f}, {NAMES.get(m, m)}  (p = {p:.2f})')

    ax.set_xlabel("basis functions")
    ax.set_ylabel("gradient wall time (s)")
    ax.grid(True, which="both", linewidth=0.3, color="#cbd5e1")
    ax.legend(fontsize=8, frameon=False)

    tx.axis("off")
    rows = [[f, NAMES.get(m, m), f'{p:.2f}']
            for (f, m), (_, _, p) in series.items()]
    t = tx.table(cellText=rows, colLabels=["functional", "method", "p"],
                 colWidths=[0.3, 0.42, 0.18], cellLoc="left",
                 loc="upper center", bbox=[0, 0.42, 1, 0.48])
    t.auto_set_font_size(False)
    t.set_fontsize(8.5)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("#cbd5e1")
        cell.set_linewidth(0.4)
        if r == 0:
            cell.set_facecolor("#1e293b")
            cell.set_text_props(color="white", fontweight="bold")

    tx.text(0.0, 0.34,
            "Cost fitted as proportional to nao to the p, over every basis of\n"
            "the run. One auxiliary basis serves all of them, so naux does not\n"
            "grow with nao: the exponent of the resolution of the identity is\n"
            "in the orbital dimension alone, and the ratio between the methods\n"
            "widens here partly because the fitting set is held still.",
            fontsize=8, color="#334155", va="top", linespacing=1.6)

    pdf.savefig(fig, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    parser.add_argument("--out", default=None)
    parser.add_argument("--basis", default=None,
                        help="render only the rows of this basis")
    parser.add_argument("--no-pdf", action="store_true")
    args = parser.parse_args()

    if args.out:
        out = Path(args.out)
    elif args.basis:
        stem = Path(args.path).with_suffix("")
        out = Path(f"{stem}_{args.basis}.md")
    else:
        out = Path(args.path).with_suffix(".md")

    out.write_text(render(args.path, args.basis) + "\n")
    print(f"written to {out}")

    if not args.no_pdf:
        print(f"written to {render_pdf(args.path, out.with_suffix('.pdf'), args.basis)}")


if __name__ == "__main__":
    main()
