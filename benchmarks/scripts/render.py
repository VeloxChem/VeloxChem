"""A run's records as the table BENCHMARKS.md carries.

    python render.py ../data/scf/2026-09-15_m4max_caffeine.json [--functional HF]

Writes the table beside its data file, as <name>.md, so the numbers and the table
made from them travel together and a table is never hand-typed.

Speedups are computed here and not stored: they are read against the `full` row of
the same molecule, basis and functional in the same file, so a re-measured column
can never leave a stale ratio behind it.
"""
import argparse
import json
from pathlib import Path

NAMES = {
    "full": "full four-centre",
    "ri_jk_conventional": "RI-JK veloxchem",
    "ri_jk_simd": "RI-JK simd",
    "ri_jk_simd_direct": "RI-JK simd",
}

MODES = {"in_memory": "in memory", "direct": "direct"}


def method_name(row):
    """What the row ran, mode included.

    The mode is read back from the driver rather than from what was asked for:
    `automatic` chooses, and a table which said only "simd" could not tell a run
    that held the B vectors from one that swept the integrals every build.
    """
    name = NAMES.get(row["method"], row["method"])
    mode = MODES.get(row.get("ri_mode"))
    return f"{name}, {mode}" if mode else name


def provenance_line(run):
    machine = run["machine"]
    config = run["config"]
    dirty = ", working tree dirty" if run.get("dirty") else ""
    return (f"Measured at `{run['commit']}`{dirty} on {machine['name']} "
            f"({machine['cpu']}, {machine['cores']} cores), "
            f"{config['mpi_ranks']} rank of {config['omp_threads']} threads, "
            f"veloxchem {run['veloxchem']}, {run['date'][:10]}.")


HEADER = ["basis", "nao", "fitting set", "naux", "method", "wall", "B vectors",
          "2e build", "XC", "rest", "iters", "energy", "build x", "whole x"]


def blocks(doc, functional=None):
    """The rows of each functional's table, formatted, as one place.

    Both the markdown and the pdf are rendered from this, so they cannot come to
    disagree about what a cell says.
    """
    rows = [r for r in doc["rows"] if r]

    for want in sorted({r["functional"] for r in rows},
                       key=lambda f: (f != "HF", f)):
        if functional and want != functional:
            continue
        here = [r for r in rows if r["functional"] == want]
        cells = []

        for basis in dict.fromkeys(r["basis"] for r in here):
            block = [r for r in here if r["basis"] == basis]
            full = next((r for r in block if r["method"] == "full"), None)
            # NOTE: the fitting set belongs to the block, not to its first row --
            # that row is the four-centre build, which has no auxiliary basis.
            fitted = next((r for r in block if r["aux_basis"]), None)

            for at, row in enumerate(block):
                build_x = whole_x = "--"
                if full and row["method"] != "full" and row["fock_2e_mean"]:
                    build_x = f"{full['fock_2e_mean'] / row['fock_2e_mean']:.2f}"
                    whole_x = f"{full['wall'] / row['wall']:.2f}"

                head = ([basis, f"{row['nao']}",
                         fitted["aux_basis"] if fitted else "--",
                         f"{fitted['naux']}" if fitted else "--"]
                        if at == 0 else ["", "", "", ""])

                cells.append(head + [
                    method_name(row),
                    f"{row['wall']:.2f}",
                    f"{row.get('ri_setup', 0.0):.2f}" if row.get("aux_basis")
                    else "--",
                    f"{row['fock_2e_total']:.2f}",
                    f"{row['fock_xc_total']:.2f}", f"{row['remainder']:.2f}",
                    f"{row['iterations']}", f"{row['energy']:.8f}",
                    build_x, whole_x])

        yield want, cells


def render(path, functional=None):
    doc = json.loads(Path(path).read_text())

    out = [provenance_line(doc["run"]), ""]

    for want, cells in blocks(doc, functional):
        out += [f"#### {want}", "",
                "| " + " | ".join(HEADER) + " |",
                "| " + " | ".join(["---"] + ["---:"] * (len(HEADER) - 1)) + " |"]
        out += ["| " + " | ".join(row) + " |" for row in cells]
        out.append("")

    return "\n".join(out)


def render_pdf(path, out):
    """The same tables as a pdf, one page per functional."""
    import matplotlib
    matplotlib.use("pdf")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    doc = json.loads(Path(path).read_text())
    line = provenance_line(doc["run"]).replace("`", "")

    with PdfPages(out) as pdf:
        for want, cells in blocks(doc):
            height = 1.6 + 0.20 * len(cells)
            fig, ax = plt.subplots(figsize=(13.5, height))
            ax.axis("off")

            ax.text(0.0, 1.0, f"{doc['suite'].upper()}: {want}",
                    transform=ax.transAxes, fontsize=13, fontweight="bold")
            ax.text(0.0, 1.0 - 0.9 / height, line, transform=ax.transAxes,
                    fontsize=7.5, color="#475569")

            # a band per basis, so the three methods of one row read as a group
            shades, band = [], 0
            for row in cells:
                if row[0]:
                    band += 1
                shades.append(["#ffffff" if band % 2 else "#f1f5f9"] * len(HEADER))

            # NOTE: explicit widths. Left to itself the table divides the row
            # evenly and the two columns carrying words -- the fitting set and
            # the method -- are clipped mid-name.
            widths = [9, 4, 14, 5, 16, 6, 7, 7, 5, 6, 5, 12, 6, 6]
            widths = [w / sum(widths) for w in widths]

            table = ax.table(cellText=cells, colLabels=HEADER, cellColours=shades,
                             colWidths=widths, cellLoc="right", loc="upper center",
                             bbox=[0, 0, 1, 1.0 - 1.35 / height])
            table.auto_set_font_size(False)
            table.set_fontsize(6.6)

            for (r, c), cell in table.get_celld().items():
                cell.set_edgecolor("#cbd5e1")
                cell.set_linewidth(0.4)
                if r == 0:
                    cell.set_facecolor("#1e293b")
                    cell.set_text_props(color="white", fontweight="bold")
                elif c in (12, 13) and cells[r - 1][c] != "--":
                    cell.set_text_props(fontweight="bold")
                if c in (0, 2, 4):
                    cell.set_text_props(ha="left")
                    cell._text.set_x(0.03)

            pdf.savefig(fig, bbox_inches="tight", pad_inches=0.3)
            plt.close(fig)

    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    parser.add_argument("--functional", default=None)
    parser.add_argument("--out", default=None,
                        help="where to write it; the data file's name with .md "
                             "by default, so a table sits beside the run it came "
                             "from. '-' prints instead.")
    parser.add_argument("--no-pdf", action="store_true",
                        help="the markdown alone")
    args = parser.parse_args()

    text = render(args.path, args.functional)

    if args.out == "-":
        print(text)
        return

    out = Path(args.out) if args.out else Path(args.path).with_suffix(".md")
    out.write_text(text + "\n")
    print(f"written to {out}")

    if not args.no_pdf:
        print(f"written to {render_pdf(args.path, out.with_suffix('.pdf'))}")


if __name__ == "__main__":
    main()
