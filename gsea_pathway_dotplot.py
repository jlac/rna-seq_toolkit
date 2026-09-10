#!/usr/bin/env python3
"""
gsea_pathway_dotplot.py

Build a publication-ready dot plot comparing a user-selected set of pathways
across an arbitrary number of GSEA result tables.

  * dot COLOR  = NES        (red = enriched/up, blue = depleted/down)
  * dot SIZE   = adjusted p (bigger = more significant, log10 scaled)

Accepts clusterProfiler / fgsea / GSEA-desktop (GSEApy) style tables in
CSV, TSV, TXT, XLSX or XLS format, and auto-detects the relevant columns.

--------------------------------------------------------------------------
QUICK START
--------------------------------------------------------------------------
  # labelled inputs (label=path), pathways from a text file
  python gsea_pathway_dotplot.py \
      -i "ADP-heptose=gsea_adp.csv" \
      -i "ROSAH vs healthy donor=gsea_rosah.csv" \
      -p pathways.txt \
      -o custom_dotplot_ADP-heptose_ROSAH.pdf

  # inline pathway list, Excel workbook sheets as separate comparisons
  python gsea_pathway_dotplot.py \
      -i "Day 3=results.xlsx::D3_GSEA" -i "Day 7=results.xlsx::D7_GSEA" \
      -p "HALLMARK_INFLAMMATORY_RESPONSE,HALLMARK_TNFA_SIGNALING_VIA_NFKB" \
      -o dotplot.pdf

  # no pathway list: take the union of the top 15 hits per table
  python gsea_pathway_dotplot.py -i a.tsv -i b.tsv --top-n 15 -o dotplot.pdf

--------------------------------------------------------------------------
PATHWAY LIST FILE
--------------------------------------------------------------------------
One pathway per line. Blank lines and lines starting with '#' are ignored.
An optional tab (or ' | ') separated second field overrides the display label:

    HALLMARK_TNFA_SIGNALING_VIA_NFKB<TAB>TNFa signaling via NF-kB
    HALLMARK_INTERFERON_ALPHA_RESPONSE

Matching is fuzzy-normalised (case, underscores, spaces, punctuation and
common collection prefixes such as HALLMARK_ / GOBP_ / REACTOME_ / KEGG_ are
ignored), and is attempted against both the ID and the Description column.
--------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D

# --------------------------------------------------------------------------
# column synonyms
# --------------------------------------------------------------------------
ID_COLS = ["ID", "pathway", "Pathway", "PATHWAY", "NAME", "Name", "Term",
           "term", "TERM", "ont", "GeneSet", "gene_set", "geneset",
           "Description", "description", "DESCRIPTION"]

DESC_COLS = ["Description", "description", "DESCRIPTION", "NAME", "Name",
             "Term", "term", "pathway", "Pathway", "ID"]

NES_COLS = ["NES", "nes", "normalizedEnrichmentScore",
            "normalized_enrichment_score", "Normalized Enrichment Score"]

ES_COLS = ["ES", "es", "enrichmentScore", "enrichment_score",
           "Enrichment Score"]

PADJ_COLS = ["p.adjust", "padj", "p_adjust", "padjust", "p.adj", "FDR q-val",
             "fdr", "FDR", "FDR.q.val", "qvalue", "qvalues", "q.value",
             "qval", "adj.P.Val", "adjusted_pvalue", "Adjusted P-value",
             "p.adjusted", "FDR p-value"]

PVAL_COLS = ["pvalue", "pval", "p.value", "P.Value", "p_value", "PValue",
             "NOM p-val", "NOM.p.val", "pvalues"]

SETSIZE_COLS = ["setSize", "size", "SIZE", "set_size", "Size", "Gene Set Size"]

COLLECTION_PREFIXES = (
    "HALLMARK_", "GOBP_", "GOCC_", "GOMF_", "GO_", "REACTOME_", "KEGG_",
    "KEGG_MEDICUS_", "WP_", "BIOCARTA_", "PID_", "HP_", "MSIGDB_",
    "CP_", "C2_", "C5_", "C7_",
)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def eprint(*a):
    print(*a, file=sys.stderr)


def norm_key(s: str) -> str:
    """Normalise a pathway name for fuzzy matching."""
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return ""
    s = str(s).strip().upper()
    for pref in COLLECTION_PREFIXES:
        if s.startswith(pref):
            s = s[len(pref):]
            break
    s = re.sub(r"[^A-Z0-9]+", "", s)
    return s


def pretty_label(s: str, style: str = "upper") -> str:
    """Human-readable pathway label."""
    s = str(s).strip()
    stripped = s.upper()
    for pref in COLLECTION_PREFIXES:
        if stripped.startswith(pref):
            s = s[len(pref):]
            break
    s = s.replace("_", " ")
    s = re.sub(r"\s+", " ", s).strip()
    if style == "upper":
        return s.upper()
    if style == "title":
        return s.title()
    if style == "sentence":
        return s[:1].upper() + s[1:].lower() if s else s
    return s  # raw


def pick_col(df: pd.DataFrame, candidates, override=None, required=False,
             what="column", src=""):
    if override:
        if override in df.columns:
            return override
        raise SystemExit(
            f"[error] column '{override}' not found in {src}. "
            f"Available: {list(df.columns)}")
    lower = {str(c).strip().lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        if cand.strip().lower() in lower:
            return lower[cand.strip().lower()]
    if required:
        raise SystemExit(
            f"[error] could not find a {what} column in {src}. "
            f"Available: {list(df.columns)}. Use the matching --*-col option.")
    return None


def read_table(spec: str) -> pd.DataFrame:
    """Read csv/tsv/txt/xlsx. Excel sheet given as 'file.xlsx::SheetName'."""
    sheet = None
    path = spec
    if "::" in spec:
        path, sheet = spec.split("::", 1)
    if not os.path.exists(path):
        raise SystemExit(f"[error] input file not found: {path}")

    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls", ".xlsm"):
        df = pd.read_excel(path, sheet_name=sheet if sheet is not None else 0)
    else:
        # sniff the delimiter (handles , \t ; and whitespace)
        df = pd.read_csv(path, sep=None, engine="python", comment=None)
    df.columns = [str(c).strip() for c in df.columns]
    # GSEA-desktop tables sometimes carry a stray unnamed index column
    df = df.loc[:, ~df.columns.str.match(r"^Unnamed: \d+$")]
    return df


def parse_inputs(items):
    """Turn 'Label=path' / 'path' strings into (label, spec) pairs."""
    out = []
    for it in items:
        if "=" in it and not os.path.exists(it):
            label, spec = it.split("=", 1)
            label = label.strip()
        else:
            spec = it
            base = os.path.basename(spec.split("::")[0])
            label = os.path.splitext(base)[0]
            if "::" in it:
                label = it.split("::", 1)[1]
        out.append((label, spec.strip()))
    return out


def load_pathways(arg):
    """Return list of (query, display_override_or_None)."""
    if arg is None:
        return None
    entries = []
    if os.path.exists(arg):
        with open(arg) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                if "\t" in line:
                    q, lab = line.split("\t", 1)
                elif " | " in line:
                    q, lab = line.split(" | ", 1)
                else:
                    q, lab = line, None
                entries.append((q.strip(), lab.strip() if lab else None))
    else:
        for q in arg.split(","):
            if q.strip():
                entries.append((q.strip(), None))
    return entries


def fmt_p(p: float) -> str:
    """Format a p-value for the size legend (0.05 -> '0.05', 1e-5 -> '10^-5')."""
    if p >= 1e-3:
        s = f"{p:g}"
        return s
    e = int(round(np.log10(p)))
    return rf"$10^{{{e}}}$"


# --------------------------------------------------------------------------
# core
# --------------------------------------------------------------------------
def build_matrix(inputs, pathway_entries, args):
    """Return long-format dataframe of plotted values + ordered pathway list."""
    tables = {}
    for label, spec in inputs:
        df = read_table(spec)
        src = spec

        id_col = pick_col(df, ID_COLS, args.id_col, True, "pathway ID", src)
        desc_col = pick_col(df, DESC_COLS, args.desc_col, False,
                            "description", src) or id_col
        nes_col = pick_col(df, NES_COLS, args.nes_col, False, "NES", src)
        if nes_col is None:
            nes_col = pick_col(df, ES_COLS, None, True, "NES/ES", src)
            eprint(f"[warn] {label}: no NES column; using '{nes_col}'.")
        padj_col = pick_col(df, PADJ_COLS, args.padj_col, False,
                            "adjusted p-value", src)
        if padj_col is None:
            padj_col = pick_col(df, PVAL_COLS, None, True,
                                "adjusted p-value / p-value", src)
            eprint(f"[warn] {label}: no adjusted p column; using "
                   f"nominal '{padj_col}' for dot size.")
        size_col = pick_col(df, SETSIZE_COLS, None, False)

        sub = pd.DataFrame({
            "id": df[id_col].astype(str),
            "desc": df[desc_col].astype(str),
            "NES": pd.to_numeric(df[nes_col], errors="coerce"),
            "padj": pd.to_numeric(df[padj_col], errors="coerce"),
        })
        sub["setSize"] = (pd.to_numeric(df[size_col], errors="coerce")
                          if size_col else np.nan)
        sub["key_id"] = sub["id"].map(norm_key)
        sub["key_desc"] = sub["desc"].map(norm_key)
        tables[label] = sub
        eprint(f"[info] {label}: {len(sub)} rows from {src} "
               f"(ID='{id_col}', NES='{nes_col}', p='{padj_col}')")

    # ---- decide which pathways to plot -----------------------------------
    if pathway_entries is None:
        pool = []
        for label, sub in tables.items():
            s = sub.dropna(subset=["padj", "NES"]).sort_values("padj")
            pool.extend(s["id"].head(args.top_n).tolist())
        seen, pathway_entries = set(), []
        for p in pool:
            if norm_key(p) not in seen:
                seen.add(norm_key(p))
                pathway_entries.append((p, None))
        eprint(f"[info] no --pathways given; using union of top {args.top_n} "
               f"per table -> {len(pathway_entries)} pathways")

    # ---- match & assemble -------------------------------------------------
    records, ordered, missing_all = [], [], []
    for query, override in pathway_entries:
        key = norm_key(query)
        display = override if override else pretty_label(query, args.label_style)
        found_any = False
        canonical = None

        for label, sub in tables.items():
            hit = sub[(sub["key_id"] == key) | (sub["key_desc"] == key)]
            if hit.empty and args.partial_match:
                hit = sub[sub["key_id"].str.contains(key, regex=False, na=False)
                          | sub["key_desc"].str.contains(key, regex=False,
                                                         na=False)]
            if hit.empty:
                continue
            if len(hit) > 1:
                hit = hit.iloc[[hit["padj"].fillna(1).values.argmin()]]
            row = hit.iloc[0]
            if pd.isna(row["NES"]) or pd.isna(row["padj"]):
                continue
            found_any = True
            if canonical is None and override is None:
                canonical = pretty_label(
                    row["desc"] if args.use_description else row["id"],
                    args.label_style)
            records.append({
                "pathway": display, "comparison": label,
                "NES": float(row["NES"]), "padj": float(row["padj"]),
                "setSize": row["setSize"], "source_id": row["id"],
            })

        if not found_any:
            missing_all.append(query)
            if args.drop_missing:
                continue
        if override is None and canonical:
            for r in records:
                if r["pathway"] == display:
                    r["pathway"] = canonical
            display = canonical
        ordered.append(display)

    if missing_all:
        eprint("[warn] not found in any table: " + "; ".join(missing_all))
    if not records:
        raise SystemExit("[error] nothing to plot — no pathway matched. "
                         "Try --partial-match, or check your pathway names.")

    long = pd.DataFrame(records)

    # significance filter
    if args.max_padj is not None:
        keep = long["padj"] <= args.max_padj
        eprint(f"[info] dropping {int((~keep).sum())} dots with "
               f"padj > {args.max_padj}")
        long = long[keep]
    if args.min_nes is not None:
        long = long[long["NES"].abs() >= args.min_nes]
    if long.empty:
        raise SystemExit("[error] all dots filtered out; relax --max-padj/--min-nes.")

    # drop rows that ended up with no dots at all
    ordered = [p for p in ordered if p in set(long["pathway"])] \
        if args.drop_empty_rows else ordered
    ordered = list(dict.fromkeys(ordered))

    # ---- ordering ---------------------------------------------------------
    if args.order == "nes":
        mean_nes = long.groupby("pathway")["NES"].mean()
        ordered = list(mean_nes.reindex(ordered).sort_values(
            ascending=False).index)
    elif args.order == "padj":
        best = long.groupby("pathway")["padj"].min()
        ordered = list(best.reindex(ordered).sort_values().index)
    elif args.order == "alpha":
        ordered = sorted(ordered)
    elif args.order == "cluster":
        mat = long.pivot_table(index="pathway", columns="comparison",
                               values="NES").reindex(ordered).fillna(0)
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list
            from scipy.spatial.distance import pdist
            if len(mat) > 2:
                Z = linkage(pdist(mat.values), method="average")
                ordered = [mat.index[i] for i in leaves_list(Z)]
        except ImportError:
            eprint("[warn] scipy not available; keeping input order.")

    if args.reverse_rows:
        ordered = ordered[::-1]

    comparisons = [lab for lab, _ in inputs]
    return long, ordered, comparisons


def make_plot(long, pathways, comparisons, args):
    n_row, n_col = len(pathways), len(comparisons)

    # ---- size mapping: area ~ -log10(padj), clipped ----------------------
    floor_p = 10.0 ** (-args.p_cap)
    long = long.copy()
    long["padj_clip"] = long["padj"].clip(lower=floor_p, upper=1.0)
    long["nlp"] = -np.log10(long["padj_clip"])

    def size_of(nlp):
        frac = np.clip(nlp / args.p_cap, 0.0, 1.0) ** args.size_gamma
        return args.min_size + frac * (args.max_size - args.min_size)

    long["size"] = size_of(long["nlp"].values)

    # ---- colour mapping ---------------------------------------------------
    if args.nes_limit is not None:
        vmax = args.nes_limit
    else:
        vmax = float(np.nanmax(np.abs(long["NES"].values)))
        vmax = max(1.0, np.ceil(vmax * 2) / 2)
    cmap = plt.get_cmap(args.cmap)
    norm = Normalize(vmin=-vmax, vmax=vmax)

    # ---- figure geometry (inches) ----------------------------------------
    max_lab = max((len(p) for p in pathways), default=10)
    left = args.left_margin if args.left_margin else \
        min(4.6, 0.35 + 0.088 * max_lab)
    max_x = max((len(c) for c in comparisons), default=10)
    bottom = args.bottom_margin if args.bottom_margin else \
        min(3.2, 0.45 + 0.062 * max_x)
    right_panel = 1.85
    top = 0.35
    pw = args.col_width * n_col
    ph = args.row_height * n_row
    W = left + pw + right_panel
    H = top + ph + bottom
    if args.figsize:
        W, H = args.figsize

    fig = plt.figure(figsize=(W, H), dpi=args.dpi)
    ax = fig.add_axes([left / W, bottom / H, pw / W, ph / H])

    # ---- grid + dots ------------------------------------------------------
    x_of = {c: i for i, c in enumerate(comparisons)}
    y_of = {p: i for i, p in enumerate(pathways)}

    ax.set_xlim(-0.5, n_col - 0.5)
    ax.set_ylim(n_row - 0.5, -0.5)          # first pathway at the top
    ax.set_xticks(range(n_col))
    ax.set_yticks(range(n_row))
    ax.set_xticklabels(comparisons, rotation=args.x_rotation,
                       ha="right" if args.x_rotation else "center",
                       rotation_mode="anchor", fontsize=args.tick_fontsize)
    ax.set_yticklabels(pathways, fontsize=args.tick_fontsize)
    ax.grid(True, which="major", color="0.75", linestyle=":", linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_linewidth(1.0)
    ax.tick_params(length=3, width=1.0)

    sub = long[long["pathway"].isin(y_of) & long["comparison"].isin(x_of)]
    ax.scatter(
        [x_of[c] for c in sub["comparison"]],
        [y_of[p] for p in sub["pathway"]],
        s=sub["size"], c=sub["NES"], cmap=cmap, norm=norm,
        edgecolors=args.edge_color, linewidths=args.edge_width,
        zorder=3, clip_on=False,
    )

    # mark cells with no result
    if args.mark_missing:
        present = set(zip(sub["pathway"], sub["comparison"]))
        for p in pathways:
            for c in comparisons:
                if (p, c) not in present:
                    ax.plot(x_of[c], y_of[p], marker="x", color="0.7",
                            markersize=5, markeredgewidth=1.2, zorder=2)

    if args.title:
        ax.set_title(args.title, fontsize=args.title_fontsize, pad=10)

    # ---- NES colourbar (top right) ---------------------------------------
    cb_h = min(ph * 0.45, 2.0)
    cb_ax = fig.add_axes([(left + pw + 0.28) / W,
                          (bottom + ph - cb_h) / H,
                          0.16 / W, cb_h / H])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cb_ax)
    cb.outline.set_linewidth(0.8)
    cb.ax.tick_params(labelsize=args.legend_fontsize, length=3)
    cb.ax.set_title(args.color_label, fontsize=args.legend_fontsize + 1,
                    pad=8, loc="left")

    # ---- p.adjust size legend (below colourbar) --------------------------
    legend_ps = [p for p in args.legend_p if p >= floor_p * 0.999]
    handles = []
    for p in legend_ps:
        nlp = -np.log10(min(max(p, floor_p), 1.0))
        handles.append(Line2D([], [], linestyle="none", marker="o",
                              markerfacecolor="0.75", markeredgecolor="0.25",
                              markeredgewidth=args.edge_width,
                              markersize=np.sqrt(size_of(nlp)),
                              label=fmt_p(p)))
    leg_y = (bottom + ph - cb_h - 0.35) / H
    leg = fig.legend(handles=handles, labels=[h.get_label() for h in handles],
                     loc="upper left",
                     bbox_to_anchor=((left + pw + 0.20) / W, leg_y),
                     frameon=False, title=args.size_label,
                     labelspacing=1.15, handletextpad=1.0, borderpad=0.0,
                     fontsize=args.legend_fontsize)
    leg.get_title().set_fontsize(args.legend_fontsize + 1)
    leg._legend_box.align = "left"

    # ---- write ------------------------------------------------------------
    out = args.output
    fig.savefig(out, bbox_inches="tight", transparent=args.transparent)
    eprint(f"[done] wrote {out}")
    for extra in args.also:
        stem = os.path.splitext(out)[0]
        p2 = f"{stem}.{extra.lstrip('.')}"
        fig.savefig(p2, bbox_inches="tight", transparent=args.transparent)
        eprint(f"[done] wrote {p2}")
    plt.close(fig)


# --------------------------------------------------------------------------
def get_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Dot plot of selected pathways across multiple GSEA "
                    "result tables (colour = NES, size = adjusted p).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("QUICK START")[1] if "QUICK START" in __doc__ else None)

    io = ap.add_argument_group("input / output")
    io.add_argument("-i", "--input", action="append", required=True,
                    metavar="[LABEL=]FILE[::SHEET]",
                    help="GSEA table; repeat once per comparison. Column "
                         "order in the plot follows the order given.")
    io.add_argument("-p", "--pathways", metavar="FILE|LIST",
                    help="File with one pathway per line (optional TAB-"
                         "separated display label), or a comma-separated list. "
                         "Omit to use --top-n.")
    io.add_argument("--top-n", type=int, default=10,
                    help="If no --pathways: union of top N by padj per table "
                         "(default 10).")
    io.add_argument("-o", "--output", default="gsea_dotplot.pdf",
                    help="Output file (.pdf/.png/.svg/.tiff). Default "
                         "gsea_dotplot.pdf")
    io.add_argument("--also", nargs="*", default=[], metavar="EXT",
                    help="Additional formats to write, e.g. --also png svg")
    io.add_argument("--save-table", metavar="TSV",
                    help="Write the plotted values to a TSV.")

    cols = ap.add_argument_group("column overrides (auto-detected by default)")
    cols.add_argument("--id-col")
    cols.add_argument("--desc-col")
    cols.add_argument("--nes-col")
    cols.add_argument("--padj-col")

    sel = ap.add_argument_group("matching / filtering / ordering")
    sel.add_argument("--partial-match", action="store_true",
                     help="Fall back to substring matching of pathway names.")
    sel.add_argument("--use-description", action="store_true",
                     help="Label rows with the Description column instead of ID.")
    sel.add_argument("--max-padj", type=float, default=None,
                     help="Hide dots with padj above this (default: show all).")
    sel.add_argument("--min-nes", type=float, default=None,
                     help="Hide dots with |NES| below this.")
    sel.add_argument("--drop-missing", action="store_true",
                     help="Drop pathways absent from every table.")
    sel.add_argument("--drop-empty-rows", action="store_true",
                     help="Drop rows left with no dots after filtering.")
    sel.add_argument("--mark-missing", action="store_true",
                     help="Draw a grey x where a pathway is absent.")
    sel.add_argument("--order", default="input",
                     choices=["input", "nes", "padj", "alpha", "cluster"],
                     help="Row order (default: as listed in --pathways).")
    sel.add_argument("--reverse-rows", action="store_true")

    app = ap.add_argument_group("appearance")
    app.add_argument("--cmap", default="RdBu_r",
                     help="Diverging colormap (default RdBu_r: blue=down, "
                          "red=up).")
    app.add_argument("--nes-limit", type=float, default=None,
                     help="Symmetric NES colour limit (default: from data).")
    app.add_argument("--p-cap", type=float, default=10.0,
                     help="-log10(padj) at which dot size saturates "
                          "(default 10 = 1e-10).")
    app.add_argument("--min-size", type=float, default=25.0,
                     help="Marker area for the least significant dot.")
    app.add_argument("--max-size", type=float, default=420.0,
                     help="Marker area for the most significant dot.")
    app.add_argument("--size-gamma", type=float, default=1.0,
                     help="Exponent on the size scale (<1 boosts small dots).")
    app.add_argument("--legend-p", type=float, nargs="*",
                     default=[0.5, 0.05, 0.001, 1e-5, 1e-10],
                     help="p values shown in the size legend.")
    app.add_argument("--label-style", default="upper",
                     choices=["upper", "title", "sentence", "raw"],
                     help="How pathway names are rendered (default UPPER).")
    app.add_argument("--color-label", default="NES")
    app.add_argument("--size-label", default="p.adjust")
    app.add_argument("--title", default=None)
    app.add_argument("--edge-color", default="black")
    app.add_argument("--edge-width", type=float, default=0.8)
    app.add_argument("--col-width", type=float, default=0.95,
                     help="Inches per comparison column.")
    app.add_argument("--row-height", type=float, default=0.72,
                     help="Inches per pathway row.")
    app.add_argument("--figsize", type=float, nargs=2, default=None,
                     metavar=("W", "H"), help="Override figure size (inches).")
    app.add_argument("--left-margin", type=float, default=None)
    app.add_argument("--bottom-margin", type=float, default=None)
    app.add_argument("--x-rotation", type=float, default=45)
    app.add_argument("--tick-fontsize", type=float, default=13)
    app.add_argument("--legend-fontsize", type=float, default=11)
    app.add_argument("--title-fontsize", type=float, default=14)
    app.add_argument("--font", default=None,
                     help="Font family, e.g. Arial or Helvetica.")
    app.add_argument("--dpi", type=int, default=300)
    app.add_argument("--transparent", action="store_true")

    return ap.parse_args(argv)


def main(argv=None):
    args = get_args(argv)
    if args.font:
        plt.rcParams["font.family"] = args.font
    plt.rcParams["pdf.fonttype"] = 42     # editable text in Illustrator
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"

    inputs = parse_inputs(args.input)
    entries = load_pathways(args.pathways)
    long, pathways, comparisons = build_matrix(inputs, entries, args)

    if args.save_table:
        (long.sort_values(["pathway", "comparison"])
             .to_csv(args.save_table, sep="\t", index=False))
        eprint(f"[done] wrote {args.save_table}")

    make_plot(long, pathways, comparisons, args)


if __name__ == "__main__":
    main()
