"""
Usage
-----
$ python blast_project.py \
    --blast Project1_BLASToutput_2024.txt \
    --fasta wildolive_GCF_002742605.1_protein.faa \
    --outdir results

This single script performs all three required tasks:
1. Generates two hit‑lists with the requested filters and prints comparison statistics.
2. Builds a heat‑map (identity and E‑value) for the shared hits.
3. Extracts the shared hit sequences into a new FASTA file.
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
from typing import Iterable, Set
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ──────────────────────────────────────────────────────────────────────────────
#                                   CONSTANTS
# ──────────────────────────────────────────────────────────────────────────────
IDENTITY_THRESHOLD = 35.0    
EVALUE_THRESHOLD = 1e-7         

# Column names expected in the BLAST tab‑separated file
BLAST_COLUMNS = [
    "query_ID",
    "subject_ID",
    "%identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "q_start",
    "q_end",
    "s_start",
    "s_end",
    "evalue",
    "bit_score",
]

# ──────────────────────────────────────────────────────────────────────────────
#                                 CORE LOGIC
# ──────────────────────────────────────────────────────────────────────────────

def load_blast_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", names=BLAST_COLUMNS, header=0)
    df["%identity"] = pd.to_numeric(df["%identity"], errors="coerce")
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    return df


def make_lists(df: pd.DataFrame, outdir: Path) -> tuple[Set[str], Set[str]]:
    list1 = df.loc[df["%identity"] >= IDENTITY_THRESHOLD, "subject_ID"].unique()
    list2 = df.loc[df["evalue"] <= EVALUE_THRESHOLD, "subject_ID"].unique()
    (outdir / "list1_identity_ge35.txt").write_text("\n".join(list1) + "\n")
    (outdir / "list2_evalue_le1e-7.txt").write_text("\n".join(list2) + "\n")
    return set(list1), set(list2)


def compare_lists(list1: Set[str], list2: Set[str]) -> None:
    common = list1 & list2
    only_a = list1 - list2
    only_b = list2 - list1

    print("\n─── List comparison ───────────────────────────────────────────────")
    print(f"Hits fulfilling (a) only         : {len(only_a):3d}")
    print(f"Hits fulfilling (b) only         : {len(only_b):3d}")
    print(f"Hits fulfilling both criteria   : {len(common):3d}\n")

    if len(list1) > len(list2):
        print("Criterion (a) (identity ≥ 35 %) retrieves more candidates more than (b).")
    elif len(list1) < len(list2):
        print("Criterion (b) (E‑value ≤ 1 × 10⁻⁷) retrieves more candidates more than (a).")
    else:
        print("Both criteria retrieve the same number of candidates.")


def plot_heatmap(df: pd.DataFrame, shared_ids: Iterable[str], out_png: Path) -> None:
    subset = df[df["subject_ID"].isin(shared_ids)].copy()
    subset = subset.drop_duplicates("subject_ID")
    subset.sort_values("%identity", ascending=False, inplace=True)

    identity = subset["%identity"].to_numpy().reshape(-1, 1)
    evalue = subset["evalue"].to_numpy().reshape(-1, 1)

    fig, axes = plt.subplots(
        ncols=2,
        figsize=(4, 0.35 * len(subset) + 1),
        gridspec_kw={"width_ratios": [1, 1]},
    )

    im0 = axes[0].imshow(identity, aspect="auto", cmap="viridis")
    im1 = axes[1].imshow(evalue, aspect="auto", cmap="magma_r")

    # heatmap
    axes[0].set_yticks(np.arange(len(subset)))
    axes[0].set_yticklabels(subset["subject_ID"], fontsize=7)
    axes[1].set_yticks([])
    axes[0].set_xticks([0])
    axes[1].set_xticks([0])
    axes[0].set_xticklabels(["% identity"], rotation=90)
    axes[1].set_xticklabels(["E‑value"], rotation=90)
    cbar_kw = {"fraction": 0.046, "pad": 0.04}
    fig.colorbar(im0, ax=axes[0], **cbar_kw, label="% identity")
    fig.colorbar(im1, ax=axes[1], **cbar_kw, label="E‑value")
    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
#                       FASTA SEQUENCE EXTRACTION
# ──────────────────────────────────────────────────────────────────────────────

def fasta_iterator(path: Path) -> Iterable[tuple[str, str]]:
    with path.open() as fh:
        header = None
        seq_lines: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield header, "".join(seq_lines)


def extract_sequences(fasta_path: Path, wanted_ids: Set[str], out_fasta: Path) -> None:
    with out_fasta.open("w") as out:
        for header, seq in fasta_iterator(fasta_path):
            primary_id = header.split()[0]
            if primary_id in wanted_ids:
                out.write(f">{header}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i : i + 60] + "\n")

# ──────────────────────────────────────────────────────────────────────────────
#                                   MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="...",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("results"))
    args = parser.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)

    # ── Task 1 ────────────────────────────────────────────────────────────
    blast_df = load_blast_table(args.blast)
    list1, list2 = make_lists(blast_df, args.outdir)    
    compare_lists(list1, list2)

    shared = list1 & list2
    if not shared:
        print("No shared hits – nothing more to do.")
        sys.exit()

    # ── Task 2 ────────────────────────────────────────────────────────────
    heatmap_png = args.outdir / "shared_hits_heatmap.png"
    plot_heatmap(blast_df, shared, heatmap_png)
    print(f"Heat‑map saved as: {heatmap_png}")

    # ── Task 3 ────────────────────────────────────────────────────────────
    fasta_out = args.outdir / "shared_hits_sequences.fasta"
    extract_sequences(args.fasta, shared, fasta_out)
    print(f"FASTA of shared hit sequences written to: {fasta_out}\n")


if __name__ == "__main__":
    main()
