#!/usr/bin/env python3
from pathlib import Path
import subprocess
import click
import pandas as pd

COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand", "stitle"
]


def run_blast(fasta, blastdb, out, threads, evalue, max_targets):
    cmd = [
        "blastn",
        "-task", "blastn",
        "-query", str(fasta),
        "-db", str(blastdb),
        "-outfmt", "6 " + " ".join(COLS),
        "-out", str(out),
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_targets),
        "-dust", "no",
        "-num_threads", str(threads),
    ]
    subprocess.run(cmd, check=True)
    if out.stat().st_size == 0:
        return pd.DataFrame(columns=COLS)
    return pd.read_csv(out, sep='\t', header=None, names=COLS)


@click.command()
@click.option('--probes', required=True, help='probe_qc.tsv from the QC step')
@click.option('--human_blastdb', multiple=True, required=True, help='Human RNA BLAST database (repeatable)')
@click.option('--output_dir', required=True, help='Directory to write outputs')
@click.option('--threads', default=8, show_default=True)
@click.option('--evalue', default='1e-5', show_default=True)
@click.option('--max_targets', default=100, show_default=True)
@click.option('--pident_min', default=90.0, show_default=True)
@click.option('--qcov_min', default=0.9, show_default=True)
def main(probes, human_blastdb, output_dir, threads, evalue, max_targets, pident_min, qcov_min):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta = output_dir / "human_blast_query.fa"
    raw_out = output_dir / "human_blast_results.tsv"
    filtered_out = output_dir / "human_blast_filtered.tsv"
    summary_out = output_dir / "human_blast_summary.tsv"

    probe_df = pd.read_csv(probes, sep='\t')

    with open(fasta, "w") as fh:
        for probe_id, seq in zip(probe_df["probe_id"], probe_df["hyb_region_full"]):
            fh.write(f">{probe_id}\n{seq}\n")

    lengths = dict(zip(probe_df["probe_id"], probe_df["hyb_region_full"].str.len()))

    frames = []
    for blastdb in human_blastdb:
        name = Path(blastdb).name
        tmp_out = output_dir / f"_human_{name}.tsv"
        df_db = run_blast(fasta, blastdb, tmp_out, threads, evalue, max_targets)
        df_db.insert(0, "db", name)
        frames.append(df_db)
        tmp_out.unlink()

    df = pd.concat(frames, ignore_index=True)
    df["qcov"] = df["length"] / df["qseqid"].map(lengths)
    df = df.sort_values(["qseqid", "bitscore"], ascending=[True, False])
    df.to_csv(raw_out, sep='\t', index=False)

    df_pass = df[(df["pident"] >= pident_min) & (df["qcov"] >= qcov_min)].copy()
    df_pass.to_csv(filtered_out, sep='\t', index=False)

    summary = probe_df[["probe_id"]].copy()
    if df_pass.empty:
        summary["human_n_hits"] = 0
        summary["human_hit"] = False
    else:
        best = df_pass.groupby("qseqid", as_index=False).first()
        best = best[["qseqid", "db", "sseqid", "stitle", "pident", "qcov", "sstrand"]]
        best.columns = ["probe_id", "human_best_db", "human_best_subject", "human_best_title",
                        "human_best_pident", "human_best_qcov", "human_best_strand"]

        counts = df_pass.groupby("qseqid").size().reset_index()
        counts.columns = ["probe_id", "human_n_hits"]

        summary = summary.merge(best, on="probe_id", how="left")
        summary = summary.merge(counts, on="probe_id", how="left")
        summary["human_n_hits"] = summary["human_n_hits"].fillna(0).astype(int)
        summary["human_hit"] = summary["human_n_hits"] > 0

    summary.to_csv(summary_out, sep='\t', index=False)

    n_hit = int(summary["human_hit"].sum())
    print(f"{len(df)} alignments, {len(df_pass)} passing "
          f"pident>={pident_min} & qcov>={qcov_min}, "
          f"{n_hit} of {len(probe_df)} probes with a passing hit")


if __name__ == '__main__':
    main()
