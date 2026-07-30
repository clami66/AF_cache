#!/usr/bin/env python3
"""
Batch ipSAE scoring for the AF_cache pipeline.

This is a thin batch wrapper around ipsae.py (Roland Dunbrack, Fox Chase Cancer
Center; MIT license). If you use the ipSAE scores produced here, please cite:

    Dunbrack RL. ipSAE: a scoring metric for the confidence of AlphaFold and
    other structure prediction models in protein-protein interactions.
    bioRxiv 2025.02.10.637595. https://doi.org/10.1101/2025.02.10.637595

Runs ipsae.py on every model of every prediction in an AF2 (AF_data_multimer/)
or AF3 (AF3_data_multimer/) output directory and writes a single TSV of ipSAE
scores per pair. No intermediate files are kept: ipsae.py is run inside a
temporary directory and all of its side-files (.txt / _byres.txt / .pml) are
discarded once the score has been read.

Output columns:
    * if every pair has exactly one model:  pair    ipsae
    * otherwise (multiple models per pair):  pair    min_ipsae    max_ipsae    avg_ipsae

AF2 layout (one dir per pair):
    AF_data_multimer/<pair>/
        ranking_debug.json
        <pair>--ranked_0.pdb ...
        result_<model_name>.pkl ...

AF3 layout (one dir per pair):
    AF3_data_multimer/<pair>/
        seed-<s>_sample-<i>/
            <pair>_seed-<s>_sample-<i>_model.cif
            <pair>_seed-<s>_sample-<i>_confidences.json
            <pair>_seed-<s>_sample-<i>_summary_confidences.json
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from statistics import mean


def pae_dist_strings(pae_cutoff: float, dist_cutoff: float):
    pae_str = str(int(pae_cutoff))
    if pae_cutoff < 10:
        pae_str = "0" + pae_str
    dist_str = str(int(dist_cutoff))
    if dist_cutoff < 10:
        dist_str = "0" + dist_str
    return pae_str, dist_str


def detect_backend(run_dir: Path) -> str:
    """Guess whether a prediction dir is AF2 or AF3 output."""
    if (run_dir / "ranking_debug.json").exists():
        return "af2"
    if any(run_dir.glob("seed-*_sample-*")) or any(run_dir.glob("*_model.cif")):
        return "af3"
    return "unknown"


def get_af2_models(run_dir: Path):
    """Return list of (structure_pdb, pae_pkl, extra_files) in ranking order."""
    with open(run_dir / "ranking_debug.json") as f:
        order = json.load(f)["order"]
    models = []
    for rank, model_name in enumerate(order):
        # AF2 PDBs are prefixed with the pair id: e.g. <pair>--ranked_0.pdb
        pdbs = sorted(run_dir.glob(f"*ranked_{rank}.pdb"))
        pkl = run_dir / f"result_{model_name}.pkl"
        if not pkl.exists():
            alt = sorted(run_dir.glob(f"*result_{model_name}.pkl"))
            pkl = alt[0] if alt else pkl
        if pdbs and pkl.exists():
            models.append((pdbs[0], pkl, []))
    return models


def get_af3_models(run_dir: Path):
    """Return list of (structure_cif, pae_json, extra_files) for each sample."""
    models = []
    sample_dirs = sorted(d for d in run_dir.glob("seed-*_sample-*") if d.is_dir())
    search_dirs = sample_dirs if sample_dirs else [run_dir]
    for d in search_dirs:
        cifs = sorted(d.glob("*_model.cif"))
        paes = sorted(p for p in d.glob("*_confidences.json")
                      if "summary" not in p.name)
        summaries = sorted(d.glob("*_summary_confidences.json"))
        if cifs and paes:
            # ipsae.py derives the summary path from the pae path, so it must be
            # staged alongside in the temp dir.
            models.append((cifs[0], paes[0], summaries[:1]))
    return models


def parse_ipsae_from_txt(txt_path: Path) -> float:
    """Extract the (best) ipSAE from the 'max' rows of an ipsae .txt file."""
    values = []
    with open(txt_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 6 and parts[4] == "max":
                try:
                    values.append(float(parts[5]))
                except ValueError:
                    pass
    if not values:
        raise ValueError(f"No 'max' row found in {txt_path}")
    return max(values)


def run_ipsae_for_model(ipsae_script, structure, pae, extra, pae_cutoff,
                        dist_cutoff, python_exe=None) -> float:
    """Run ipsae.py in a temp dir so side-files never land next to the inputs."""
    interpreter = python_exe if python_exe else sys.executable
    tmpdir = Path(tempfile.mkdtemp())
    try:
        tmp_struct = tmpdir / structure.name
        tmp_pae = tmpdir / pae.name
        tmp_struct.symlink_to(structure.resolve())
        tmp_pae.symlink_to(pae.resolve())
        for extra_file in extra:
            (tmpdir / extra_file.name).symlink_to(extra_file.resolve())

        result = subprocess.run(
            [interpreter, str(ipsae_script), str(tmp_pae), str(tmp_struct),
             str(pae_cutoff), str(dist_cutoff)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr.decode().strip())

        pae_str, dist_str = pae_dist_strings(pae_cutoff, dist_cutoff)
        ext = ".pdb" if tmp_struct.suffix == ".pdb" else ".cif"
        stem = str(tmp_struct).replace(ext, f"_{pae_str}_{dist_str}")
        return parse_ipsae_from_txt(Path(stem + ".txt"))
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _worker(task):
    return run_ipsae_for_model(
        Path(task["ipsae_script"]), Path(task["structure"]), Path(task["pae"]),
        [Path(e) for e in task["extra"]],
        task["pae_cutoff"], task["dist_cutoff"], python_exe=task["python_exe"],
    )


def main():
    parser = argparse.ArgumentParser(description="Batch ipSAE scoring for AF_cache.")
    parser.add_argument("--predictions-dir", required=True, type=Path,
                        help="Directory whose subdirectories are per-pair AF outputs")
    parser.add_argument("--ipsae-script", required=True, type=Path,
                        help="Path to ipsae.py")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV file")
    parser.add_argument("--backend", choices=["af2", "af3", "auto"], default="auto",
                        help="Prediction backend (default: auto-detect per pair)")
    parser.add_argument("--python", type=str, default=None,
                        help="Python interpreter for ipsae.py (default: same as this script)")
    parser.add_argument("--pae-cutoff", type=float, default=10.0)
    parser.add_argument("--dist-cutoff", type=float, default=10.0)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    args = parser.parse_args()

    run_dirs = sorted(d for d in args.predictions_dir.iterdir() if d.is_dir())
    if not run_dirs:
        sys.exit(f"No prediction subdirectories found in {args.predictions_dir}")

    # Build a flat task list across all pairs / models.
    tasks = []
    per_pair_expected = defaultdict(int)
    for run_dir in run_dirs:
        pair = run_dir.name
        backend = args.backend if args.backend != "auto" else detect_backend(run_dir)
        if backend == "af2":
            models = get_af2_models(run_dir)
        elif backend == "af3":
            models = get_af3_models(run_dir)
        else:
            print(f"SKIP {pair}: could not determine backend / no models found")
            continue
        if not models:
            print(f"SKIP {pair}: no scorable models found")
            continue
        for idx, (structure, pae, extra) in enumerate(models):
            per_pair_expected[pair] += 1
            tasks.append({
                "pair": pair, "model_idx": idx,
                "ipsae_script": str(args.ipsae_script),
                "structure": str(structure), "pae": str(pae),
                "extra": [str(e) for e in extra],
                "pae_cutoff": args.pae_cutoff, "dist_cutoff": args.dist_cutoff,
                "python_exe": args.python,
            })

    print(f"Submitting {len(tasks)} models across {len(per_pair_expected)} pairs "
          f"using {args.workers} workers...")

    results = defaultdict(list)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(_worker, t): t for t in tasks}
        for future in as_completed(futures):
            t = futures[future]
            label = f"{t['pair']}/model_{t['model_idx']}"
            try:
                val = future.result()
                results[t["pair"]].append(val)
                print(f"  done {label}  ipsae={val:.4f}")
            except Exception as e:
                print(f"  WARN {label}: {e}")

    # A pair gets multi-model columns only if it actually has >1 scored model.
    multi_model = any(len(v) > 1 for v in results.values())

    rows = []
    for run_dir in run_dirs:
        pair = run_dir.name
        vals = results.get(pair)
        if not vals:
            continue
        if multi_model:
            rows.append({
                "pair": pair,
                "min_ipsae": round(min(vals), 6),
                "max_ipsae": round(max(vals), 6),
                "avg_ipsae": round(mean(vals), 6),
            })
        else:
            rows.append({"pair": pair, "ipsae": round(vals[0], 6)})

    fieldnames = (["pair", "min_ipsae", "max_ipsae", "avg_ipsae"]
                  if multi_model else ["pair", "ipsae"])
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} pairs to {args.output}")


if __name__ == "__main__":
    main()
