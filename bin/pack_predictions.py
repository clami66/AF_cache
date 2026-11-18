#!/usr/bin/env python3
"""
Pack the per-pair prediction directories of one AF_cache job chunk into a
single zip archive.

run_af2_jobs / run_af3_jobs write one directory per predicted complex, each
holding many files. On cluster file systems the number of files, not their
size, is often the limit, so with --pack_outputs the pipeline stores a chunk's
predictions in one zip file instead:

    AF_data_multimer/<pair>/...   ->   AF_data_multimer/<chunk>.zip
                                           <pair>/...

Paths inside the archive are the same as on disk. Symlinks (e.g. AF2's
ranked_<n>.pdb) are stored as copies of the file they point to. The pair
directories are only removed once the archive has been written and every
member CRC-checked, so a failed pack never loses predictions.

Read the archives with unzip, or Python's zipfile module.
"""

import argparse
import os
import shutil
import sys
import zipfile
from pathlib import Path


def pack(pred_dir: Path, archive: Path, compresslevel: int) -> list:
    """Write every pair directory in pred_dir to archive; return the dirs packed."""
    pair_dirs = sorted(d for d in pred_dir.iterdir() if d.is_dir())
    if not pair_dirs:
        sys.exit(f"No prediction directories found in {pred_dir}")

    partial = archive.with_name(archive.name + ".partial")
    with zipfile.ZipFile(partial, "w", zipfile.ZIP_DEFLATED,
                         compresslevel=compresslevel) as zf:
        for pair_dir in pair_dirs:
            for root, dirs, files in os.walk(pair_dir):
                dirs.sort()
                rel_root = Path(root).relative_to(pred_dir).as_posix()
                if not dirs and not files:
                    # keep empty directories (e.g. AF2's msas/)
                    zf.writestr(rel_root + "/", "")
                for name in sorted(files):
                    path = Path(root, name)
                    if not path.exists():
                        print(f"WARN skipping broken symlink {path}")
                        continue
                    zf.write(path, f"{rel_root}/{name}")

    with zipfile.ZipFile(partial) as zf:
        bad = zf.testzip()
    if bad is not None:
        sys.exit(f"Archive check failed on {bad}; keeping {pred_dir} as it is")

    partial.replace(archive)
    return pair_dirs


def main():
    parser = argparse.ArgumentParser(
        description="Pack AF_cache prediction directories into one zip archive.")
    parser.add_argument("pred_dir", type=Path,
                        help="Directory whose subdirectories are per-pair AF outputs")
    parser.add_argument("archive", type=Path, help="Zip file to write")
    # Level 1 compresses the large result pickles nearly as well as the
    # default and is several times faster, which matters on a GPU node.
    parser.add_argument("--compresslevel", type=int, default=1)
    parser.add_argument("--keep", action="store_true",
                        help="Keep the pair directories after packing")
    args = parser.parse_args()

    pair_dirs = pack(args.pred_dir, args.archive, args.compresslevel)
    if not args.keep:
        for d in pair_dirs:
            shutil.rmtree(d)
    print(f"Packed {len(pair_dirs)} predictions into {args.archive}")


if __name__ == "__main__":
    main()
