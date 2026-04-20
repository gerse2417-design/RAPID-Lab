"""
Step 3 — AlphaFold2 / AF2-multimer structure validation (standalone).

Reads the Step 1 backbone PDB plus the Step 2 `design.fasta`, runs AF2 on each
sequence, and writes per-design PDBs and an `mpnn_results.csv` summary.

Invoked as a subprocess by run_make_ari_v4.py (Docker + host orchestrators).

Requires colabdesign on PYTHONPATH and AF2 params at `--param_dir`.

Chain 처리:
  --contigs로부터 binder/partial/fixbb 프로토콜을 자동 감지하고, 그에 맞는
  af_model을 빌드해 multi-chain PDB(target+binder 등)도 정상 처리한다.
"""
import argparse
import csv
import os
import sys

import numpy as np

from colabdesign.af import mk_af_model

# 같은 디렉토리(/app/backends)의 헬퍼
from _chain_helpers import detect_protocol


def read_fasta(path):
    """Yield (header, seq, score) tuples. Header = id before first space."""
    entries = []
    header, score, seq_buf = None, 0.0, []

    def _flush():
        if header is not None:
            entries.append((header, "".join(seq_buf), score))

    with open(path) as f:
        for raw in f:
            line = raw.rstrip()
            if line.startswith(">"):
                _flush()
                parts = line[1:].split()
                header = parts[0]
                score = 0.0
                for tok in parts[1:]:
                    if tok.startswith("score="):
                        try:
                            score = float(tok.split("=", 1)[1])
                        except ValueError:
                            pass
                seq_buf = []
            elif line:
                seq_buf.append(line)
        _flush()
    return entries


def main():
    p = argparse.ArgumentParser(description="AF2 / AF2-multimer validation")
    p.add_argument("--pdb", required=True,
                   help="Step 1 backbone PDB (used as fixbb template).")
    p.add_argument("--fasta", required=True,
                   help="Step 2 design.fasta (headers design{m}_n{n}).")
    p.add_argument("--loc", required=True,
                   help="Output directory (typically {output_dir}/af2/).")
    p.add_argument("--contigs", required=True,
                   help="Same contig string passed to RFdiffusion (post Step 0 "
                        "renumber). Used for protocol/chain detection.")
    p.add_argument("--num_recycles", type=int, default=3)
    p.add_argument("--use_multimer", action="store_true")
    p.add_argument("--initial_guess", action="store_true",
                   help="mk_af_model의 initial_guess 플래그.")
    p.add_argument("--copies", type=int, default=1,
                   help="Homooligomer copies (partial/fixbb 프로토콜).")
    p.add_argument("--rm_aa", default="C",
                   help="MPNN과 동일 인자 — chain helper 호환용 (AF2에서는 영향 없음).")
    p.add_argument("--param_dir", required=True,
                   help="Path to directory containing params/params_*.npz.")
    p.add_argument("--num_designs", type=int, default=1)
    args = p.parse_args()

    os.makedirs(args.loc, exist_ok=True)
    all_pdb_dir = os.path.join(args.loc, "all_pdb")
    os.makedirs(all_pdb_dir, exist_ok=True)

    # ── 프로토콜 감지 + af_model 빌드 ─────────────────────────────────────
    protocol, af_kwargs, prep_flags, fixed_pos = detect_protocol(
        args.contigs,
        rm_aa=args.rm_aa,
        copies=args.copies,
        use_multimer=args.use_multimer,
        initial_guess=args.initial_guess,
        param_dir=args.param_dir,
    )
    print(f"[run_af2] protocol={protocol} prep_flags_keys={list(prep_flags.keys())}")

    af = mk_af_model(**af_kwargs)

    all_entries = read_fasta(args.fasta)
    if not all_entries:
        print(f"[run_af2] ERROR: no sequences in {args.fasta}", file=sys.stderr)
        sys.exit(1)

    # CSV 컬럼: 프로토콜에 따라 i_ptm/i_pae 노출 여부 결정.
    if protocol == "binder":
        af_terms = ["plddt", "i_ptm", "ptm", "i_pae", "rmsd"]
    elif args.use_multimer or args.copies > 1:
        af_terms = ["plddt", "ptm", "i_ptm", "pae", "i_pae", "rmsd"]
    else:
        af_terms = ["plddt", "ptm", "pae", "rmsd"]
    fieldnames = ["design", "n", "mpnn"] + af_terms + ["seq"]
    rows = []

    for m in range(args.num_designs):
        if args.num_designs > 1 and args.pdb.endswith("_0.pdb"):
            pdb_m = args.pdb[:-len("_0.pdb")] + f"_{m}.pdb"
        else:
            pdb_m = args.pdb

        if not os.path.exists(pdb_m):
            print(f"[run_af2] skip missing design{m}: {pdb_m}", file=sys.stderr)
            continue

        af.prep_inputs(pdb_filename=pdb_m, **prep_flags)
        if protocol == "partial":
            p_idx = np.where(fixed_pos)[0]
            af.opt["fix_pos"] = p_idx[p_idx < af._len]

        entries_m = [
            (h, s, sc) for (h, s, sc) in all_entries
            if h.startswith(f"design{m}_")
        ]
        if not entries_m:
            print(f"[run_af2] design{m}: no FASTA entries, skipping")
            continue

        best_rmsd, best_n = float("inf"), None
        for header, seq, mpnn_score in entries_m:
            try:
                n = int(header.split("_n")[-1])
            except ValueError:
                n = 0
            sub_seq = seq[-af._len:]

            af.predict(seq=sub_seq, num_recycles=args.num_recycles, verbose=False)
            log = af.aux["log"]

            row = {
                "design": m,
                "n": n,
                "mpnn": mpnn_score,
                "seq": sub_seq,
            }
            for k in af_terms:
                row[k] = float(log.get(k, 0.0))
            rows.append(row)

            af.save_current_pdb(
                os.path.join(all_pdb_dir, f"design{m}_n{n}.pdb")
            )
            if row["rmsd"] < best_rmsd:
                best_rmsd, best_n = row["rmsd"], n

        af.save_pdb(os.path.join(args.loc, f"best_design{m}.pdb"))
        print(
            f"[run_af2] design{m}: {len(entries_m)} seqs, "
            f"best n={best_n} rmsd={best_rmsd:.2f}"
        )

    csv_path = os.path.join(args.loc, "mpnn_results.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"[run_af2] wrote {len(rows)} predictions → {csv_path}")
    if len(rows) == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
