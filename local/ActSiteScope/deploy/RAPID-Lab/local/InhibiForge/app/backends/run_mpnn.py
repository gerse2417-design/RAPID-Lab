"""
Step 2 — ProteinMPNN sequence design (standalone).

Reads backbone PDB(s) from Step 1 (RFdiffusion) output, samples sequences with
ColabDesign's ProteinMPNN wrapper, and writes them to a FASTA file that Step 3
(run_af2.py) will consume.

Invoked as a subprocess by run_make_ari_v4.py (Docker + host orchestrators).

Requires colabdesign on PYTHONPATH — caller sets
`export PYTHONPATH=$PYTHONPATH:{COLABDESIGN_DIR}` before launching this script.

Chain 처리:
  --contigs로부터 binder/partial/fixbb 프로토콜을 자동 감지하고,
  af_model을 그에 맞게 prep한 뒤 mpnn.get_af_inputs(af_model)로 chain/fix
  정보를 MPNN에 전달한다 (designability_test.py 패턴 동일).
"""
import argparse
import os
import sys

import numpy as np

from colabdesign.af import mk_af_model
from colabdesign.mpnn import mk_mpnn_model

# 같은 디렉토리(/app/backends)의 헬퍼 — sys.path[0]가 스크립트 디렉토리이므로 import 가능
from _chain_helpers import detect_protocol


def main():
    p = argparse.ArgumentParser(description="ProteinMPNN sequence design")
    p.add_argument("--pdb", required=True,
                   help="Path to Step 1 backbone PDB. For num_designs>1, the "
                        "script substitutes _0.pdb → _{m}.pdb for each design.")
    p.add_argument("--loc", required=True,
                   help="Output directory (typically {output_dir}/mpnn/).")
    p.add_argument("--contigs", required=True,
                   help="Same contig string passed to RFdiffusion (post Step 0 "
                        "renumber). Used for protocol/chain detection.")
    p.add_argument("--param_dir", required=True,
                   help="AF2 params directory (mk_af_model 요구).")
    p.add_argument("--num_seqs", type=int, default=50,
                   help="Total sequences to sample per design.")
    p.add_argument("--num_designs", type=int, default=1,
                   help="Number of Step 1 designs to iterate over.")
    p.add_argument("--batch", type=int, default=8,
                   help="MPNN batch size per sampling iteration.")
    p.add_argument("--mpnn_sampling_temp", type=float, default=0.1)
    p.add_argument("--use_soluble", action="store_true",
                   help="Use soluble MPNN weights instead of original.")
    p.add_argument("--use_multimer", action="store_true",
                   help="AF2 multimer mode (binder protocol에 권장).")
    p.add_argument("--initial_guess", action="store_true",
                   help="mk_af_model의 initial_guess 플래그.")
    p.add_argument("--copies", type=int, default=1,
                   help="Homooligomer copies (partial/fixbb 프로토콜).")
    p.add_argument("--rm_aa", default="C",
                   help="Amino acids to exclude from sampling (comma-separated).")
    args = p.parse_args()

    os.makedirs(args.loc, exist_ok=True)
    fasta_path = os.path.join(args.loc, "design.fasta")

    # ── 프로토콜 감지 + af_model 빌드 ─────────────────────────────────────
    protocol, af_kwargs, prep_flags, fixed_pos = detect_protocol(
        args.contigs,
        rm_aa=args.rm_aa,
        copies=args.copies,
        use_multimer=args.use_multimer,
        initial_guess=args.initial_guess,
        param_dir=args.param_dir,
    )
    print(f"[run_mpnn] protocol={protocol} prep_flags_keys={list(prep_flags.keys())}")

    af_model = mk_af_model(**af_kwargs)
    mpnn = mk_mpnn_model(weights="soluble" if args.use_soluble else "original")

    total_written = 0
    with open(fasta_path, "w") as f:
        for m in range(args.num_designs):
            if args.num_designs > 1 and args.pdb.endswith("_0.pdb"):
                pdb_m = args.pdb[:-len("_0.pdb")] + f"_{m}.pdb"
            else:
                pdb_m = args.pdb

            if not os.path.exists(pdb_m):
                print(f"[run_mpnn] skip missing design{m}: {pdb_m}",
                      file=sys.stderr)
                continue

            # ── chain/fix 정보를 af_model에 prep ──────────────────────────
            af_model.prep_inputs(pdb_filename=pdb_m, **prep_flags)
            if protocol == "partial":
                p_idx = np.where(fixed_pos)[0]
                af_model.opt["fix_pos"] = p_idx[p_idx < af_model._len]

            # ── af_model의 chain/fix 정보를 MPNN에 전달 ───────────────────
            mpnn.get_af_inputs(af_model)

            num_batches = max(1, args.num_seqs // max(1, args.batch))
            out = mpnn.sample(
                num=num_batches,
                batch=args.batch,
                temperature=args.mpnn_sampling_temp,
            )

            seqs = out["seq"]
            scores = out["score"]
            for n in range(len(seqs)):
                seq_clean = str(seqs[n]).replace("/", "")
                score_val = float(scores[n])
                f.write(f">design{m}_n{n} score={score_val:.4f}\n{seq_clean}\n")
            total_written += len(seqs)
            print(f"[run_mpnn] design{m}: {len(seqs)} sequences from {pdb_m}")

    print(f"[run_mpnn] wrote {total_written} sequences to {fasta_path}")
    if total_written == 0:
        print("[run_mpnn] ERROR: no sequences written", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
