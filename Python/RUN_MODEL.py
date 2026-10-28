#!/usr/bin/env python3.12
"""
Inference only: read Rust's input file, sample seeds with the trained model, write a seed file.
No training happens here (see train_model.py).

  python3.12 RUN_MODEL.py --input input.txt --output seeds.txt [--n 16] [--temp 1.0] [--ckpt PATH] [--seed 0]

Input  (blocks, any number of targets):      Output (n blocks per target):
  Sequence : SSSSNNNNSSSSNN...                 Sequence : GCGCAAAAGCGCAA...
  Structure : ((((....))))..                   Structure : ((((....))))..

Exit codes: 0 ok, 1 some targets failed (others still written), 2 fatal (no checkpoint / no output).
"""
import argparse
import random
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message="enable_nested_tensor")   # keep stderr clean for Rust

import torch
from ML_MODEL import load, sample_seeds

HERE = Path(__file__).resolve().parent
DEFAULT_CKPT = HERE / "checkpoints" / "best.pt"

COMP = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
NUCLEOTIDES = 'AUGC'
SPECIAL_CASE = 'NKS'


def comp(nuc):
    if nuc not in COMP:
        raise ValueError(f"Invalid nucleotide: {nuc!r}")
    return COMP[nuc]


def complement(seq, reverse=True):
    """Complement of an RNA string. Reverse by default, so it pairs antiparallel."""
    c = ''.join(comp(n) for n in seq)
    return c[::-1] if reverse else c


def insert_placeholder(seq, reverse=True, rng=random):
    """
    Replace placeholder letters with random sequences.
      - Uppercase placeholder (X, Q, ...): random sequence, same everywhere.
      - Lowercase placeholder (x, q, ...): complement of that sequence.
    Placeholder lengths are taken from the longest run of the uppercase
    or lowercase letter in the input.
    """
    seq_chars = list(seq)

    
    placeholders = []
    for ch in seq:
        key = ch.upper()
        if key not in NUCLEOTIDES and key not in placeholders and key not in SPECIAL_CASE:
            placeholders.append(key)

    
    lengths = {p: 0 for p in placeholders}
    i = 0
    while i < len(seq):
        key = seq[i].upper()
        if key in lengths:
            j = i
            while j < len(seq) and seq[j].upper() == key and seq[j].isupper() == seq[i].isupper():
                j += 1
            lengths[key] = max(lengths[key], j - i)
            i = j
        else:
            i += 1

    
    placeholder_seq = {
        p: ''.join(rng.choice(NUCLEOTIDES) for _ in range(lengths[p]))
        for p in placeholders
    }

    
    out = []
    i = 0
    while i < len(seq):
        ch = seq[i]
        key = ch.upper()
        if key in placeholder_seq:
            j = i
            while j < len(seq) and seq[j] == ch:
                j += 1
            run_len = j - i
            base = placeholder_seq[key]
            if ch.isupper():
                out.append(base[:run_len])
            else:
                # Complement the portion this run covers
                out.append(complement(base, reverse=reverse)[:run_len]
                           if not reverse else
                           complement(base[:run_len], reverse=True))
            i = j
        else:
            out.append(ch)
            i += 1

    return ''.join(out)

    

def read_targets(path):
    """Parse 'Sequence : <constraint>' / 'Structure : <dot-bracket>' blocks -> list of (db, cons)."""
    targets, cons, db = [], None, None
    with open(path) as f:
        for line in f:
            if ":" not in line:
                continue
            key, val = line.split(":", 1)
            key, val = key.strip().lower(), val.strip()
            if key == "sequence":
                cons = val
            elif key == "structure":
                db = val
            if cons and db:
                db = db.replace("[", ".").replace("]", ".").replace("{", "(").replace("}", ")")
                cons = insert_placeholder(cons)
                targets.append((db, cons))
                cons = db = None
    return targets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--n", type=int, default=16, help="seeds per target")
    ap.add_argument("--temp", type=float, default=1.0)
    ap.add_argument("--ckpt", default=str(DEFAULT_CKPT))
    ap.add_argument("--seed", type=int, default=None, help="set for reproducible seeds")
    a = ap.parse_args()

    if a.seed is not None:
        random.seed(a.seed)
        torch.manual_seed(a.seed)
    if not Path(a.ckpt).exists():
        print(f"ERROR: checkpoint not found: {a.ckpt} (run train_model.py pretrain)", file=sys.stderr)
        sys.exit(2)

    targets = read_targets(a.input)
    if not targets:
        print(f"ERROR: no targets parsed from {a.input}", file=sys.stderr)
        sys.exit(2)

    model = load(a.ckpt)          # loaded once, reused for every target
    failed, written = 0, 0

    
    with open(a.output, "w") as out:
        for db, cons in targets:
            try:
                seeds = sample_seeds(model, db, cons, n=a.n, temp=a.temp)
            except ValueError as e:            # length mismatch, bad characters, infeasible constraint
                print(f"WARNING: skipped target ({e}): {db}", file=sys.stderr)
                failed += 1
                continue
            for s in seeds:
                out.write(f"Sequence : {s}\nStructure : {db}\n\n")
                written += 1

    print(f"wrote {written} seeds for {len(targets) - failed}/{len(targets)} targets to {a.output}")
    sys.exit(1 if failed else (0 if written else 2))


if __name__ == "__main__":
    main()