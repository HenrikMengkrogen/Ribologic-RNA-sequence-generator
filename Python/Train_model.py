#!/usr/bin/env python3.12
"""
Training driver for ML_MODEL.py.

  python train_model.py pretrain [--n 20000] [--epochs 30]
      synthetic data -> train with validation / early stopping -> checkpoints/best.pt + buffer.jsonl

  python train_model.py update [--epochs 2]
      reads Rust's training_data_positive/negative.txt -> one self-training round
      -> checkpoints/round_XX.pt + best.pt + buffer.jsonl; consumed Rust files are archived

Expects ML_MODEL.py (your renamed seed_model.py, including read_design_file) in the same folder.
"""
import argparse
import json
import random
import time
from pathlib import Path

import torch
from ML_MODEL import (SeedModel, synth_pairs, train, loss_fn, collate,
                      save, load, update_model, read_design_file)

# Define this first.
HERE = Path(__file__).resolve().parent

# Then derive paths from HERE.
CKPT_DIR = HERE / "checkpoints"
BEST = CKPT_DIR / "best.pt"
BUFFER = CKPT_DIR / "buffer.jsonl"

RUST_DIR = HERE / "Training_data"
RUST_FILES = [
    RUST_DIR / "training_data_positive.txt",
    RUST_DIR / "training_data_negative.txt",
]
ARCHIVE = RUST_DIR / "archive"

# Create required output directories.
CKPT_DIR.mkdir(parents=True, exist_ok=True)
RUST_DIR.mkdir(parents=True, exist_ok=True)
# ------------------------------------------------------------ data io
def save_data(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for st, seq, cons in data:
            f.write(json.dumps({"struct": st, "seq": seq, "cons": cons}) + "\n")


def load_data(path):
    with open(path) as f:
        rows = [json.loads(l) for l in f if l.strip()]
    return [(r["struct"], r["seq"], r["cons"]) for r in rows]


# ------------------------------------------------------------ helpers
@torch.no_grad()
def val_loss(model, val, bs=64):
    """Mean validation loss. Constraints are re-randomised each call, so expect some noise."""
    model.eval()
    dev = next(model.parameters()).device
    losses = [loss_fn(model, [t.to(dev) for t in collate(val[k:k + bs])]).item()
              for k in range(0, len(val), bs)]
    return sum(losses) / len(losses)


def next_round_number():
    return len(list(CKPT_DIR.glob("round_*.pt")))


# ------------------------------------------------------------ modes
def pretrain(args):
    random.seed(args.seed)
    torch.manual_seed(args.seed)
    data = synth_pairs(args.n)
    random.shuffle(data)
    val, tr = data[:args.n_val], data[args.n_val:]

    model = SeedModel()
    best, bad = float("inf"), 0
    for epoch in range(args.epochs):
        train(model, tr, epochs=1)
        vl = val_loss(model, val)
        print(f"epoch {epoch}: val loss {vl:.4f}", flush=True)
        if vl < best - 1e-3:
            best, bad = vl, 0
            save(model, BEST)
        else:
            bad += 1
            if epoch >= 5 and bad >= args.patience:
                print("early stop")
                break
    save_data(BUFFER, tr)                      # replay base for later update rounds
    print(f"best val loss {best:.4f}; saved {BEST} and {BUFFER} ({len(tr)} examples)")


def update(args):
    if not BEST.exists():
        raise SystemExit("no checkpoint found: run `pretrain` first")
    items = []
    for p in RUST_FILES:
        items += read_design_file(str(p))
    if not items:
        raise SystemExit("no new designs from Rust, nothing to do")

    model = load(BEST)
    buffer = load_data(BUFFER) if BUFFER.exists() else synth_pairs(5000)
    print(f"{len(items)} new designs, buffer had {len(buffer)} examples")

    buffer = update_model(model, items, buffer, epochs=args.epochs)

    r = next_round_number()
    save(model, CKPT_DIR / f"round_{r:02d}.pt")   # versioned, so a bad round can be rolled back
    save(model, BEST)
    save_data(BUFFER, buffer)

    # archive consumed Rust files so the next round only sees new designs
    ARCHIVE.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    for p in RUST_FILES:
        if p.exists():
            p.rename(ARCHIVE / f"{p.stem}_{stamp}{p.suffix}")
    print(f"saved round_{r:02d}.pt, buffer now {len(buffer)} examples")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    p = sub.add_parser("pretrain")
    p.add_argument("--n", type=int, default=20000)
    p.add_argument("--n_val", type=int, default=1000)
    p.add_argument("--epochs", type=int, default=30)
    p.add_argument("--patience", type=int, default=4)
    p.add_argument("--seed", type=int, default=0)
    u = sub.add_parser("update")
    u.add_argument("--epochs", type=int, default=2)
    args = ap.parse_args()
    pretrain(args) if args.mode == "pretrain" else update(args)