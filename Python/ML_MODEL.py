#!/usr/bin/env python3.12


"""
Constraint-aware ML seed module for RNA inverse folding.

Input per target:  dot-bracket  +  constraint string of the same length
  A C G U : known base (copied through unchanged)
  N       : any base
  S       : G or C          (use on both sides of a pair -> GC / CG pair)
  K       : G or U          (use on both sides of a pair -> GU / UG wobble pair)
  (other IUPAC codes W R Y M B D H V are also accepted)

Decoding is masked, so every sampled seed satisfies the constraints and every base pair is a
valid pair that is compatible with both constraint characters. The model only chooses among
the options that remain.

Requires: pip install torch ViennaRNA   (import RNA; used only to label training data)
"""
import random

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Categorical
import RNA

NUC = "ACGU"
PAIRS = ["AU", "UA", "GC", "CG", "GU", "UG"]
SYM = {".": 0, "(": 1, ")": 2}
LUT = torch.full((4, 4), -1, dtype=torch.long)          # (nuc_i, nuc_j) -> pair class
for k, p in enumerate(PAIRS):
    LUT[NUC.index(p[0]), NUC.index(p[1])] = k

IUPAC = {"A": "A", "C": "C", "G": "G", "U": "U", "N": "ACGU", "S": "CG", "K": "GU",
         "W": "AU", "R": "AG", "Y": "CU", "M": "AC", "B": "CGU", "D": "AGU", "H": "ACU", "V": "ACG"}
CODES = list(IUPAC)
CODE_IDX = {c: i for i, c in enumerate(CODES)}
ALLOWED = torch.tensor([[n in IUPAC[c] for n in NUC] for c in CODES])       # (15, 4) bool


# ------------------------------------------------------------ constraints
def normalize_cons(cons, n):
    if cons is None:
        return "N" * n
    cons = cons.upper().replace("T", "U")
    if len(cons) != n:
        raise ValueError(f"constraint length {len(cons)} != structure length {n}")
    bad = set(cons) - set(IUPAC)
    if bad:
        raise ValueError(f"unsupported constraint characters: {bad}")
    return cons


def check_feasible(db, cons):
    """Every base pair must have at least one valid pair type allowed by both constraint characters."""
    for i, j in enumerate(pair_table(db)):
        if j > i and not any(a in IUPAC[cons[i]] and b in IUPAC[cons[j]] for a, b in PAIRS):
            raise ValueError(f"constraint {cons[i]}-{cons[j]} at positions {i},{j} allows no valid base pair")


def satisfies(seq, cons):
    return all(c in IUPAC[k] for c, k in zip(seq, cons))


def random_constraint(db, seq):
    """Training augmentation: derive a constraint string from a known (structure, sequence) pair."""
    n = len(db)
    if random.random() < 0.2:
        return "N" * n                                  # keep an unconstrained mode
    pt, p_known = pair_table(db), random.uniform(0.0, 0.6)
    out = ["N"] * n
    for i, j in enumerate(pt):
        if j == -1:
            if random.random() < p_known:
                out[i] = seq[i]
            elif random.random() < 0.15 and seq[i] in "GC":
                out[i] = "S"
            elif random.random() < 0.15 and seq[i] in "GU":
                out[i] = "K"
        elif j > i:
            pair, r = seq[i] + seq[j], random.random()
            if r < p_known:
                out[i], out[j] = seq[i], seq[j]
            elif r < p_known + 0.15:                    # one side known, partner free
                k = random.choice([i, j])
                out[k] = seq[k]
            elif random.random() < 0.6:                 # symbolic pair type
                if pair in ("GC", "CG"):
                    out[i] = out[j] = "S"
                elif pair in ("GU", "UG"):
                    out[i] = out[j] = "K"
    return "".join(out)


# ------------------------------------------------------------ structure features
def pair_table(db):
    stack, pt = [], [-1] * len(db)
    for i, c in enumerate(db):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            pt[i], pt[j] = j, i
    return pt


def loop_types(db):
    """0 external, 1 hairpin, 2 bulge/interior, 3 multiloop, 4 stem, 5 stem end."""
    pt, n = pair_table(db), len(db)
    children, parent, stack = {}, [-1] * n, []
    for i, c in enumerate(db):
        if c == "(":
            children[i] = 0
            if stack:
                children[stack[-1]] += 1
            stack.append(i)
        elif c == ")":
            stack.pop()
        else:
            parent[i] = stack[-1] if stack else -1
    out = []
    for i, c in enumerate(db):
        if c == ".":
            p = parent[i]
            out.append(0 if p < 0 else min(children[p], 2) + 1)
        else:
            o = i if c == "(" else pt[i]
            stem_end = not (o + 1 < n and pt[o + 1] == pt[o] - 1)
            out.append(5 if stem_end else 4)
    return out


def collate(items):
    """items: list of (db, seq_or_None, cons_or_None) -> padded tensors.
    cons=None with a known seq -> random constraint (training); with no seq -> all N."""
    B, L = len(items), max(len(d) for d, _, _ in items)
    sym = torch.zeros(B, L, dtype=torch.long)
    loop = torch.zeros(B, L, dtype=torch.long)
    pt = torch.full((B, L), -1, dtype=torch.long)
    seq = torch.zeros(B, L, dtype=torch.long)
    cons_t = torch.full((B, L), CODE_IDX["N"], dtype=torch.long)
    pad = torch.ones(B, L, dtype=torch.bool)
    for b, (db, s, c) in enumerate(items):
        n = len(db)
        if c is None:
            c = random_constraint(db, s) if s is not None else "N" * n
        sym[b, :n] = torch.tensor([SYM[ch] for ch in db])
        loop[b, :n] = torch.tensor(loop_types(db))
        pt[b, :n] = torch.tensor(pair_table(db))
        cons_t[b, :n] = torch.tensor([CODE_IDX[ch] for ch in c])
        pad[b, :n] = False
        if s is not None:
            seq[b, :n] = torch.tensor([NUC.index(ch) for ch in s])
    return sym, loop, pt, pad, seq, cons_t


# ------------------------------------------------------------ model
class SeedModel(nn.Module):
    """Transformer with partner-biased attention and constraint-masked output heads."""

    def __init__(self, d=128, heads=4, layers=4, maxlen=512):
        super().__init__()
        self.cfg = dict(d=d, heads=heads, layers=layers, maxlen=maxlen)
        self.heads = heads
        self.sym, self.loop = nn.Embedding(3, d), nn.Embedding(6, d)
        self.cons, self.pos = nn.Embedding(len(CODES), d), nn.Embedding(maxlen, d)
        layer = nn.TransformerEncoderLayer(d, heads, 4 * d, 0.1, batch_first=True, norm_first=True)
        self.enc = nn.TransformerEncoder(layer, layers)
        self.pbias = nn.Parameter(torch.tensor(2.0))
        self.nuc_head = nn.Linear(d, 4)
        self.pair_head = nn.Linear(d, 6)
        self.register_buffer("allowed_tab", ALLOWED, persistent=False)
        self.register_buffer("pair_a", torch.tensor([NUC.index(p[0]) for p in PAIRS]), persistent=False)
        self.register_buffer("pair_b", torch.tensor([NUC.index(p[1]) for p in PAIRS]), persistent=False)

    def forward(self, sym, loop, pt, pad, cons):
        B, L = sym.shape
        x = (self.sym(sym) + self.loop(loop) + self.cons(cons)
             + self.pos(torch.arange(L, device=sym.device))[None])
        bias = torch.zeros(B, L, L, device=sym.device)
        b, i = (pt >= 0).nonzero(as_tuple=True)
        bias[b, i, pt[b, i]] = 1.0
        mask = (bias * self.pbias).masked_fill(pad[:, None, :], -1e9).repeat_interleave(self.heads, 0)
        h = self.enc(x, mask=mask)
        pidx = pt.clamp(min=0)
        h_partner = h.gather(1, pidx[..., None].expand_as(h))

        allowed = self.allowed_tab[cons]                                        # (B, L, 4)
        allowed_partner = allowed.gather(1, pidx[..., None].expand(-1, -1, 4))
        pair_ok = allowed[..., self.pair_a] & allowed_partner[..., self.pair_b]  # (B, L, 6)

        nl = self.nuc_head(h).masked_fill(~allowed, -1e9)
        pl = self.pair_head(h + h_partner).masked_fill(~pair_ok, -1e9)
        return nl, pl, allowed, pair_ok


def loss_fn(model, batch):
    sym, loop, pt, pad, seq, cons = batch
    nl, pl, allowed, pair_ok = model(sym, loop, pt, pad, cons)
    tgt = LUT.to(seq.device)[seq, seq.gather(1, pt.clamp(min=0))]
    # only positions with a real choice (>1 allowed option) carry learning signal
    unpaired = (pt < 0) & ~pad & (allowed.sum(-1) > 1)
    opening = (sym == 1) & (tgt >= 0) & (pair_ok.sum(-1) > 1)
    loss = nl.sum() * 0.0
    if unpaired.any():
        loss = loss + F.cross_entropy(nl[unpaired], seq[unpaired])
    if opening.any():
        loss = loss + F.cross_entropy(pl[opening], tgt[opening])
    return loss


@torch.no_grad()
def sample_seeds(model, db, cons=None, n=16, temp=1.0):
    """n seeds for one target. Constraints are enforced by masking, pairs are sampled jointly."""
    model.eval()
    cons = normalize_cons(cons, len(db))
    check_feasible(db, cons)
    dev = next(model.parameters()).device
    sym, loop, pt, pad, _, cons_t = [t.to(dev) for t in collate([(db, None, cons)] * n)]
    nl, pl, _, _ = model(sym, loop, pt, pad, cons_t)
    nuc = Categorical(logits=nl / temp).sample().cpu()
    pr = Categorical(logits=pl / temp).sample().cpu()
    pt_l, seeds = pair_table(db), []
    for k in range(n):
        s = [NUC[nuc[k, i]] for i in range(len(db))]
        for i, j in enumerate(pt_l):
            if j > i:
                s[i], s[j] = PAIRS[pr[k, i]]
        s = "".join(s)
        assert satisfies(s, cons), "decoded seed violates constraint (bug)"
        seeds.append(s)
    return seeds


# ------------------------------------------------------------ data
def random_seq(n):
    gc = random.uniform(0.3, 0.7)
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(random.choices(NUC, w, k=n))


def synth_pairs(n_samples, lo=20, hi=100):
    """(MFE structure, sequence, None) triples. The None means a fresh random constraint is
    derived every time the example is batched. TODO: broaden structure coverage."""
    out = []
    for _ in range(n_samples):
        s = random_seq(random.randint(lo, hi))
        out.append((RNA.fold(s)[0], s, None))
    return out


def relabel(items):
    """Hindsight relabeling. items: list of (seq, target_db_or_None, cons_or_None).
    Each sequence is paired with the structure it ACTUALLY folds into. The original constraint
    is kept only if the design solved its target; otherwise a random constraint is derived."""
    out = []
    for seq, target, cons in items:
        st = RNA.fold(seq)[0]
        out.append((st, seq, cons if (cons and target == st) else None))
    return out


def train(model, data, epochs=3, bs=32, lr=3e-4):
    dev = next(model.parameters()).device
    opt = torch.optim.AdamW(model.parameters(), lr=lr)
    model.train()
    data = list(data)
    for _ in range(epochs):
        random.shuffle(data)
        for k in range(0, len(data), bs):
            opt.zero_grad()
            batch = [t.to(dev) for t in collate(data[k:k + bs])]
            loss_fn(model, batch).backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()


def update_model(model, items, buffer, fresh_frac=0.3, epochs=2):
    """One self-training step. items: list of (seq, target_db, cons) from your Rust generator
    (target_db / cons may be None). Includes failed designs: they are relabeled, not discarded."""
    buffer.extend(relabel(items))
    fresh = synth_pairs(int(fresh_frac * len(buffer)))      # guards against drift / collapse
    train(model, buffer + fresh, epochs=epochs)
    return buffer


def save(model, path):
    torch.save({"cfg": model.cfg, "state": model.state_dict()}, path)


def load(path):
    ck = torch.load(path, map_location="cpu")
    m = SeedModel(**ck["cfg"])
    m.load_state_dict(ck["state"])
    return m

def read_input_file(path):
    """Parse blocks of 'Sequence : <constraint>' and 'Structure : <dot-bracket>'."""
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
                targets.append((db, cons))
                cons = db = None
    return targets


def read_design_file(path):
    """Rust training-data file -> list of (seq, target_db, None) for update_model."""
    items, seq, st = [], None, None
    try:
        with open(path) as f:
            for line in f:
                if ":" not in line:
                    continue
                key, val = line.split(":", 1)
                key, val = key.strip().lower(), val.strip()
                if key == "sequence":
                    seq = val.upper().replace("T", "U")
                elif key == "structure":
                    st = val
                if seq and st:
                    items.append((seq, st, None))
                    seq = st = None
    except FileNotFoundError:
        pass
    return items

def make_seed_file(model, in_path, out_path, n=16, temp=1.0):
    with open(out_path, "w") as f:
        for db, cons in read_input_file(in_path):
            for s in sample_seeds(model, db, cons, n=n, temp=temp):
                f.write(f"Sequence  : {s}\nStructure : {db}\n\n")