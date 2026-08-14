"""kfm extend -- add your own measurements to a released KFM model.

Works on BOTH released v2 models. The layout is read from the bundle's
MANIFEST.json, so you point it at a bundle and it does the right thing:

  LigASeqLigB   potency      [ ligand A 1038 | sequence 480 | ligand B 1038 ] = 2556
  SeqALigSeqB   selectivity  [ sequence A 480 | ligand 1038 | sequence B 480 ] = 1998

Neither side shares data. You fit trees on your machine and they are concatenated
into the released forest; our corpus never reaches you and yours never reaches us.
What crosses the boundary is the feature recipe, which both bundles already ship.

    python kfm_extend.py --base ./LigASeqLigB_v2_potency     --data mine.csv --out ./mine
    python kfm_extend.py --base ./SeqALigSeqB_v2_selectivity --data mine.csv --out ./mine

CSV, one comparison per row, header required. Unknown columns are ignored.

POTENCY (two ligands, one target)
    smiles_a, smiles_b                 the two ligands             REQUIRED
    gene  OR  sequence                 the target                  REQUIRED

SELECTIVITY (one ligand, two targets)
    smiles                             the ligand                  REQUIRED
    gene_a OR sequence_a               target A                    REQUIRED
    gene_b OR sequence_b               target B                    REQUIRED

OUTCOME, either model, one of:
    pic50_a, pic50_b                   larger = more potent
    relation_a, relation_b             '=' (default), '>' or '<'
  or
    winner                             'A' or 'B'

`relation` describes the POTENCY: '>' means at least this potent. Rows are used
only when the two intervals are disjoint. Two identical exact values are a tie,
kept and entered both ways so the model returns about 0.5.

Every usable row is entered twice with the label reversed: ligands exchanged
for potency, sequences exchanged for selectivity. Do not supply both orders.

Read KFM_EXTEND_SPEC.md before relying on this. Two things it says that matter:
your added trees vote everywhere, including where you hold no data, so expect a
local gain and a small global cost; and the released operating point does
not survive extension and must be re-measured on your own held-out data.
"""
import argparse
import copy
import csv
import hashlib
import json
import os
import shutil
import sys
import time

import numpy as np

POINT, CEIL, FLOOR = 0, 1, 2          # exact, "at most", "at least"
LIG_DIMS, SEQ_DIMS = 1038, 480
ELEMENTS = ["C", "N", "O", "S", "F", "Cl", "Br", "I", "P"]


def ligand_matrix(smiles_list):
    """Morgan COUNT fingerprint r=2/1024, then the 14 descriptors, in order.

    Implemented here rather than imported, because the two released bundles
    expose different predict.py APIs and this tool must serve both.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    F = np.zeros((len(smiles_list), LIG_DIMS), dtype=np.float32)
    for i, s in enumerate(smiles_list):
        m = Chem.MolFromSmiles(s)
        fp = np.asarray(gen.GetCountFingerprintAsNumPy(m), dtype=np.float32)
        c = {e: 0 for e in ELEMENTS}
        for at in m.GetAtoms():
            if at.GetSymbol() in c:
                c[at.GetSymbol()] += 1
        F[i] = np.concatenate([fp, np.array(
            [Descriptors.MolWt(m), m.GetNumHeavyAtoms(), m.GetNumBonds(),
             Descriptors.NumRotatableBonds(m), Descriptors.RingCount(m)]
            + [c[e] for e in ELEMENTS], dtype=np.float32)])
    return F


def parseable(smiles_set):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    return {s for s in smiles_set if Chem.MolFromSmiles(s) is not None}


class Bundle:
    """Reads geometry and lookups from a released bundle, either model."""

    def __init__(self, path, load_model=True, kind=None):
        """load_model=False reads only the encoders, sequence vectors and gene map,
        which are ESM2 embeddings of public UniProt sequences. kfm_buildnew uses
        that path so no released weights are touched."""
        import joblib
        self.path = path
        self.manifest = json.load(open(os.path.join(path, "MANIFEST.json")))
        fo = str(self.manifest.get("feature_order", ""))
        self.model_file = self.manifest.get("model_file") or self._find_model()
        if load_model:
            self.model = joblib.load(os.path.join(path, self.model_file))
            w = self.model.n_features_in_
        else:
            self.model = None
            w = ({"potency": 2 * LIG_DIMS + SEQ_DIMS,
                  "selectivity": 2 * SEQ_DIMS + LIG_DIMS}[kind] if kind else
                 (2 * LIG_DIMS + SEQ_DIMS if fo.strip().startswith("[ ligand")
                  else 2 * SEQ_DIMS + LIG_DIMS))
        if fo.strip().startswith("[ ligand") or w == 2 * LIG_DIMS + SEQ_DIMS:
            self.kind, self.width = "LSL", 2 * LIG_DIMS + SEQ_DIMS
        elif fo.strip().startswith("[ sequence") or w == 2 * SEQ_DIMS + LIG_DIMS:
            self.kind, self.width = "SLS", 2 * SEQ_DIMS + LIG_DIMS
        else:
            raise SystemExit(f"cannot identify this model: feature_order={fo!r}, "
                             f"n_features_in_={w}")
        if w != self.width:
            raise SystemExit(f"manifest says {self.kind} ({self.width}) but the model "
                             f"has {w} features")
        self.vectors = np.load(os.path.join(
            path, self._first(["sequence_vectors.npz", "kinase_vectors.npz"])))["vectors"]
        j = json.load(open(os.path.join(
            path, self._first(["sequence_index.json", "kinase_index.json"]))))
        self.seq_key_to_row = j["seq_key_to_row"]
        self.gene_to_seq_key = j.get("gene_to_seq_key") or json.load(
            open(os.path.join(path, "gene_map.json")))["gene_to_seq_key"]
        self.c1 = int(np.where(self.model.classes_ == 1)[0][0]) if self.model else 1
        hp = self.manifest.get("hyperparameters") or {}
        self.leaf = hp.get("min_samples_leaf") or self.manifest.get("min_samples_leaf") or 20

    def _find_model(self):
        c = [f for f in os.listdir(self.path) if f.endswith(".joblib")]
        if len(c) != 1:
            raise SystemExit(f"cannot identify the model file in {self.path}: {c}")
        return c[0]

    def _first(self, names):
        for n in names:
            if os.path.exists(os.path.join(self.path, n)):
                return n
        raise SystemExit(f"none of {names} found in {self.path}")

    def seq_vector(self, kind, value, emb):
        if kind == "gene":
            return self.vectors[self.seq_key_to_row[self.gene_to_seq_key[value]]]
        s = emb.normalise(value) if emb else "".join(value.split()).upper()
        k = hashlib.md5(s.encode()).hexdigest()
        if k in self.seq_key_to_row:
            return self.vectors[self.seq_key_to_row[k]]
        if emb is None:
            raise SystemExit("A sequence in your file is not in the bundle; embedding "
                             "needs torch and transformers.")
        return emb.embed_sequence(s)

    def write_row(self, X, i, first, middle, last):
        if self.kind == "LSL":
            X[i, :LIG_DIMS] = first
            X[i, LIG_DIMS:LIG_DIMS + SEQ_DIMS] = middle
            X[i, LIG_DIMS + SEQ_DIMS:] = last
        else:
            X[i, :SEQ_DIMS] = first
            X[i, SEQ_DIMS:SEQ_DIMS + LIG_DIMS] = middle
            X[i, SEQ_DIMS + LIG_DIMS:] = last


def log(m):
    print(m, flush=True)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


# ── outcome logic, shared by both models ────────────────────────────────────
def _relation(v):
    v = (v or "=").strip()
    if v in ("", "=", "=="):
        return POINT
    if v in (">", ">=", "=>"):
        return FLOOR                       # at least this potent
    if v in ("<", "<=", "=<"):
        return CEIL                        # at most this potent
    raise ValueError(f"relation {v!r} not understood; use =, > or <")


def _decide(pa, ka, pb, kb):
    """(usable, a_wins, tie). The two intervals must be disjoint."""
    if ka == POINT and kb == POINT:
        return (True, True, True) if pa == pb else (True, pa > pb, False)
    if ka == POINT and kb == CEIL:
        return (pa > pb), True, False
    if ka == CEIL and kb == POINT:
        return (pb > pa), False, False
    if ka == POINT and kb == FLOOR:
        return (pa < pb), False, False
    if ka == FLOOR and kb == POINT:
        return (pb < pa), True, False
    if ka == FLOOR and kb == CEIL:
        return (pa > pb), True, False
    if ka == CEIL and kb == FLOOR:
        return (pb > pa), False, False
    return False, False, False             # two bounds the same way


MEASURE_CAP_DEFAULT = 5000


def _is_measurement_file(cols, kind):
    """One measurement per row: a ligand, a target, a value. The natural shape of
    a knowledgebase or ChEMBL extract, and what users actually have.

    Detected only when the paired columns are ABSENT, so an existing paired file
    is never reinterpreted.
    """
    if "smiles" not in cols or not ({"pic50", "pvalue", "pactivity"} & cols):
        return False
    if not ({"gene", "sequence"} & cols):
        return False
    paired = ({"smiles_a", "smiles_b"} & cols) or ({"gene_a", "gene_b",
                                                    "sequence_a", "sequence_b"} & cols)
    return not paired


def pair_measurements(rows, kind, cap=MEASURE_CAP_DEFAULT, seed=0, log=lambda s: None):
    """Turn one-measurement-per-row into the paired rows the models are fitted on.

    Potency  (LSL): group by TARGET, pair the ligands measured on it.
    Selectivity (SLS): group by LIGAND, pair the targets it was measured against.

    Same input file either way; only the grouping differs. Callability is NOT
    decided here -- the pairs are emitted with both values and both relations and
    the existing _decide() rule applies downstream, so censored readings behave
    exactly as they do for a hand-built file.
    """
    import itertools
    import random
    import statistics

    val_key = None
    for k in ("pic50", "pvalue", "pactivity"):
        if any(k in {c.lower().strip() for c in r} for r in rows[:1]):
            val_key = k
            break
    tgt_key = "gene" if any("gene" in {c.lower().strip() for c in r} for r in rows[:1]) else "sequence"

    # Collapse duplicates by MEDIAN, never by most potent: taking the best value
    # reliably selects unit errors rather than real potency.
    cell, rel = {}, {}
    bad = 0
    for r in rows:
        g = {k.lower().strip(): (v or "").strip() for k, v in r.items()}
        smi, tgt, raw = g.get("smiles"), g.get(tgt_key), g.get(val_key)
        if not (smi and tgt and raw):
            bad += 1
            continue
        try:
            v = float(raw)
        except ValueError:
            bad += 1
            continue
        cell.setdefault((tgt, smi), []).append(v)
        rel[(tgt, smi)] = g.get("relation", "=") or "="
    if bad:
        log(f"  dropped, incomplete or unparsable        {bad:,}")
    flat = {k: statistics.median(v) for k, v in cell.items()}
    dupes = sum(len(v) - 1 for v in cell.values())
    if dupes:
        log(f"  duplicate measurements merged by median  {dupes:,}")

    groups = {}
    for (tgt, smi), v in flat.items():
        key = tgt if kind == "LSL" else smi          # target for potency, ligand for selectivity
        groups.setdefault(key, []).append((smi, tgt, v))

    rng = random.Random(seed)
    out, capped = [], 0
    for key, members in groups.items():
        if len(members) < 2:
            continue
        pairs = list(itertools.combinations(range(len(members)), 2))
        if len(pairs) > cap:                          # one deep target must not become the model
            capped += 1
            rng.shuffle(pairs)
            pairs = pairs[:cap]
        for i, j in pairs:
            if rng.random() < 0.5:                    # randomise which member is A, or the
                i, j = j, i                           # forest learns to read column position
            (sa, ta, va), (sb, tb, vb) = members[i], members[j]
            # str(): the paired path downstream reads every cell as text, exactly
            # as csv.DictReader hands it over for a hand-built file.
            row = {"pic50_a": str(va), "pic50_b": str(vb),
                   "relation_a": rel[(ta, sa)], "relation_b": rel[(tb, sb)]}
            if kind == "LSL":
                row.update({"smiles_a": sa, "smiles_b": sb, tgt_key: key})
            else:
                row.update({"smiles": key,
                            tgt_key + "_a": ta, tgt_key + "_b": tb})
            out.append(row)

    # A measurement is only usable if something else shares its group: another
    # ligand on the same target for potency, another target for the same ligand
    # for selectivity. Report what fell out, because a file that looks large can
    # be almost entirely unpairable for one of the two layouts and the user has
    # no other way to see it.
    paired_by = "target" if kind == "LSL" else "ligand"
    orphan_groups = [k for k, m in groups.items() if len(m) < 2]
    orphan_rows = sum(len(groups[k]) for k in orphan_groups)
    used_rows = len(flat) - orphan_rows
    log(f"  measurements usable                      {used_rows:,} of {len(flat):,}")
    if orphan_rows:
        log(f"  UNUSED, no second measurement to pair with against the same "
            f"{'target' if kind == 'LSL' else 'ligand'}:")
        log(f"    measurements                           {orphan_rows:,} "
            f"({100.0*orphan_rows/len(flat):.1f}%)")
        log(f"    {paired_by}s with only one measurement{'':<3} {len(orphan_groups):,}")
    log(f"  {paired_by}s with 2+ measurements{'':<14} "
        f"{sum(1 for m in groups.values() if len(m) > 1):,}")
    if capped:
        log(f"  groups capped at {cap:,} pairs{'':<12} {capped:,}")
    log(f"  comparisons generated                    {len(out):,}")
    if not out:
        log(f"  Nothing pairable. For {'potency' if kind == 'LSL' else 'selectivity'} "
            f"you need at least one {paired_by} carrying two measurements.")
    return out


def read_csv(path, b, have_emb, measure_cap=MEASURE_CAP_DEFAULT, measure_seed=0):
    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8-sig")))
    if not rows:
        raise SystemExit(f"{path} has no rows")
    cols = {c.lower().strip() for c in rows[0]}
    if _is_measurement_file(cols, b.kind):
        log(f"\n{path} is one measurement per row. Pairing them "
            f"{'by target' if b.kind == 'LSL' else 'by ligand'}:")
        rows = pair_measurements(rows, b.kind, cap=measure_cap,
                                 seed=measure_seed, log=log)
        if not rows:
            raise SystemExit(
                f"{path} produced no comparisons. Every target needs at least two "
                "measured ligands (potency), or every ligand at least two measured "
                "targets (selectivity).")
        cols = {c.lower().strip() for c in rows[0]}
    if b.kind == "LSL":
        for c in ("smiles_a", "smiles_b"):
            if c not in cols:
                raise SystemExit(f"{path} needs {c!r} for a potency model")
        if not ({"gene", "sequence"} & cols):
            raise SystemExit(f"{path} needs a 'gene' or 'sequence' column")
    else:
        if "smiles" not in cols:
            raise SystemExit(f"{path} needs 'smiles' for a selectivity model")
        if not ({"gene_a", "sequence_a"} & cols) or not ({"gene_b", "sequence_b"} & cols):
            raise SystemExit(f"{path} needs gene_a/sequence_a and gene_b/sequence_b")
    by_value = {"pic50_a", "pic50_b"} <= cols
    if not by_value and "winner" not in cols:
        raise SystemExit(f"{path} needs pic50_a and pic50_b, or a 'winner' column")

    emb = None
    if cols & {"sequence", "sequence_a", "sequence_b"}:
        if not have_emb:
            raise SystemExit("This file supplies raw sequences. Embedding needs torch "
                             "and transformers:  pip install torch transformers")
        sys.path.insert(0, b.path)
        import embed as emb

    out, drop = [], {}
    def bump(k):
        drop[k] = drop.get(k, 0) + 1

    def target(g, sfx):
        if g.get("sequence" + sfx):
            return ("seq", g["sequence" + sfx])
        if g.get("gene" + sfx):
            return ("gene", g["gene" + sfx])
        return None

    for r in rows:
        g = {k.lower().strip(): (v or "").strip() for k, v in r.items()}
        if b.kind == "LSL":
            la, lb = g.get("smiles_a"), g.get("smiles_b")
            ta = tb = target(g, "")
        else:
            la = lb = g.get("smiles")
            ta, tb = target(g, "_a"), target(g, "_b")
        if not la or not lb:
            bump("blank SMILES"); continue
        if ta is None or tb is None:
            bump("no target"); continue
        if any(t[0] == "gene" and t[1] not in b.gene_to_seq_key for t in (ta, tb)):
            bump("gene not in bundle"); continue
        if by_value and g.get("pic50_a") and g.get("pic50_b"):
            try:
                pa, pb = float(g["pic50_a"]), float(g["pic50_b"])
                ka, kb = _relation(g.get("relation_a")), _relation(g.get("relation_b"))
            except ValueError as e:
                bump(f"unparseable value ({e})"); continue
            usable, a_wins, tie = _decide(pa, ka, pb, kb)
            if not usable:
                bump("readings do not order (two bounds the same way, or overlapping)")
                continue
        else:
            w = g.get("winner", "").upper()
            if w not in ("A", "B"):
                bump("winner is not A or B"); continue
            a_wins, tie = (w == "A"), False
        out.append((la, lb, ta, tb, 1 if a_wins else 0, tie))

    # A partner's file will contain structures RDKit cannot read. The bundles'
    # predict.py refuses to drop those silently, which is right when scoring and
    # wrong when importing, so they are filtered and counted here instead.
    smis = {s for r in out for s in (r[0], r[1])}
    good = parseable(smis)
    if smis - good:
        keep = [r for r in out if r[0] in good and r[1] in good]
        drop["SMILES RDKit could not parse"] = len(out) - len(keep)
        out = keep
    return (out, sorted({s for r in out for s in (r[0], r[1])}),
            list({t for r in out for t in (r[2], r[3])}), drop, emb)


def build(usable, smis, tgts, b, emb):
    L = ligand_matrix(smis); li = {s: i for i, s in enumerate(smis)}
    V = {t: b.seq_vector(t[0], t[1], emb) for t in tgts}
    n = len(usable)
    X = np.empty((2 * n, b.width), dtype=np.float32)
    y = np.empty(2 * n, dtype=np.int8)
    for i, (la, lb, ta, tb, lab, tie) in enumerate(usable):
        first, mid, last = ((L[li[la]], V[ta], L[li[lb]]) if b.kind == "LSL"
                            else (V[ta], L[li[la]], V[tb]))
        b.write_row(X, i, first, mid, last)
        b.write_row(X, n + i, last, mid, first)
        y[i] = lab
        y[n + i] = lab if tie else 1 - lab
    return X, y


def score(model, c1, usable, smis, tgts, b, emb):
    """Both orders averaged, which is how the shipped predictors behave."""
    rows = [r for r in usable if not r[5]]
    if not rows:
        return float("nan"), 0
    L = ligand_matrix(smis); li = {s: i for i, s in enumerate(smis)}
    V = {t: b.seq_vector(t[0], t[1], emb) for t in tgts}
    n = len(rows)
    X = np.empty((2 * n, b.width), dtype=np.float32)
    for i, (la, lb, ta, tb, lab, _t) in enumerate(rows):
        first, mid, last = ((L[li[la]], V[ta], L[li[lb]]) if b.kind == "LSL"
                            else (V[ta], L[li[la]], V[tb]))
        b.write_row(X, i, first, mid, last)
        b.write_row(X, n + i, last, mid, first)
    p = model.predict_proba(X)[:, c1]
    avg = 0.5 * (p[:n] + (1.0 - p[n:]))
    y = np.array([r[4] for r in rows])
    return float(((avg > 0.5).astype(int) == y).mean()), n


def design(rows, smis, tgts, b, emb):
    """Both orders stacked: first n rows forward, next n reversed."""
    L = ligand_matrix(smis); li = {s: i for i, s in enumerate(smis)}
    V = {t: b.seq_vector(t[0], t[1], emb) for t in tgts}
    n = len(rows)
    X = np.empty((2 * n, b.width), dtype=np.float32)
    for i, (la, lb, ta, tb, lab, _t) in enumerate(rows):
        first, mid, last = ((L[li[la]], V[ta], L[li[lb]]) if b.kind == "LSL"
                            else (V[ta], L[li[la]], V[tb]))
        b.write_row(X, i, first, mid, last)
        b.write_row(X, n + i, last, mid, first)
    return X, np.array([r[4] for r in rows]), n


def sweep_curve(base, theirs, c1, X, y, n, points):
    """Accuracy of base+N of theirs, for every N, from ONE pass over the trees.

    A forest's vote at N trees is the mean over its first N, and the first N
    trees of a larger fit are identical to an N-tree fit at the same seed
    (verified). So one fit and one cumulative pass gives the whole curve.
    """
    pb = base.predict_proba(X)[:, c1]
    cum = np.zeros(X.shape[0], dtype=np.float64)
    out, nb = {}, base.n_estimators
    nxt = sorted(points)
    for k, t in enumerate(theirs.estimators_, 1):
        cum += t.predict_proba(X)[:, c1]
        if nxt and k == nxt[0]:
            pm = (nb * pb + cum) / (nb + k)          # merged vote
            avg = 0.5 * (pm[:n] + (1.0 - pm[n:]))    # both orders
            out[k] = float(((avg > 0.5).astype(int) == y).mean())
            nxt.pop(0)
    avg0 = 0.5 * (pb[:n] + (1.0 - pb[n:]))
    return float(((avg0 > 0.5).astype(int) == y).mean()), out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Add your data to a released KFM model.")
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--base", default=None,
                    help="bundle to extend. Defaults to the installed bundle for "
                         "--model. May itself be an extended bundle.")
    ap.add_argument("--model", choices=["potency", "selectivity"], default=None,
                    help="which installed model to extend; shorthand for --base")
    ap.add_argument("--trees", type=int, default=20)
    ap.add_argument("--holdout", default=None)
    ap.add_argument("--pairs-per-group", type=int, default=MEASURE_CAP_DEFAULT,
                    metavar="N",
                    help="only when --data is one measurement per row: the most\n"
                         "comparisons to draw from any single group, so one deep\n"
                         "target cannot become the model. Default %(default)s.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-samples-leaf", type=int, default=None,
                    help="defaults to the base model's own value")
    ap.add_argument("--name", default=None, help="name for the merged model")
    ap.add_argument("--label", default=None, help="a note recorded in the extension log")
    ap.add_argument("--sweep", action="store_true",
                    help="measure accuracy against tree count on YOUR holdout, then "
                         "build at the chosen count. Requires --holdout. One fit; the "
                         "curve comes from evaluating prefixes of the same trees.")
    ap.add_argument("--sweep-points", default="5,10,20,40,80,150,300")
    ap.add_argument("--holdout-breadth", default=None,
                    help="comparisons on targets your data does NOT cover, used to "
                         "measure what extension costs. Defaults to the reference set "
                         "shipped in the bundle, which spans the whole panel, so you "
                         "normally do not need to supply one. Pass 'none' to skip it.")
    ap.add_argument("--max-breadth-loss", type=float, default=0.005,
                    help="the recommendation is the largest tree count whose loss on "
                         "--holdout-breadth stays within this")
    ap.add_argument("--allow-small", action="store_true",
                    help="permit fewer than 1,000 usable rows. For demonstration only; "
                         "a handful of rows produces trees that are noise and they vote "
                         "on every prediction.")
    a = ap.parse_args(argv)

    if a.base is None:
        from . import bundles as _b
        want = a.model or "potency"
        if not _b.is_present(want):
            raise SystemExit(
                f"The {want} model is not installed, so there is nothing to extend.\n"
                f"  Get it with:  ./install.sh {want}\n"
                f"  Or point at a bundle you already have:  --base /path/to/bundle")
        a.base = _b.bundle_dir(want)

    import joblib, sklearn
    from sklearn.ensemble import RandomForestClassifier

    b = Bundle(a.base)
    log(f"base bundle : {a.base}")
    log(f"model       : {b.manifest.get('name', b.model_file)}")
    log(f"kind        : {b.kind}   {b.manifest.get('feature_order','')}")
    log(f"base forest : {b.model.n_estimators} trees, {b.model.n_features_in_} features")

    env = b.manifest.get("environment") or {}
    pinned = env.get("scikit_learn") or env.get("sklearn")
    if pinned and str(pinned).split(".")[:2] != sklearn.__version__.split(".")[:2]:
        raise SystemExit(f"scikit-learn {sklearn.__version__} but this bundle was built "
                         f"with {pinned}. Trees from different minor versions must not "
                         f"be merged.")
    try:
        import torch, transformers  # noqa: F401
        have_emb = True
    except ImportError:
        have_emb = False

    usable, smis, tgts, drop, emb = read_csv(a.data, b, have_emb, a.pairs_per_group, a.seed)
    total = len(usable) + sum(drop.values())
    log(f"\nyour data   : {a.data}")
    log(f"  rows read                 {total:,}")
    log(f"  usable comparisons        {len(usable):,}")
    for k, v in sorted(drop.items(), key=lambda x: -x[1]):
        log(f"  dropped, {k:<52} {v:,}")
    ties = sum(1 for r in usable if r[5])
    if ties:
        log(f"  ties kept as 0.5          {ties:,}")
    log(f"  distinct targets          {len(tgts):,}")
    log(f"  distinct ligands          {len(smis):,}")
    if not usable:
        raise SystemExit("No usable rows. See the drop reasons above.")
    if len(usable) < 1000 and not a.allow_small:
        raise SystemExit(f"Only {len(usable):,} usable comparisons. Fewer than 1,000 "
                         "cannot support a useful tree; the result would be noise voting "
                         "on every prediction. Pass --allow-small to override for a "
                         "format demonstration.")
    if len(usable) < 1000:
        log("  WARNING: --allow-small. These trees are fitted on too little data to be "
            "meaningful and will degrade the model. Demonstration only.")

    points = sorted({int(x) for x in a.sweep_points.split(",") if x.strip()})
    if a.sweep and a.holdout_breadth is None:
        _ref_name = ("breadth_reference_potency.csv" if b.kind == "LSL"
                     else "breadth_reference_selectivity.csv")
        ref = os.path.join(a.base, _ref_name)
        if not os.path.isfile(ref):
            # Bundles downloaded before the breadth sets shipped do not carry
            # them. The package copy is the same ChEMBL-derived file.
            _pkg = os.path.join(os.path.dirname(os.path.dirname(
                os.path.abspath(__file__))), "examples", _ref_name)
            if os.path.isfile(_pkg):
                ref = _pkg
        if os.path.exists(ref):
            a.holdout_breadth = ref
            log(f"breadth check: using the bundle's reference set, {os.path.basename(ref)}")
        else:
            log("breadth check: no reference set in this bundle and none supplied")
    elif a.sweep and str(a.holdout_breadth).lower() == "none":
        a.holdout_breadth = None
        log("breadth check: skipped at your request")
    if a.sweep and not a.holdout:
        raise SystemExit("--sweep needs --holdout: the whole point is to measure the "
                         "trade-off on your data rather than guess it.")
    n_fit = max(points) if a.sweep else a.trees

    X, y = build(usable, smis, tgts, b, emb)
    log(f"\ndesign matrix {X.shape[0]:,} x {X.shape[1]:,} "
        f"({X.nbytes / 2**30:.2f} GB) after the swap")
    leaf = a.min_samples_leaf if a.min_samples_leaf is not None else b.leaf
    t0 = time.time()
    yours = RandomForestClassifier(n_estimators=n_fit, min_samples_leaf=leaf,
                                   max_features="sqrt", n_jobs=-1,
                                   random_state=a.seed).fit(X, y)
    log(f"fitted {n_fit} trees (min_samples_leaf={leaf}) in {time.time()-t0:.0f}s")
    del X, y

    chosen = a.trees
    if a.sweep:
        log("\n" + "=" * 74)
        log("SWEEP: accuracy against tree count, measured on your own holdout")
        log("One fit; the curve is prefixes of the same trees, which is exact.")
        log("=" * 74)
        h, hs, ht, _hd, hemb = read_csv(a.holdout, b, have_emb, a.pairs_per_group, a.seed)
        hemb = hemb or emb
        Xh, yh, nh = design([r for r in h if not r[5]], hs, ht, b, hemb)
        base_acc, curve = sweep_curve(b.model, yours, b.c1, Xh, yh, nh, points)
        del Xh
        bcurve, base_b = None, None
        if a.holdout_breadth:
            hb, hbs, hbt, _x, hbe = read_csv(a.holdout_breadth, b, have_emb, a.pairs_per_group, a.seed)
            hbe = hbe or emb
            Xb, yb, nbn = design([r for r in hb if not r[5]], hbs, hbt, b, hbe)
            base_b, bcurve = sweep_curve(b.model, yours, b.c1, Xb, yb, nbn, points)
            del Xb
        log(f"\nyour holdout      : {nh:,} comparisons, released model {base_acc:.4f}")
        if bcurve is not None:
            log(f"breadth holdout   : targets your data does not cover, "
                f"released model {base_b:.4f}")
        log("")
        hdr = f"{'--trees':>8} {'share':>7} {'yours':>9} {'gain':>8}"
        if bcurve is not None:
            hdr += f" {'breadth':>9} {'cost':>8}"
        log(hdr)
        rec = None
        for k in points:
            share = k / (b.model.n_estimators + k)
            line = (f"{k:>8} {share:>6.1%} {curve[k]:>9.4f} "
                    f"{curve[k]-base_acc:>+8.4f}")
            if bcurve is not None:
                loss = base_b - bcurve[k]
                line += f" {bcurve[k]:>9.4f} {-loss:>+8.4f}"
                if loss <= a.max_breadth_loss:
                    rec = k
            log(line)
        if bcurve is None:
            log("\nNo recommendation: without --holdout-breadth this sweep measured "
                "only what you gain, never what you give up on the rest of the panel. "
                "Supply one, or choose deliberately.")
            chosen = a.trees
            log(f"Building at --trees {chosen}.")
        else:
            chosen = rec or min(points)
            log(f"\nRecommended: {chosen} trees, the largest whose breadth loss stays "
                f"within {a.max_breadth_loss:.3f}. Building at {chosen}.")
        log("=" * 74)
        yours.estimators_ = yours.estimators_[:chosen]
        yours.n_estimators = chosen

    if b.model.n_features_in_ != yours.n_features_in_:
        raise SystemExit("feature width mismatch; refusing to merge")
    if list(b.model.classes_) != list(yours.classes_):
        raise SystemExit("class order mismatch; refusing to merge")

    merged = copy.copy(b.model)
    merged.estimators_ = list(b.model.estimators_) + list(yours.estimators_)
    merged.n_estimators = len(merged.estimators_)
    share = chosen / merged.n_estimators
    log(f"\nmerged forest: {merged.n_estimators} trees, your share {share:.1%}")

    os.makedirs(a.out, exist_ok=True)
    for f in sorted(os.listdir(a.base)):
        s = os.path.join(a.base, f)
        if os.path.isfile(s) and not f.endswith(".joblib") and f != "MANIFEST.json":
            shutil.copy2(s, os.path.join(a.out, f))
    shutil.copy2(os.path.abspath(__file__), os.path.join(a.out, "kfm_extend.py"))
    base_name = str(b.manifest.get("name", os.path.splitext(b.model_file)[0]))
    new_name = a.name or (base_name if base_name.endswith("_extended")
                          else base_name.replace(" ", "_") + "_extended")
    # Keep the BASE model's filename so each bundle's own prediction utility
    # loads the extended model with no code change. The new name is carried in
    # the manifest and in the directory name, not in the filename.
    dest = os.path.join(a.out, b.model_file)
    joblib.dump(merged, dest, compress=3)
    log(f"wrote {dest} ({os.path.getsize(dest)/2**30:.2f} GB)")
    log(f"  filename kept as {b.model_file} so the bundle's own predict.py loads it")

    man = dict(b.manifest)
    man["name"] = new_name
    man["model_file"] = b.model_file
    man["derived_from"] = {"name": base_name,
                           "model_sha256": b.manifest.get("model_sha256")}
    man.setdefault("extensions", []).append({
        "when": time.strftime("%Y-%m-%dT%H:%M:%S"), "label": a.label,
        "model_kind": b.kind, "trees_added": chosen,
        "chosen_by_sweep": bool(a.sweep), "rows_read": total,
        "usable_comparisons": len(usable), "ties_kept": ties, "dropped": drop,
        "distinct_targets": len(tgts), "distinct_ligands": len(smis),
        "seed": a.seed, "min_samples_leaf": leaf,
        "n_estimators_after": merged.n_estimators,
        "contributor_vote_share": round(share, 4),
        "base_model_sha256": b.manifest.get("model_sha256"),
        "scikit_learn": sklearn.__version__,
        "allow_small_override": bool(len(usable) < 1000),
    })
    man["model_sha256"] = _sha256(dest)
    WITHDRAWN = ("This model has been extended with data we did not measure. The "
                 "released figures do not describe it. Measure it on your own "
                 "held-out comparisons before quoting any number, and re-measure "
                 "the prediction-strength operating point on your own data.")
    # Withdraw any block that reports performance. The two bundles name these
    # differently, so match on the key name rather than assuming one schema.
    hits = [k for k in list(man)
            if any(t in k.lower() for t in ("performance", "accuracy", "strength",
                                            "metric", "validation", "result",
                                            "holdout", "benchmark"))]
    for key in hits:
        man[key] = {"WITHDRAWN": WITHDRAWN, "released_figures_for_reference": man[key]}
    man["EXTENDED_MODEL_WARNING"] = WITHDRAWN
    if hits:
        log(f"  withdrew {len(hits)} performance block(s): {', '.join(hits)}")
    json.dump(man, open(os.path.join(a.out, "MANIFEST.json"), "w"), indent=2)
    log("wrote MANIFEST.json with the extension log; released figures withdrawn")

    refs = os.path.join(a.out, "reference_predictions.json")
    if os.path.exists(refs):
        os.replace(refs, refs + ".base")
        log("reference_predictions.json renamed to .base: added trees change every "
            "prediction, so the released reference values no longer apply")

    if a.holdout:
        log(f"\nholdout: {a.holdout}")
        h, hs, ht, hd, hemb = read_csv(a.holdout, b, have_emb, a.pairs_per_group, a.seed)
        hemb = hemb or emb
        before, n = score(b.model, b.c1, h, hs, ht, b, hemb)
        c1m = int(np.where(merged.classes_ == 1)[0][0])
        after, _ = score(merged, c1m, h, hs, ht, b, hemb)
        log(f"  {n:,} scorable comparisons")
        log(f"  released model {before:.4f}")
        log(f"  merged model   {after:.4f}   {after-before:+.4f}")
        if after < before:
            log("  NOTE: your trees vote on every prediction, including targets you hold "
                "no data for. Try fewer trees, or broader data.")

    log(f"\nDone. New model {new_name} written to {a.out}.")
    log(f"The base bundle at {a.base} was opened read-only and is unchanged.")
    log("Your data was never written to the bundle; only trees fitted from it.")


if __name__ == "__main__":
    main()
