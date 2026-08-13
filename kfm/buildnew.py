"""kfm buildnew -- fit a model on YOUR data alone, in the KFM input format.

This trains from scratch. No KFM weights are loaded, merged, or consulted. What
you get is your model, fitted on your comparisons, with the same feature layout
and the same interface as the released models, so it drops into the same
prediction utility and the same downstream code.

    python kfm_buildnew.py --layout potency     --data mine.csv --out ./my_model
    python kfm_buildnew.py --layout selectivity --data mine.csv --out ./my_model

CSV is the same format kfm_extend uses.

POTENCY (two ligands, one target)
    smiles_a, smiles_b                 the two ligands             REQUIRED
    gene  OR  sequence                 the target                  REQUIRED

SELECTIVITY (one ligand, two targets)
    smiles                             the ligand                  REQUIRED
    gene_a OR sequence_a               target A                    REQUIRED
    gene_b OR sequence_b               target B                    REQUIRED

OUTCOME, either layout, one of:
    pic50_a, pic50_b                   larger = more potent
    relation_a, relation_b             '=' (default), '>' or '<'
  or
    winner                             'A' or 'B'

`relation` describes the POTENCY: '>' means at least this potent. That is the
opposite of ChEMBL and most assay exports, where the qualifier sits on the
concentration, so IC50 > 1000 nM becomes pic50 6.0 with relation '<'.

WHY YOU MIGHT WANT THIS INSTEAD OF kfm extend
---------------------------------------------
Measured on a contributor holding 119,660 comparisons across five targets:

    model                              their 5 targets   10 other targets
    released KFM, 300 trees                     0.7027             0.7346
    kfm extend, +20 of their trees              0.7131             0.7326
    kfm buildnew, 300 trees, their data alone   0.7765             0.6192

Above roughly 10,000 comparisons, a model fitted on your own data beats the
released model on the targets that data covers, and beats an extended model too.
What you give up is everywhere else: a from-scratch model does not improve off
your own targets no matter how much of your own data you add. Coverage is what
the released model provides; depth on your targets is what your data provides.

Use kfm buildnew when you have deep data and care about those targets.
Use kfm extend when you want to keep coverage across the whole panel.
Run both and compare on your own holdout, which --holdout makes cheap.

SEQUENCE VECTORS
----------------
--recipe points at a released bundle and reads ONLY its sequence vectors, gene
map and encoders, never its weights. Those vectors are ESM2 embeddings of public
UniProt sequences and are recomputable from public inputs. If you would rather
not use a bundle at all, supply raw sequences in your CSV and pass
--recipe with a bundle for the ESM2 code, or embed them yourself with embed.py.
"""
import argparse
import importlib.util
import json
import os
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
try:
    from . import extend as kx
except ImportError:                                          # standalone in a bundle
    _spec = importlib.util.spec_from_file_location(
        "kfm_extend", os.path.join(HERE, "kfm_extend.py"))
    kx = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(kx)

DEFAULT_LEAF = {"potency": 20, "selectivity": 8}
KIND = {"potency": "LSL", "selectivity": "SLS"}
SUPPORT = ["predict.py", "embed.py", "sequence_vectors.npz", "kinase_vectors.npz",
           "sequence_index.json", "kinase_index.json", "gene_map.json",
           "targets.csv", "targets.json", "requirements.txt", "embed_reference.json",
           "requirements-embed.txt"]


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Fit a KFM-format model on your own data alone.")
    ap.add_argument("--layout", required=True, choices=["potency", "selectivity"])
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--recipe", default=None,
                    help="a released bundle to read sequence vectors and encoders "
                         "from. Its weights are never loaded.")
    ap.add_argument("--trees", type=int, default=300)
    ap.add_argument("--min-samples-leaf", type=int, default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--holdout", default=None)
    ap.add_argument("--compare-to", default=None,
                    help="a released bundle to score on the same holdout, so you can "
                         "see what your model gains and what it gives up")
    ap.add_argument("--name", default=None)
    ap.add_argument("--allow-small", action="store_true")
    a = ap.parse_args(argv)

    import joblib, sklearn
    from sklearn.ensemble import RandomForestClassifier

    recipe = a.recipe
    if recipe is None:
        from . import bundles as _b
        if not _b.is_present(a.layout):
            raise SystemExit(
                f"No {a.layout} bundle installed to read the recipe from.\n"
                f"  Get it with:  ./install.sh {a.layout}\n"
                f"  Or point at one you already have:  --recipe /path/to/bundle\n"
                f"Only the sequence vectors and encoders are read; the weights are not.")
        recipe = _b.bundle_dir(a.layout)
    if not os.path.isdir(recipe):
        raise SystemExit(f"--recipe not found: {recipe}")

    b = kx.Bundle(recipe, load_model=False, kind=a.layout)
    b.kind = KIND[a.layout]
    b.width = {"LSL": 2 * kx.LIG_DIMS + kx.SEQ_DIMS,
               "SLS": 2 * kx.SEQ_DIMS + kx.LIG_DIMS}[b.kind]
    kx.log(f"layout      : {a.layout}  ({b.kind}, {b.width} features)")
    kx.log(f"recipe from : {recipe}")
    kx.log("             sequence vectors and encoders only; no weights are loaded")

    try:
        import torch, transformers  # noqa: F401
        have_emb = True
    except ImportError:
        have_emb = False

    usable, smis, tgts, drop, emb = kx.read_csv(a.data, b, have_emb)
    total = len(usable) + sum(drop.values())
    kx.log(f"\nyour data   : {a.data}")
    kx.log(f"  rows read                 {total:,}")
    kx.log(f"  usable comparisons        {len(usable):,}")
    for k, v in sorted(drop.items(), key=lambda x: -x[1]):
        kx.log(f"  dropped, {k:<52} {v:,}")
    ties = sum(1 for r in usable if r[5])
    if ties:
        kx.log(f"  ties kept as 0.5          {ties:,}")
    kx.log(f"  distinct targets          {len(tgts):,}")
    kx.log(f"  distinct ligands          {len(smis):,}")
    if not usable:
        raise SystemExit("No usable rows. See the drop reasons above.")
    if len(usable) < 1000 and not a.allow_small:
        raise SystemExit(f"Only {len(usable):,} usable comparisons. A model fitted on "
                         "fewer than 1,000 is not worth running. Pass --allow-small to "
                         "override for a format demonstration.")
    if len(usable) < 10000:
        kx.log(f"\n  NOTE: {len(usable):,} comparisons. Below roughly 10,000 the released "
               "model is likely to beat this one even on your own targets. Pass "
               "--compare-to to check rather than assume.")

    X, y = kx.build(usable, smis, tgts, b, emb)
    kx.log(f"\ndesign matrix {X.shape[0]:,} x {X.shape[1]:,} "
           f"({X.nbytes / 2**30:.2f} GB) after the swap")
    leaf = a.min_samples_leaf if a.min_samples_leaf is not None else DEFAULT_LEAF[a.layout]
    t0 = time.time()
    rf = RandomForestClassifier(n_estimators=a.trees, min_samples_leaf=leaf,
                                max_features="sqrt", n_jobs=-1,
                                random_state=a.seed).fit(X, y)
    kx.log(f"fitted {a.trees} trees (min_samples_leaf={leaf}) in {time.time()-t0:.0f}s")
    del X, y

    os.makedirs(a.out, exist_ok=True)
    for f in SUPPORT:
        src = os.path.join(recipe, f)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(a.out, f))
    shutil.copy2(os.path.abspath(__file__), os.path.join(a.out, "kfm_buildnew.py"))
    shutil.copy2(kx.__file__, os.path.join(a.out, "kfm_extend.py"))
    name = a.name or f"my_{a.layout}_model"
    dest = os.path.join(a.out, b.model_file)     # the base filename, so predict.py loads it
    joblib.dump(rf, dest, compress=3)
    kx.log(f"\nwrote {dest} ({os.path.getsize(dest)/2**30:.2f} GB)")
    kx.log(f"  filename kept as {b.model_file} so the bundle's predict.py loads it")

    json.dump({
        "name": name,
        "model_file": b.model_file,
        "built_by": "kfm_buildnew",
        "PROVENANCE": ("Fitted from scratch on the contributor's own data. No KFM "
                       "weights were loaded, merged or consulted. The sequence vectors "
                       "and encoders were read from a released bundle; those are ESM2 "
                       "embeddings of public UniProt sequences and are recomputable "
                       "from public inputs."),
        "layout": a.layout,
        "feature_order": (b.manifest.get("feature_order")
                          if b.manifest.get("feature_order") else None),
        "total_dims": b.width,
        "encoders": {"ligand": "Morgan count fingerprint radius 2, 1024 bits, then 14 "
                               "descriptors: MolWt, HeavyAtomCount, NumBonds, "
                               "NumRotatableBonds, RingCount, C, N, O, S, F, Cl, Br, I, P",
                     "sequence": "ESM2-t12-35M, mean over residues, first 480 dimensions"},
        "output": "two values summing to 1; class 1 means A wins",
        "augmentation": "every comparison entered twice with the label reversed",
        "training": {"source_file": os.path.basename(a.data),
                     "rows_read": total, "usable_comparisons": len(usable),
                     "ties_kept": ties, "dropped": drop,
                     "distinct_targets": len(tgts), "distinct_ligands": len(smis),
                     "rows_after_swap": 2 * len(usable)},
        "hyperparameters": {"n_estimators": a.trees, "min_samples_leaf": leaf,
                            "max_features": "sqrt", "random_state": a.seed},
        "model_sha256": kx._sha256(dest),
        "environment": {"python": sys.version.split()[0],
                        "scikit_learn": sklearn.__version__},
        "MEASURED_PERFORMANCE": ("None. This model has never been evaluated by us. "
                                 "Measure it on your own held-out comparisons before "
                                 "quoting any number. Nothing about the released "
                                 "models' accuracy applies to it."),
        "COVERAGE_WARNING": ("A model fitted only on your data is weak on targets your "
                             "data does not cover, and more of your own data does not "
                             "fix that. Measured on a five-target contributor set: "
                             "0.6027 on unrelated targets at 2,000 comparisons and "
                             "0.6192 at 119,660, against 0.7346 for the released model."),
    }, open(os.path.join(a.out, "MANIFEST.json"), "w"), indent=2)
    kx.log("wrote MANIFEST.json")

    if a.holdout:
        kx.log(f"\nholdout: {a.holdout}")
        h, hs, ht, hd, hemb = kx.read_csv(a.holdout, b, have_emb)
        hemb = hemb or emb
        c1 = int(np.where(rf.classes_ == 1)[0][0])
        mine, n = kx.score(rf, c1, h, hs, ht, b, hemb)
        kx.log(f"  {n:,} scorable comparisons")
        kx.log(f"  your model      {mine:.4f}")
        if a.compare_to:
            rb = kx.Bundle(a.compare_to)
            theirs, _ = kx.score(rb.model, rb.c1, h, hs, ht, rb, hemb)
            kx.log(f"  released model  {theirs:.4f}   yours {mine-theirs:+.4f}")
            kx.log("  Remember this holdout covers the targets your data covers. It "
                   "says nothing about the rest of the panel, where the released model "
                   "is the stronger one.")

    kx.log(f"\nDone. {name} written to {a.out}.")
    kx.log("No KFM weights are in this model. It is yours.")


if __name__ == "__main__":
    main()
