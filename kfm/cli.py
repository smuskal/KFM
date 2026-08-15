"""Command line for the Kinase Foundation Model v2 models.

    kfm v1           --target ABL1 --ligand "<SMILES> <name>"   # predicted pIC50
    kfm potency      --target ABL1 --ligand "<SMILES> <name>" --ligand "<SMILES> <name>"
    kfm selectivity  --target ABL1 --target GSK3B --ligand "<SMILES> <name>"
    kfm targets potency
    kfm download potency --accept-licence

Both models answer a COMPARISON and return a probability. Neither returns a
potency, an affinity, or a statement that a compound is active. Read the
ordering, and use the confidence to decide which parts of it to act on.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

from . import bundles, io as kio

CUTOFF = 0.70

# Measured accuracy per confidence band, from each model's own report. These are
# NOT interchangeable: the two were evaluated on different test sets, and the
# same confidence means something different on each.
BANDS = {
    "potency":     [(0.90, "87.2%"), (0.80, "85.0%"), (0.70, "83.2%"),
                    (0.60, "74.6%"), (0.00, "58.4%")],
    "selectivity": [(0.90, "99.3%"), (0.80, "96.4%"), (0.70, "88.4%"),
                    (0.60, "74.4%"), (0.00, "57.9%")],
}


def band_accuracy(model: str, conf: float) -> str:
    for lo, acc in BANDS[model]:
        if conf >= lo:
            return acc
    return BANDS[model][-1][1]


def _emit(args, payload, text: str) -> None:
    if args.json:
        json.dump(payload, sys.stdout, indent=1)
        sys.stdout.write("\n")
    else:
        print(text)


# ---------------------------------------------------------------------------
# potency: N compounds against ONE kinase
# ---------------------------------------------------------------------------
def cmd_potency(args) -> int:
    # Arguments are checked BEFORE the model is loaded, so a typo answers
    # immediately and says what is wrong -- rather than spending a multi-gigabyte
    # load first, or telling someone to download a model when the real problem is
    # that they passed one compound.
    ligands = kio.read_smiles(args.ligand)
    if len(ligands) < 2:
        raise SystemExit("This model compares compounds, so it needs at least two "
                         "(--ligand twice, or --ligand @file.smi).")
    if not args.target or len(args.target) != 1:
        raise SystemExit("Give exactly one --target for potency. To compare two "
                         "kinases, use the selectivity command.")

    P = bundles.load_predict("potency")
    m = P.load(bundles.bundle_dir("potency"))

    target = args.target[0]
    kw = {}
    if target.startswith("@") or len(target) > 30:
        kw["sequence"] = kio.read_sequence(target)
    else:
        kw["gene"] = target

    smis = [l["smiles"] for l in ligands]
    # rank() scores every unordered pair in BOTH ligand orders and averages
    # them. That averaging is the model's own and is not optional: the two
    # orders disagree on the winner for 10.07% of pairs, on which single-order
    # accuracy is exactly 0.5000 -- a coin flip decided by argument order.
    ranked = P.rank(m, smis, **kw)

    by_smiles = {l["smiles"]: l for l in ligands}
    rows, payload = [], []
    for i, (smi, score, strength) in enumerate(ranked):
        name = by_smiles.get(smi, {}).get("name") or f"#{i + 1}"
        # Column order and content deliberately match the website's results
        # table: #, Compound, SMILES, Score, Confidence. One run described two
        # ways should not look like two different results. Band accuracy is
        # still in --json, and the bands themselves are printed below the table.
        rows.append([i + 1, name, smi, f"{score:.3f}", f"{strength:.2f}"])
        payload.append({"rank": i + 1, "name": name, "smiles": smi,
                        "score": float(score), "confidence": float(strength),
                        "band_accuracy": band_accuracy("potency", strength)})

    text = (f"Potency ranking against {target}\n\n"
            + kio.fmt_table(
                ["#", "Compound", "SMILES", "Score", "Confidence"], rows,
                ["<", "<", "<", ">", ">"])
            + "\n\nSCORE orders the table: how often the model expects this compound to"
              "\nbe the more potent one, averaged over every comparison against the"
              "\nothers you supplied. 1.00 beat them all, 0.00 lost to them all, 0.50"
              "\nbroke even. Change the other compounds and it moves - it is a position"
              "\nwithin your set, not a property of the molecule."
              "\n\nCONFIDENCE is how decisive those comparisons were, 0.50 being a coin"
              "\nflip. It does NOT say the compound is good: one the model confidently"
              "\nranks LAST also has a high confidence. Run with --json for the"
              "\nmeasured accuracy of each confidence band on held-out ChEMBL data."
              "\n\nThis model compares compounds against ONE kinase. For selectivity"
              "\nbetween kinases, use: kfm selectivity")
    _emit(args, {"model": "LigASeqLigB_v2_potency", "target": target,
                 "ranking": payload}, text)
    return 0


# ---------------------------------------------------------------------------
# selectivity: compounds across 2..5 kinases
# ---------------------------------------------------------------------------
def cmd_selectivity(args) -> int:
    # Arguments before the model, for the same reason as potency above.
    ligands = kio.read_smiles(args.ligand)
    if not ligands:
        raise SystemExit("Give at least one --ligand.")
    targets = list(dict.fromkeys(args.target or []))
    if len(targets) < 2:
        raise SystemExit("This model compares kinases, so it needs at least two "
                         "(--target ABL1 --target GSK3B).")
    if len(targets) > 5:
        raise SystemExit(
            "At most 5 targets. Every pair is scored, so 5 kinases is already 10 "
            "comparisons per compound; beyond that the output stops being readable.")

    P = bundles.load_predict("selectivity")
    m = P.SeqALigSeqB(bundles.bundle_dir("selectivity"))

    import itertools
    pairs = list(itertools.combinations(targets, 2))

    # One predict_proba over every (compound, pair) row. Ligand featurisation is
    # computed once per compound and reused across pairs.
    triples, index = [], []
    for li, lig in enumerate(ligands):
        for a, b in pairs:
            triples.append((a, lig["smiles"], b))
            index.append((li, a, b))
    probs = m.compare_many_genes(triples)

    p = {}
    for (li, a, b), row in zip(index, probs):
        p[(li, a, b)] = float(row[1])     # index 1 = P(a prefers the ligand)
        p[(li, b, a)] = float(row[0])

    rows, payload = [], []
    for li, lig in enumerate(ligands):
        # A kinase's score is the mean probability it beats each of the others.
        # Across a row these average to exactly 0.5 by construction, so 0.5 means
        # "no preference" -- it is a midpoint, not a low value.
        scores = {t: sum(p[(li, t, o)] for o in targets if o != t) / (len(targets) - 1)
                  for t in targets}
        top = max(scores, key=scores.get)
        # Same columns as the website: Compound, SMILES, then one column per
        # kinase. The page carries "prefers" in the heat-map colour rather than
        # a column, and the winner is simply the largest number in the row.
        rows.append([lig.get("name") or f"#{li + 1}", lig["smiles"]]
                    + [f"{scores[t]:.2f}" for t in targets])
        payload.append({"name": lig.get("name"), "smiles": lig["smiles"],
                        "scores": scores, "prefers": top,
                        "confidence": scores[top],
                        "band_accuracy": band_accuracy("selectivity", scores[top]),
                        "pairs": {f"{a}|{b}": p[(li, a, b)] for a, b in pairs}})

    # Most decisive compound first, matching the web page and the emailed report.
    order = sorted(range(len(payload)), key=lambda i: -payload[i]["confidence"])
    rows = [rows[i] for i in order]
    payload = [payload[i] for i in order]

    text = (f"Selectivity across {', '.join(targets)}\n"
            f"{len(pairs)} comparison{'s' if len(pairs) != 1 else ''} per compound\n\n"
            + kio.fmt_table(
                ["Compound", "SMILES"] + targets, rows,
                ["<", "<"] + [">"] * len(targets))
            + "\n\nEach number is how strongly that compound prefers that kinase over"
              "\nthe others listed. 0.50 means NO PREFERENCE - the signal is distance"
              "\nfrom 0.50 in either direction. The scores are relative to the kinases"
              "\nyou chose: adding another moves them all."
              "\n\nA PREFERRED SIDE IS NOT ACTIVITY. The probabilities sum to 1, so one"
              "\nkinase always wins - including for a compound that binds none of them."
              "\n\nThe highest number in a row is the kinase that compound prefers."
              "\nRun with --json for the per-pair probabilities and the measured"
              "\naccuracy of each confidence band.")
    _emit(args, {"model": "SeqALigSeqB_v2_selectivity", "targets": targets,
                 "results": payload}, text)
    return 0



# ---------------------------------------------------------------------------
# v1: predicted pIC50 for compounds against ONE kinase
#
# This is a REGRESSION, not a comparison, which is the whole difference from the
# v2 models. It returns a number on the pIC50 scale - and that number is only
# meaningful WITHIN one kinase. The v1 scorer puts each target on its own scale,
# so subtracting two targets' scores reports scale differences as if they were
# selectivity, which is exactly the mistake the v2 selectivity model exists to
# avoid. Per-target accuracy ranges from 0.09 to 0.83 across the panel.
# ---------------------------------------------------------------------------
def cmd_v1(args) -> int:
    ligands = kio.read_smiles(args.ligand)
    if not ligands:
        raise SystemExit("Give at least one --ligand.")
    if not args.target or len(args.target) != 1:
        raise SystemExit("Give exactly one --target. Version 1 scores compounds "
                         "against a single kinase; its scores are NOT comparable "
                         "between kinases.")

    P = bundles.load_predict("v1")
    m = P.load(bundles.bundle_dir("v1"))

    target = args.target[0]
    kw = {}
    if target.startswith("@") or len(target) > 30:
        kw["sequence"] = kio.read_sequence(target)
    else:
        kw["gene"] = target

    smis = [l["smiles"] for l in ligands]
    ranked = P.rank(m, smis, **kw)

    by_smiles = {l["smiles"]: l for l in ligands}
    rows, payload = [], []
    for i, (smi, score) in enumerate(ranked):
        name = by_smiles.get(smi, {}).get("name") or f"#{i + 1}"
        rows.append([i + 1, name, f"{score:.2f}", smi])
        payload.append({"rank": i + 1, "name": name, "smiles": smi,
                        "predicted_pIC50": float(score)})

    text = (f"Version 1 - predicted pIC50 against {target}\n\n"
            + kio.fmt_table(["#", "Compound", "pIC50", "SMILES"], rows,
                            ["<", "<", ">", "<"])
            + "\n\nThese are predicted pIC50 values on THIS kinase's own scale."
              "\nDo NOT compare them between kinases: version 1 scores each target"
              "\nseparately, so a difference between two targets reports scale"
              "\ndifferences as if they were selectivity. Use the version 2"
              "\nselectivity model for that question. Per-target accuracy across the"
              "\npanel ranges from 0.09 to 0.83.")
    _emit(args, {"model": "kfm_v1_potency", "target": target,
                 "ranking": payload}, text)
    return 0

def cmd_targets(args) -> int:
    P = bundles.load_predict(args.model)
    if args.model == "v1":
        m = P.load(bundles.bundle_dir("v1"))
        names = sorted(m["gene_to_seq_key"])
    elif args.model == "potency":
        m = P.load(bundles.bundle_dir("potency"))
        names = sorted(m["gene_to_seq_key"])
    else:
        m = P.SeqALigSeqB(bundles.bundle_dir("selectivity"))
        names = m.genes()
    if args.search:
        q = args.search.upper()
        names = [n for n in names if q in n.upper()]
    print("\n".join(names))
    print(f"\n{len(names)} target{'s' if len(names) != 1 else ''}", file=sys.stderr)
    return 0


def cmd_extend(args) -> int:
    """Hand off to the extension tool, which owns its own flags.

    REMAINDER rather than re-declaring every flag here: the two would drift, and
    a flag that silently stops being passed through is exactly the kind of thing
    that produces a merged model nobody can account for.
    """
    from . import extend
    return extend.main(["--help"]) or 0


def cmd_buildnew(args) -> int:
    from . import buildnew
    return buildnew.main(["--help"]) or 0


def cmd_where(args) -> int:
    """Where the weights are, and whether they are there.

    Exists because "I ran the download, where did it go?" is the first question
    anyone asks, and the honest answer depends on two environment variables.
    Printing the resolved path beats describing the rule.
    """
    import os
    root = bundles.cache_root()
    print(f"\nModels directory: {root}")
    if os.environ.get("KFM_HOME"):
        print("                  (set by KFM_HOME)")
    elif os.path.isdir(root):
        print("                  (found by walking up from the current directory)")
    else:
        print("                  (not created yet; `kfm download` will put it here)")
    print()
    for name in ("v1", "potency", "selectivity"):
        d = bundles.bundle_dir(name)
        here = bundles.is_present(name)
        override = os.environ.get("KFM_BUNDLE_" + name.upper())
        mark = "present" if here else "not downloaded"
        print(f"  {name:<12} {mark:<15} {d}")
        if override:
            print(f"  {'':<12} {'':<15} (from KFM_BUNDLE_{name.upper()})")
    ram = bundles.total_ram_gb()
    print("\nMemory. These forests expand about 7x from their file size when"
          "\nloaded, so what matters is RAM, not disk:")
    for name in ("v1", "potency", "selectivity"):
        need = bundles.MODELS[name].get("ram_gb")
        verdict = ""
        if ram and need:
            verdict = "ok on this machine" if ram >= need * 1.15 else "TOO SMALL on this machine"
        print("  %-12s needs ~%4.1f GB RAM   %s" % (name, need, verdict))
    if ram:
        print("  this machine has %.0f GB" % ram)

    print("\nA run reads these from disk. Nothing is fetched unless you run"
          "\n`kfm download`, so predictions work with no network at all.\n")
    return 0


def cmd_download(args) -> int:
    bundles.download(args.model, accept_licence=args.accept_licence)
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="kfm", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--ligand", "-l", action="append", metavar="SMILES[ NAME]|@FILE",
                       help="A SMILES, optionally followed by a name; or @file.smi / "
                            "@file.csv. Repeatable.")
        p.add_argument("--target", "-t", action="append", metavar="GENE|@FILE",
                       help="A gene symbol, or @file containing a sequence. Repeatable.")
        p.add_argument("--json", action="store_true", help="Machine-readable output.")

    p = sub.add_parser("potency", help="Rank compounds against one kinase")
    common(p); p.set_defaults(fn=cmd_potency)

    p = sub.add_parser("selectivity", help="Compare 2-5 kinases for each compound")
    common(p); p.set_defaults(fn=cmd_selectivity)

    p = sub.add_parser("v1", help="Version 1: predicted pIC50 against one kinase")
    common(p); p.set_defaults(fn=cmd_v1)

    p = sub.add_parser("targets", help="List the kinases a model covers")
    p.add_argument("model", choices=["v1", "potency", "selectivity"])
    p.add_argument("--search", "-s", help="Filter, case-insensitive substring")
    p.set_defaults(fn=cmd_targets)

    p = sub.add_parser("download", help="Fetch a model bundle")
    p.add_argument("model", choices=["v1", "potency", "selectivity"])
    p.add_argument("--accept-licence", "--accept-license", action="store_true",
                   dest="accept_licence",
                   help="Confirm you accept LICENSE-MODELS.txt (research and "
                        "evaluation use only).")
    p.set_defaults(fn=cmd_download)

    p = sub.add_parser(
        "extend",
        help="Add your own measurements to either v2 model, on your machine",
        description="Fit trees on your own data and concatenate them into the "
                    "released forest. Works on both v2 models: "
                    "--model potency or --model selectivity. Neither side shares "
                    "data: your CSV never leaves your machine and our corpus "
                    "never reaches you. What crosses is the feature recipe, which "
                    "both bundles already ship. See docs/KFM_EXTEND_SPEC.md.")
    # Listed so it appears in `kfm --help`; the real parsing happens in main()
    # before this parser runs, because these flags belong to the extend tool.
    p.set_defaults(fn=cmd_extend)

    p = sub.add_parser(
        "buildnew",
        help="Fit a model on ONLY your data, with no KFM weights involved",
        description="Fits a model on your comparisons alone. No KFM weights are "
                    "loaded, merged or consulted -- only the published feature "
                    "recipe. The result has the same layout and interface as a "
                    "released model. Use it when you have deep data on the "
                    "targets you care about; use `kfm extend` when you want to "
                    "keep the whole panel. See docs/EXTEND.md.")
    p.set_defaults(fn=cmd_buildnew)

    p = sub.add_parser("where",
                       help="Print where the models are cached on this machine")
    p.set_defaults(fn=cmd_where)
    return ap


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    # `extend` owns its own flags, so it is dispatched BEFORE the main parser
    # sees them. argparse.REMAINDER cannot do this: it only collects tokens
    # after the first non-flag one, so a leading --data is consumed by the
    # parent parser and the subcommand never gets it. Re-declaring every flag
    # here would work until the two copies drift, which is worse.
    if argv and argv[0] == "extend":
        from . import extend
        return extend.main(argv[1:]) or 0
    if argv and argv[0] == "buildnew":
        from . import buildnew
        return buildnew.main(argv[1:]) or 0

    args = build_parser().parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
