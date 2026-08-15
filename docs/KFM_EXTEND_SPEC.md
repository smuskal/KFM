# kfm extend - specification

**Works on both released v2 models.** The tool reads the layout from the bundle's
`MANIFEST.json` and adapts: `LigASeqLigB` potency is `[ ligand A 1038 | sequence 480 |
ligand B 1038 ]` = 2556, `SeqALigSeqB` selectivity is `[ sequence A 480 | ligand 1038 |
sequence B 480 ]` = 1998. It also picks up each model's own `min_samples_leaf` (20 and 8).

Lets a third party add their own measurements to a released KFM model **without either
side sharing data**. They fit trees on their own machine and concatenate them into the
released forest. Our corpus never leaves us; theirs never leaves them.

---

## 1. Why this works

A random forest is an ensemble of independently grown trees. Prediction is a vote average
over all of them. Trees fitted on different data can therefore be merged by concatenating
`estimators_`, which is exactly how the released model was built: 300 trees pooled from 15
independent 20-tree fits.

What crosses the boundary is the **feature recipe**, not the corpus. The partner needs the
column layout, the ligand encoder and the ESM2 recipe - all of which already ship in
`predict.py` and `embed.py`.

## 2. What is preserved, and what is not

**Preserved exactly.** The released trees are untouched. Verified by object identity, not
by assertion: the first 300 estimators of a merged forest are the same objects. Nothing is
refitted and nothing is overwritten.

**Not preserved.** Every prediction moves, because the output is a vote average. In our
reference run, adding 20 trees to 300 shifted 100% of predictions by a mean of 0.012 in
probability. The evidence is preserved; the conclusion is not. Any tool that implies
otherwise is lying.

**Influence is the tree count and nothing else.** N added to 300 gives the partner
`N/(300+N)` of the vote. That is the whole dial.

**The 0.70 operating point does not survive.** It was measured on the released model.
A merged model needs its own threshold measured on the partner's held-out data. `extend`
refuses to copy our number forward.

**Antisymmetry does survive**, because `compare()` averages both ligand orders
structurally rather than relying on the forest being symmetric.

## 3. The cost, measured

A partner holding data on five targets fitted 20 trees on 400,000 of their pairs:

| | released, 300 trees | merged, 320 trees | change |
|---|---|---|---|
| the partner's five targets | 0.7148 | 0.7513 | **+0.037** |
| everywhere else | 0.8841 | 0.8826 | −0.0015 |

Their trees vote everywhere, including where they hold nothing. Concentrated data buys a
large local gain and pays a small global tax. Both numbers should be shown to the user
rather than only the first.

## 4. Input format

**CSV, UTF-8, one comparison per row, header required.** CSV because every ELN, LIMS and
spreadsheet exports it. Column order does not matter; unknown columns are ignored so a
partner can pass an export straight through.

### Required: the two ligands

| column | meaning |
|---|---|
| `smiles_a` | ligand A, any RDKit-parseable SMILES |
| `smiles_b` | ligand B |

### Required: the target, one of two ways

| column | meaning |
|---|---|
| `gene` | a gene symbol in the released bundle (500 available, see `targets.csv`) |
| `sequence` | a raw amino-acid sequence, for a target not in the bundle |

`sequence` takes precedence if both are present. Sequences are embedded with `embed.py`
and therefore require `torch` and `transformers`. **This is how a partner extends the
protein axis**, not just the chemistry.

### Required: the outcome, one of two ways

**Preferred - give the measurements and let the tool decide:**

| column | meaning |
|---|---|
| `pic50_a`, `pic50_b` | pIC50 (or pKi/pKd) on any consistent scale, larger = more potent |
| `relation_a`, `relation_b` | optional; `=` (default), `>` or `<` |

`relation` describes the **potency value**, not a concentration: `>` means "at least this
potent". This is the reading most people expect, and it is deliberately the opposite of the
KKB source table's convention so partners never meet our internal quirk.

**It is also the opposite of ChEMBL's, and of most assay exports, so conversion inverts.**
A record of `IC50 > 1000 nM` is a statement about concentration and means the compound is
*weak*. Converted to pIC50 that is 6.0 with our relation `<`, at most this potent. A
contributor who copies the qualifier across without inverting will label every weak
compound as potent. The shipped examples include converted ChEMBL bounds so the
transformation is visible rather than described:

    IC50  > 1000 nM   ->  pic50 6.0, relation '<'
    IC50  <   10 nM   ->  pic50 8.0, relation '>'
    Ki    > 50000 nM  ->  pic50 4.3, relation '<'

Rows are used only when the two intervals are disjoint and the winner is decidable:

- point vs point - usable when the values differ
- point vs a `<` bound - usable when the point is above the bound
- point vs a `>` bound - usable when the point is below the bound
- `>` vs `<` - usable when the `>` value exceeds the `<` value
- two bounds in the same direction - never usable, dropped and counted

Two identical exact values are a **tie**: kept, and entered with the same label in both
orders so the model learns to return about 0.5 rather than inventing a winner.

**Alternative - give the answer directly:**

| column | meaning |
|---|---|
| `winner` | `A` or `B` |

Use this when the partner has a ranking but no numbers. Ties are not expressible this way.

### Example

    smiles_a,gene,smiles_b,pic50_a,pic50_b,relation_a,relation_b
    CC(C)n1nc(...)c2c(N)ncnc12,ABL1,COc1cc(...)c(Cl)cc1Cl,4.50,8.96,=,=
    CCN(CC)...,EGFR,c1ccc2...,9.10,5.00,>,=

## 5. Swap augmentation

Every usable row is entered **twice**, `[A|seq|B]` with its label and `[B|seq|A]` with the
label reversed, exactly as the released model was trained. The partner does not supply
both orders and must not; the tool does it.

## 6. Command

    python extend.py --data mydata.csv --trees 20 --out ./my_bundle

| flag | default | meaning |
|---|---|---|
| `--data` | required | the partner's CSV |
| `--out` | required | directory for the merged bundle |
| `--base` | `.` | the bundle to extend; may itself be an extended bundle |
| `--trees` | 20 | how many trees to fit. This is the influence dial |
| `--holdout` | none | a second CSV, same format, scored before and after |
| `--seed` | 0 | |
| `--min-samples-leaf` | 20 | matches the release; raising it shrinks the artifact |
| `--name` | `<base>_extended` | name for the merged model, e.g. `LigASeqLigB_v2_potency_acme` |
| `--label` | none | a human label for this extension, recorded in the log |

The merged model is a **new file with its own name in a new directory**. The base bundle is
opened read-only. `predict.py` reads the model filename from `MANIFEST.json`, so an
extended bundle loads under its own name with no code change.

Extending an already-extended bundle is supported and is the normal path as a partner's
data grows: run it again against the previous output and the provenance log accumulates.

## 7. Guardrails, all fatal

A merged pickle that loads but scores nonsense is the worst failure available here, so
each of these stops the run rather than warning:

1. feature width mismatch against the base model
2. class order mismatch
3. scikit-learn major.minor differing from the base bundle's manifest
4. fewer than 1,000 usable rows after filtering, which cannot support a useful tree
5. every row unusable - reports which rule dropped them

## 8. What the merged bundle contains

    <out>/
      LigASeqLigB_v2_potency.joblib     merged forest
      MANIFEST.json                     base manifest + extension log
      predict.py, embed.py              copied unchanged
      kinase_vectors.npz, kinase_index.json, targets.csv, targets.json
      requirements.txt
      reference_predictions.json        REGENERATED for the merged model
      reference_cases.csv               copied unchanged
      verify_export.py                  copied unchanged

`reference_predictions.json` must be regenerated, because adding trees changes every
prediction and the released reference values would otherwise fail the self-test for the
wrong reason. The base model's SHA256 is recorded so provenance survives the regeneration.

The extension log records, per extension: when, how many trees, how many rows, how many
were dropped and why, the seed, the partner's label if supplied, and the resulting tree
count and partner vote share. A merged model is auditable without access to anyone's data.

## 9. Licence

The released weights are a derived work of the KKB, a commercial database. A merged forest
is a derived work of those weights. Extension must be named explicitly in the licence
rather than left to inference; a supported workflow the terms arguably prohibit is worse
than no workflow. Academic and non-profit use is intended to be unrestricted here; the
question is industrial use and redistribution of a merged artifact.

---

## 10. Validated 2026-08-12

Run against the released bundle with a simulated partner holding 120,000 comparisons on
five targets (EGFR, ERBB2, KDR, MAPK14, BTK).

| check | result |
|---|---|
| released model byte-identical after the run | sha256 `96e3b448…` unchanged, mtime unchanged |
| ingestion | 120,000 rows read, 119,647 usable, 353 dropped for unparseable SMILES, 411 ties kept |
| fit | 20 trees on 239,294 rows (2.28 GB) in 2 s |
| merge | 320 trees, contributor share 6.2% |
| released trees preserved | first 300 estimators are the same objects |
| merged bundle self-test | ALL CHECKS PASSED, max difference 4.72e-10 |
| merged bundle loads under its own name | `LigASeqLigB_v2_potency_acme`, 320 trees |
| partner holdout, 5,955 comparisons | 0.6947 released, **0.7053** merged, +0.0106 |
| chained second extension | 320 to 330 trees, both entries in the extension log |
| accuracies withdrawn in the merged manifest | `measured_performance` and `prediction_strength` both marked WITHDRAWN |

Guardrails were exercised and all fired: a 1-row file was refused; two same-direction
bounds were dropped as unorderable; identical exact values were kept as a tie; unparseable
SMILES were dropped and counted rather than raising.

One bug was found by validating rather than assuming. `predict.ligand_features` raises on
any unparseable SMILES, which is correct when scoring and wrong when importing a partner's
file, where a few bad structures are certain. `extend` now filters and counts them.

`example_extend_data.csv` in the bundle is a 200-row file in the exact input format.

---

## 11. Both models, validated 2026-08-12

`kfm_extend.py` is one tool for both. The CSV schema differs only in how the comparison is
laid out:

| | potency (LigASeqLigB) | selectivity (SeqALigSeqB) |
|---|---|---|
| ligand columns | `smiles_a`, `smiles_b` | `smiles` |
| target columns | `gene` or `sequence` | `gene_a`/`sequence_a` and `gene_b`/`sequence_b` |
| swap augmentation | ligands exchanged | sequences exchanged |
| row width | 2,556 | 1,998 |
| base `min_samples_leaf` | 20 | 8 |

**The merged model keeps the base model's filename** - `LigASeqLigB_v2_potency.joblib` and
`model.joblib` respectively - so each bundle's own `predict.py` loads it with no code
change. The new name lives in `MANIFEST.json` and in the directory name. This was found by
testing: the selectivity `predict.py` hardcodes `model.joblib`, so renaming the artifact
broke it.

Validated end to end with ten ChEMBL comparisons per model:

| check | potency | selectivity |
|---|---|---|
| layout detected from the manifest | LSL, 2,556 | SLS, 1,998 |
| merged | 300 to 305 trees | 200 to 205 trees |
| loads in the bundle's own predict.py | yes | yes |
| ten ChEMBL comparisons, released model | 10/10 | 10/10 |
| ten ChEMBL comparisons, extended model | 10/10 | 10/10 |
| base model byte-identical afterwards | yes | yes |
| performance blocks withdrawn | 2 | via `EXTENDED_MODEL_WARNING` |

Performance blocks are matched by key name rather than by a fixed schema, because the two
bundles name them differently, and an `EXTENDED_MODEL_WARNING` is always written at the top
level so no manifest can silently carry released figures forward.

`--allow-small` permits fewer than 1,000 usable rows for format demonstrations, with a loud
warning. The shipped ten-row examples need it. Real use does not.

Example files, both derived from ChEMBL rather than the KKB, so they carry no licence
encumbrance and can ship publicly:

    examples/example_potency_chembl_10.csv
    examples/example_selectivity_chembl_10.csv

---

## 12. `kfm buildnew` - the companion utility

`kfm_buildnew.py` fits a model on the contributor's data **alone**. No KFM weights are
loaded, merged or consulted. The output has the same feature layout and the same interface
as a released model, so it drops into the same `predict.py` and the same downstream code.
It exists because the measurement below says `extend` is the wrong tool for a
data-rich contributor, and pretending otherwise would be dishonest.

### Why it exists

Measured on a contributor holding 119,660 comparisons across five targets:

| model | their 5 targets | 10 other targets |
|---|---|---|
| released KFM, 300 trees | 0.7027 | 0.7346 |
| `kfm extend`, +20 of their trees | 0.7131 | 0.7326 |
| **`kfm buildnew`, 300 trees, their data alone** | **0.7761** | 0.6192 |

And by data volume, from-scratch against the released model on the contributor's own
targets:

| their comparisons | their own model | released model |
|---|---|---|
| 2,000 | 0.6797 | 0.7027 |
| 10,000 | 0.7165 | 0.7027 |
| 40,000 | 0.7462 | 0.7027 |
| 119,660 | 0.7765 | 0.7027 |

**The crossover is near 10,000 comparisons.** Above it, a contributor gets a better model on
their own targets by ignoring the released weights entirely. What they give up is
everywhere else, and more of their own data does not fix it: 0.6027 on unrelated targets at
2,000 comparisons and 0.6192 at 119,660, against 0.7346 for the released model.

Coverage is what the released model provides. Depth on your targets is what your data
provides. `buildnew` warns below 10,000 comparisons and `--compare-to` scores a released
bundle on the same holdout so the choice is measured rather than argued.

### Command

    python kfm_buildnew.py --layout potency|selectivity --data mine.csv --out ./my_model

| flag | default | meaning |
|---|---|---|
| `--layout` | required | `potency` or `selectivity`; there is no base model to read it from |
| `--recipe` | the matching released bundle | read sequence vectors, gene map and encoders from here. **Weights are never loaded.** |
| `--trees` | 300 | |
| `--min-samples-leaf` | 20 potency, 8 selectivity | matches the released settings |
| `--holdout`, `--compare-to` | none | score your model, and optionally a released one, on the same file |

Input CSV is identical to `kfm extend`, including qualifiers.

### On the recipe dependency

`--recipe` reads only `sequence_vectors.npz`, the gene map and the encoders. Those vectors
are ESM2 embeddings of public UniProt sequences and are recomputable from public inputs;
supplying raw sequences in the CSV avoids the gene map entirely. The released `.joblib` is
never opened.

### Validated 2026-08-12

| check | potency | selectivity |
|---|---|---|
| fitted | 300 trees on 119,660 comparisons in 29 s | 20 trees on the 16-row example |
| artifact | 0.07 GB | small |
| loads in the bundle's own `predict.py` | yes | yes |
| scores the ChEMBL examples | 10/10 | 7/10 |
| released model byte-identical afterwards | yes | yes |
| holdout, 7,974 comparisons | **0.7761** against 0.7027 released, +0.0735 | - |

The selectivity demo scoring 7/10 against the released model's 10/10 is the expected and
correct result for a model fitted on sixteen rows, and is left in the record rather than
tuned away.

The manifest carries no performance figures at all, only `MEASURED_PERFORMANCE: None` and a
`COVERAGE_WARNING` with the off-target numbers above, because we have never evaluated a
contributor's model and nothing about the released accuracies applies to it.

---

## 13. `--sweep` - measure the tree count instead of guessing it

Section 3 gives our measured trade-off curve. It is the right shape but it is our
contributor's curve, not yours. `--sweep` produces yours.

    python kfm_extend.py --base ./LigASeqLigB_v2_potency --data mine.csv --out ./mine \
        --sweep --holdout mine_holdout.csv --holdout-breadth other_targets.csv

| flag | default | meaning |
|---|---|---|
| `--sweep` | off | measure accuracy against tree count, then build at the chosen count |
| `--sweep-points` | `5,10,20,40,80,150,300` | the counts to evaluate |
| `--holdout` | required with `--sweep` | comparisons on the targets your data covers |
| `--holdout-breadth` | the bundle's shipped reference set | comparisons on targets your data does **not** cover. You normally do not supply this; pass `none` to skip the check |
| `--max-breadth-loss` | 0.005 | the recommendation is the largest count whose breadth loss stays within this |

### It is exact, and it costs one fit

A forest's vote at N trees is the mean over its first N, and **the first N trees of a
larger fit are identical to an N-tree fit at the same seed** - verified by comparing tree
structure and thresholds, not assumed. So `--sweep` fits once at `max(sweep_points)`,
accumulates per-tree votes in a single pass, and reads every point off the running mean.
Seven points cost one fit rather than seven.

### Example output

    your holdout      : 7,974 comparisons, released model 0.7027
    breadth holdout   : targets your data does not cover, released model 0.7346

     --trees   share     yours     gain   breadth     cost
           5   1.6%    0.7047  +0.0020    0.7347  +0.0001
          10   3.2%    0.7075  +0.0049    0.7331  -0.0015
          20   6.2%    0.7129  +0.0103    0.7319  -0.0028
          40  11.8%    0.7193  +0.0167    0.7290  -0.0056
          80  21.1%    0.7353  +0.0326    0.7256  -0.0090
         150  33.3%    0.7504  +0.0478    0.7150  -0.0196
         300  50.0%    0.7612  +0.0586    0.6957  -0.0389

    Recommended: 20 trees, the largest whose breadth loss stays within 0.005.

The chosen count is then used to build, so a sweep is one command rather than a
diagnostic followed by a rerun. The count and the fact that a sweep chose it are recorded
in the extension log.

### The breadth check is automatic

Each bundle ships a reference set of ChEMBL comparisons spanning the whole panel -
`breadth_reference_potency.csv`, 3,744 comparisons over 468 targets, and
`breadth_reference_selectivity.csv`, 4,000 over 334. `--sweep` uses it by default, so the
normal invocation is just:

    python kfm_extend.py --base ./LigASeqLigB_v2_potency --data mine.csv --out ./mine \
        --sweep --holdout mine_holdout.csv

There are only two other states. Pass `--holdout-breadth none` and the table prints gains
only, no recommendation is made, and the build falls back to `--trees` - because a sweep
that sees what you gain and never what you give up is the thing that makes people over-set
the count. And `--sweep` with no `--holdout` at all stops with an explanation, since there
is nothing to measure against.

The reference sets are ChEMBL-derived, so they ship publicly with no licence encumbrance.

### Validated 2026-08-12

One fit of 300 trees in 29 s produced all seven points; the built model's independent
holdout check returned 0.7129, matching the sweep's 0.7129 exactly. The refusal paths were
exercised: no breadth holdout gives a table and no recommendation, and no holdout at all
exits with an explanation.
