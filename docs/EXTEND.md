# Add-ons: `kfm extend` and `kfm buildnew`

Two utilities for using your own data with the version 2 models. Both run
entirely on your machine. No data is transmitted.

| | `kfm extend` | `kfm buildnew` |
|---|---|---|
| what it does | merges trees fitted on your data into a released model | fits a model on your data alone |
| KFM weights | loaded and merged | none loaded |
| coverage | the whole 500-kinase panel | your targets only |
| licence | `LICENSE-MODELS.txt` §3(d) | no KFM encumbrance |

Measured on a contributor with 119,660 comparisons across five targets:

| model | their 5 targets | 10 other targets |
|---|---|---|
| released KFM | 0.7027 | 0.7346 |
| `kfm extend`, +20 trees | 0.7131 | 0.7326 |
| `kfm buildnew` | 0.7761 | 0.6192 |

`buildnew` overtakes the released model on the contributor's own targets at
roughly 10,000 comparisons: 0.6797 at 2,000, 0.7165 at 10,000, 0.7462 at 40,000,
0.7765 at 119,660. Off-target accuracy stays near 0.61 regardless of volume.

---

# CSV format

Shared by both tools. UTF-8, header required. Column order does not matter.
Unknown columns are ignored.

**Two accepted shapes.** Supply whichever you have; the tools detect which it is.

| Shape | One row is | You need |
|---|---|---|
| **Measurements** | one ligand, one target, one value | `smiles`, `gene` or `sequence`, `pic50` |
| **Comparisons** | both members already paired | the columns in the tables below |

A measurement file is the natural export from a knowledgebase or ChEMBL, and it
is usually what you have. The same measurement file builds either model: for
**potency** the tools pair the ligands measured on each target, and for
**selectivity** they pair the targets each ligand was measured against. Nothing
else about the file changes.

```
smiles,gene,pic50,relation
CC(C)n1nc(...)c2c(N)ncnc12,ABL1,8.0,=
COc1cc(...)c(Cl)cc1Cl,ABL1,5.0,=
```

Pairing is reported as it happens: duplicate measurements of the same ligand on
the same target are collapsed by **median** (never by most potent, which
reliably selects unit errors), which member is A is **randomised** so the forest
cannot score by reading column position, and any group offering more than
`--pairs-per-group` comparisons is sampled down to it (default 5,000) so a single
deeply screened target cannot become the model. Censored readings are handled by
the same disjoint-interval rule used for a hand-built file.

If you already hold comparisons, supply them directly and nothing is generated.

**Potency**

| column | meaning |
|---|---|
| `smiles_a`, `smiles_b` | the two ligands |
| `gene` or `sequence` | the target |

**Selectivity**

| column | meaning |
|---|---|
| `smiles` | the ligand |
| `gene_a` or `sequence_a` | the first target |
| `gene_b` or `sequence_b` | the second target |

**Outcome, either form**

| column | meaning |
|---|---|
| `pic50_a`, `pic50_b` | measurements, larger = more potent |
| `relation_a`, `relation_b` | optional: `=` (default), `>` or `<` |
| `winner` | `A` or `B`, instead of potencies |

Raw sequences are embedded with the bundle's ESM2 recipe and require
`pip install -r requirements-sequences.txt`. Gene symbols do not.

## Qualifiers

`relation` describes the **potency**:

- `>` means at least this potent
- `<` means at most this potent

ChEMBL and most assay exports place the qualifier on the **concentration**.
Conversion inverts it:

| source record | our columns |
|---|---|
| `IC50 > 1000 nM` | `pic50` 6.0, `relation` `<` |
| `IC50 < 10 nM` | `pic50` 8.0, `relation` `>` |
| `Ki > 50000 nM` | `pic50` 4.3, `relation` `<` |

Copying the qualifier across without inverting it labels weak compounds as
potent, and produces no error. ChEMBL holds 306,190 censored kinase records with
no pChEMBL value.

## Which rows are used

A pair is used only when the two intervals are disjoint:

| case | used |
|---|---|
| exact against exact, values differ | yes |
| exact above a `<` bound | yes |
| exact below a `>` bound | yes |
| a `>` bound above a `<` bound | yes |
| two bounds in the same direction | no |
| exact against exact, values identical | yes, as a tie |

Ties are entered with the same label in both orders. Every usable row is entered
twice, once in each order with the label reversed; do not supply both orders.
Unparseable SMILES are dropped and counted.

## Shipped examples

| file | rows | usable |
|---|---|---|
| `examples/extend_potency_example.csv` | 18 | 16 (2 same-direction bounds) |
| `examples/extend_selectivity_example.csv` | 16 | 16 |
| `examples/measurements_example.csv` | 60 measurements | builds **either** model: paired by target for potency, by ligand for selectivity |

Both are ChEMBL-derived. Both require `--allow-small`.

---

# `kfm extend`

```bash
./kfm.sh extend --model potency     --data mine.csv --out ./mine --trees 20
./kfm.sh extend --model selectivity --data mine.csv --out ./mine --trees 20
```

## `--trees`

`--trees N` fits N trees on your data and merges them into the 300-tree released
forest, giving your data `N / (300 + N)` of the vote on every prediction,
including targets your data does not cover.

Measured: 120,000 comparisons on five targets. Released baselines 0.7027 on those
targets, 0.7346 on ten others.

| `--trees` | share | their targets | other targets |
|---|---|---|---|
| 5 | 1.6% | 0.7049 (+0.002) | 0.7346 (0.000) |
| 20 | 6.2% | 0.7124 (+0.010) | 0.7326 (−0.002) |
| 60 | 16.7% | 0.7302 (+0.028) | 0.7272 (−0.007) |
| 150 | 33.3% | 0.7502 (+0.048) | 0.7151 (−0.020) |
| 300 | 50.0% | 0.7602 (+0.058) | 0.6966 (−0.038) |

Artifact size at 300 added trees: 0.83 GB, from 0.76 GB.

## `--sweep`

Measures the above on your own data, then builds at the recommended count.
Requires `--holdout`.

```bash
./kfm.sh extend --model potency --data mine.csv --out ./mine \
    --sweep --holdout mine_holdout.csv
```

```
your holdout      : 7,974 comparisons, released model 0.7027
breadth holdout   : targets your data does not cover, released model 0.7188

 --trees   share     yours     gain   breadth     cost
       5   1.6%    0.7047  +0.0020    0.7202  +0.0013
      10   3.2%    0.7075  +0.0049    0.7194  +0.0005
      20   6.2%    0.7129  +0.0103    0.7191  +0.0003
      40  11.8%    0.7193  +0.0167    0.7180  -0.0008
      80  21.1%    0.7353  +0.0326    0.7151  -0.0038
     150  33.3%    0.7504  +0.0478    0.7116  -0.0073
     300  50.0%    0.7612  +0.0586    0.6960  -0.0228

Recommended: 80 trees, the largest whose breadth loss stays within 0.005.
```

`--holdout` is your comparisons on targets you cover. `--holdout-breadth` is
comparisons on targets you do not; it defaults to the reference set shipped in
the bundle - `breadth_reference_potency.csv`, 3,744 comparisons over 468 targets,
or `breadth_reference_selectivity.csv`, 4,000 over 334, both ChEMBL-derived.
`--holdout-breadth none` skips the breadth check.

The sweep fits once at the largest point and reads each row from a cumulative
per-tree pass. The first N trees of a larger fit are identical to an N-tree fit
at the same seed.

### What the breadth set is for, and when to replace it

Adding your trees buys accuracy on your targets and costs accuracy everywhere
else. Nothing in your own data can measure that cost, because by definition your
data does not cover the targets it falls on. The breadth set is the control group
for that second half: a fixed set of comparisons on targets you are *not*
contributing, scored identically at every tree count, so the `cost` column means
something. Without it the sweep could report `gain` alone and would recommend the
largest tree count every time.

The shipped sets are ChEMBL-derived and are held out from the released models -
the released potency model scores 0.7188 on `breadth_reference_potency.csv`,
its generalisation accuracy rather than a memorisation signature.

**They are plain CSVs in the bundle, and you can and sometimes should replace
them.** The format is exactly the contribution format described above, so any
file you could pass to `--data` you can pass to `--holdout-breadth`. Two reasons
to do so:

- **Relevance.** If the panel you care about is narrower than the kinome - a
  family, a therapeutic area, a set of anti-targets you must not lose - build a
  breadth set over *those* targets. A recommendation tuned to average loss across
  468 targets is not tuned to the twenty you actually ship against.
- **Resolution.** Accuracy measured on a finite set carries noise, and the
  default `--max-breadth-loss` of 0.005 is a fine threshold. On the shipped
  potency set, the paired difference between nested tree counts carries roughly
  ±0.004 at the counts where the recommendation is usually made, so the
  recommendation is being read at about one standard error. Roughly quadrupling
  the number of comparisons halves that noise. If you are choosing between
  adjacent sweep rows whose `cost` figures differ by less than about 0.005, that
  difference is not resolvable on the shipped set - either enlarge the set or
  treat those rows as tied and take the smaller tree count.

Passing more comparisons costs only scoring time, which is linear and small
relative to the fit.

## Switches

| switch | default | meaning |
|---|---|---|
| `--data` | required | your CSV |
| `--out` | required | directory for the new model |
| `--model` | `potency` | which installed model to extend |
| `--base` | from `--model` | explicit bundle path; may itself be extended |
| `--trees` | 20 | trees fitted on your data |
| `--sweep` | off | measure, then build at the recommended count. Requires `--holdout` |
| `--sweep-points` | `5,10,20,40,80,150,300` | counts to evaluate |
| `--holdout` | none | your comparisons on targets you cover |
| `--holdout-breadth` | the bundle's reference set | comparisons on targets you do not cover; `none` skips |
| `--max-breadth-loss` | 0.005 | recommendation is the largest count within this |
| `--name` | `<base>_extended` | name for the merged model |
| `--label` | none | note recorded in the extension log |
| `--pairs-per-group` | 5000 | measurement input only: most comparisons drawn from any one group |
| `--seed` | 0 | also seeds pair sampling and A/B randomisation |
| `--min-samples-leaf` | the base model's | 20 potency, 8 selectivity |
| `--allow-small` | off | permit fewer than 1,000 usable rows |

## Output

Written to a new directory. The base bundle is opened read-only and is
byte-identical afterwards.

The merged model keeps the base bundle's model filename, so the bundle's own
`predict.py` loads it unchanged. `--name` and the output directory carry the new
identity, recorded in `MANIFEST.json`.

```bash
KFM_BUNDLE_POTENCY=./mine ./kfm.sh potency -t ABL1 -l @compounds.smi
```

The manifest withdraws every performance block and adds
`EXTENDED_MODEL_WARNING`. `reference_predictions.json` is renamed to `.base`.
The base model's SHA256 is recorded.

Released accuracy figures and the 0.70 operating point do not apply to an
extended model. Adding 20 trees to 300 changes 100% of predictions, by a mean of
0.012 in probability. The released trees themselves are unchanged, verified by
object identity. Re-measure accuracy and operating point on your own held-out
data before quoting either.

The manifest carries an extension log: date, tree count, rows read, rows dropped
and why, seed, label, resulting tree count and contributor share. Extending an
extended bundle is supported; the log accumulates.

---

# `kfm buildnew`

Fits a model on your data alone. No KFM weights are loaded, merged or consulted.
Output has the same feature layout and interface as a released model and loads in
the same `predict.py`.

```bash
./kfm.sh buildnew --layout potency --data mine.csv --out ./my_model --trees 300
```

`--recipe` reads `sequence_vectors.npz`, the gene map and the encoders only. The
weights file is never opened. Those vectors are ESM2 embeddings of public UniProt
sequences.

Below 10,000 comparisons the tool warns that the released model is likely to
score higher on your own targets. `--compare-to` measures that:

```bash
./kfm.sh buildnew --layout potency --data mine.csv --out ./mine \
    --holdout mine_holdout.csv --compare-to ./kfm-models/LigASeqLigB_v2_potency
```

The manifest carries `MEASURED_PERFORMANCE: None` and a `COVERAGE_WARNING`, and
no accuracy figures. The output keeps the base bundle's model filename.

Fitted on the 16-row selectivity example, the result scores 7 of 10 on cases the
released model scores 10 of 10.

## Switches

| switch | default | meaning |
|---|---|---|
| `--layout` | required | `potency` or `selectivity` |
| `--data` | required | your CSV |
| `--out` | required | directory for the new model |
| `--recipe` | the matching installed bundle | sequence vectors, gene map and encoders only |
| `--trees` | 300 | |
| `--holdout` | none | score your model on this file |
| `--compare-to` | none | also score a released bundle on the same file |
| `--min-samples-leaf` | 20 potency, 8 selectivity | |
| `--pairs-per-group` | 5000 | measurement input only: most comparisons drawn from any one group |
| `--name`, `--seed`, `--allow-small` | | as for `extend` |

---

# Guardrails

Five conditions stop the run:

1. feature width mismatch against the base model
2. class order mismatch
3. scikit-learn minor version differing from the base bundle's manifest
4. fewer than 1,000 usable rows
5. no usable rows, with a breakdown of which rule dropped them

`--allow-small` overrides condition 4.

---

# Licence

**`kfm extend`.** A merged forest is a derivative work of the released weights,
which are a derivative work of the Kinase Knowledgebase.

- Academic and non-profit research and teaching: unrestricted.
- Commercial use: requires a licence, as the released weights do.
- Redistribution of a merged model: requires written permission, including
  between affiliates.
- Your input data remains yours.

`LICENSE-MODELS.txt` §3(d).

**`kfm buildnew`.** A model fitted on your own data with a published recipe is
not a derivative work of the weights. No KFM licence encumbrance applies.

Distilling a released model's outputs into a new model is a separate act and
remains prohibited (`LICENSE-MODELS.txt` §4).

---

Specification and validation record: [`KFM_EXTEND_SPEC.md`](KFM_EXTEND_SPEC.md).
