# Using KFM on a different target family

KFM was trained and validated on kinases, but nothing in the method is specific
to them. The featuriser, the pairwise formulation and the forest work on any
protein family for which you have sequences and paired activity data.

This has been done end to end at least once, for GPCRs, using a ChEMBL pull of
1.09M activity rows across 405 receptors.

## What KFM supplies, and what you supply

| | |
|---|---|
| **KFM supplies** | The feature recipe (Morgan count fingerprint r=2/1024 plus 14 RDKit descriptors; ESM2 `t12_35M_UR50D` mean-pooled sequence embedding), the pairwise layouts, pairing, the A/B label-reversal swap, the censored-value logic, and the forest. None of it assumes a kinase. |
| **You supply** | A CSV of measurements, and the judgement about which measurements belong in the same comparison. |

## 1. Supply measurements. KFM pairs them.

One row per measurement is enough:

```
smiles,gene,pic50,relation
CC(C)n1nc(...)c2c(N)ncnc12,ADRB2,8.0,=
COc1cc(...)c(Cl)cc1Cl,ADRB2,5.0,=
```

Use `sequence` in place of `gene` for a non-kinase family — see section 3.

The **same file builds either model**. For potency the tools pair the ligands
measured on each target; for selectivity they pair the targets each ligand was
measured against. A measurement is used by whichever pairings its partners
allow, and one measurement often feeds both models.

Both tools report what they could not use:

```
  measurements usable                      29 of 60
  UNUSED, no second measurement to pair with against the same ligand:
    measurements                           31 (51.7%)
```

Read that number. A file that is rich for potency is often thin for
selectivity — selectivity needs the *same compound measured on two targets*,
which most screening data does not contain in quantity.

KFM also collapses duplicate measurements by median, randomises which member of
each pair is A, and caps any one group at `--pairs-per-group` (default 5,000) so
a deeply screened target cannot become the model. Callability of censored
readings is decided by the disjoint-interval rule.

**What KFM cannot decide for you.** These are judgements about your data, and
getting them wrong is silent:

- **Match the endpoint within a pair.** Filter to one assay type before you
  start. Comparing a K<sub>i</sub> against an IC<sub>50</sub> asks the model to
  learn the offset between two assay conventions and report it as potency.
- **Derive pActivity from value and unit yourself.** Do not trust
  `pchembl_value`. In one GPCR case it was present *only* on erroneous rows.
- **Key targets by accession, not gene symbol.** Symbol-to-accession is
  many-to-many and the errors are silent: in one GPCR roster, gene `CHRM5` was
  paired with `P08913`, which is ADRA2A.
- **Drop percentage endpoints.** Inhibition and E<sub>max</sub> are not
  concentrations and must not be log-transformed.
- **Bound plausible pActivity for your family.** A range of (-1, 13] suits
  GPCRs. A family with genuinely sub-picomolar binders needs a higher ceiling or
  real data will be discarded as unit errors.

## 2. Training

```bash
kfm buildnew --layout potency     --data my_measurements.csv --out ./my_potency
kfm buildnew --layout selectivity --data my_measurements.csv --out ./my_selectivity
```

`buildnew` fits on your data alone; no KFM weights are loaded or merged. Use
`kfm extend` only if you want to add trees to a released kinase model, which is
not what you want for a different family. Both tools accept the same input.

Defaults match the released models: 300 trees, minimum leaf 20 for potency and 8
for selectivity.

## 3. Use `sequence`, not `gene`

Gene names resolve only through a bundle's lookup table, and the released
bundles are kinase bundles. For another family, supply the raw amino-acid
sequence in `sequence` (potency) or `sequence_a` / `sequence_b` (selectivity).
KFM embeds it with the same ESM2 recipe, which needs the optional dependencies:

```bash
pip install -r requirements-sequences.txt
```

Sequence lookups are keyed by a hash of the normalised sequence, so a sequence
from another family misses the bundled index and is embedded fresh. A hit can
only occur for a byte-identical sequence, which is the same protein — the
sequence path cannot return a wrong target.

`buildnew` will not copy the recipe's gene and vector lookup tables into your
model unless the recipe actually covers the targets you fitted on, and it tells
you when it skips them. Score the model with the same `sequence` columns you
trained with.

## 4. Measure it honestly

**Hold out compounds, not comparisons.** A compound that appears in a training
pair has already shown the model its structure and its potency. Meeting it again
in a test pair is not a held-out prediction.

**Report where on the novelty scale your number sits.** Following Mattsson and
Walters (bioRxiv 2026.06.29.735309), partition your test set by each compound's
maximum Tanimoto to the compounds you fitted on, on a *binary* Morgan r=2
*2,048-bit* fingerprint — deliberately not the count-based 1,024-bit vector the
model consumes:

| Tier | Max Tanimoto |
|---|---|
| Novel | < 0.35 |
| Distant | 0.35 – 0.50 |
| Related | 0.50 – 0.70 |
| Familiar | 0.70 – 1.00 |
| Fingerprint identity | = 1.00 |

Expect a steep gradient. In the GPCR build, potency accuracy ran from 0.681 at
fingerprint identity down to 0.441 — below chance — where both ligands were
novel. A pooled figure hides that completely.

Splitting on exact structure does **not** remove the identity tier: stereoisomers,
salts and tautomers with different InChIKeys collapse onto the same fingerprint.
In the GPCR build, 3,979 of 43,508 held-out ligands still sat at Tanimoto 1.0
despite a compound-disjoint split.

**Run the sequence ablation. The answer is family-specific.** Re-score with the
sequence block zeroed, and separately with the sequence vectors *permuted* —
each target consistently given another target's embedding, so target identity
survives as a lookup key while its correspondence to biology is destroyed. The
gap between those two says whether the model is using the embedding as biology
or as a per-target index.

| | full | zeroed | permuted |
|---|---|---|---|
| Kinase, released potency model, held-out ChEMBL | 0.7188 | 0.6140 (−0.105) | 0.6033 (−0.116) |
| GPCR build, temporally split | 0.607 | 0.553 (−0.054) | 0.601 (−0.007) |

On the GPCR model permuting costs almost nothing: the embedding is acting as an
arbitrary per-target index rather than as transferable biology, which is the
proteochemometric failure mode Mattsson and Walters describe. On the released
kinase model the opposite holds — permuting costs *more* than zeroing, which is
what you see when a model has learned real sequence structure and is then handed
misleading vectors. Neither result predicts the other. Measure your own.

*Kinase figures: 3,720 held-out ChEMBL comparisons across 465 genes, mean of five
random derangements, standard deviation 0.005.*

## 5. Checklist

- [ ] One assay endpoint throughout
- [ ] pActivity derived from value and unit, not taken from `pchembl_value`
- [ ] Targets keyed by accession
- [ ] Percentage endpoints dropped
- [ ] The UNUSED count read, and acceptable for the model you are building
- [ ] No compound in both the training and held-out splits
- [ ] Accuracy reported per novelty tier, not only pooled
- [ ] Accuracy compared against zeroed and permuted sequence baselines
- [ ] The released operating point re-measured on your own held-out data; it
      does not carry over
