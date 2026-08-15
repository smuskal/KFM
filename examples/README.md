# Example test sets

Two small files, drawn from the Eidogen-Sertanty Kinase Knowledgebase, with the
**right answer encoded in each compound's name**. Run them and you can see for
yourself whether the models order things correctly, rather than taking a headline
accuracy figure on trust.

These are held-out examples in the sense that matters here - you can check the
model's output against a measurement without having to look anything up.

| File | Use with | Truth encoded in the name |
|---|---|---|
| `ABL1_potency.smi` | `kfm v1`, `kfm potency` | measured pIC50 against ABL1, and its rank |
| `selectivity_MTOR_PIK3CA_PIK3CG.smi` | `kfm selectivity` | measured pIC50 against each of the three kinases, and which one wins |

---

## Potency - works for both versions

Five ABL1 compounds spanning **8 log units**, from pIC50 12.10 down to 4.00.

```bash
python -m kfm potency -t ABL1 -l @examples/ABL1_potency.smi
```

```
#  Compound                     SMILES                                                                                 Score  Confidence
-  ---------------------------  -------------------------------------------------------------------------------------  -----  ----------
1  ABL1_rank1_KKB_pIC50_12.100  COc1ccc(cc1)S(=O)(=O)N(CC(C)C)CC(O)C(Cc2ccccc2)NC(=O)C3CN(C(=O)O3)c4cccc(c4)C(=O)C     0.641        0.64
2  ABL1_rank2_KKB_pIC50_10.140  COc1ccc(cc1F)S(=O)(=O)N(CC(C)C)CC(O)C(Cc2ccccc2)NC(=O)C3CN(C(=O)O3)c4cccc(c4)C(F)(F)F  0.635        0.63
3  ABL1_rank3_KKB_pIC50_8.050   COc1ccc(cc1)S(=O)(=O)N(CC(C)C)CC(O)C(CCCCC=C)NC(=O)c2cccc(O)c2OCCCC=C                  0.566        0.62
4  ABL1_rank4_KKB_pIC50_6.050   CN1CCC(C(O)C1)c2c(O)cc(O)c3C(=O)C=C(Oc23)c4ccccc4Cl                                    0.443        0.67
5  ABL1_rank5_KKB_pIC50_4.000   CC1=CN2CC(=O)NN=C2C=C1                                                                 0.215        0.78
```

**The order is exactly right - 1, 2, 3, 4, 5.**

Notice the confidences, though: most of these comparisons sit in the 0.60–0.70
band, where the model is right about three quarters of the time. Getting all five
in order here is a good result, not a guaranteed one. The one pair it is most
confident about (0.78) is the weakest compound against the rest - which is the
easy call.

The same file works with version 1, which returns a value rather than a
comparison:

```bash
python -m kfm v1 -t ABL1 -l @examples/ABL1_potency.smi
```

```
#  Compound                     pIC50
-  ---------------------------  -----
1  ABL1_rank1_KKB_pIC50_12.100   9.23
2  ABL1_rank2_KKB_pIC50_10.140   8.91
3  ABL1_rank3_KKB_pIC50_8.050    7.79
4  ABL1_rank4_KKB_pIC50_6.050    5.74
5  ABL1_rank5_KKB_pIC50_4.000    4.36
```

Also correctly ordered, and worth comparing against the true values in the names:
the middle of the range is close (7.79 predicted against 8.05 measured; 5.74
against 6.05), but **the top end is compressed** - 9.23 predicted for a compound
measured at 12.10. That is characteristic of a regression fitted on a distribution
where extreme potencies are rare, and it is one reason to prefer the version 2
ordering when the question is "which of these is better" rather than "how potent
is this".

---

## Selectivity

Three compounds, each measured against **MTOR, PIK3CA and PIK3CG**.

```bash
python -m kfm selectivity -t MTOR -t PIK3CA -t PIK3CG \
    -l @examples/selectivity_MTOR_PIK3CA_PIK3CG.smi
```

```
Compound                                            SMILES                                                                        MTOR  PIK3CA  PIK3CG
--------------------------------------------------  ----------------------------------------------------------------------------  ----  ------  ------
cpd2_MTOR_9.22_PIK3CA_6.91_PIK3CG_6.70_best_MTOR    CN(C)c1cccc(c1)c2ccnc3c2c(\C=C\4/Oc5cc(O)cc(O)c5C4=O)cn3C                     0.94    0.47    0.09
cpd1_MTOR_3.31_PIK3CA_8.46_PIK3CG_7.65_best_PIK3CA  CC(C)N1CCCN(CC1)c2nc(nc(n2)c3ccc(NC(=O)Nc4ccc(cc4)C(=O)N5CCNCC5)cc3)N6CCOCC6  0.30    0.90    0.31
cpd3_MTOR_7.30_PIK3CA_6.69_PIK3CG_5.80_best_MTOR    FC(F)(F)c1ccc(NC(=O)Nc2ccc(cc2)c3nc(nc(n3)N4CCOCC4)N5CCOCC5)nc1               0.81    0.62    0.07
```

**Three out of three.** In every row the largest number sits under the kinase
named in that compound's `best_` label.

Two details worth reading properly:

- The scores track the *size* of the gap, not just its direction. cpd2 wins on
  MTOR by more than two log units and scores 0.94; cpd3 wins by 0.6 log units and
  scores 0.81. The model is less certain when the real difference is smaller,
  which is the behaviour you want.
- **0.50 means no preference.** cpd1's MTOR score of 0.30 says it *disprefers*
  MTOR - correct, it is measured 5 log units weaker there.

---

## What these examples cannot tell you

Eight compounds is a demonstration, not a validation. The published accuracies
come from 1.8 million (potency) and 3.1 million (selectivity) held-out ChEMBL
comparisons, and they are the numbers to quote - see
<https://kinasefoundationmodel.com/v2/limitations.html>.

In particular, these examples were chosen to span wide potency gaps, which is
where the models are strongest. **Compounds that are close in potency come back
near 0.50 and are declined rather than guessed** - that is correct behaviour, and
you will see much more of it on a real series of analogues than you see here.

## File format

Tab- or space-separated, one compound per line, `SMILES` then an optional name.
Lines starting with `#` are comments. The same format `kfm` accepts anywhere via
`-l @file`; `.csv` with a `smiles` column works too.
