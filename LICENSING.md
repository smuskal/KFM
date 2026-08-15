# Licensing - how this repository is split, and why

**You do not need a second GitHub repository.** One repository can carry two
licenses as long as it is unambiguous which file is under which. That is the
normal arrangement for open-source code that ships restricted model weights
(Meta's Llama, Stability's models and BioNeMo all do a version of this), and it
is what this repository does.

## The split

| What | License | File |
|---|---|---|
| All source code, CLI tools, docs, examples | **Apache License 2.0** | `LICENSE` |
| Trained model artefacts - the `.joblib` forests, the `*_vectors.npz` caches, their index and gene-map files | **Research and evaluation only** | `LICENSE-MODELS.txt` |

The boundary is the trained artefact. Anything a computer *ran* is Apache-2.0.
Anything that *came out of training on the KKB* is under the model license.

## Why not Apache-2.0 for the weights too

The weights are a derived work of the Kinase Knowledgebase, which is the
commercial asset. Apache-2.0 is irrevocable and permits commercial use and
redistribution, so releasing the weights under it would give away, permanently,
the thing customers pay for - and it could not be walked back afterwards, because
every copy already distributed keeps its grant.

Apache-2.0 is exactly right for the *code*: it is permissive, it is compatible
with every dependency here, and unlike MIT it carries an explicit patent grant,
which matters when a company publishes a method.

## What you need to do on GitHub

You already created https://github.com/smuskal/KFM with Apache-2.0, which is the
correct choice for the repository's primary license. GitHub reads the root
`LICENSE` file to label the repo, and that stays Apache-2.0. To add the model
license, all you need is:

1. Keep `LICENSE` exactly as GitHub created it (Apache-2.0).
2. Add `LICENSE-MODELS.txt` (in this directory) at the repository root.
3. Add `LICENSING.md` (this file) at the root.
4. Make sure the `README.md` states the split near the top - ours does.
5. Put a `NOTICE` file at the root. Apache-2.0 has a specific mechanism for
   this: a `NOTICE` file travels with redistributions, so it is the correct
   place to record that the repository contains separately-licensed artefacts.

GitHub will show the repository as "Apache-2.0" in its sidebar. That label
describes the repository's primary license, not every file in it, and the
`LICENSE-MODELS.txt` plus the README statement are what make the weights'
terms binding on anyone who downloads them. This is standard and defensible.

## The one thing to decide before publishing weights

Nothing above matters if the weights are not actually in the repository. Model
files here run to gigabytes and **GitHub rejects any file over 100 MB**, so the
forests cannot be committed directly even if you want them public. In practice
that gives you a useful choice:

- **Recommended: do not put the weights on GitHub at all.** Ship the code, the
  CLI and the documentation; have the CLI download the weights on first use from
  a location you control, behind a click-through of `LICENSE-MODELS.txt`. You
  keep the ability to see who took a copy, and to change the terms later.
- If you do want them publicly downloadable, use GitHub Releases (2 GB per file)
  or a public GCS bucket, and link from the README. The potency bundle is 780 MB
  and fits; the selectivity forest is 11.7 GB raw, 3.0 GB recompressed, and
  needs GCS either way.

**Recommendation: publish the code now, and keep the weights gated.** It gives
you the credibility of an open method without giving away the database, and it
is far easier to loosen later than to tighten.
