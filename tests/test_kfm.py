"""Tests for the KFM command line tools.

Two tiers, and the split matters:

  * Tests that need NO model bundle - argument parsing, file reading, the
    confidence bands, the download licence gate. These always run.

  * Tests that DO need a bundle, marked `needs_bundle`. They are skipped with a
    clear reason when the weights are not installed, and they check the numbers
    that must never move: the worked examples from the two published reports.

Run everything:            python -m pytest tests/ -v
Run only what always runs: python -m pytest tests/ -v -m "not needs_bundle"

To include the model tests, install the weights or point at a bundle you have:

    export KFM_BUNDLE_POTENCY=/path/to/LigASeqLigB_v2_potency
    export KFM_BUNDLE_SELECTIVITY=/path/to/SeqALigSeqB_v2_selectivity
"""
import json
import os
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kfm import bundles, io as kio          # noqa: E402
from kfm.cli import band_accuracy           # noqa: E402

# The worked examples from the two reports. These are the canaries: if either
# moves, the installed library versions are not the pinned ones, or the bundle
# is not the published one. Neither is a thing to shrug at.
BOSUTINIB = "COc1cc(Nc2c(cnc3cc(OCCCN4CCN(C)CC4)c(OC)cc23)C#N)c(Cl)cc1Cl"
PP1_TYPE = "CC(C)n1nc(c2cccnc2)c3c(N)ncnc13"
DASATINIB = "Cc1nc(Nc2ncc(s2)C(=O)Nc2c(C)cccc2Cl)cc(n1)N1CCN(CCO)CC1"

needs_potency = pytest.mark.skipif(
    not bundles.is_present("potency"),
    reason="potency bundle not installed (kfm download potency --accept-licence)")
needs_selectivity = pytest.mark.skipif(
    not bundles.is_present("selectivity"),
    reason="selectivity bundle not installed")


def _is_released(model):
    """Is the installed bundle THE released one, rather than one someone built.

    The worked-example values below belong to a specific released model. A model
    built by `buildnew`, or any future release, returns different numbers by
    construction, so asserting them there fails a healthy bundle. Those tests are
    skipped instead; `test_reproduces_its_own_reference_predictions` covers every
    bundle, whichever it is.
    """
    if not bundles.is_present(model):
        return False
    try:
        man = json.load(open(os.path.join(bundles.bundle_dir(model),
                                          "MANIFEST.json")))
    except Exception:                                        # noqa: BLE001
        return False
    return man.get("built_by") != "kfm_buildnew" and not man.get("extensions")


released_potency = pytest.mark.skipif(
    not _is_released("potency"),
    reason="not the released potency bundle; its own canary is tested instead")
released_selectivity = pytest.mark.skipif(
    not _is_released("selectivity"),
    reason="not the released selectivity bundle; its own canary is tested instead")


@pytest.mark.needs_bundle
@pytest.mark.parametrize("model", ["potency", "selectivity"])
def test_reproduces_its_own_reference_predictions(model):
    """Every bundle must reproduce the answers it shipped with.

    This is the check that survives a new model. It reads the bundle's own
    reference_cases.csv and reference_predictions.json rather than any literal,
    so it is valid for the released models, for a version 3, and for a model a
    user built themselves.
    """
    if not bundles.is_present(model):
        pytest.skip(f"{model} bundle not installed")
    d = bundles.bundle_dir(model)
    if not os.path.exists(os.path.join(d, "reference_predictions.json")):
        pytest.skip(f"{model} bundle predates the canary")
    r = subprocess.run([sys.executable, "-m", "kfm", "verify", model],
                       capture_output=True, text=True,
                       cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    assert r.returncode == 0, r.stdout + r.stderr


# ---------------------------------------------------------------------------
# Input parsing - no bundle needed
# ---------------------------------------------------------------------------
class TestReadSmiles:
    def test_bare_smiles_has_no_name(self):
        assert kio.read_smiles(["CCO"]) == [{"smiles": "CCO", "name": None}]

    def test_name_after_whitespace(self):
        assert kio.read_smiles(["CCO ethanol"]) == [
            {"smiles": "CCO", "name": "ethanol"}]

    def test_name_may_contain_spaces(self):
        # "compound 12a" is one name, not a name and a stray token.
        assert kio.read_smiles(["CCO compound 12a"])[0]["name"] == "compound 12a"

    def test_smi_file(self, tmp_path):
        f = tmp_path / "in.smi"
        f.write_text("# a comment\nCCO ethanol\nCCCO propanol\n\n")
        got = kio.read_smiles(["@" + str(f)])
        assert got == [{"smiles": "CCO", "name": "ethanol"},
                       {"smiles": "CCCO", "name": "propanol"}]

    def test_csv_file(self, tmp_path):
        f = tmp_path / "in.csv"
        f.write_text("name,smiles\nethanol,CCO\npropanol,CCCO\n")
        got = kio.read_smiles(["@" + str(f)])
        assert [g["smiles"] for g in got] == ["CCO", "CCCO"]
        assert [g["name"] for g in got] == ["ethanol", "propanol"]

    def test_csv_without_smiles_column_is_a_clear_error(self, tmp_path):
        f = tmp_path / "bad.csv"
        f.write_text("name,structure_id\nethanol,1\n")
        with pytest.raises(SystemExit) as e:
            kio.read_smiles(["@" + str(f)])
        assert "SMILES column" in str(e.value)

    def test_missing_file_names_the_path(self, tmp_path):
        with pytest.raises(SystemExit) as e:
            kio.read_smiles(["@" + str(tmp_path / "nope.smi")])
        assert "nope.smi" in str(e.value)

    def test_hash_inside_a_smiles_is_not_a_comment(self, tmp_path):
        # '#' is a triple bond in SMILES. Treating it as a comment mid-line
        # would silently corrupt the structure, so it only counts at line start.
        f = tmp_path / "in.smi"
        f.write_text("CC#N acetonitrile\n")
        assert kio.read_smiles(["@" + str(f)])[0]["smiles"] == "CC#N"


class TestReadSequence:
    def test_uppercases_and_strips_non_letters(self):
        assert kio.read_sequence(" ac de\nfg ") == "ACDEFG"

    def test_strips_fasta_header(self, tmp_path):
        f = tmp_path / "s.fasta"
        f.write_text(">sp|P00519|ABL1_HUMAN\nMLEICLKLVG\nSSKQ\n")
        assert kio.read_sequence("@" + str(f)) == "MLEICLKLVGSSKQ"

    def test_lowercase_is_uppercased(self):
        # Not cosmetic: the ESM2 tokenizer has no lowercase vocabulary, so a
        # lowercase sequence collapses to a single <unk> and scores as garbage
        # with no error at all.
        assert kio.read_sequence("mleicl") == "MLEICL"


# ---------------------------------------------------------------------------
# Confidence bands - no bundle needed
# ---------------------------------------------------------------------------
class TestBands:
    def test_potency_bands(self):
        assert band_accuracy("potency", 0.95) == "89.3%"
        assert band_accuracy("potency", 0.85) == "90.0%"
        assert band_accuracy("potency", 0.55) == "60.3%"

    def test_selectivity_bands(self):
        assert band_accuracy("selectivity", 0.95) == "99.3%"
        assert band_accuracy("selectivity", 0.85) == "96.4%"
        assert band_accuracy("selectivity", 0.55) == "57.9%"

    def test_the_two_models_are_not_interchangeable(self):
        # Measured on different test sets. Copying one model's bands to the
        # other would overstate potency by eleven points in this band.
        assert band_accuracy("potency", 0.85) != band_accuracy("selectivity", 0.85)

    def test_boundaries_land_in_the_higher_band(self):
        assert band_accuracy("potency", 0.70) == "87.0%"
        assert band_accuracy("potency", 0.699) == "77.6%"


class TestBandsBelongToTheModelThatEarnedThem:
    """A band is a measurement of one model. It must not be lent to another.

    The table in cli.py is keyed by the COMMAND, which is right only while the
    bundle behind that command is the released one. Pointed at a model built by
    `kfm buildnew` or grown by `kfm extend`, the CLI used to quote the released
    model's accuracy for a model nobody has ever evaluated, contradicting that
    bundle's own manifest.
    """

    def test_a_released_bundle_is_unchanged(self):
        """No manifest, or one without bands, keeps the published table."""
        assert band_accuracy("potency", 0.85, None) == "90.0%"
        assert band_accuracy("potency", 0.85, {"name": "released"}) == "90.0%"
        assert band_accuracy("selectivity", 0.85, None) == "96.4%"

    def test_a_model_you_built_reports_no_accuracy(self):
        assert band_accuracy("potency", 0.85, {"built_by": "kfm_buildnew"}) is None

    def test_an_extended_model_reports_no_accuracy(self):
        assert band_accuracy("potency", 0.95,
                             {"EXTENDED_MODEL_WARNING": "figures withdrawn"}) is None
        assert band_accuracy("potency", 0.95,
                             {"extensions": [{"trees_added": 5}]}) is None

    def test_a_bundle_carrying_its_own_bands_uses_them(self):
        """How version 3 will ship its numbers: in the manifest, not in code."""
        m = {"confidence_bands": [[0.90, "91.1%"], [0.70, "80.0%"], [0.00, "55.0%"]]}
        assert band_accuracy("potency", 0.95, m) == "91.1%"
        assert band_accuracy("potency", 0.75, m) == "80.0%"
        assert band_accuracy("potency", 0.10, m) == "55.0%"


# ---------------------------------------------------------------------------
# Table formatting - no bundle needed
# ---------------------------------------------------------------------------
class TestTable:
    def test_nothing_is_truncated(self):
        long_smiles = "C" * 120
        out = kio.fmt_table(["SMILES"], [[long_smiles]])
        assert long_smiles in out
        assert "…" not in out and "..." not in out

    def test_columns_line_up(self):
        out = kio.fmt_table(["a", "b"], [["x", "yy"], ["zzz", "w"]])
        lines = out.splitlines()
        assert len(set(len(l) for l in lines)) == 1


# ---------------------------------------------------------------------------
# Bundle resolution and the licence gate - no bundle needed
# ---------------------------------------------------------------------------
class TestBundles:
    def test_env_override_wins(self, monkeypatch, tmp_path):
        monkeypatch.setenv("KFM_BUNDLE_POTENCY", str(tmp_path))
        assert bundles.bundle_dir("potency") == str(tmp_path)

    def test_models_land_in_the_working_directory(self, monkeypatch, tmp_path):
        """Weights go to ./kfm-models, where the user can see them.

        They used to default to ~/.cache/kfm. That hid several gigabytes in a
        directory people do not look in and cannot easily find again, so the
        default moved to the directory the download is run from.
        """
        monkeypatch.delenv("KFM_BUNDLE_POTENCY", raising=False)
        monkeypatch.delenv("KFM_HOME", raising=False)
        monkeypatch.chdir(tmp_path)
        d = bundles.bundle_dir("potency")
        assert d == str(tmp_path / "kfm-models" / "LigASeqLigB_v2_potency")

    def test_a_subdirectory_still_finds_the_models(self, monkeypatch, tmp_path):
        """Resolution walks up, the way git finds .git.

        Without this, running a command from a subdirectory of the project would
        silently miss models already downloaded at its root and offer to fetch
        several gigabytes again.
        """
        monkeypatch.delenv("KFM_BUNDLE_POTENCY", raising=False)
        monkeypatch.delenv("KFM_HOME", raising=False)
        (tmp_path / "kfm-models").mkdir()
        deep = tmp_path / "a" / "b"
        deep.mkdir(parents=True)
        monkeypatch.chdir(deep)
        assert bundles.cache_root() == str(tmp_path / "kfm-models")

    def test_the_models_directory_is_gitignored(self):
        """The weights cannot be committed, which is what the old
        outside-the-repo default was really protecting against. Now that they
        land in the working directory, .gitignore is what enforces it."""
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        with open(os.path.join(repo, ".gitignore")) as fh:
            assert "kfm-models/" in fh.read().split()

    def test_download_refuses_without_licence_acceptance(self):
        with pytest.raises(SystemExit) as e:
            bundles.download("potency", accept_licence=False)
        msg = str(e.value)
        assert "licensed separately" in msg and "--accept-licence" in msg

    def test_missing_bundle_says_how_to_get_it(self, monkeypatch, tmp_path):
        monkeypatch.setenv("KFM_BUNDLE_POTENCY", str(tmp_path / "absent"))
        with pytest.raises(SystemExit) as e:
            bundles.load_predict("potency")
        assert "kfm download potency" in str(e.value)


# ---------------------------------------------------------------------------
# The published numbers - these need the weights
# ---------------------------------------------------------------------------
def run_cli(*args):
    """Invoke the CLI the way a user does, and return parsed JSON."""
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = subprocess.run([sys.executable, "-m", "kfm", *args, "--json"],
                         capture_output=True, text=True, cwd=repo, timeout=900)
    assert out.returncode == 0, out.stderr[-2000:]
    return json.loads(out.stdout)


@needs_potency
@pytest.mark.needs_bundle
class TestPotencyNumbers:
    @released_potency
    def test_worked_example_from_the_report(self):
        """ABL1, bosutinib vs the PP1-type compound.

        The report prints 0.846 and 0.151 - its SINGLE-ORDER figures, which sum
        to 0.997. rank() averages both ligand orders, and (0.846 + (1-0.151))/2
        = 0.8475. So 0.847 IS the report's number correctly aggregated. If this
        ever reads 0.846, the averaging has been lost.
        """
        d = run_cli("potency", "-t", "ABL1",
                    "-l", f"{BOSUTINIB} bosutinib", "-l", f"{PP1_TYPE} PP1-type")
        top = d["ranking"][0]
        assert top["name"] == "bosutinib"
        assert round(top["score"], 3) == 0.847
        assert round(d["ranking"][1]["score"], 3) == 0.153

    def test_scores_of_a_pair_sum_to_one(self):
        d = run_cli("potency", "-t", "ABL1",
                    "-l", f"{BOSUTINIB} a", "-l", f"{PP1_TYPE} b")
        assert abs(sum(r["score"] for r in d["ranking"]) - 1.0) < 1e-9

    def test_identical_compounds_score_exactly_a_half(self):
        # The strongest single sanity check available: a molecule against itself
        # cannot be more potent than itself.
        d = run_cli("potency", "-t", "ABL1",
                    "-l", f"{PP1_TYPE} one", "-l", f"{PP1_TYPE} two")
        for r in d["ranking"]:
            assert abs(r["score"] - 0.5) < 1e-9

    def test_argument_order_does_not_change_the_answer(self):
        a = run_cli("potency", "-t", "ABL1",
                    "-l", f"{BOSUTINIB} bosutinib", "-l", f"{PP1_TYPE} PP1-type")
        b = run_cli("potency", "-t", "ABL1",
                    "-l", f"{PP1_TYPE} PP1-type", "-l", f"{BOSUTINIB} bosutinib")
        sa = {r["name"]: round(r["score"], 6) for r in a["ranking"]}
        sb = {r["name"]: round(r["score"], 6) for r in b["ranking"]}
        assert sa == sb


@needs_selectivity
@pytest.mark.needs_bundle
class TestSelectivityNumbers:
    @released_selectivity
    def test_worked_example_from_the_report(self):
        """ABL1 / dasatinib / GSK3B must give 0.9737 and 0.0263."""
        d = run_cli("selectivity", "-t", "ABL1", "-t", "GSK3B",
                    "-l", f"{DASATINIB} dasatinib")
        r = d["results"][0]
        assert r["prefers"] == "ABL1"
        assert round(r["scores"]["ABL1"], 4) == 0.9737
        assert round(r["scores"]["GSK3B"], 4) == 0.0263

    def test_scores_average_to_a_half_across_the_panel(self):
        """The property the heat map's diverging colour scale depends on.

        A compound's scores across the selected kinases average to exactly 0.50,
        which is what makes 0.50 a real "no preference" midpoint rather than a
        low value.
        """
        d = run_cli("selectivity", "-t", "ABL1", "-t", "SRC", "-t", "GSK3B",
                    "-t", "EGFR", "-l", f"{DASATINIB} dasatinib")
        scores = d["results"][0]["scores"]
        assert abs(sum(scores.values()) / len(scores) - 0.5) < 1e-9

    def test_known_drugs_land_on_their_known_targets(self):
        gefitinib = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
        d = run_cli("selectivity", "-t", "ABL1", "-t", "SRC", "-t", "GSK3B",
                    "-t", "EGFR",
                    "-l", f"{DASATINIB} dasatinib", "-l", f"{gefitinib} gefitinib")
        got = {r["name"]: r["prefers"] for r in d["results"]}
        assert got["gefitinib"] == "EGFR"
        assert got["dasatinib"] in ("ABL1", "SRC")   # a dual ABL/SRC inhibitor

    def test_at_most_five_targets(self):
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        out = subprocess.run(
            [sys.executable, "-m", "kfm", "selectivity",
             *sum([["-t", g] for g in
                   ("ABL1", "SRC", "GSK3B", "EGFR", "CDK2", "MAPK1")], []),
             "-l", "CCO"],
            capture_output=True, text=True, cwd=repo, timeout=300)
        assert out.returncode != 0
        assert "At most 5 targets" in (out.stdout + out.stderr)


class TestCliContract:
    """Argument handling that must fail clearly rather than confusingly."""

    def _run(self, *args):
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        return subprocess.run([sys.executable, "-m", "kfm", *args],
                              capture_output=True, text=True, cwd=repo, timeout=300)

    def test_potency_needs_two_compounds(self):
        out = self._run("potency", "-t", "ABL1", "-l", "CCO")
        assert out.returncode != 0
        assert "at least two" in (out.stdout + out.stderr)

    def test_selectivity_needs_two_targets(self):
        out = self._run("selectivity", "-t", "ABL1", "-l", "CCO")
        assert out.returncode != 0
        assert "at least two" in (out.stdout + out.stderr)

    def test_help_works_without_any_bundle(self):
        out = self._run("--help")
        assert out.returncode == 0
        assert "potency" in out.stdout and "selectivity" in out.stdout


class TestEnvironmentFilesAgree:
    """environment.yml repeats requirements.txt's pins, so they must match.

    The repetition is forced: conda copies the pip block to a temp file, so a
    relative `-r requirements.txt` in environment.yml fails outright. Repeating
    the versions is the reliable form, and this test is what stops the two
    copies drifting apart -- a drift that would hand conda users a different
    scikit-learn from pip users and change the numbers with no warning.
    """

    @staticmethod
    def _pins(lines):
        out = {}
        for raw in lines:
            line = raw.strip().lstrip("-").strip()
            if not line or line.startswith("#"):
                continue
            # Compare every constraint, not only "==". rdkit is deliberately a
            # floor rather than a pin, and an earlier version of this helper
            # skipped anything without "==" -- so the one dependency most likely
            # to drift between the two files was the one it did not check.
            for op in ("==", ">=", "~=", ">"):
                if op in line:
                    name, version = line.split(op, 1)
                    out[name.strip()] = op + version.strip()
                    break
        return out

    def test_pins_match(self):
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        req = self._pins(open(os.path.join(repo, "requirements.txt")))
        env = self._pins(open(os.path.join(repo, "environment.yml")))
        for pkg, version in req.items():
            assert pkg in env, f"{pkg} is in requirements.txt but not environment.yml"
            assert env[pkg] == version, (
                f"{pkg} pinned to {version} in requirements.txt but "
                f"{env[pkg]} in environment.yml")

    def test_environment_yml_has_no_relative_requirements_include(self):
        """A `-r requirements.txt` DIRECTIVE breaks `conda env create`.

        Conda copies the pip block to a temp file, so the relative path resolves
        against /tmp. Comments mentioning it are fine - only real directives are
        checked, which is why this looks at stripped, non-comment lines rather
        than searching the whole file.
        """
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        for raw in open(os.path.join(repo, "environment.yml")):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            directive = line.lstrip("-").strip()
            assert not directive.startswith("-r "), (
                f"environment.yml has a relative requirements include: {line}")


@needs_potency
@pytest.mark.needs_bundle
class TestShippedExamples:
    """The examples in examples/ must keep producing the documented result.

    They are the first thing a new user runs, and a wrong answer there is worse
    than no example at all. Ground truth is encoded in each compound's name, so
    these assertions read it back out of the file rather than hard-coding it.
    """

    def test_potency_example_orders_by_measured_pic50(self):
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        d = run_cli("potency", "-t", "ABL1",
                    "-l", "@" + os.path.join(repo, "examples", "ABL1_potency.smi"))
        got = [r["name"] for r in d["ranking"]]
        want = sorted(got, key=lambda n: -float(n.rsplit("_", 1)[1]))
        assert got == want, f"expected {want}, got {got}"

    def test_example_file_names_encode_their_truth(self):
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        for line in open(os.path.join(repo, "examples", "ABL1_potency.smi")):
            if not line.strip() or line.startswith("#"):
                continue
            name = line.split("\t")[1].strip()
            float(name.rsplit("_", 1)[1])          # parses, or the file is wrong


class TestMeasurementInput:
    """One measurement per row is paired by the tools, not by the user.

    The same file must build either model: potency pairs the ligands measured on
    each target, selectivity pairs the targets each ligand was measured against.
    """

    def _rows(self, layout):
        import kfm.extend as kx
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        src = os.path.join(repo, "examples", "measurements_example.csv")
        import csv as _csv
        rows = list(_csv.DictReader(open(src)))
        kind = "LSL" if layout == "potency" else "SLS"
        return kx.pair_measurements(rows, kind, cap=5000, seed=0)

    def test_potency_pairs_by_target(self):
        out = self._rows("potency")
        assert out, "no comparisons generated"
        for r in out:
            assert "smiles_a" in r and "smiles_b" in r and "gene" in r
            assert r["smiles_a"] != r["smiles_b"]

    def test_selectivity_pairs_by_ligand(self):
        out = self._rows("selectivity")
        assert out, "no comparisons generated"
        for r in out:
            assert "smiles" in r and "gene_a" in r and "gene_b" in r
            assert r["gene_a"] != r["gene_b"]

    def test_a_wins_is_not_column_position(self):
        """If A always held the winner the forest would score by reading the slot."""
        out = self._rows("potency")
        a_wins = sum(1 for r in out if float(r["pic50_a"]) > float(r["pic50_b"]))
        frac = a_wins / len(out)
        assert 0.3 < frac < 0.7, f"A-wins fraction {frac:.2f} is not balanced"

    def test_paired_files_are_not_reinterpreted(self):
        """A file that already holds comparisons must be left alone."""
        import kfm.extend as kx
        cols = {"smiles_a", "smiles_b", "gene", "pic50_a", "pic50_b"}
        assert not kx._is_measurement_file(cols, "LSL")
        assert kx._is_measurement_file({"smiles", "gene", "pic50"}, "LSL")


class TestChunkPlanning:
    """The arithmetic that decides whether a fit has to be chunked.

    No bundle needed, so these always run. The two reference points are the
    released models' own design matrices, and they are what makes this testable
    rather than a guess.
    """

    def test_matrix_size_matches_the_released_selectivity_model(self):
        from kfm.buildnew import matrix_gb
        # SeqALigSeqB v2: 4,340,117 comparisons, 1,998 features, and its spec
        # records the resulting matrix as 64.6 GB.
        assert round(matrix_gb(4_340_117, 1998), 1) == 64.6

    def test_a_fit_that_fits_is_not_chunked(self):
        from kfm.buildnew import plan_chunks
        assert plan_chunks(5_000, 2556, budget_gb=32, trees=300) == 1

    def test_a_fit_that_does_not_fit_is_split(self):
        from kfm.buildnew import plan_chunks
        # The released potency pool: 22.5M comparisons at 2,556 features is
        # 428 GB, so a 31 GB budget cannot take it in one pass.
        assert plan_chunks(22_500_000, 2556, budget_gb=31, trees=300) == 14

    def test_never_more_chunks_than_trees(self):
        from kfm.buildnew import plan_chunks
        assert plan_chunks(10**9, 2556, budget_gb=1, trees=8) == 8

    def test_every_chunk_gets_a_share_of_every_group(self):
        """The sampling rule, and it is not cosmetic.

        A uniform split hands one chunk all of a rare target and the others
        none, so most trees never see it. Round-robin within each group keeps
        every target in every chunk.
        """
        from kfm.buildnew import stratified_chunks
        usable = ([("a", "b", ("gene", "EGFR"), None, 1, False)] * 100
                  + [("c", "d", ("gene", "CSK"), None, 1, False)] * 6)
        parts = stratified_chunks(usable, "LSL", 3, seed=0)
        assert len(parts) == 3
        for p in parts:
            genes = {usable[i][2][1] for i in p}
            assert genes == {"EGFR", "CSK"}, f"a chunk missed a target: {genes}"
        assert sum(len(p) for p in parts) == len(usable), "rows were lost or duplicated"


def run_tool(*args, timeout=1800):
    """Invoke `extend` or `buildnew` the way a user does.

    Not run_cli: those two own their own flags and are dispatched before the
    main parser, so they take no --json and print a human log. Returns the
    combined output so a test can assert on what the tool told the user.
    """
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = subprocess.run([sys.executable, "-m", "kfm", *args],
                         capture_output=True, text=True, cwd=repo, timeout=timeout)
    assert out.returncode == 0, (out.stdout + out.stderr)[-3000:]
    return out.stdout + out.stderr


def _example(name):
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(repo, "examples", name)


@needs_potency
@pytest.mark.needs_bundle
class TestBuildNewEndToEnd:
    """`kfm buildnew` fits on the user's data alone and writes a usable bundle.

    The pairing tests above check the arithmetic in isolation. These run the
    tool, which is the only way to catch a break between pairing, featurising,
    fitting and writing -- the seam where the Aug 2026 measurement-input change
    landed. Everything is written under tmp_path, so nothing touches ./kfm-models.
    """

    def test_fits_a_potency_model_from_one_measurement_per_row(self, tmp_path):
        out = tmp_path / "my_potency"
        log = run_tool("buildnew", "--layout", "potency",
                       "--data", _example("measurements_example.csv"),
                       "--out", str(out), "--trees", "8", "--allow-small")
        man = json.load(open(out / "MANIFEST.json"))
        assert man["layout"] == "potency"
        assert any(f.endswith(".joblib") for f in os.listdir(out)), \
            f"no model written; tool said:\n{log[-1500:]}"

    def test_chunked_and_pooled_gives_the_trees_asked_for(self, tmp_path):
        """--trees is the total across all chunks, not per chunk.

        Getting this wrong multiplies the forest by the chunk count, which is
        the kind of error that shows up as a mysteriously enormous model file
        rather than as a failure.
        """
        out = tmp_path / "pooled"
        run_tool("buildnew", "--layout", "potency",
                 "--data", _example("measurements_example.csv"),
                 "--out", str(out), "--trees", "12", "--chunks", "3",
                 "--allow-small")
        import joblib
        man = json.load(open(out / "MANIFEST.json"))
        assert man["training"]["chunks"] == 3
        assert man["hyperparameters"]["n_estimators"] == 12
        forest = joblib.load(out / man["model_file"])
        assert forest.n_estimators == 12
        assert len(forest.estimators_) == 12

    def test_chunks_are_checkpointed_and_reused(self, tmp_path):
        """An interrupted chunked fit must cost one chunk, not the whole run."""
        out = tmp_path / "ckpt"
        run_tool("buildnew", "--layout", "potency",
                 "--data", _example("measurements_example.csv"),
                 "--out", str(out), "--trees", "9", "--chunks", "3",
                 "--allow-small")
        chunks = sorted((out / "chunks").glob("chunk_*.joblib"))
        assert len(chunks) == 3, [c.name for c in chunks]
        stamps = {c.name: c.stat().st_mtime_ns for c in chunks}
        log = run_tool("buildnew", "--layout", "potency",
                       "--data", _example("measurements_example.csv"),
                       "--out", str(out), "--trees", "9", "--chunks", "3",
                       "--allow-small")
        assert "reusing" in log
        for c in sorted((out / "chunks").glob("chunk_*.joblib")):
            assert c.stat().st_mtime_ns == stamps[c.name], f"{c.name} was refitted"

    @needs_selectivity
    def test_the_same_file_builds_the_selectivity_layout(self, tmp_path):
        """The point of one-measurement-per-row: one file, either model.

        Potency pairs the ligands measured on each target; selectivity pairs the
        targets each ligand was measured against. If this ever needs a different
        input file from the test above, that property has been lost.

        Guarded on the selectivity bundle as well as potency: with no --recipe,
        buildnew reads the recipe from the bundle matching --layout, so this
        needs the selectivity bundle installed even though it loads no weights.
        """
        out = tmp_path / "my_selectivity"
        run_tool("buildnew", "--layout", "selectivity",
                 "--data", _example("measurements_example.csv"),
                 "--out", str(out), "--trees", "8", "--allow-small")
        assert json.load(open(out / "MANIFEST.json"))["layout"] == "selectivity"

    def test_no_kfm_weights_are_claimed(self, tmp_path):
        """buildnew reads a bundle for the feature recipe and nothing else.

        `derived_from` is extend's record of a merge. Its presence here would
        mean released weights had reached a model advertised as free of them,
        which is a licensing claim as much as a technical one.
        """
        out = tmp_path / "clean"
        run_tool("buildnew", "--layout", "potency",
                 "--data", _example("measurements_example.csv"),
                 "--out", str(out), "--trees", "8", "--allow-small")
        assert "derived_from" not in json.load(open(out / "MANIFEST.json"))


@needs_potency
@pytest.mark.needs_bundle
class TestExtendEndToEnd:
    """`kfm extend` merges trees fitted on the user's data into a released model.

    Deliberately no assertion that the extended model still returns 0.847: it
    must not. Added trees change every prediction, which is exactly why the tool
    withdraws the released figures. Asserting the published number here would
    encode the opposite of what the tool guarantees.
    """

    @staticmethod
    def _extension(out):
        return json.load(open(os.path.join(out, "MANIFEST.json")))["extensions"][-1]

    def test_merges_and_records_what_it_did(self, tmp_path):
        out = tmp_path / "extended"
        run_tool("extend", "--model", "potency",
                 "--data", _example("extend_potency_example.csv"),
                 "--out", str(out), "--trees", "5", "--allow-small",
                 "--label", "test suite")
        rec = self._extension(out)
        assert rec["trees_added"] == 5
        assert rec["label"] == "test suite"
        assert rec["model_kind"] == "LSL"
        assert rec["usable_comparisons"] == 16, \
            "the shipped example is documented as 18 rows, 16 usable"

    def test_the_forest_actually_grew_by_the_trees_added(self, tmp_path):
        """A merge that silently added nothing would still write a manifest."""
        out = tmp_path / "extended"
        run_tool("extend", "--model", "potency",
                 "--data", _example("extend_potency_example.csv"),
                 "--out", str(out), "--trees", "5", "--allow-small")
        import joblib
        man = json.load(open(out / "MANIFEST.json"))
        merged = joblib.load(out / man["model_file"])
        assert merged.n_estimators == self._extension(out)["n_estimators_after"]
        assert len(merged.estimators_) == merged.n_estimators

    def test_released_figures_are_withdrawn(self, tmp_path):
        """The extended model is not the model the reports describe."""
        out = tmp_path / "extended"
        run_tool("extend", "--model", "potency",
                 "--data", _example("extend_potency_example.csv"),
                 "--out", str(out), "--trees", "5", "--allow-small")
        man = json.load(open(out / "MANIFEST.json"))
        assert "EXTENDED_MODEL_WARNING" in man
        assert man["derived_from"]["model_sha256"] != man["model_sha256"]
        assert not os.path.exists(out / "reference_predictions.json"), \
            "released reference predictions must not survive a merge"

    def test_the_model_filename_is_kept_so_predict_py_still_loads_it(self, tmp_path):
        out = tmp_path / "extended"
        run_tool("extend", "--model", "potency",
                 "--data", _example("extend_potency_example.csv"),
                 "--out", str(out), "--trees", "5", "--allow-small")
        # Resolve the base filename the way the tool does. A bundle manifest is
        # not required to carry `model_file`; Bundle falls back to finding the
        # single .joblib, and a test that assumed the key was always present
        # failed against a bundle that omits it.
        import kfm.extend as kx
        base_file = kx.Bundle(bundles.bundle_dir("potency"),
                              load_model=False, kind="potency").model_file
        assert json.load(open(out / "MANIFEST.json"))["model_file"] == base_file
        assert os.path.exists(out / base_file)


@needs_selectivity
@pytest.mark.needs_bundle
class TestShippedSelectivityExample:
    def test_preferred_kinase_matches_the_measurement(self):
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        f = os.path.join(repo, "examples", "selectivity_MTOR_PIK3CA_PIK3CG.smi")
        d = run_cli("selectivity", "-t", "MTOR", "-t", "PIK3CA", "-t", "PIK3CG",
                    "-l", "@" + f)
        for r in d["results"]:
            expected = r["name"].rsplit("_best_", 1)[1]
            assert r["prefers"] == expected, (
                f"{r['name']}: model said {r['prefers']}, measurement says {expected}")
