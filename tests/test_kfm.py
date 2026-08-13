"""Tests for the KFM command line tools.

Two tiers, and the split matters:

  * Tests that need NO model bundle — argument parsing, file reading, the
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


# ---------------------------------------------------------------------------
# Input parsing — no bundle needed
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
# Confidence bands — no bundle needed
# ---------------------------------------------------------------------------
class TestBands:
    def test_potency_bands(self):
        assert band_accuracy("potency", 0.95) == "87.2%"
        assert band_accuracy("potency", 0.85) == "85.0%"
        assert band_accuracy("potency", 0.55) == "58.4%"

    def test_selectivity_bands(self):
        assert band_accuracy("selectivity", 0.95) == "99.3%"
        assert band_accuracy("selectivity", 0.85) == "96.4%"
        assert band_accuracy("selectivity", 0.55) == "57.9%"

    def test_the_two_models_are_not_interchangeable(self):
        # Measured on different test sets. Copying one model's bands to the
        # other would overstate potency by eleven points in this band.
        assert band_accuracy("potency", 0.85) != band_accuracy("selectivity", 0.85)

    def test_boundaries_land_in_the_higher_band(self):
        assert band_accuracy("potency", 0.70) == "83.2%"
        assert band_accuracy("potency", 0.699) == "74.6%"


# ---------------------------------------------------------------------------
# Table formatting — no bundle needed
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
# Bundle resolution and the licence gate — no bundle needed
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
# The published numbers — these need the weights
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
    def test_worked_example_from_the_report(self):
        """ABL1, bosutinib vs the PP1-type compound.

        The report prints 0.846 and 0.151 — its SINGLE-ORDER figures, which sum
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
        against /tmp. Comments mentioning it are fine — only real directives are
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
