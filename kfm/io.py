"""Reading compounds and sequences from arguments or files.

One rule throughout: **a name is never invented and never truncated.** If the
input gave a name, that name comes back in the output; if it did not, the row is
labelled by position. Matching a result to the thing you submitted is the whole
job of these functions.
"""
from __future__ import annotations

import csv
import os
import re


def read_smiles(values) -> list[dict]:
    """Turn --ligand arguments into [{smiles, name}].

    Each value is either a SMILES, "SMILES name", or @file. A file may be
    .smi/.txt (whitespace or tab separated, one per line) or .csv (a header
    naming a smiles column, and optionally a name column).
    """
    out: list[dict] = []
    for v in values or []:
        if v.startswith("@"):
            out.extend(_read_file(v[1:]))
        else:
            out.append(_split_one(v))
    return out


def _split_one(text: str) -> dict:
    parts = re.split(r"[\s\t,]+", text.strip(), maxsplit=1)
    return {"smiles": parts[0],
            "name": parts[1].strip() if len(parts) > 1 and parts[1].strip() else None}


def _read_file(path: str) -> list[dict]:
    path = os.path.expanduser(path)
    if not os.path.isfile(path):
        raise SystemExit(f"No such file: {path}")
    rows: list[dict] = []
    if path.lower().endswith(".csv"):
        with open(path, newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh)
            cols = {c.lower().strip(): c for c in (rdr.fieldnames or [])}
            smi_col = next((cols[c] for c in ("smiles", "smile", "structure")
                            if c in cols), None)
            name_col = next((cols[c] for c in ("name", "id", "compound", "title")
                             if c in cols), None)
            if not smi_col:
                raise SystemExit(
                    f"{path}: no SMILES column found. Expected a header naming "
                    f"one of: smiles, smile, structure.")
            for rec in rdr:
                smi = (rec.get(smi_col) or "").strip()
                if smi:
                    rows.append({"smiles": smi,
                                 "name": (rec.get(name_col) or "").strip() or None
                                 if name_col else None})
        return rows
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            # '#' starts a comment, but only at the start of a line: it is a
            # legal SMILES character elsewhere and stripping it mid-string would
            # silently corrupt a structure.
            if not line or line.startswith("#"):
                continue
            rows.append(_split_one(line))
    return rows


def read_sequence(value: str) -> str:
    """A sequence from a literal, or @file for FASTA/plain text.

    Normalisation happens in the model's own predict module, which is where it
    must happen, but FASTA headers are stripped here because a '>' line is a
    file-format artefact rather than part of the sequence.
    """
    if value.startswith("@"):
        path = os.path.expanduser(value[1:])
        if not os.path.isfile(path):
            raise SystemExit(f"No such file: {path}")
        with open(path, encoding="utf-8") as fh:
            value = "".join(l for l in fh if not l.startswith(">"))
    return re.sub(r"[^A-Za-z]", "", value).upper()


def fmt_table(headers, rows, aligns=None) -> str:
    """A fixed-width table that lines up in a terminal.

    Nothing is truncated. A long SMILES makes a wide table; that is preferable
    to a table that cannot be pasted back into anything.
    """
    cols = len(headers)
    aligns = aligns or ["<"] * cols
    widths = [len(str(h)) for h in headers]
    for r in rows:
        for i in range(cols):
            widths[i] = max(widths[i], len(str(r[i])))
    line = "  ".join(f"{str(headers[i]):{aligns[i]}{widths[i]}}" for i in range(cols))
    out = [line, "  ".join("-" * widths[i] for i in range(cols))]
    for r in rows:
        out.append("  ".join(
            f"{str(r[i]):{aligns[i]}{widths[i]}}" for i in range(cols)))
    return "\n".join(out)
