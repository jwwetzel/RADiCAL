#!/usr/bin/env python3
"""
gen_authors.py — RADiCAL authorship / collaboration pipeline.

SOURCE OF TRUTH is per test-beam era, hand-edited:
    data/<era>/metadata/authors.yaml     (one per campaign — 2023, 2024, ...)

This tool reads every era file and GENERATES (never hand-edit the outputs):
    collaboration/people.yaml                          aggregate roster; years_active COMPUTED
    data/<era>/metadata/generated/authors.tex          elsarticle \\author + \\affiliation block
    data/<era>/metadata/generated/zenodo_creators.json creators for that era's Zenodo deposit
    CITATION.cff                                        repo-level (from --citation-era, default latest)

WHY per-era source + generated rollup: authorship is time-varying (people join/leave,
affiliations change, not every member authors every paper). Keeping each era self-contained
matches the RAD_YEAR era-container model; ORCID is the identity join key, so the aggregate
roster and years_active can be COMPUTED without a second hand-maintained list drifting.

The tool also CHECKS consistency: the same ORCID must carry the same name across eras
(catches typos / disambiguation errors for free).

Usage (from repo root):
    python3 tools/gen_authors.py              # generate everything
    python3 tools/gen_authors.py --check      # validate only, write nothing (CI-friendly)
    python3 tools/gen_authors.py --era 2023   # regenerate one era's artifacts
    python3 tools/gen_authors.py --citation-era 2023

Dependencies: none (stdlib only). Uses PyYAML if importable; otherwise a vendored
strict-subset YAML loader that handles the authors.yaml grammar (block maps, block
sequences, inline [lists], scalars, null, comments) and errors clearly on anything else.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"

# Which era's roster populates the repo-level CITATION.cff. Pinned to the campaign
# whose analysis this repository currently publishes — NOT "the latest era", which
# would silently flip the citation to a future, unpublished roster. Override: --citation-era.
CITATION_ERA = "2023"


# ---------------------------------------------------------------------------
# YAML loading — PyYAML if present, else a vendored strict-subset loader
# ---------------------------------------------------------------------------
def load_yaml(text: str, where: str = "<yaml>"):
    try:
        import yaml  # type: ignore
        return yaml.safe_load(text)
    except ModuleNotFoundError:
        return _mini_load(text, where)


def _strip_comment(line: str) -> str:
    """Drop a trailing '# ...' comment not inside quotes."""
    out, q = [], None
    for ch in line:
        if q:
            out.append(ch)
            if ch == q:
                q = None
        elif ch in ('"', "'"):
            q = ch
            out.append(ch)
        elif ch == "#":
            break
        else:
            out.append(ch)
    return "".join(out).rstrip()


def _scalar(tok: str, where: str):
    tok = tok.strip()
    if tok == "" or tok in ("null", "~"):
        return None
    if tok == "[]":
        return []
    if tok.startswith("[") and tok.endswith("]"):
        inner = tok[1:-1].strip()
        if not inner:
            return []
        return [_scalar(p, where) for p in _split_commas(inner)]
    if tok in ("{", "|", ">") or tok.startswith("{") or tok.startswith("&") or tok.startswith("*"):
        raise ValueError(
            f"{where}: unsupported YAML construct {tok!r} in the vendored loader. "
            f"Install PyYAML (pip3 install pyyaml) to use flow mappings / anchors / block scalars."
        )
    if (tok[0] == '"' and tok[-1] == '"') or (tok[0] == "'" and tok[-1] == "'"):
        return tok[1:-1]
    if tok in ("true", "True"):
        return True
    if tok in ("false", "False"):
        return False
    try:
        return int(tok)
    except ValueError:
        return tok  # bare string (years like "2023" are quoted in source → stay str)


def _split_commas(s: str):
    parts, buf, q = [], [], None
    for ch in s:
        if q:
            buf.append(ch)
            if ch == q:
                q = None
        elif ch in ('"', "'"):
            q = ch
            buf.append(ch)
        elif ch == ",":
            parts.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        parts.append("".join(buf))
    return [p.strip() for p in parts]


def _mini_load(text: str, where: str):
    toks = []
    for raw in text.split("\n"):
        s = _strip_comment(raw)
        if s.strip() == "":
            continue
        indent = len(s) - len(s.lstrip(" "))
        toks.append((indent, s.strip()))
    idx = [0]

    def peek():
        return toks[idx[0]] if idx[0] < len(toks) else (None, None)

    def parse(cur_indent):
        i, content = peek()
        if i is None or i < cur_indent:
            return None
        if content.startswith("- "):
            return parse_seq(i)
        return parse_map(i)

    def parse_map(indent):
        d = {}
        while True:
            i, content = peek()
            if i is None or i != indent or content.startswith("- "):
                break
            key, sep, val = content.partition(":")
            if not sep:
                raise ValueError(f"{where}: expected 'key: value', got {content!r}")
            idx[0] += 1
            key, val = key.strip(), val.strip()
            d[key] = _scalar(val, where) if val != "" else parse(indent + 1)
        return d

    def parse_seq(indent):
        lst = []
        while True:
            i, content = peek()
            if i is None or i != indent or not content.startswith("- "):
                break
            item = content[2:].strip()
            idx[0] += 1
            # "- key: value ..." → a block map whose keys align at indent+2
            if (":" in item) and not item.startswith("[") and item[0] not in ('"', "'"):
                key, _, val = item.partition(":")
                m = {key.strip(): (_scalar(val.strip(), where) if val.strip() else parse(indent + 2))}
                while True:
                    i2, c2 = peek()
                    if i2 is None or i2 != indent + 2 or c2.startswith("- "):
                        break
                    k2, sep2, v2 = c2.partition(":")
                    if not sep2:
                        raise ValueError(f"{where}: expected 'key: value', got {c2!r}")
                    idx[0] += 1
                    m[k2.strip()] = _scalar(v2.strip(), where) if v2.strip() else parse(indent + 3)
                lst.append(m)
            else:
                lst.append(_scalar(item, where))
        return lst

    result = parse(0)
    if idx[0] != len(toks):
        i, content = peek()
        raise ValueError(f"{where}: could not parse near {content!r} (indentation?)")
    return result if result is not None else {}


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
def era_files():
    out = []
    for d in sorted(DATA.glob("[0-9][0-9][0-9][0-9]")):
        f = d / "metadata" / "authors.yaml"
        if f.exists():
            out.append((d.name, f))
    return out


def _split_name(name, family, given):
    if family or given:
        return (family or ""), (given or "")
    toks = name.split()
    if len(toks) == 1:
        return toks[0], ""
    return toks[-1], " ".join(toks[:-1])  # naive: last token = family


def load_era(era, path):
    doc = load_yaml(path.read_text(encoding="utf-8"), str(path))
    doc.setdefault("era", era)
    doc.setdefault("authors", [])
    doc.setdefault("affiliations", [])
    doc.setdefault("contributors", [])
    aff = {a["id"]: a for a in doc["affiliations"]}
    for au in doc["authors"]:
        for k in au.get("affil", []) or []:
            if k not in aff:
                raise ValueError(f"{path}: author {au.get('name')!r} references unknown affiliation id {k!r}")
    return doc


# ---------------------------------------------------------------------------
# Consistency check
# ---------------------------------------------------------------------------
def check(eras):
    problems, warns = [], []
    no_orcid = []  # (era, name)
    by_orcid = {}  # orcid -> (name, era) first seen
    for era, doc in eras:
        for au in doc["authors"]:
            name, orcid = au.get("name"), au.get("orcid")
            if not name:
                problems.append(f"{era}: an author entry has no name")
                continue
            if orcid:
                if not _valid_orcid(orcid):
                    problems.append(f"{era}: {name}: malformed ORCID {orcid!r} (want 0000-0000-0000-000X)")
                elif orcid in by_orcid and by_orcid[orcid][0] != name:
                    problems.append(
                        f"ORCID {orcid} is {by_orcid[orcid][0]!r} in {by_orcid[orcid][1]} "
                        f"but {name!r} in {era} — same person must have the same name."
                    )
                else:
                    by_orcid.setdefault(orcid, (name, era))
            else:
                no_orcid.append((era, name))
            fam, giv = _split_name(name, au.get("family"), au.get("given"))
            if not (au.get("family") or au.get("given")) and (len(name.split()) > 2 or "-" in name):
                warns.append(
                    f"{era}: {name}: auto-split to family={fam!r} given={giv!r} — "
                    f"add explicit family:/given: if that is wrong (compound surname?)."
                )
    if no_orcid:
        by_era = {}
        for era, name in no_orcid:
            by_era.setdefault(era, []).append(name)
        for era in sorted(by_era):
            names = by_era[era]
            warns.append(f"{era}: {len(names)} author(s) without ORCID (the identity join key): "
                         + ", ".join(names))
    return problems, warns


def _valid_orcid(o):
    import re
    return bool(re.fullmatch(r"\d{4}-\d{4}-\d{4}-\d{3}[\dX]", o or ""))


# ---------------------------------------------------------------------------
# Emitters
# ---------------------------------------------------------------------------
def _yq(s):  # YAML-quote a scalar
    if s is None:
        return "null"
    return '"' + str(s).replace('"', '\\"') + '"'


def gen_people_yaml(eras):
    """Aggregate roster; identity keyed by ORCID (or name slug if none). years_active COMPUTED."""
    people = {}  # key -> record
    order = []
    for era, doc in eras:
        for au in doc["authors"]:
            key = au.get("orcid") or ("name:" + au["name"])
            rec = people.get(key)
            if rec is None:
                rec = {"name": au["name"], "orcid": au.get("orcid"),
                       "years": [], "affiliations": {}, "roles": {}}
                people[key] = rec
                order.append(key)
            if era not in rec["years"]:
                rec["years"].append(era)
            affs = doc_aff_names(doc, au.get("affil", []) or [])
            if affs:
                rec["affiliations"][era] = affs
            if au.get("roles"):
                rec["roles"][era] = au["roles"]
    lines = [
        "# collaboration/people.yaml — GENERATED by tools/gen_authors.py. DO NOT EDIT.",
        "# Aggregate RADiCAL roster across all test-beam eras. Source of truth is",
        "# data/<era>/metadata/authors.yaml; years_active is COMPUTED (which eras list the ORCID).",
        "# Regenerate:  python3 tools/gen_authors.py",
        "people:",
    ]
    for key in sorted(order, key=lambda k: people[k]["name"].split()[-1].lower()):
        r = people[key]
        lines.append(f"  - name: {_yq(r['name'])}")
        lines.append(f"    orcid: {_yq(r['orcid']) if r['orcid'] else 'null'}")
        lines.append(f"    years_active: [{', '.join(_yq(y) for y in sorted(r['years']))}]")
        if r["affiliations"]:
            lines.append("    affiliations_by_era:")
            for era in sorted(r["affiliations"]):
                lines.append(f'      "{era}": [{", ".join(_yq(a) for a in r["affiliations"][era])}]')
    return "\n".join(lines) + "\n"


def doc_aff_names(doc, ids):
    aff = {a["id"]: a for a in doc["affiliations"]}
    return [aff[i]["name"] for i in ids if i in aff]


def gen_authors_tex(doc):
    """elsarticle \\author + \\affiliation block — a drop-in \\input for the manuscript.
    Built with f-strings (LaTeX is full of literal % and \\, so %-formatting is unsafe)."""
    used = []
    for au in doc["authors"]:
        for k in au.get("affil", []) or []:
            if k not in used:
                used.append(k)
    letter = {k: chr(ord("a") + i) for i, k in enumerate(used)}
    aff = {a["id"]: a for a in doc["affiliations"]}
    out = [f"% GENERATED by tools/gen_authors.py from data/{doc['era']}/metadata/authors.yaml — DO NOT EDIT.",
           "% Requires \\usepackage{orcidlink} and a \\cortext in the preamble.", ""]
    for au in doc["authors"]:
        keys = ",".join(letter[k] for k in (au.get("affil", []) or []))
        name = au["name"].replace(" ", "~")
        orc = ("\\," + "\\orcidlink{" + au["orcid"] + "}") if au.get("orcid") else ""
        cor = "\\corref{cor}" if "corresponding" in (au.get("roles", []) or []) else ""
        out.append(f"\\author[{keys}]{{{name}{orc}{cor}}}")
    out.append("\\cortext[cor]{Corresponding author.}")
    out.append("")
    for k in used:
        a = aff[k]
        bits = [f"organization={{{a['name']}}}"]
        if a.get("city"):
            bits.append(f"city={{{a['city']}}}")
        if a.get("region"):
            bits.append(f"state={{{a['region']}}}")
        if a.get("country"):
            bits.append(f"country={{{a['country']}}}")
        out.append(f"\\affiliation[{letter[k]}]{{{', '.join(bits)}}}")
    return "\n".join(out) + "\n"


def gen_zenodo_creators(doc):
    aff = {a["id"]: a for a in doc["affiliations"]}
    creators = []
    for au in doc["authors"]:
        fam, giv = _split_name(au["name"], au.get("family"), au.get("given"))
        entry = {"name": f"{fam}, {giv}".strip(", ")}
        ids = au.get("affil", []) or []
        if ids and ids[0] in aff:
            entry["affiliation"] = aff[ids[0]]["name"]
        if au.get("orcid"):
            entry["orcid"] = au["orcid"]
        creators.append(entry)
    return json.dumps({"creators": creators}, indent=2, ensure_ascii=False) + "\n"


def gen_citation_cff(doc, eras_all):
    fams = []
    for au in doc["authors"]:
        fam, giv = _split_name(au["name"], au.get("family"), au.get("given"))
        line = ["  - family-names: %s" % _yq(fam), "    given-names: %s" % _yq(giv)]
        if au.get("orcid"):
            line.append("    orcid: %s" % _yq("https://orcid.org/" + au["orcid"]))
        fams.append("\n".join(line))
    return (
        "# GENERATED by tools/gen_authors.py (authors = %s era). DO NOT EDIT.\n" % doc["era"]
        + "cff-version: 1.2.0\n"
        + "message: \"If you use this software or data, please cite it as below.\"\n"
        + "title: \"RADiCAL — precision-timing shashlik calorimeter analysis\"\n"
        + "type: software\n"
        + "repository-code: \"https://github.com/jwwetzel/RADiCAL\"\n"
        + "authors:\n" + "\n".join(fams) + "\n"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Generate RADiCAL authorship artifacts from per-era authors.yaml")
    ap.add_argument("--check", action="store_true", help="validate only; write nothing")
    ap.add_argument("--era", help="regenerate a single era's artifacts (e.g. 2023)")
    ap.add_argument("--citation-era", help="which era populates CITATION.cff (default: latest era with authors)")
    args = ap.parse_args()

    files = era_files()
    if not files:
        print("No data/<era>/metadata/authors.yaml files found.", file=sys.stderr)
        return 1
    eras = [(era, load_era(era, path)) for era, path in files]

    problems, warns = check(eras)
    for w in warns:
        print("  warn: " + w)
    if problems:
        print("\nCONSISTENCY ERRORS:", file=sys.stderr)
        for p in problems:
            print("  ✗ " + p, file=sys.stderr)
        return 2
    populated = [(e, d) for e, d in eras if d["authors"]]
    print("checked %d era file(s), %d with authors; %d warning(s)." % (len(eras), len(populated), len(warns)))
    if args.check:
        return 0

    # aggregate roster
    (REPO / "collaboration").mkdir(exist_ok=True)
    (REPO / "collaboration" / "people.yaml").write_text(gen_people_yaml(populated), encoding="utf-8")
    print("wrote collaboration/people.yaml")

    # per-era artifacts
    for era, doc in eras:
        if args.era and era != args.era:
            continue
        if not doc["authors"]:
            continue
        gdir = DATA / era / "metadata" / "generated"
        gdir.mkdir(exist_ok=True)
        (gdir / "authors.tex").write_text(gen_authors_tex(doc), encoding="utf-8")
        (gdir / "zenodo_creators.json").write_text(gen_zenodo_creators(doc), encoding="utf-8")
        print("wrote data/%s/metadata/generated/{authors.tex,zenodo_creators.json}" % era)

    # repo-level CITATION.cff
    cera = args.citation_era or CITATION_ERA
    if cera and cera in dict(eras):
        cdoc = dict(eras)[cera]
        (REPO / "CITATION.cff").write_text(gen_citation_cff(cdoc, eras), encoding="utf-8")
        print("wrote CITATION.cff (authors: %s era)" % cera)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
