"""Validated sequence and alignment operations derived from the research workflow.

BLAST coordinates are 1-based and inclusive. Every alignment keeps its own ID
so scores for different locations on one subject cannot overwrite each other.
"""
from dataclasses import dataclass
from pathlib import Path
import math

DNA = set("ACGTRYSWKMBDHVN")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def validate_dna(sequence):
    value = sequence.upper()
    if not value or not set(value) <= DNA:
        raise ValueError("Expected a non-empty IUPAC DNA sequence")
    return value


def reverse_complement(sequence):
    return validate_dna(sequence).translate(COMPLEMENT)[::-1]


def read_fasta(path):
    records = {}
    name, chunks = None, []

    def flush():
        if name is not None:
            if name in records:
                raise ValueError(f"Duplicate FASTA ID: {name}")
            records[name] = validate_dna("".join(chunks))

    for number, raw in enumerate(
        Path(path).read_text(encoding="utf-8").splitlines(), 1
    ):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            parts = line[1:].split()
            if not parts:
                raise ValueError(f"Empty FASTA header at line {number}")
            name, chunks = parts[0], []
        else:
            if name is None:
                raise ValueError(f"Sequence before FASTA header at line {number}")
            chunks.append(line)
    flush()
    return records


def write_fasta(path, records):
    Path(path).write_text(
        "".join(f">{name}\n{seq}\n" for name, seq in records.items()), encoding="utf-8"
    )


def repeat_identity(left, right):
    """Global alignment identity; matches divided by aligned length.

    This compact implementation uses match=2, mismatch=-1, linear gap=-2.
    The archived research version used Biopython's affine-gap globalms.
    """
    left, right = validate_dna(left), validate_dna(right)
    n, m = len(left), len(right)
    scores = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        scores[i][0] = -2 * i
    for j in range(m + 1):
        scores[0][j] = -2 * j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            scores[i][j] = max(
                scores[i - 1][j - 1] + (2 if left[i - 1] == right[j - 1] else -1),
                scores[i - 1][j] - 2,
                scores[i][j - 1] - 2,
            )
    i, j, matches, length = n, m, 0, 0
    while i or j:
        length += 1
        if (
            i
            and j
            and scores[i][j]
            == scores[i - 1][j - 1] + (2 if left[i - 1] == right[j - 1] else -1)
        ):
            matches += left[i - 1] == right[j - 1]
            i -= 1
            j -= 1
        elif i and scores[i][j] == scores[i - 1][j] - 2:
            i -= 1
        else:
            j -= 1
    return matches / length


def select_homologs(path, coverage=0.9, identity=0.9):
    if not 0 <= coverage <= 1 or not 0 <= identity <= 1:
        raise ValueError("Protein thresholds must be fractions in [0, 1]")
    selected = {}
    for number, line in enumerate(
        Path(path).read_text(encoding="utf-8").splitlines(), 1
    ):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 13:
            raise ValueError(f"Expected 13 protein BLAST columns at line {number}")
        pident, evalue, qcovs = float(fields[2]), float(fields[10]), float(fields[12])
        if (
            not all(math.isfinite(x) for x in (pident, evalue, qcovs))
            or not 0 <= pident <= 100
            or not 0 <= qcovs <= 100
            or evalue < 0
        ):
            raise ValueError(f"Invalid protein BLAST values at line {number}")
        if evalue <= 1e-5 and qcovs / 100 >= coverage and pident / 100 >= identity:
            level = (
                "LEVEL-2"
                if evalue <= 1e-30 and qcovs >= 98 and pident >= 98
                else "LEVEL-3"
            )
            if fields[1] not in selected or level == "LEVEL-2":
                selected[fields[1]] = level
    return selected


def expand_spacers(original, homologs, database, repeat, threshold=0.8):
    if not 0 <= threshold <= 1:
        raise ValueError("Repeat threshold must be in [0, 1]")
    expanded = {name: validate_dna(seq) for name, seq in original.items()}
    seen = set(expanded.values())
    provenance = {name: {"source": "input", "level": "LEVEL-1"} for name in expanded}
    for protein, level in sorted(homologs.items()):
        if protein not in database:
            raise ValueError(f"Missing homolog metadata: {protein}")
        record = database[protein]
        forward = repeat_identity(record["repeat"], repeat)
        backward = repeat_identity(reverse_complement(record["repeat"]), repeat)
        if max(forward, backward) < threshold:
            continue
        reverse = backward > forward
        for index, seq in enumerate(record["spacers"], 1):
            value = reverse_complement(seq) if reverse else validate_dna(seq)
            if value in seen:
                continue
            name = f"{record['array_id']}_{index}"
            if name in expanded:
                raise ValueError(f"Spacer ID collision: {name}")
            seen.add(value)
            expanded[name] = value
            provenance[name] = {
                "source": protein,
                "level": level,
                "reverse_complemented": reverse,
            }
    return expanded, provenance


@dataclass(frozen=True)
class BlastHit:
    hit_id: str
    query: str
    subject: str
    pident: float
    length: int
    mismatches: int
    gaps: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float


def read_blast(path, lengths):
    hits, seen = [], set()
    for number, line in enumerate(
        Path(path).read_text(encoding="utf-8").splitlines(), 1
    ):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) != 12:
            raise ValueError(f"Expected 12 BLAST columns at line {number}")
        if f[0] not in lengths:
            raise ValueError(f"Unknown query ID: {f[0]}")
        hit = BlastHit(
            f"hit_{number:04d}",
            f[0],
            f[1],
            float(f[2]),
            *[int(x) for x in f[3:10]],
            float(f[10]),
            float(f[11]),
        )
        if not all(math.isfinite(x) for x in (hit.pident, hit.evalue, hit.bitscore)):
            raise ValueError(f"Non-finite BLAST values at line {number}")
        if (
            not 0 <= hit.pident <= 100
            or hit.evalue < 0
            or hit.bitscore < 0
            or min(hit.length, hit.qstart, hit.qend, hit.sstart, hit.send) < 1
        ):
            raise ValueError(f"Invalid BLAST values at line {number}")
        if (
            min(hit.mismatches, hit.gaps) < 0
            or max(hit.qstart, hit.qend) > lengths[hit.query]
        ):
            raise ValueError(f"Invalid query coordinates at line {number}")
        if hit.mismatches > hit.length or hit.gaps > hit.length:
            raise ValueError(f"Invalid alignment counts at line {number}")
        key = tuple(f)
        if key not in seen:
            hits.append(hit)
            seen.add(key)
    return hits


def is_significant(hit, qlen):
    start, end = sorted((hit.qstart, hit.qend))
    full = start == 1 and end == qlen and hit.mismatches <= 3
    one_missing = (
        end - start + 1 == qlen - 1
        and (start, end) in {(2, qlen), (1, qlen - 1)}
        and hit.mismatches <= 2
    )
    return hit.evalue <= 1 and hit.gaps == 0 and (full or one_missing)


def flank_sequence(hit, genome, flank=400):
    if flank < 0:
        raise ValueError("Flank length must be non-negative")
    sequence = validate_dna(genome)
    start, end = sorted((hit.sstart, hit.send))
    if end > len(sequence):
        raise ValueError("Subject coordinates exceed genome length")
    value = sequence[max(0, start - 1 - flank) : min(len(sequence), end + flank)]
    return reverse_complement(value) if hit.sstart > hit.send else value


def merge_evidence(hits, lengths, scores, threshold=0.9):
    if not 0 <= threshold <= 1:
        raise ValueError("Score threshold must be in [0, 1]")
    unknown = set(scores) - {h.hit_id for h in hits}
    if unknown:
        raise ValueError(f"Scores have unknown alignment IDs: {sorted(unknown)}")
    for key, value in scores.items():
        if (
            not isinstance(value, (int, float))
            or not math.isfinite(value)
            or not 0 <= value <= 1
        ):
            raise ValueError(f"Invalid compatibility score: {key}")
    result = []
    for hit in hits:
        strict = is_significant(hit, lengths[hit.query])
        score = scores.get(hit.hit_id)
        evidence = (
            "strict_alignment"
            if strict
            else "embedding_supported"
            if score is not None and score >= threshold
            else "unresolved"
        )
        result.append(
            {
                "hit_id": hit.hit_id,
                "query_id": hit.query,
                "subject_id": hit.subject,
                "evidence": evidence,
                "compatibility_score": score,
                "accepted": evidence != "unresolved",
            }
        )
    return result
