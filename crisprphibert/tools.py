"""Run external alignment with explicit arguments and checked exit status."""
from pathlib import Path
import subprocess


def local_blast(query, database, output, relaxed=False):
    query = Path(query).absolute()
    if not query.is_file():
        raise FileNotFoundError(query)
    output = Path(output).absolute()
    if query == output:
        raise ValueError("Output must differ from input")
    if output.exists():
        raise FileExistsError(output)
    command = [
        "blastn",
        "-query",
        str(query),
        "-db",
        str(Path(database).resolve()),
        "-out",
        str(output),
        "-outfmt",
        "6",
    ]
    if relaxed:
        command.extend(["-task", "blastn-short", "-word_size", "7", "-dust", "no"])
    subprocess.run(command, check=True, shell=False, timeout=3600)
    if not output.is_file():
        raise RuntimeError("BLAST did not create its output")
    return output
