# Public implementation change record

This version was prepared in October 2026 from the archived research implementation. Bohan Wu's research contribution is described in the README. The changes below describe preparation of the public version and its current verification suite.

| Finding in the archived source | Display treatment |
|---|---|
| Source code stored in a Word `.doc` file | Recovered the entry point and six modules; normalized non-breaking spaces and removed Word field/footer artifacts before preparing the public implementation. |
| Local BLAST wrote only to `outfile.temp` while its caller expected `outfile` | The display wrapper writes to the requested output and checks the external command's exit status. |
| Terminal exception checked start/end with `or`, could drop records or admit wider truncation | Exact one-base terminal span is required; every parsed record receives an evidence state. |
| Flank code mixed 1-based and Python coordinates, permitted negative slices and lost reverse-strand information | Coordinates are normalized explicitly, clipped at both boundaries, and reverse complement uses the original strand. |
| Fragment IDs included coordinates but merge compared scores with subject IDs | A unique alignment ID is carried through extraction and scoring, allowing multiple locations on the same subject. |
| Strict results could be skipped when either reference/query FASTA was empty | Strict alignments are retained independently of supplemental model scores. |
| Protein `qcovs` used percentage values against fractional thresholds | Converted percentages to fractions before comparison. |
| Broad exceptions could turn failures into empty predictions | Inputs and failures produce specific errors; the demo records provenance. |
| Encoder pooled padding tokens, lacked `eval()` and assumed named outputs | Added attention-mask-aware pooling, inference mode, CPU transfer and tuple-output handling. |
| Covariance inversion and scoring could be numerically unstable | Regularized covariance, pseudo-inverse, seeded GMM and a stable sigmoid. |
| Shell string construction and raw pickle dependencies | The public code uses argument arrays for commands and JSON fixtures for metadata. |
| Output comparison across operating systems | Text fixtures and generated FASTA/JSON/TSV outputs use LF consistently. A verification script checks the reference output and raw input SHA-256 hashes; GitHub Actions repeats that check on Linux. |

The compact repeat-alignment implementation uses a linear-gap global alignment; the archived source used Biopython `pairwise2.globalms` with affine-gap scoring. The display implementation is conceptually similar and has its own test coverage; it is not a byte-for-byte scientific reproduction of that alignment step.
