# Validation record

Recorded on **2026-10-11**. The [JSON record](../validation/local-run.json) supplies the counts and environment used for the chart in the README.

| Check group | Passed | Scope |
|---|---:|---|
| Sequence evidence | 34 | FASTA validation, repeat agreement, deduplication, BLAST parsing, strict-hit rules, coordinates, score joins and the offline data flow |
| BLAST submission client | 8 | Batch limits, encoded requests, RID handling, mocked waiting/ready responses and retrieval pacing |
| Reference statistics | 6 | Sample requirements, PCA/GMM scores, covariance regularization and invalid input handling |
| Total regression checks | 48 | All passed; no skipped cases in this run |

Local environment: Python **3.11.5**, NumPy **2.4.6**, scikit-learn **1.9.1**, Windows. GitHub Actions repeats the suite on Ubuntu with Python 3.11 and the declared statistics dependencies. Its current result is linked by the README badge.

## Reproducible fixture output

The sequence example completed and matched the checked-in reference summary, FASTA files, provenance, alignment table and candidate ranking. SHA-256 checks verified all **8 input fixture files**. The checked-in HTML report is a human-readable view of that output.

| Synthetic fixture | Result |
|---|---|
| Spacer processing | 3 inputs become 7 unique spacers |
| Evidence decisions | 6 records: 3 strict, 1 illustrative-score-supported, 2 unresolved |
| Candidate aggregation | 2 candidates have accepted spacer support |
| Statistical example | 18 generated reference embeddings, 4 dimensions and 2 query embeddings |

The main sequence example reads manually assigned illustrative scores. The statistical example separately runs PCA/GMM on generated numeric vectors. Neither example executes DNABERT-S weights or measures biological prediction accuracy.

## Execution scope

| Component | What was exercised |
|---|---|
| Sequence data flow | Actual Python execution with labelled synthetic fixtures |
| PCA/GMM scoring | Actual NumPy/scikit-learn execution with generated embeddings |
| NCBI upload client | Injected mock transport; no live BLAST queue submission |
| Local BLAST wrapper | Mocked subprocess; installed BLAST+ was not executed |
| DNABERT-S adapter | Included source; real-weight inference was not run |
| Complete biological pipeline | Requires genome inputs, annotation tools and appropriate databases; not executed in this validation |

## Reproduce

```bash
python -m pip install -e ".[statistics]"
python -m unittest discover -s tests -v
python -m crisprphibert --demo --outdir demo-output
python scripts/verify_demo.py --output demo-output
python scripts/score_embeddings.py --output demo-output/statistics.json
```

The output directory must be new. The first 42 sequence/client checks can run without optional libraries; the six statistical cases are skipped if NumPy or scikit-learn is absent.

本记录统计工程检查及合成示例的运行结果。真实 DNABERT-S 推理、在线 BLAST 提交和有标签的生物学评估不在本次执行范围内。
