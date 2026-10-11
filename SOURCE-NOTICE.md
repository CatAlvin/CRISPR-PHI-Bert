# Sources and attribution / 来源与归属

## CRISPRPHIBert team

CRISPRPHIBert is a team research project. This public implementation was prepared from archived source and technical documentation in October 2026. It retains the sequence-evidence workflow and adds a compact executable example and regression suite. The [change record](docs/CHANGES.md) describes differences from the archived implementation.

Bohan Wu / 吴博涵 participated in **data preprocessing, BLAST upload workflows and part of the data-flow integration**. The project was supervised by **Dr. Dan Wang / 王丹老师**. Team members are **祁琛、陈乐怡、黄子倩、吴博涵、刘昇瑞**.

本项目为团队研究成果。吴博涵参与负责数据预处理、BLAST 上传流程与部分数据流集成。公开版本根据归档源码和技术文档整理；展示版中的代码调整和回归检查属于本次整理工作。

## Software registration

The contributor supplied registration **2025SR1133680**, acquisition method **原始取得** and rights scope **全部权利**. This is a contributor-provided record; the certificate and registry were not independently verified for this release. The record is not an open-source license grant.

参与者提供的软件登记号为 **2025SR1133680**，权利取得方式为“原始取得”，权利范围为“全部权利”。本次发布未独立核验登记证书或登记系统，也未据此另行授予开源许可。

## Related research and upstream tools

| Component | Source and credit |
|---|---|
| Cas/repeat-guided spacer expansion | CRISPRinter and [PAMPHLET research](https://www.sciencedirect.com/science/article/pii/S1673852724002881). The archived CRISPRinter source credits Yishu Ma, Xintong Huang and Sijia Ma, supervised by Dan Wang. |
| DNA embeddings | [DNABERT-S official implementation](https://github.com/MAGICS-LAB/DNABERT_S), [model](https://huggingface.co/zhihan1996/DNABERT-S) and [paper](https://arxiv.org/abs/2402.08777), by its upstream authors. |
| CRISPR-Cas annotation | [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper) and [MinCED](https://github.com/ctSkennerton/minced). |
| Protein prediction | [Prodigal](https://github.com/hyattpd/Prodigal). |
| Sequence search | [NCBI BLAST API](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html) and [usage guidance](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html). |

Pretrained models, upstream software, databases and publications retain their own authorship and terms. The CRISPRinter experiments belong to that historical work; this repository reports engineering checks on its own labelled synthetic fixtures.

预训练模型、上游软件、数据库与论文保留各自作者归属和使用条款。CRISPRinter 的历史实验属于该项研究。本仓库的验证记录基于明确标注的合成示例，不将历史研究数据计为本项目或个人的实验成果。
