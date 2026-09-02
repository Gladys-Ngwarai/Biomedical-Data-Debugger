# Biomedical Data Debugger

## Artifact vs. True Biology

Biomedical datasets can contain technical artifacts that appear to be genuine biological signals.

**Biomedical Data Debugger** is an evidence-driven platform designed to help researchers investigate suspicious findings in biomedical data.

### Core Pipeline

```text
Upload biomedical data
        ↓
Quality control
        ↓
Detect suspicious signals
        ↓
Investigate possible artifacts
        ↓
Explain the evidence
        ↓
Recommend a correction
        ↓
Apply correction
        ↓
Re-analyze
        ↓
Compare before vs. after
```

### Initial Module — scRNA-seq

The first implementation focuses on single-cell RNA sequencing.

The debugger will investigate potential issues such as:

* Doublets
* Ambient RNA contamination
* Stress-related expression
* Batch effects
* Low-quality cells
* Suspicious cell populations
* Biological coherence

### Future Modules

The architecture will be designed to support additional biomedical data types:

* Whole-genome / exome sequencing
* Bulk RNA-seq
* Proteomics
* Other omics datasets

### Core Principle

> **Don't trust the biological conclusion until you debug the data.**

The system will not simply label a finding as "true" or "false."

Instead, it will evaluate the evidence supporting the finding, identify possible technical explanations, recommend appropriate corrections, and determine whether the biological conclusion remains supported after debugging.

---

## Project Status

🚧 Early development

**Current focus:** Foundation and scRNA-seq quality control.
