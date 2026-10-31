# VARAN

[![DOI](https://zenodo.org/badge/788270006.svg)](https://zenodo.org/doi/10.5281/zenodo.12806060)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-003366)](https://bioinformatics-policlinicogemelli.github.io/Varan/)

<p align="center">
<img src="docs/img/logo_VARAN.png" alt="Varan logo" style="height: 260px; width:260px;"/>
</p>

<p align="justify">
Varan turns raw variant-calling output (VCF, MAF, Illumina TSO500 <code>CombinedVariantOutput</code>/<code>MetricsOutput</code>, and more) into a fully structured, validated cancer study ready to load into <a href="https://www.cbioportal.org/">cBioPortal</a> - no manual file wrangling required. It also manages the full lifecycle of an existing study: versioning, merging new batches in, extracting or removing samples, and re-validating along the way.
</p>

📚 **[Full documentation](https://bioinformatics-policlinicogemelli.github.io/Varan/)** · 🚀 **[Quick Start](https://bioinformatics-policlinicogemelli.github.io/Varan/quickstart.html)** · 📖 **[User Guide](https://bioinformatics-policlinicogemelli.github.io/Varan/user_guide.html)**

## Why Varan

- **One command, a whole cBioPortal study.** Point Varan at your VCFs/MAFs (or a TSV of per-sample file paths) and get a validated, versioned study folder out - clinical/patient tables, mutations, CNA, structural variants, case lists and meta files all generated and cross-consistent.
- **Built for real lab workflows.** Input is a TSV listing each sample's file paths, not a pre-combined folder - because sequencing and analysis files usually can't be moved into one place. SNV, CNV and fusion oncogenicity filters are each independently configurable, so a lab can, say, keep CNV VUS while excluding SNV VUS.
- **Study lifecycle management.** Update a study with a new batch of samples, extract a subset into its own study, or remove samples - each operation produces a new, validated, version-tracked folder rather than mutating data in place.
- **TSO500-aware.** Reads TMB, MSI, fusions, HRD/Genomic Instability Score, and exon-level BRCA1/BRCA2 CNVs straight out of Illumina's `CombinedVariantOutput`/`MetricsOutput` files, alongside standard VCF/MAF input.
- **Validated, not just generated.** Every run is checked against cBioPortal's own offline validator, with an HTML report summarizing what passed, what needs review, and what's missing.
- **Reproducible by construction.** Deterministic, version-controlled output folders (`study_v1`, `study_v2`, ...) plus a per-run report documenting exactly which filters, thresholds and reference files were used.

## Quick Start

Varan ships as a Docker image, so there's nothing to install beyond Docker itself.

```bash
# See all available options
docker run --rm -it varan -h

# Build a study from raw input
docker run --rm -it \
  -v ./TEST:/test_folder \
  -v ./conf.ini:/conf.ini \
  varan -i test_folder/Input -o test_folder/Output/output_test -c mixed
```

Full walk-through with synthetic test data (raw files → validated study in under 30 minutes): see the **[Quick Start guide](https://bioinformatics-policlinicogemelli.github.io/Varan/quickstart.html)**.

Every setting - reference files, filter thresholds, OncoKB annotation, clinical header customization - lives in a single `conf.ini` file, which can be pointed at explicitly with `-C <path>` if you keep several around (e.g. per assay type or per environment).

### Workflow orchestration (experimental)

A [`Snakefile`](Snakefile) wraps Varan's four top-level operations (create/update/extract/remove) as Snakemake rules driven by [`config.yaml`](config.yaml), adding dependency tracking, per-run logs and conda-environment management on top of the existing CLI - without changing any pipeline logic. It's a starting point, not yet exercised against a real Snakemake install or cluster environment: review `config.yaml` and test before relying on it.

## Core Operations

| Operation | Flag | What it does |
|---|---|---|
| **Create** | `-i` | Build a new validated study folder from raw input files |
| **Update** | `-u` | Merge a new batch of samples into an existing study, producing a new version |
| **Extract** | `-e` | Pull a subset of samples out of a study into a new one |
| **Remove** | `-r` | Drop a subset of samples from a study, producing a new version without them |

Every operation validates its output against cBioPortal's own `validateData.py` and produces `report_VARAN.html` / `report_validate.html` summarizing the run.

### Tissue and plasma samples

`[Sample_Type] TYPE` in `conf.ini` sets the run's specimen type, `Tissue` or `Plasma` (the older `Solid`/`Liquid` still work). An optional `SPECIMEN_TYPE` column in `sample.tsv` overrides it per sample, so one study can hold tissue and plasma samples - e.g. a tissue run and a plasma run merged with `-u`. Every sample's type is written to `data_clinical_sample.txt`. Settings that differ between the two go in `[Filters.Plasma]`, `[MSI.Plasma]`, `[TMB.Plasma]`, `[SigMA.Plasma]`, `[GenePanel.Plasma]` (and `.Tissue`) sections, read before the general ones.

### Clonal hematopoiesis flags (plasma)

Plasma variants that may come from white blood cell clones rather than the tumor are flagged, never removed, in `data_mutations_extended.txt` (`CH.EVIDENCE`, `CH.REASON`, shown in cBioPortal's mutation tables) and counted per sample (`N_CH_SUSPECT`, `CH_SUSPECT_GENES`). Evidence comes from an optional matched white blood cell VCF (`wbc_path` column in `sample.tsv`), tissue samples of the same patient in the study, serial plasma samples, and the gene list in [`ch_genes.tsv`](ch_genes.tsv) - see `[CH]` in `conf.ini` and [`clonal_hematopoiesis.py`](clonal_hematopoiesis.py). Without white blood cell sequencing the origin can't be confirmed: research use only.

### Serial ctDNA

With an optional `COLLECTION_DATE` column in `sample.tsv` (and optionally `BASELINE`), the plasma samples of each patient are compared over time ([`serial_ctdna.py`](serial_ctdna.py), `[Serial]` in `conf.ini`): ctDNA level (max or mean VAF, clonal hematopoiesis and germline-like variants left out), change from baseline and molecular response (thresholds from the ctMoniTR analyses, configurable), new and lost variants, and known acquired-resistance mutations from [`resistance_mutations.tsv`](resistance_mutations.tsv) with their references. Samples are ordered in cBioPortal's patient view through `data_timeline_specimen.txt`; the dates themselves stay in the study's `service/` folder - cBioPortal only gets days from the patient's first sample. Research use only.

### Clinical timeline

An optional `treatments.tsv` (next to `patient.tsv` or `sample.tsv`, or `[Timeline] TREATMENTS` in `conf.ini`) lists each patient's treatments - agent, line and a therapy class from [`therapy_classes.tsv`](therapy_classes.tsv) - and status events (progression, RECIST response, death), with dates; `patient.tsv` can add a `DIAGNOSIS_DATE`. They become `data_timeline_treatment.txt` and `data_timeline_status.txt`, shown with the sample collections in cBioPortal's patient view, on one day 0 (the diagnosis when every patient has one, otherwise each patient's first event) - see [`clinical_timeline.py`](clinical_timeline.py). Unreadable dates, treatments ending before they start and values outside the allowed ones stop the run. With treatments, serial ctDNA compares each plasma sample with the last one before the treatment it is on, and flags a resistance mutation only if the patient received a treatment of its class. Dates stay in `service/`, as above.

### Actionability and knowledge bases

- **OncoKB actionability** (with `-k`): per sample highest therapeutic and resistance level and the alterations at the `[Actionability] ONCOKB_LEVELS` levels (default 1-2), over mutations, fusions and CNAs; the OncoKB data version is recorded and shown in the report, with a cohort summary by level and tumour type - [`actionability.py`](actionability.py).
- **CIViC** (on by default, no token): CIViC's open evidence (CC0) for every study, downloaded on first use (about 5 MB, seconds to match) and dated in the report. Kept in separate columns (`CIVIC.*` in the MAF, `CIVIC_*` per sample); levels are for any tumour type, each item names its disease - [`civic.py`](civic.py).
- **Tumour-agnostic biomarkers**: NTRK/RET fusions, BRAF V600E, MSI-H, TMB-H (tumour-agnostic approvals) kept apart from pan-tumour biomarkers (ALK, ROS1, FGFR, ERBB2, HRD), each Detected / Not detected / Not tested, with FDA and EMA approvals in [`tumor_agnostic.tsv`](tumor_agnostic.tsv) - [`tumor_agnostic.py`](tumor_agnostic.py).
- **ESCAT**: tiers from a table the lab curates (`[ESCAT] FILE`, format in [`escat_template.tsv`](escat_template.tsv)), matched to each sample's OncoTree code and its parents - [`escat.py`](escat.py).
- **Case lists** for OncoKB-actionable, CIViC A-B, tumour-agnostic biomarker, MSI-H, TMB-high, tissue and plasma samples.

Research use only.

## Documentation

| | |
|---|---|
| 🚀 [Quick Start](https://bioinformatics-policlinicogemelli.github.io/Varan/quickstart.html) | Get a study built from synthetic data in under 30 minutes |
| ⚙️ [Installation](https://bioinformatics-policlinicogemelli.github.io/Varan/installation_procedures.html) | Docker and manual setup instructions |
| 📖 [User Guide](https://bioinformatics-policlinicogemelli.github.io/Varan/user_guide.html) | Full reference: every flag, every `conf.ini` field, output structure |

## Issue Reporting

Found a bug or have a feature request? Please [open an issue on GitHub](https://github.com/bioinformatics-policlinicogemelli/Varan/issues) - include your `conf.ini` (redacted of any sensitive paths) and the relevant log from `Logs/` when reporting a run failure.

## Citation

If you use Varan in your research, please cite:

Parrillo C., Kulesko M., Persiani F., De Marco L., Petescia P., Mastrantoni L.,
Nero C., Minucci A., Giacò L.
**Varan: a tool for managing mutational data and creating cancer studies in cBioPortal.**
*NAR Genomics and Bioinformatics*, 2025.
https://doi.org/10.1093/nargab/lqaf107

<details>
<summary><strong>BibTeX</strong></summary>

```bibtex
@article{parrillo2025varan,
  title   = {Varan: a tool for managing mutational data and creating cancer studies in cBioPortal},
  author  = {Parrillo, C. and Kulesko, M. and Persiani, F. and De Marco, L. and
             Petescia, P. and Mastrantoni, L. and Nero, C. and Minucci, A. and Giacò, L.},
  journal = {NAR Genomics and Bioinformatics},
  year    = {2025},
  doi     = {10.1093/nargab/lqaf107}
}
```
</details>

## Disclaimer

⚠️ Varan is intended for research use only - not for patient treatment, diagnosis, and/or medical records.

## License

This project is distributed under the Apache License 2.0. See the [LICENSE](LICENSE) file for details.
