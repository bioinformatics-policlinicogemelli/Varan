# Guardant input into Varan — status and how to compile `sample.tsv`

## Honest status: what's actually integrated vs. not

**Integrated** (branch `stabilize`, 2026-07-19): the underlying mechanism
Guardant needs — the MSI/TMB VALUE+THR precedence system in `walk.py`
(`_resolve_biomarker_threshold`, used by `fill_from_file`) and the
"drop the column entirely if never populated anywhere" behavior. This
means Guardant no longer needs `create_Varan_input.py`'s `SUM_JSD=0/999999`
trick, or a synthesized fake `CombinedVariantOutput.tsv` at all — you can
fill `sample.tsv`'s own `MSI`/`MSI_THR`/`TMB`/`TMB_THR` columns directly.

**Not integrated**: `create_Varan_input.py`'s own file-parsing logic (reading
`.msi_call.hdr.tsv`, `.cnv_call.hdr.tsv`, `_finalmetadata.xml`, etc. and
producing `sample.tsv` automatically) was never ported into Varan as a
module or script. The original script's content isn't available anymore in
this session (it was pasted as chat text a couple of turns ago and wasn't
re-attached in the later, shortened message, so it was lost across a
context handoff) - if you still want an automated Guardant-BIP-output ->
`sample.tsv` converter shipped as part of Varan (vs. compiling `sample.tsv`
by hand per batch), resend `create_Varan_input.py` and the example BIP
output files and it can be built properly against the real file shapes,
the same way the splice-variant format was verified against real files
rather than guessed.

**Gap found while checking this**: Guardant's `.cnv_call.hdr.tsv` is a
plain TSV, not a VCF - Varan's `cnv_path` column (`get_cnv_from_folder` in
walk.py) requires a `.vcf` file. There's no converter for this either. The
template below leaves `cnv_path` blank for Guardant samples - CNV from
ctDNA is lower-confidence anyway given typical tumor fraction in liquid
biopsy, so deprioritizing it isn't purely a shortcut, but flag if you
actually need Guardant CNV calls in the study and it can be prioritized.

## How to compile `sample.tsv` for Guardant samples today

See `Templates/sample_guardant_example.tsv` for a filled-in example (fake
IDs, not real patient data). Column by column:

| Column | Guardant source | Notes |
|---|---|---|
| `SAMPLE_ID` | `finalmetadata.xml` AccessionId / CustomReportFields.sampleID | |
| `PATIENT_ID` | `finalmetadata.xml` SubjectId | |
| `ONCOTREE_CODE` | `finalmetadata.xml` Diagnosis / PrimarySourceTissue | mapped by hand to an OncoTree code, same as any other pipeline |
| `snv_path` | Guardant's own `.vcf` | ingested directly, same SNV pipeline as Illumina |
| `cnv_path` | — | leave blank, see gap above |
| `comb_path` | — | leave blank - no need to synthesize a fake CombinedVariantOutput.tsv anymore |
| `MSI` | — | leave blank - Guardant's MSI call (`finalmetadata.xml` Specimen/MicrosatelliteInstability, e.g. "MS-stable") is text, not a number, and `MSI` is a NUMBER-typed clinical column - putting text there breaks cBioPortal validation |
| `TMB` | — | leave blank - Guardant360 CDx doesn't appear to report TMB (no TMB field anywhere in `finalmetadata.xml`; the CDx label is DNA alteration + MSI, not TMB - bTMB is a GuardantOMNI/research-only metric on a different product). Left blank for every sample in a run, this column (and `TMB_THR`) will be dropped entirely from `data_clinical_sample.txt` rather than filled with "NA" |
| `MSI_THR` | `finalmetadata.xml` Specimen/MicrosatelliteInstability | translate Guardant's own call to Varan's vocabulary: "MS-stable" -> `Stable`, otherwise `Unstable`. This is trusted as-is (not recomputed from conf.ini) since `MSI` is blank and `MSI_THR` is filled - the "value absent but THR present, keep it" case |
| `TMB_THR` | — | leave blank, same reasoning as `TMB` |

### "Golden" extra fields (BRCA Reversion, AutoQC metrics, etc.)

Already works with no new code needed: `write_clinical_sample()` keeps any
extra column present in `sample.tsv` beyond the ten required ones and
carries it straight through into `data_clinical_sample.txt`. So a
`BRCA_REVERSION` or `GERMLINE_PATHOGENIC_COUNT` column can just be added to
`sample.tsv` today - Illumina-sourced samples would leave it blank, Guardant
ones would fill it in. What's *not* built yet is the "drop the column
entirely if nobody ever fills it" behavior for these - that only exists for
MSI/TMB right now. Worth generalizing later if a run with no Guardant
samples starts shipping a `data_clinical_sample.txt` full of blank golden
columns, but not urgent.
