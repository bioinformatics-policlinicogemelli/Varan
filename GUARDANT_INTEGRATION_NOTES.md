# Guardant input into Varan — status and how to compile `sample.tsv`

## Status: `create_Varan_input.py` reviewed and corrected against real data

You sent the real script plus a real example run
(`C:\Users\mikik\Desktop\Progetti\Varan\es. output run\`, per-sample folder
`20250630\` + the run-level files). Read it properly rather than integrating
it blind, per your instruction. The corrected version lives at
`create_Varan_input.py` in this branch's repo root - it's still a standalone
script outside Varan itself ("extra Varan"), not merged into Varan's own
codebase, since you said you'd rather have your own script compile
`sample.tsv` than do anything by hand.

**Real bugs found and fixed, each verified by re-running the corrected
function against the real example VCF/TSV files (not just read, actually
executed):**

1. **Germline contamination (the serious one).** `process_vcf()`'s
   SNV-writing branch wrote every non-structural VCF row to the output SNV
   VCF regardless of the `FILTER` column. In the real example VCF: 225 of
   343 non-structural rows are `FILTER=FAIL` with `GT=1/1` (homozygous) and
   a dbSNP rsID - classic germline common-SNP signatures - versus 118 real
   `FILTER=PASS` somatic calls (`GT=0/1`, low VAF, COSMIC IDs). Germline
   variants were outnumbering real somatic calls roughly 2:1 flowing
   unfiltered into Varan's SNV pipeline. Fixed: only `FILTER=PASS` rows are
   now written. Verified: 110 real point-mutation rows survive (down from
   225 total non-structural rows).
2. **A dead filter that never fired.** `load_cnv_tsv_ordered()`'s
   `if cn_value == 2.0: continue` was meant to drop copy-neutral genes, but
   `copy_number` is a continuous value (2.07, 1.84, 3.17, ...) that's
   essentially never exactly `2.0` - it excluded nothing in the real file.
   Harmless in practice only because `process_vcf()` already has a second,
   correct filter further down (gates PASS/FAIL on the real `call` column:
   0 = no significant call, 1/2 = deletion/amplification) - removed the
   dead line rather than "fixing" a redundant check.
3. **Fusion breakend rows nearly got misclassified as SNVs.** While fixing
   `is_struct`'s over-broad `"SVTYPE=" in info_str` check (see #4), a first
   pass narrowed it to `SVTYPE=CNV/DUP/DEL` specifically - which, on the
   real file, caused 4 `SVTYPE=BND` (fusion breakend) rows to fall through
   into the *SNV* branch instead (their own `FILTER=PASS` meant they
   weren't caught by fix #1 either). Caught this by re-running against the
   real file and comparing row counts before shipping it (114 written vs.
   110 expected) - added an explicit `SVTYPE=BND` skip. Verified: 0 BND
   rows leak into the SNV output now.
4. **Fragile gene\<->copy-number pairing (latent, not yet triggered).**
   `is_struct`'s original `"SVTYPE=" in info_str` check also matched
   `SVTYPE=BND`, meaning BND rows previously consumed a slot in
   `ordered_cnv_data` via the same sequential `cnv_index` counter used for
   real CNV rows. In the one real file checked this was harmless (all 4 BND
   rows happen to come *after* all 18 real CNV-track rows in the VCF), but
   it's an order-dependent design: if a future file interleaves BND rows
   among CNV rows, gene<->copy-number pairing would silently misattribute
   copy numbers to the wrong genes. Narrowed the match to
   `SVTYPE=CNV/DUP/DEL` (see #3's fix) removes this risk going forward.
5. **MSI: real numeric score was being thrown away.** `.msi_call.hdr.tsv`
   has an actual `msi_score` (27 in the real example) plus `msi_status`
   (`MSI-H` in the example) - the old code discarded both in favor of the
   binary `SUM_JSD=0/999999` encoding for Varan's now-removed workaround.
   Fixed: `MSI` gets the real numeric score, `MSI_THR` gets `msi_status`
   translated to Varan's own `Stable`/`Unstable` vocabulary. **Caveat**: the
   real example only had an `MSI-H` case - the exact stable-case string was
   never seen, so the status match checks substrings (`"MSI-H"`/`"MSS"`/
   `"MSI-L"`/`"STABLE"`/`"UNSTABLE"`, case-insensitive) rather than one
   hardcoded compound string, and leaves `MSI_THR` blank (does not guess)
   for any status string it doesn't recognize - flag the first real
   MSS/MSI-L example you get so this can be confirmed.
6. **Fusions no longer go through a synthesized `CombinedVariantOutput.tsv`
   at all - architectural fix, not just a bug fix.** The old
   `create_combined_file()` wrote a fake `[MSI]`/`[Fusions]` file per
   sample and pointed `comb_path` at it. Checked what that actually does in
   `walk.py`: `transform_input()` only populates the `CombinedOutput` temp
   folder for samples whose `comb_path` is non-empty, but
   `_walk_setup()`'s check for "does this run have CombinedOutput" is
   **folder-existence, not per-sample** - so as soon as *any* sample in a
   batch had a non-empty `comb_path` (needed for its fusions), the *whole
   run's* MSI/TMB handling would silently switch from the new
   sample.tsv-driven VALUE+THR precedence path (`fill_from_file`) to the
   native-TSO500-CombinedOutput path (`fill_from_combined`) - defeating the
   very mechanism built to fix MSI for Guardant. Fixed by routing fusions
   through Varan's *other*, independent fusion-ingestion path instead:
   a plain `FUSIONS/*.tsv` file in Varan's own `data_sv.txt` column shape
   (`Sample_Id`/`SV_Status`/`Class`/`Site1_Hugo_Symbol`/`Site2_Hugo_Symbol`/
   `Normal_Paired_End_Read_Count`/`Event_Info`/`RNA_Support`), read by
   `fill_fusion_from_temp()`, which never touches `CombinedOutput` at all.
   `comb_path` is now always blank. See
   `Templates/sample_guardant_example_fusions.tsv` for the real shape
   (built from the real example's 2 confirmed fusions, ROS1-CD74 and
   RET-NCOA4). Pass it as the third path to `varan.py -i sample.tsv
   patient.tsv fusions.tsv`.

**One more thing found, not yet fixed - needs your call:**
`fill_fusion_from_temp()` in `walk.py` applies its own **hardcoded**
`min_read_count = 15` cutoff, independent of `conf.ini`'s
`THRESHOLD_FUSION` (used by the CombinedOutput fusion path) and independent
of Guardant's own `call=1` confidence flag. In the real example, RET-NCOA4
is a real Guardant-confirmed fusion (`call=1`) with only 10 supporting
molecules - it would be **silently dropped** by this hardcoded threshold
even though Guardant's own algorithm already validated it. Whether 15 is
the right cutoff for Guardant's molecule-count scale (which isn't
necessarily comparable to Illumina's read-count scale) is a real open
question, not something to guess at - flagging rather than changing it
blindly.

## Remaining gap: CNV format

Guardant's `.cnv_call.hdr.tsv` (via the corrected `process_vcf()`) now
writes a real per-gene CNV VCF (`{sample}.cnv.vcf`, verified against the
real file: correctly identifies MET/MYC/ERBB2 as the 3 amplified genes out
of 18 tracked). `cnv_path` in the template now points at this - the
earlier "leave it blank" guidance was wrong; the VCF conversion turned out
to already work once the two structural-variant bugs above were fixed. Not
a gap anymore.

## How to compile `sample.tsv` for Guardant samples

Run the corrected `create_Varan_input.py` - it produces `sample.tsv` and a
companion `{run_id}_fusions.tsv` directly; see
`Templates/sample_guardant_example.tsv` and
`Templates/sample_guardant_example_fusions.tsv` for what the output looks
like (built from the real example sample `20250630_30ng_22`, plus a second
row using the redacted `52084217_00102139` for reference - fake paths, not
real patient data). Column by column:

| Column | Guardant source | Notes |
|---|---|---|
| `SAMPLE_ID` | `finalmetadata.xml` AccessionId / `_metadata.xml` AccessionId | |
| `PATIENT_ID` | `finalmetadata.xml` SubjectId (or AccessionId again for the RUO/QCI-only path, which has no SubjectId) | |
| `ONCOTREE_CODE` | `finalmetadata.xml` Diagnosis, mapped via `dict.csv` | unchanged from the original script's logic, not reviewed further this pass |
| `snv_path` | corrected script's `{sample}.snv.vcf` | now PASS-only, no germline contamination |
| `cnv_path` | corrected script's `{sample}.cnv.vcf` | now correctly gene-attributed, see gap section above |
| `comb_path` | — | always blank now - see fix #6 |
| `MSI` | `.msi_call.hdr.tsv` `msi_score` | real number now, not a placeholder |
| `TMB` | — | blank - Guardant360 CDx has no TMB field anywhere in either metadata XML variant checked. Left blank for every sample with no CombinedOutput fallback, the column disappears entirely from `data_clinical_sample.txt` rather than showing "NA" |
| `MSI_THR` | `.msi_call.hdr.tsv` `msi_status`, mapped to Stable/Unstable | kept as-is by Varan's VALUE+THR precedence (both VALUE and THR now populated - this is the "both present, trust THR, ignore conf.ini" case, logged as a warning by design) |
| `TMB_THR` | — | blank, same reasoning as `TMB` |

Fusions: separate `{run_id}_fusions.tsv`, passed as the third `-i` argument
to `varan.py`. See the `min_read_count=15` caveat above before trusting it
blindly for a real batch.

### "Golden" extra fields (BRCA Reversion, AutoQC metrics, etc.)

Already works with no new code needed: `write_clinical_sample()` keeps any
extra column present in `sample.tsv` beyond the ten required ones and
carries it straight through into `data_clinical_sample.txt`. So a
`BRCA_REVERSION` or `GERMLINE_PATHOGENIC_COUNT` column can just be added to
`sample.tsv` today - Illumina-sourced samples would leave it blank, Guardant
ones would fill it in. What's *not* built yet is the "drop the column
entirely if nobody ever fills it" behavior for these - that only exists for
MSI/TMB right now. Not yet wired into `create_Varan_input.py` either - the
real board-summary/AutoQC file (BRCA Reversion, germline pathogenic count,
etc.) wasn't among the files re-sent this turn; send it when ready and this
can be extended the same way MSI was.
