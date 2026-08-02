# SigMA integration — feasibility and cost/benefit analysis

## Status: research/analysis only, no pipeline wiring

This is a feasibility report, not a wiring change. `walk_folder()`, `varan.py`,
and the `Snakefile` are untouched. What *did* change this round: two draft
files (`sigma_filter.py`, `sigma_cancer_type_map.py`) got corrected in place
where research turned up a concrete, verifiable error — see "Fixes made this
round" at the end. Everything else here is intended to let a human decide
whether/how to proceed, per the branch's existing draft-then-decide pattern
(same style as `GUARDANT_INTEGRATION_NOTES.md` on `stabilize-guardant`).

Sources actually read for this pass (not assumed): the SigMA GitHub repo
(README, `R/run.R`, `R/list_tumor_types.R` lookup attempt), its rendered
Roxygen docs on rdrr.io (`run()` man page), its wiki (Parameter choices,
Quick start, FAQs, MMRD input file format), the original Gulhan et al. 2019
*Nature Genetics* paper ("Detecting the mutational signature of homologous
recombination deficiency in clinical samples"), and the follow-up clinical
validation paper Zhao et al. 2022 *Clinical Cancer Research*
("Mutational Signature 3 Detected from Clinical Panel Sequencing Is
Associated with Responses to Olaparib in Breast and Ovarian Cancers",
PMC9623231) which is the closest published analogue to Varan's actual
use case (MSK-IMPACT-style panel, tumor-only, oncology clinic). Read
2026-08-02.

---

## 1. What SigMA actually is

SigMA ("Signature Multivariate Analysis", parklab/Park lab, Harvard) is an R
package that infers a tumor's mutational-signature profile from its SNV
trinucleotide-context spectrum and, as its flagship use case, produces a
score for whether **COSMIC Signature 3** (associated with homologous
recombination deficiency / HRD, historically linked to BRCA1/2 and other
HR-pathway gene defects) is present — specifically tuned to work with the
low SNV counts that targeted gene panels produce, where classic
signature-decomposition methods (built for WGS-scale mutation counts) fail.

Core method (from the paper and `R/run.R`): for each candidate tumor
subtype cluster, SigMA computes (1) a **likelihood** that the sample's
96-dimensional trinucleotide spectrum matches that cluster's known
signature composition, (2) a **cosine similarity** between the sample's
spectrum and the cluster's reference spectrum, and (3) an **NNLS
(non-negative least squares) exposure** of Signature 3 fit directly to the
sample. These three feature families are then combined by a **pre-trained
gradient boosting model (GBM)** per tumor type into one final score. This
is a classifier applied at inference time, not a model that trains itself
per run — so per-sample cost should be dominated by evaluating a small
pre-fit model against a 96-length vector, not by any heavy optimization.
(See §4, runtime — this is inference from the method's design, not a
documented benchmark; nothing published states run time directly.)

Published performance (Gulhan et al. 2019): 74% sensitivity at a 10% false
positive rate for detecting Signature-3-positive tumors from low-SNV-count
panel data — the number the SigMA project itself cites as its headline
result. Applied to 878 breast tumors from the real MSK-IMPACT panel, SigMA
called 202 (23%) likely Sig3+, 121 (14%) at a stricter cutoff. Downstream
biological validation in that paper: SigMA-HRD-positive cell lines
responded to PARP inhibitors, and SigMA-HRD-positive ovarian cancer
patients had longer overall survival on platinum-based therapy — i.e. the
signal has real, published clinical-actionability backing, not just a
statistical association.

---

## 2. What SigMA outputs — full enumeration, sample-level vs. variant-level

**Everything SigMA's `run()` produces is sample-level (one row per tumor
sample / per input MAF), not variant-level.** Individual mutations are
consumed to build the aggregate 96-context spectrum and then discarded —
SigMA does not annotate or return per-variant fields. This matters for
Varan: it would be a new *sample-level* column set (like MSI/TMB already
are in `data_clinical_sample.txt`), not a new MAF/`data_mutations.txt`
column.

Confirmed output columns (from `R/run.R` source + rdrr.io `run()` docs +
FAQ):

| Column (pattern) | Sample- or variant-level | Meaning |
|---|---|---|
| `total_snvs` | sample | SNV count actually used (post make_matrix filtering to SNPs) |
| `Signature_3_c{1..N}_ml` | sample | Per-cluster likelihood the spectrum matches a Sig3-positive tumor subtype |
| `Signature_3_ml` | sample | Sum of the above — aggregate Sig3 likelihood |
| `exp_sig3` | sample | NNLS-fit exposure (raw count) of Signature 3 |
| `rat_sig3` | sample | `exp_sig3 / total_snvs` — Sig3 exposure as a fraction of the sample's mutation burden |
| `Signature_3_mva` | sample | **The** final score: GBM-combined likelihood + cosine similarity + NNLS exposure. This is the number a significance/positivity threshold would be applied to. |
| `pass_mva` / `pass_mva_strict` | sample | Binary Sig3-positive/negative calls, from SigMA's own built-in cutoffs on `Signature_3_mva` (`do_assign=True`) — a permissive and a strict variant |
| `categ` | sample | A categorical Sig3 call, only present with `lite_format=True` (condensed output mode) |
| `Signature_msi_ml`, `Signature_pole_msi` | sample | Only when `check_msi=True` — see next section, this is **not** the same thing as clinical MSI status |

`conf.ini`'s existing `CHECK_MSI` flag maps directly and correctly to
`run(check_msi=...)`. What it actually does, confirmed from `R/run.R`:
it adds two extra likelihood-matching steps ("median_catalog" and
"decompose" against "average" and "cosmic_tissue" catalogs) whose results
get a `_msi` column suffix, and the resulting `Signature_msi_ml` /
`Signature_pole_msi` columns are used internally (per the code, filtered
at `> 0.99`) to **flag samples whose mutational spectrum looks
MSI-high or POLE-mutant** — because an MSI-high or POLE-ultramutated
spectrum can spuriously resemble the Signature-3 spectrum and confound the
HRD call. **This is a confounder-detection safety check for the Sig3 call
itself, not a clinical MSI-H/MSS caller.** Varan already has its own,
independent, orthogonal MSI mechanism (`msisensor`-based,
`[MSI] THRESHOLD_SITES`/`THRESHOLD_MSI`/`THRESHOLD_MSI_LIQUID` in
`conf.ini`, feeding `data_clinical_sample.txt`'s `MSI`/`MSI_THR` columns).
The two must not be conflated in a future GUI — `CHECK_MSI` in `[SigMA]`
is "should SigMA sanity-check its own Sig3 call against an
MSI-like/POLE-like spectrum", not "compute this sample's MSI status".

There is also a separate, more elaborate **MMRD (mismatch-repair
deficiency) / POLE classification module** documented on its own wiki page
("MMRD input file format"), which takes a materially different, larger
input (the 96 SBS columns *plus* `nins`/`ndel` indel counts, repeat-region-
overlapping indel counts `nmsi_ins`/`nmsi_del`, and an optional MSIsensor
score column) computed via `bedtools` or a built-in fallback. This is **not**
what `CHECK_MSI` in `run()` triggers — it is a distinct, heavier workflow
that isn't part of the drafted `[SigMA]` conf.ini section at all. Flagging
this so nobody assumes `CHECK_MSI=True` already gets Varan MMRD/POLE
detection — it doesn't; that would be a separate, later scope decision.

### Does SigMA cover LOH / HRD / copy-number analysis? No — verified, not inferred.

This needed a definitive answer and got one, from three independent angles:

1. **Source code.** `make_matrix()` (the function that turns a MAF/VCF into
   SigMA's input) filters explicitly to `Variant_Type == "SNP"` and further
   validates that both the reference and alternate allele are exactly one
   base long — i.e. it structurally cannot ingest a copy-number segment,
   a B-allele-frequency track, or even an indel. There is no code path in
   SigMA that reads copy number, LOH, or structural variant data at all.
2. **The original paper.** SigMA is presented as inferring Signature 3 *from
   SNVs*; HRDetect (a different, comparator method cited in the paper) is
   mentioned as a separate mutational-signature-based predictor, not as
   something SigMA incorporates.
3. **The clinical validation paper (Zhao et al. 2022, PMC9623231)**, which
   is the most directly relevant source since it's literally SigMA applied
   to MSK-IMPACT-style panel data: the authors *separately* computed a
   **GIS (genomic instability score)** via the **scarHRD** algorithm from
   WES data — an entirely different, LOH/copy-number-based method — and
   compared it against SigMA's substitution-based Sig3 call. They found
   **78% concordance**, explicitly framing this as the two methods
   capturing *different* aspects of HR deficiency (mutational scarring vs.
   allelic imbalance), not as SigMA producing the same answer by another
   route.

**Conclusion: SigMA is strictly SNV/trinucleotide-spectrum-based. It has no
LOH, no allele-specific copy number, and no HRD/genomic-instability score
in its own scope, and its own published validation explicitly treats
LOH-based HRD scoring (scarHRD/myChoice-style) as a separate, complementary
analysis — not something SigMA computes or subsumes.** If Varan ever wants
an LOH/HRD-score output, that is a genuinely separate development effort
(a CNV/B-allele-frequency-based tool, needs either SNP-array-quality allele
frequencies or WGS — Varan's existing `data_cna_hg19.*` pipeline computes
discrete copy-number calls and log2 fold-change from panel data, but does
not ingest B-allele frequency and does not compute any LOH/genomic-
instability metric today). This is a real gap to flag plainly, not paper
over: **"SigMA integration" does not mean "HRD score integration."** If HRD
scoring by LOH is actually the clinical goal, SigMA alone doesn't get there.

---

## 3. Inputs required, cross-checked against what Varan already has

**Genome build.** `make_matrix()` supports both hg19 and hg38 via
`BSgenome.Hsapiens.UCSC.hg19`/`hg38` R packages, selected with a
`ref_genome_name` parameter ("hg19" or "hg38") or a passed-in BSgenome
object. Varan's own pipeline is entirely hg19/GRCh37 today —
`conf.ini`'s `VEP_DATA = test_folder/vep_cache/homo_sapiens_vep_111_GRCh37`,
and the CNA pipeline's own file naming (`data_cna_hg19.seg`,
`data_cna_hg19.seg.fc.txt`) confirms this explicitly. **Compatible**: use
`ref_genome_name="hg19"`. No build mismatch risk as long as this stays
consistent with the rest of the pipeline (it currently is).

**MAF/VCF column requirements.** `make_matrix(file_type="maf")` requires:
`Chromosome`, `Start_position`, `End_position`, `Reference_Allele`,
`Tumor_Seq_Allele1`, `Tumor_Seq_Allele2`, `Tumor_Sample_Barcode`, and
(used for its own SNP-only filtering) `Variant_Type`. These are all
standard MAF-spec columns that `vcf2maf.pl` — which Varan already runs via
`conf.ini`'s `VCF2MAF` path — produces by definition. **Not independently
re-verified against an actual Varan-produced `.maf` file in this pass**
(no example `.maf` is committed in-repo to check column names against
directly) — flagged as a "confirm on first real run" item, but there is no
structural reason to expect a gap here since it's the same vcf2maf.pl
output Varan's own `filter_clinvar.py`/`sigma_filter.py` already read.

**SNVs only, no indels.** `make_matrix()` filters to `Variant_Type=="SNP"`
and single-base ref/alt itself. `prepare_sigma_maf()` in `sigma_filter.py`
does **not** duplicate that filter — confirmed by reading the source, this
is intentional and correct (no need to filter twice), and the module
docstring has been corrected this round to say so explicitly rather than
pointing at a `project_future_implementations.md` file that doesn't exist
in this repo.

**Minimum mutation count.** SigMA's own `snv_cutoff` parameter (default 5,
tunable per tumor type/platform internally — e.g. the source uses 4 for
prostate, 3 for osteosarcoma/pancreatic-neuroendocrine on panel data, and
higher cutoffs for WES/WGS) is the built-in floor. The clinical validation
paper is directly informative here since it used a near-identical setting
to Varan's: on a real MSK-IMPACT-like cohort they used **≥5 SNVs**, and
77% of real panel samples cleared that bar; on a 360-gene panel they used
**≥4 SNVs**, met by 75% of samples. **Practical implication for Varan: a
meaningful minority (roughly a quarter) of real panel samples will not
have enough SNVs for a reliable call and will come back negative/
uninformative by construction, not due to any integration bug.** This is
inherent to the method on panel-scale data, not something to "fix."

**Matched normal vs. tumor-only.** SigMA is designed for and validated on
tumor-only panel data — `make_matrix()` only requires
`Tumor_Sample_Barcode`, no normal-sample handling exists in the source at
all. This is exactly the SigMA use case the paper is built around
(getting a usable HRD-adjacent signal without a matched normal). It lines
up with why `sigma_filter.py`'s `prepare_sigma_maf()` does its own
population-AF + VAF-exclude-band germline scrub before handing data to
SigMA — that scrub is compensating for the lack of a matched normal, and
is correctly scoped to do exactly that (nothing to fix there).

**tumor_type / data platform vs. what Varan has available.** Varan's
`sample.tsv`/`walk.py` already carries a mandatory per-sample
`ONCOTREE_CODE` column (validated in `_walk_setup`'s required-columns
check) and one global `SAMPLE_TYPE` (`Solid`/`Liquid`, `[Sample_Type]`
in `conf.ini`) — see §6 for how well `sigma_cancer_type_map.py`'s
OncoTree→tumor_type bridge holds up. **`SAMPLE_TYPE` itself is not
consumed by SigMA at all** — there is no liquid-biopsy/cfDNA-specific
SigMA model. This is worth flagging as an open question rather than
silently assuming it's fine: SigMA's published validation (Gulhan 2019,
Zhao 2022) is entirely on tissue/FFPE-derived tumor content; a cfDNA
sample's VAF spectrum (variable, often lower, tumor-fraction-dependent)
was never part of that training/validation population, and Varan
explicitly supports and is currently configured for `[Sample_Type] TYPE =
Liquid`. Whether SigMA's output means anything reliable on a liquid
sample is genuinely unresolved by anything published — this needs a human
decision (e.g. gate SigMA to solid samples only, or run it but label
liquid-sample output as exploratory/unvalidated).

---

## 4. Runtime

**No documented runtime number exists** — not in the README, not in any
wiki page (Quick start, Parameter choices, FAQs, Installation), not in the
paper, and no GitHub issue found discussing performance/runtime/memory
despite a direct search. The closest thing to a number is qualitative,
from the clinical validation paper: it states only that "the
implementation of Sig3 calculation adds minimal compute cost and time"
relative to alternative HRD methods (like scarHRD/GIS from WES) — no
figure attached.

Given that, this is an **inference from the algorithm's design**, not a
sourced fact, and should be labeled as such to whoever decides how to
proceed: `run()` at inference time does a likelihood match, a cosine
similarity, an NNLS fit, and a pre-trained GBM prediction, all against
small precomputed reference matrices (on the order of tens of signature
vectors × 96 trinucleotide contexts) and a 96-length input vector per
sample. There's no MCMC, no per-run model training, no simulation loop —
the GBM classifier is already fit; only prediction happens per sample.
That points to **per-sample compute in the range of low seconds**, with
process/package-load overhead (R startup, loading the pre-built reference
`.rda` model objects) likely dominating actual data-driven compute time —
similar in character to Varan's existing per-sample `vcf2maf.pl`/
`MafAnnotator.py` subprocess calls, which are also dominated by fixed
per-invocation overhead rather than input size. **This needs an actual
benchmark run before anyone commits to a runtime budget** — recommend
timing one real sample locally (`Rscript` calling `run()` once) before
deciding always-on vs. opt-in at scale, since "should be fast" is not the
same as "is fast."

### Recommendation: optional flag, not an always-on pipeline step

Given the above, a conf.ini-gated `RUN_SIGMA = True/False` toggle (future
GUI checkbox) is the right call over wiring SigMA unconditionally into
every run, for several independent reasons, not just runtime uncertainty:

1. **Runtime is unmeasured.** Until benchmarked, defaulting to always-on
   risks silently adding unknown per-sample latency to every single run,
   including runs where nobody asked for or will look at a Sig3 score.
2. **The tumor_type mapping is not fully validated** (§6) — running it
   unconditionally on every sample, including tumor types whose mapping is
   still a documentation-derived guess, risks either silent wrong answers
   (mapped to the generic "other" model when a better model might exist)
   or an outright SigMA-side crash for the two confirmed
   panel-data-unsafe types (medullo, ewing — see §6).
3. **A meaningful fraction of samples (~25%, per §3) won't have enough
   SNVs for a real answer anyway** — for those, an always-on step just
   produces a negative-by-construction result, which reads to a clinician
   as "no HRD signature" when it may really mean "not enough data," a
   distinction the report/GUI would need to make very clear.
4. **LOH/HRD coverage gap (§2)** means Sig3 alone answers only part of
   "is this tumor HRD," not the whole clinical question — appropriate to
   treat as an optional, exploratory signal until/unless it's paired with
   an LOH-based method, not baked in as a default clinical field.
5. **This is exactly the pattern Varan already uses everywhere else.**
   OncoKB annotation (`-k` / `[OncoKB] ONCOKB`), the MSI caller, TMB, and
   fusion filtering are all already conf.ini/flag-gated, off by default or
   explicitly opted into per run — SigMA fits the same "gate behind an
   explicit flag" architecture already established, not a new pattern.

**Parallelization is a non-issue either way.** Per-sample SigMA runs are
embarrassingly parallel across samples (each sample's MAF → spectrum → run()
is fully independent), exactly the same shape as Varan's existing
per-sample `vcf2maf.pl`/`MafAnnotator.py` subprocess loop in `walk.py`. The
`Snakefile`'s existing DAG (`walk_setup → {walk_cnv, walk_snv, walk_fusion,
walk_clinical}`, all independent, all depending only on `walk_setup`'s
output) is the natural home for a future `walk_sigma` stage: it would
depend on `walk_snv`'s per-sample annotated MAF output and nothing else,
slotting in as a fifth parallel rule rather than needing to be threaded
inline into the SNV stage itself. That's a strong argument for the
flag-gated, separate-stage shape being not just safer but also the
*easier* integration, when the time comes.

---

## 5. Every SigMA parameter a user should be able to configure — proposed `[SigMA]` conf.ini surface

The existing draft has `VAF_MIN`, `VAF_EXCLUDE_BANDS`, `DATA_PLATFORM`,
`CHECK_MSI`. Missing, per the parameter audit above:

```ini
[SigMA]
; Master toggle - off by default until runtime is benchmarked and the
; tumor_type mapping (sigma_cancer_type_map.py) is confirmed against a
; live OncoTree lookup and a live SigMA install. Future GUI: single
; checkbox, "Run mutational signature (SigMA) analysis".
RUN_SIGMA = False

; --- existing, unchanged ---
VAF_MIN = 0.0
VAF_EXCLUDE_BANDS = [[0.45, 0.55], [0.90, 1.0]]
DATA_PLATFORM = msk
; SigMA's own MSI/POLE confounder safety-check on the Sig3 call itself -
; NOT Varan's clinical MSI status (see [MSI] section). Do not present
; these as the same knob in a future GUI.
CHECK_MSI = False

; --- new: SigMA's own snv_cutoff, currently hardcoded to run()'s default
; of 5 if wired naively. Verified real-world panel usage varies this
; (4 for a 360-gene panel, 5 for MSK-IMPACT-like, per the clinical
; validation paper) - expose it so a lab can tune for their own panel size
; rather than inherit an unrelated lab's default.
SNV_CUTOFF = 5

; --- new: manual override of the automatic ONCOTREE_CODE -> SigMA
; tumor_type mapping (sigma_cancer_type_map.py). Empty = use the automatic
; per-sample mapping. Set to force every sample in a run to one tumor_type
; (e.g. for a single-cancer-type cohort/validation run), bypassing the
; still-unvalidated OncoTree bridge entirely.
TUMOR_TYPE_OVERRIDE =

; --- new: exposes get_sigma_tumor_type()'s existing fallback_to_other
; argument, currently hardcoded True in sigma_cancer_type_map.py. If
; False, samples with an OncoTree code not in ONCOTREE_TO_SIGMA are
; excluded from SigMA entirely instead of run through the generic "other"
; pan-cancer model - a real clinical-judgment call, not a code default.
FALLBACK_TO_OTHER = True

; --- new: FAQ-documented requirement - "MVA models are trained using
; COSMIC catalog v2 and that version should be used to get accurate
; predictions." run()'s do_mva path (the Signature_3_mva score, the whole
; point of running SigMA here) needs this pinned, not silently defaulted
; to whatever SigMA's own default is (which may drift to v3 for
; non-MVA-path uses). Pinning explicitly here rather than relying on
; SigMA's internal default avoids a silent behavior change on a future
; SigMA version bump.
COSMIC_VERSION = v2

; --- new: condensed vs. full output. Full output includes every
; Signature_3_c{N}_ml per-cluster likelihood column (noisy for a clinical
; report table); lite_format=True keeps just the summary categ/score
; columns. Recommend True for anything surfaced to a clinician-facing
; report, False if the full per-cluster detail is wanted for QC/debugging.
LITE_FORMAT = True

; --- new: what "Signature-3 positive" means for Varan's own report/GUI.
; "default" = use SigMA's own built-in pass_mva/pass_mva_strict cutoffs on
; Signature_3_mva unmodified (recommended starting point - these are the
; cutoffs published/validated in the paper). A numeric value overrides
; with a lab-chosen cutoff on the same Signature_3_mva column instead -
; per SigMA's own Parameter-choices wiki guidance: "If the values disagree
; the cutoffs...need to be optimized, or a new model needs to be tuned" -
; i.e. this is an explicitly-supported thing to tune, not a hack.
SIGNATURE3_POSITIVE_THRESHOLD = default
```

Naming/typing notes for the eventual GUI: `RUN_SIGMA`/`FALLBACK_TO_OTHER`
are booleans matching Varan's existing `check_bool()` convention
(`True`/`true`/`T`/`False`/`false`/`F`/empty). `SNV_CUTOFF` is a plain
integer. `TUMOR_TYPE_OVERRIDE` is a free-text OncoTree-style string, empty
by default. `SIGNATURE3_POSITIVE_THRESHOLD` is either the literal string
`"default"` or a float — same "sentinel-string-or-number" shape a GUI
would render as a dropdown ("Use SigMA default" / "Custom") plus a numeric
field that only enables when "Custom" is picked.

---

## 6. `sigma_cancer_type_map.py` audit

**Verdict: not yet trustworthy enough to ship as-is, for a materially
different reason than the original draft flagged.** The original draft's
stated concern was "unverified child OncoTree codes" (real, still open —
see below). Actually reading SigMA's source turned up something more
concrete: **the draft's tumor_type vocabulary was itself incomplete and,
worse, didn't distinguish which tumor_type values are safe to use with
panel data (Varan's `DATA_PLATFORM = msk`) at all.**

**Corrected/verified facts (this round), already reflected in updated
docstrings in `sigma_cancer_type_map.py`:**

- The officially documented `tumor_type` vocabulary (rdrr.io's rendered
  `run()` man page) is **17 values**, not 13:
  `bladder, bone_other, breast, crc, eso, gbm, lung, lymph, medullo, osteo,
  ovary, panc_ad, panc_en, prost, stomach, thy, uterus`. The round-1
  draft's list of 13 (`+ "other"` as a 14th) omitted `crc`, `gbm`, `lung`,
  `lymph`, `thy` — real, valid values that just weren't found before.
  Note `"ewing"` and `"other"` are each independently confirmed valid
  *elsewhere* in the R source (a `stop()` message and an
  `all_catalogs[['other_multisig']]` reference respectively) but do
  **not** appear in that same rendered man-page list — an inconsistency
  inside SigMA's own documentation, not something introduced by this
  analysis.
- **Concrete, source-verified risk found in the draft mapping:**
  `R/run.R`'s own model-availability check
  (`gbm_models[[data]]`, quoted in full in the module docstring now)
  states plainly that `medullo` and `ewing` only have trained
  `do_mva=True` classifiers for **whole-exome/WGS data**, not for
  **panel data** — which is exactly Varan's `DATA_PLATFORM = msk`
  setting. Calling SigMA with `tumor_type="medullo"` or `"ewing"`,
  `data="msk"`, `do_mva=True` (the intended clinical mode) would very
  likely raise SigMA's own `stop()` error. The panel-safe (`msk`) model
  list is only explicitly named, in that same message, as 10 values:
  `eso, osteo, ovary, panc_ad, panc_en, prost, stomach, uterus, breast,
  bladder` — exactly matching the draft mapping's non-bone/non-ewing/
  non-medullo coverage. **Fixed this round**: `MDB`/`ES` are left mapped
  (they are valid tumor_type strings) but a new
  `DO_MVA_UNSAFE_FOR_PANEL_DATA` constant and inline comments now flag
  them explicitly, so a future wiring pass doesn't call `do_mva=True` on
  them blind. This is a real bug class avoided, not just a note added.
- `crc`/`gbm`/`lung`/`lymph`/`thy` — common oncology-panel tumor types
  (colorectal, glioblastoma, lung, lymphoma, thyroid) with confirmed
  SigMA models for *some* platform — were **deliberately not added** to
  the mapping this round. Nothing in the public source confirms whether
  their model exists for panel (`msk`) data specifically, as opposed to
  WES/WGS-only (the same failure mode just found for medullo/ewing).
  Guessing them in risks reintroducing exactly the bug just fixed.
  Confirming this needs a live R session
  (`library(SigMA); names(gbm_models[["msk"]])`) rather than another
  round of doc-reading — flagged as the top open item below. Until then,
  OncoTree codes for these tumor types (very likely a large share of any
  real oncology panel's sample mix — lung and colorectal especially)
  correctly fall through to `tumor_type="other"` via the existing
  `fallback_to_other` mechanism, rather than being silently mis-mapped.
- **Child-level OncoTree code verification remains genuinely open.**
  `oncotree.mskcc.org`'s API returned 403 to automated fetches again this
  round (tried the tumor-types endpoint directly, and a search for a
  mirrored static OncoTree data file) — same block the round-1 draft
  author hit. The individual codes (`BRCANOS`, `HGSOC`, `OCS`, `EOV`,
  `PRSCC`, `UTUC`, `TSTAD`, `DSTAD`, `UMEC`, `ESCC`, `PANET`, `CHDM`,
  `MDB`, etc.) are consistent with commonly-seen OncoTree/cBioPortal
  nomenclature from general domain familiarity, but that is **not** the
  same as a confirmed live lookup, and is explicitly called out as such
  now in the module docstring. **This is the one item from the original
  draft's own stated concern that is still unresolved** — a human with
  either OncoTree portal access or a downloaded OncoTree release file
  needs to check these against the real ontology before any clinical use.

---

## 7. Open questions needing a human decision

1. **Benchmark SigMA's actual runtime** (one real sample, timed) before
   deciding on a batch-size/timeout budget for the future `walk_sigma`
   stage. Nothing published gives a number; §4's estimate is inference
   from the algorithm's shape, not a fact.
2. **Confirm `crc`/`gbm`/`lung`/`lymph`/`thy` panel-data (`msk`) model
   availability** via a live SigMA R install
   (`names(gbm_models[["msk"]])`) before adding them to
   `ONCOTREE_TO_SIGMA` — lung and colorectal in particular are likely a
   large share of real samples and currently fall back to `"other"`.
3. **Confirm the remaining child-level OncoTree codes** in
   `sigma_cancer_type_map.py` against a live `oncotree.mskcc.org` query or
   a downloaded OncoTree release (API access blocked to this analysis both
   rounds).
4. **Decide the liquid/cfDNA question** (§3): gate `RUN_SIGMA` to solid
   samples only, or run it on liquid samples too but label the output as
   exploratory/unvalidated in the report — SigMA's published validation
   population is entirely tissue-derived.
5. **Decide the medullo/ewing handling** (§6): for those two OncoTree
   codes on panel data, choose one of (a) force `do_mva=False` and only
   surface the raw likelihood/exposure features, not a binary call,
   (b) route them to `tumor_type="other"` instead, or (c) exclude them
   from SigMA entirely. Not decided here — `DO_MVA_UNSAFE_FOR_PANEL_DATA`
   in `sigma_cancer_type_map.py` exists so this choice is made explicitly,
   not accidentally.
6. **Decide whether an LOH/HRD-score capability is actually wanted**
   (§2) — if yes, that is a separate, new development effort (a
   CNV/B-allele-frequency-based method), not something SigMA provides or
   this branch's scope covers.
7. **Confirm `SIGNATURE3_POSITIVE_THRESHOLD`'s default policy** — start
   with SigMA's own built-in `pass_mva`/`pass_mva_strict` cutoffs
   (recommended, since those are what's actually published/validated), or
   pre-commit to a lab-specific override immediately.

---

## Fixes made this round (docs/analysis only, no runtime behavior change)

- `sigma_filter.py`: corrected a dangling reference to a
  `project_future_implementations.md` file that does not exist anywhere in
  this repo (checked — no such file, on this branch or `stabilize`);
  docstring now points at this report and explicitly confirms, from
  reading SigMA's actual `make_matrix()` source, that the SNV-only
  filtering it relies on (not duplicated in `prepare_sigma_maf()`) is
  correctly handled downstream by SigMA itself.
- `sigma_cancer_type_map.py`: corrected the tumor_type vocabulary (17
  documented values, not 13; noted the `ewing`/`other` documentation
  inconsistency inside SigMA itself); added a source-verified,
  concrete-crash-risk warning (`DO_MVA_UNSAFE_FOR_PANEL_DATA`) for
  `medullo`/`ewing` on panel (`msk`) data with `do_mva=True`; documented
  why `crc`/`gbm`/`lung`/`lymph`/`thy` were deliberately left unmapped
  rather than guessed in.

Neither change alters any code path currently exercised by the pipeline —
both files remain unwired (no caller anywhere in `walk.py`/`varan.py`/the
`Snakefile`). Docstring/comment corrections and one new unused constant
only.
