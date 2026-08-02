# Multi-vendor input generation for Varan — architecture and status

## Round 3: conf.ini-only selection, real pipeline names, sample-type unification, and samplesheet-mode input

This round changes four things about round 2's design, based on direct
user feedback after using it. **This section is the current, authoritative
reference for vendor selection and dispatch** - round 2's section below is
kept for its still-accurate design rationale (scratch-folder reuse, the
copy-out-of-scratch mechanics, the "adapter finds nothing" sanity check)
but its `--pipeline` CLI flag and `"native"` pipeline name no longer exist;
read this section first.

### 1. Vendor selection is now conf.ini-only - no CLI flag

Round 2 added a `--pipeline` CLI flag alongside conf.ini's `[Vendor]
PIPELINE`, with the flag taking precedence. Removed: `varan.py` no longer
has a `--pipeline` argument at all. `resolve_pipeline()` now just reads
`get_config().get("Vendor", "PIPELINE", fallback="").strip()` directly -
conf.ini is the *only* place vendor selection lives. This matches the
original framing of conf.ini's `[Vendor] PIPELINE` as "meant for a future
GUI to render as a dropdown of known, registered options" - a single
source of truth is simpler for that GUI to drive than a value that can be
set two different ways with a precedence rule to explain.

### 2. Real pipeline names: `illumina_solid` / `illumina_liquid`, not `native`

The user explicitly rejected `"native"` as too vague and asked for the
actual, precise name(s) already used for Varan's own built-in ingestion
path - investigated rather than guessed:

- `README.md` (lines 12, 22): *"Varan turns raw variant-calling output
  (VCF, MAF, **Illumina TSO500** `CombinedVariantOutput`/`MetricsOutput`,
  and more)..."* and *"**TSO500-aware.** Reads TMB, MSI, fusions, HRD/
  Genomic Instability Score, and exon-level BRCA1/BRCA2 CNVs straight out
  of **Illumina's** `CombinedVariantOutput`/`MetricsOutput` files..."*
- `walk.py` already uses exactly this vocabulary in its own code comments,
  independent of anything written for this feature: line ~935's *"native
  **TSO500** CombinedOutput ingestion"*, and `fill_from_file()`'s
  docstring (~line 1795): *"...forced **non-Illumina pipelines (e.g.
  Guardant)** that already know their own threshold to encode it via a
  fake VALUE..."* - i.e. this codebase already calls Guardant
  "non-Illumina" and treats "Illumina"/"TSO500" as the term for its own
  native path, well before this multivendor effort started.
- `tsv.py` independently cross-references "Illumina's own TSO500 v2.2"
  documentation for its CombinedVariantOutput.tsv parsing.

So `"illumina"` is not a guess - it's this codebase's own pre-existing
term for "the built-in path, no adapter needed". The reason for landing on
**two** values (`illumina_solid`/`illumina_liquid`) rather than one plain
`illumina` comes from investigating `[Sample_Type] TYPE` (see point 3
below): Illumina genuinely ships two distinct, separately-named TSO500
assay variants - a solid-tumor/tissue panel and a ctDNA/liquid-biopsy
panel - and Varan's own pre-existing `SAMPLE_TYPE` global already branches
its CombinedOutput MSI handling on exactly that distinction (`walk.py`
~line 1885: `if SAMPLE_TYPE == "LIQUID": ... elif SAMPLE_TYPE == "SOLID":
...`, feeding different MSI-stability thresholds - `THRESHOLD_MSI_LIQUID`
vs. the general `THRESHOLD_MSI` eval string). Naming the two pipeline
values after the two real assay variants, rather than one generic
`illumina`, is what lets `PIPELINE` alone determine sample type
symmetrically for every registered option (see point 3) instead of only
for vendor adapters like Guardant.

`""` (blank/absent) remains a third, distinct state - deliberately *not*
defaulted to either `illumina_solid` or `illumina_liquid` - meaning
"unspecified/legacy": no adapter runs, and `[Sample_Type] TYPE` is used
exactly as literally configured, with zero reconciliation applied. This
preserves byte-identical behavior for every conf.ini that predates this
feature (including one with no `[Vendor]` section at all) or that simply
doesn't want the extra consistency check. Choosing `illumina_solid` or
`illumina_liquid` explicitly is a *strict superset* of blank's behavior
(same "no adapter" input handling) plus the sample-type assertion from
point 3 - so existing deployments can adopt the new, more precise names
at their own pace without anything breaking if they don't.

If this specific naming choice (`illumina_solid`/`illumina_liquid`) turns
out not to match what the team actually calls these two variants
day-to-day, renaming is a small, contained change: `varan.py`'s
`ILLUMINA_PIPELINE_SAMPLE_TYPES` dict (one line per value), its two
mentions in conf.ini/Templates/conf.ini's comments, and this doc section -
nothing else depends on the exact strings.

**`origin/liquid_samples` branch investigated per the user's suggestion**
(commit `5f82cbf`, *"add function to get sample_type from report"*, and
`418fb43`, *"add different behaviour for solid and liquid samples"*):
this branch is where `[Sample_Type] TYPE` and `walk.py`'s `SAMPLE_TYPE`-
gated MSI logic were originally introduced (`418fb43`) - now already
present in `stabilize` via other, earlier merges, so nothing new to port
from there. `5f82cbf`'s "get sample_type from report" function
(`write_report.py`'s `extract_sample_type_from_html()`) turned out to be
unrelated prior art for this specific question: it recovers a *previous
Varan run's own* recorded sample type by regex-scraping that run's
already-generated HTML report (for `update`/`extract` operations merging
into an existing study), not a way to derive sample type from raw vendor
input. Interesting precedent for "don't make the user re-supply
information Varan already has/can infer", but not directly reusable here.

### 3. `[Vendor] PIPELINE` now unifies with `[Sample_Type] TYPE` where safe

Implemented in `varan.py`'s new `implied_sample_type()`/
`reconcile_sample_type()`, called once, immediately after `pipeline` is
resolved and validated, and *before* `walk` is imported (this ordering is
required: `walk.py` reads `[Sample_Type] TYPE` into a module-level
`SAMPLE_TYPE` global at its own import time, from the same shared
`config_loader` singleton `reconcile_sample_type()` mutates in place -
anything decided after that import would be invisible to `walk.py`).

- `illumina_solid` → implies `Solid`; `illumina_liquid` → implies
  `Liquid` (see point 2).
- A vendor adapter can declare its own fixed sample type as a
  module-level `SAMPLE_TYPE` attribute - `vendor_adapters/guardant.py`
  now sets `SAMPLE_TYPE = "Liquid"`, since Guardant360 CDx is ctDNA
  (blood plasma) only, with no solid-tumor/tissue variant. This lives in
  the vendor's own module (like `NAME`), not in a lookup table in
  `varan.py` - each vendor's own facts belong with that vendor.
- Blank/legacy `PIPELINE`, or a hypothetical future vendor that supports
  *both* sample types (and so declines to declare a fixed `SAMPLE_TYPE`),
  implies nothing - `reconcile_sample_type()` is a no-op, and
  `[Sample_Type] TYPE` stays a fully independent, manually-set value, same
  as before this feature existed.

**Reconciliation rule**, when a pipeline *does* imply a sample type:
- `[Sample_Type] TYPE` blank → auto-filled from the implied type (this is
  the actual redundancy elimination the user asked about: a deployment
  using `PIPELINE=guardant` no longer needs to *also* set `TYPE=Liquid`
  separately).
- `[Sample_Type] TYPE` already set and it agrees (case-insensitively) →
  no-op, already consistent.
- `[Sample_Type] TYPE` already set and it *disagrees* → raises a clear
  `ValueError` naming both conflicting values, rather than silently
  trusting one over the other. Silently keeping conf.ini's literal value
  would let a real mismatch (e.g. a copy-pasted `TYPE=Solid` alongside
  `PIPELINE=guardant`) through unnoticed; silently overwriting it would
  hide whatever led to that value being set in the first place. An
  explicit, actionable error is the only option that doesn't guess.

This is intentionally *not* a full redesign of how sample type is
configured (e.g. removing `[Sample_Type]` entirely, or moving all vendor
facts into a bigger vendor-capabilities system) - that felt like more
restructuring than proportionate to what was asked. What's implemented is
a narrow, well-motivated unification for the cases where it's unambiguous
(a pipeline/vendor that is provably single-sample-type), while leaving
`[Sample_Type] TYPE` fully independent everywhere else, exactly as before.

### 4. Samplesheet (`selection.tsv`) mode: flexible columns, explicit per-file paths, graceful per-file skipping

`vendor_adapters.guardant.run()` already had a `selection` parameter
(round 1) alongside `folder` - a TSV batch mode, previously hardcoded to
exactly two columns (`sample_id`, `s3_path_run`), auto-discovering each
sample's files from its own run folder the same way whole-folder mode
does. This is now extended, in `vendor_adapters/guardant.py`, to match how
users actually want to run this: against a *subset* of samples from a
batch, not always a whole run folder, with per-sample flexibility for
files that don't follow the usual convention.

**Column schema** (only `sample_id` is required; every other column is
optional and independently recognized - see `_SELECTION_PATH_COLUMNS` and
`run()`'s docstring in `guardant.py`):

| Column | Required | Purpose |
|---|---|---|
| `sample_id` | Yes | Row is skipped (logged) if blank. |
| `s3_path_run` | No | A run folder to auto-discover this sample's files from (round 1's original column, still fully backward compatible) - optional now: a row supplying every file it needs via the explicit columns below doesn't need one. |
| `xml_path` (or `metadata_path`/`finalmetadata_path`) | No | Explicit path to this sample's metadata XML - local path or `s3://...` - overrides auto-discovery for this file only. |
| `vcf_path` (or `snv_vcf_path`) | No | Explicit path to this sample's somatic VCF. |
| `cnv_path` (or `cnv_call_path`) | No | Explicit path to this sample's `.cnv_call.hdr.tsv`. |
| `msi_path` (or `msi_call_path`) | No | Explicit path to this sample's `.msi_call.hdr.tsv`. |
| `fusion_path` (or `fus_path`/`fusion_call_path`) | No | Explicit path to this sample's `.fusion_call.hdr.tsv`. |

Any column not in this table (e.g. a `notes` or `batch_name` column a lab
wants to keep for its own bookkeeping) is simply ignored, not an error -
`csv.DictReader` reads by header name, so unrecognized columns cost
nothing. This is the "flexible columns" the user asked for: not a rigid,
positional schema, and forgiving of extra columns.

**Explicit path resolution** (`process_single_sample()`'s new
`explicit_paths` parameter): for each of the 5 file types, an explicit
path (if given) always wins over folder auto-discovery. A local path is
used as-is (never deleted - the cleanup step only removes files this
function itself downloaded into its own temp directory, tracked
separately, so a samplesheet's own files are never touched). An
`s3://...` explicit path is downloaded the same way an auto-discovered
file would be.

**Graceful per-file skipping**: `xml` and `vcf` are the only two mandatory
file types - if neither an explicit path nor `s3_path_run` auto-discovery
resolves one of those two for a given sample, that *whole sample* is
skipped (clearly logged), same as round 1's behavior for a missing VCF.
`msi`/`cnv`/`fus` are each independently optional - already true before
this round (the rest of the function already tolerated any of them being
absent) - what's new is that this graceful-skip path is now reachable
through the explicit-path/no-folder samplesheet mode too, with a clear,
per-file-type log message (`"[sample] Optional file type 'X' not
available (...) - skipping this data type for this sample."`) rather than
silently proceeding or crashing. This directly covers the scenarios the
user described - "no CNV data for this sample" and "the CNV/MSI files
don't combine well for this sample" - by simply leaving that one column
blank for that one row.

A missing **required** file (e.g. an explicit `vcf_path` that doesn't
actually exist) skips only that one sample, with a clear log message,
and processing continues with the rest of the batch - verified with a
synthetic 3-row selection.tsv (one auto-discovery row, one all-explicit
row with several optional columns blank, one row with a bad `vcf_path`)
that produced exactly the two valid samples and skipped the third.

**Wired into `varan.py`'s dispatch**: `run_vendor_adapter()` now decides
`folder=` vs. `selection=` by checking whether `-i`'s (reinterpreted) raw
input is an existing local *file* (`Path(raw_input).is_file()`) - if so,
it's a samplesheet, passed as `selection=`; otherwise (an S3 URI, a local
directory, or anything else that isn't an existing local file) it's
passed as `folder=`. This mirrors the *exact same* is_file()/is_dir()
distinction `-i`'s argument already uses elsewhere in this codebase (see
`walk._walk_setup()`) - it is a "file vs. folder" structural check, not a
vendor-identity guess (vendor identity is still 100% explicit, via
`[Vendor] PIPELINE`).

### Minimum inputs needed to start (concrete answer)

**Whole-folder mode** (`[Vendor] PIPELINE = guardant`, `-i <run_folder>`):
every sample under that folder needs, at minimum, its own somatic **VCF**
and **metadata XML** (`_finalmetadata.xml`/`_metadata.xml`) discoverable
by Guardant's filename/suffix convention (sample ID substring + suffix
match - see `find_file_in_list()`). `.cnv_call.hdr.tsv`,
`.msi_call.hdr.tsv`, and `.fusion_call.hdr.tsv` are all optional per
sample - if absent, that sample is still processed, just without CNV
gene-attribution / MSI value / fusion calls respectively.

**Samplesheet mode** (`-i <selection.tsv>`): per row, `sample_id` is
required; then either `s3_path_run` (to auto-discover the same 2
mandatory + 3 optional files as whole-folder mode) or the explicit
`xml_path`+`vcf_path` columns (the 2 mandatory ones) must resolve to real
files, or that row is skipped. `cnv_path`/`msi_path`/`fusion_path` are
optional in every case.

So, corrected from the "4 files" framing in earlier notes (which omitted
the fusion report): **5 file types total, 2 mandatory (VCF + metadata
XML), 3 optional (CNV, MSI, fusion)** - this was already true of round 1's
code (the optionality of CNV/MSI/fusion was never new), this round mainly
makes that optionality reachable through samplesheet rows too, with
clearer logging, and corrects the doc's own miscount.

### Verification performed this round (and its limits)

Same core limitation as rounds 1-2: no real vendor example data on this
machine. What was actually exercised, via the same mocked local-folder
harness (S3 helpers monkeypatched to a local directory):

- `resolve_pipeline()` against a blank conf.ini (→ `""`) and a
  `[Vendor] PIPELINE = ...`-set one.
- `reconcile_sample_type()` against 5 scenarios: blank pipeline (no-op),
  `illumina_liquid` with blank `[Sample_Type] TYPE` (auto-filled to
  `Liquid`), `guardant` with a matching pre-set `TYPE=Liquid` (no-op,
  consistent), `guardant` with a conflicting `TYPE=Solid` (raised the
  expected `ValueError`), and an unrecognized pipeline name (`
  implied_sample_type()` returns `None` gracefully rather than crashing -
  `varan.py`'s own `known_pipelines` check catches an actually-invalid
  value separately, before reconciliation ever runs).
- The real repo's own `conf.ini` (`PIPELINE` blank, `[Sample_Type]
  TYPE = Liquid`) run through `resolve_pipeline()`/`reconcile_sample_type()`
  end-to-end confirmed byte-identical behavior: pipeline resolves to `""`,
  `TYPE` stays `Liquid`, untouched.
- A synthetic 3-row `selection.tsv` (auto-discovery row, all-explicit-
  local-paths row with 3 blank optional columns, bad-explicit-path row)
  run through `guardant.run(selection=...)` directly: produced exactly 2
  samples, skipped the third with a clear message, correctly left the
  samplesheet's own local files on disk afterward (not deleted by the
  temp-file cleanup step).
- `run_vendor_adapter()`'s `folder=`/`selection=` auto-dispatch, run
  against both a directory and a file path, correctly chose each mode.
- `python varan.py --help` confirmed no `--pipeline` flag remains, and
  `-i`'s help text reflects the new conf.ini-only wording.

Not verified (same as before): a real Guardant S3 run through the real
`aws` CLI, and the downstream VEP/vcf2maf/OncoKB pipeline stages.

## Round 2: `varan.py -i` runs a vendor's raw input directly

The first round (below) built `vendor_adapters/` but left it a manual,
two-step workflow: run `create_Varan_input.py --vendor guardant ...`
yourself, then feed its output to `varan.py -i sample.tsv patient.tsv
fusions.tsv`. This round wires the adapter directly into `varan.py`, so
pointing it at a vendor's raw input and running it "just works", the same
way the existing native (Illumina/CombinedOutput) path already does.

### Selection mechanism: explicit, not auto-detected

> **Superseded by round 3, point 1**: the `--pipeline` CLI flag described
> below was removed - conf.ini's `[Vendor] PIPELINE` is now the *only*
> selection mechanism, no precedence rule needed. The reasoning for
> explicit-over-auto-detected selection (this subsection's second half)
> still stands unchanged.

A new `--pipeline` CLI flag on `varan.py` chooses which vendor produced
`-i`'s input:

- `--pipeline native` (the default) — today's existing, completely
  unchanged behavior. `-i` points at Varan's own sample.tsv/patient.tsv/
  [fusion.tsv] files, or an already-shaped SNV/CNV/CombinedOutput input
  folder.
- `--pipeline guardant` (or any other registered vendor name — see
  `vendor_adapters/__init__.py`'s `ADAPTERS`) — `-i`'s first path is
  instead that vendor's raw run input, converted automatically before the
  rest of the pipeline runs.

Precedence: an explicit `--pipeline` on the CLI always wins; otherwise
`conf.ini`'s new `[Vendor] PIPELINE` key is used; otherwise the default is
`"native"`. This mirrors how `-C/--config` already overrides the default
`conf.ini` path elsewhere in this codebase — see `varan.resolve_pipeline()`.
`[Vendor] PIPELINE` being blank/absent by default, with a fixed set of
known values (`native` + whatever's registered in `ADAPTERS`), is meant to
let a future GUI render it as a dropdown, per the user's own framing of
that option.

**Deliberately not auto-detected from `-i`'s content.** Which vendor
produced a given raw input is a consequential classification — it picks
an entire parsing code path — so it's made an explicit choice rather than
an implicit guess, matching the reasoning already applied elsewhere in
this same integration effort (the SigMA branch's medullo/ewing `do_mva`
handling: `12d217d` forced that decision explicitly rather than leaving it
an implicit per-run choice). A lightweight *supplementary* sanity check
does exist: if the selected adapter's `run()` finds no usable data in the
given input, `varan.py` raises a clear, actionable error suggesting the
wrong `--pipeline` may have been selected, rather than either silently
proceeding or failing with a confusing downstream crash — but this never
substitutes for explicit selection as the actual dispatch mechanism, it
only helps catch a human picking the wrong one.

**Naming deviation, explained (superseded by round 3, point 2).** The
user's own suggested spelling was `--pipeline illumina/guardant`. This
implementation used `"native"` instead of `"illumina"` for the
default/non-vendor value, reasoning that the existing default path isn't
inherently Illumina-specific. **Round 3 reversed this call**: the user
rejected `"native"` outright and asked for real investigation rather than
a guess, which turned up that this codebase's own existing code comments
already call the non-adapter path "Illumina"/"TSO500" (see round 3, point
2, for the actual `README.md`/`walk.py` citations) - so `"native"` wasn't
just imprecise, it was inventing a term this codebase doesn't otherwise
use. Left here for the historical record of why the *first* naming
attempt was made and why it didn't hold up; `illumina_solid`/
`illumina_liquid` are the real, current values.

### How the raw-folder-to-canonical dispatch works

> **Partially superseded by round 3, point 4**: step 3 below originally
> always called `adapter.run(folder=...)`. It now checks whether `-i`'s
> argument is a local file (→ `run(selection=...)`, samplesheet mode) or
> not (→ `run(folder=...)`, unchanged) - see round 3 for the details.
> Steps 1/2/4/5/6 and the scratch-folder mechanics are otherwise still
> accurate, just substitute "a registered vendor name (not blank/
> illumina_solid/illumina_liquid)" wherever this says `"native"`.

No new input shape was invented at the `vendor_adapters` layer — adapters
still only take `run(folder=..., selection=..., **kwargs)`, exactly as
round 1 left them. All the new logic lives in `varan.py`:

1. `varan.py`'s `__main__` block resolves `pipeline` (see above). If it's
   not `"native"` (and this isn't an update/extract/remove run, which
   don't use `-i` the same way), it calls the new
   `run_vendor_adapter(pipeline, varan_input, output_folder)`.
2. That function creates a scratch folder via the *same*
   `create_random_name_folder()`/`clear_scratch()` helpers `walk.py`
   already uses for VEP's temp files (`<output_folder>/scratch/<random>/`)
   — reused as-is, not reimplemented. `create_random_name_folder()`
   creates all needed parent directories (`mkdir(parents=True)`), so this
   works even though the versioned output folder (`<output_folder>_v1`)
   doesn't exist yet at this point in the flow — that versioned folder is
   a *sibling* path (`Path(output_folder).parent / "..._v1"`, see
   `walk.create_folder()`), never touched by this step.
3. `varan_input[0]` (the vendor's raw input — e.g. an S3 run folder for
   Guardant) is passed straight through to `adapter.run(folder=...)`,
   with `report_base_dir`/`vcf_base_dir`/`temp_local_dir` defaulted to
   subfolders of that scratch dir (so a normal run leaves nothing behind
   in any shared/production location) — overridable per-vendor via a new
   `[Vendor.<name>]` conf.ini section (e.g. `[Vendor.guardant]` /
   `dict_path = ...`) whose keys are passed straight through as `run()`
   keyword arguments.
4. If `run()` returns `None` (no data found), `run_vendor_adapter` clears
   the scratch folder and raises a `ValueError` with the "wrong vendor
   selected?" hint mentioned above — caught by `varan.py`'s existing
   top-level `except ValueError` handler, same as any other input error.
5. Otherwise, the returned `report_path`/`fusion_table_path` become the
   new `varan_input = [sample_tsv, patient_tsv_passthrough, fusion_tsv]`,
   and execution falls straight through into the *same* `varan(...)` call
   the native path already used — from this point on there is no
   difference between a vendor-mode run and a hand-built `-i sample.tsv
   patient.tsv fusions.tsv` run. `transform_input()`/`check_folders()`
   (inside `_walk_setup()`) copy the referenced VCF/fusion files out of
   scratch into the run's real output folder during this same call, so
   nothing downstream depends on scratch surviving past it.
6. On success, `varan.py` removes the scratch folder after `varan(...)`
   returns. On failure, it's deliberately left in place for debugging —
   matching the existing precedent in `walk._walk_process_snv`, whose own
   VEP scratch folder is likewise only cleared on the success path.

A third `-i` element (a fusion file) is not accepted in vendor mode — the
adapter always generates its own — and is ignored with a logged warning
if one is passed, rather than silently doing something unexpected with it.

### Backward compatibility

`varan()` (the pipeline function itself) was not modified at all — every
line of the vendor dispatch logic lives in `varan.py`'s `__main__` block,
strictly before the existing `varan(...)` call. Leaving `[Vendor]
PIPELINE` blank/absent in conf.ini (there's no CLI flag as of round 3 -
see above) reaches `varan(...)` with `varan_input` completely untouched
from what argparse parsed — identical to every prior behavior, including
anyone still using the manual two-step Guardant flow from round 1 (that
flow still works unchanged: it just produces a sample.tsv/fusions.tsv pair
you can pass to `-i` yourself, exactly as before, with `[Vendor] PIPELINE`
staying blank the whole time).

### Known limitation: not wired into `walk_stage.py`/the Snakemake DAG

`walk_stage.py` (used by `Snakefile`'s per-stage `walk_setup`/`cnv`/`snv`/
`fusion`/`clinical` rules) calls `walk._walk_setup()` directly, bypassing
`varan.py`'s `__main__` entirely — so `--pipeline` has no effect there
yet. This was left alone deliberately: it's a separate, lower-level entry
point explicitly documented as carrying "no pipeline logic" of its own,
and wiring vendor dispatch into it wasn't part of what was asked this
round (the Snakemake-DAG migration itself is explicitly out of scope,
being handled on a separate branch). Anyone driving the Snakemake path
with a vendor input still needs to run `create_Varan_input.py`
(or `vendor_adapters.<name>.run(...)`) by hand first, same as round 1.

### Verification performed this round (and its limits)

Same core limitation as round 1: no real vendor example data is available
on this machine. What was actually verified, using the mocked local-folder
harness built for this round (S3 helpers monkeypatched to operate on a
local directory instead of shelling out to `aws s3`):

- `resolve_pipeline()`'s precedence (CLI override > conf.ini > "native"
  default) was exercised directly against real `ConfigParser` instances.
- The full `run_vendor_adapter()` dispatch was run end-to-end against a
  local mock "S3" folder built from the same synthetic Guardant fixtures
  used in round 1 (one sample's `_finalmetadata.xml`/`.vcf`/
  `.msi_call.hdr.tsv`/`.cnv_call.hdr.tsv`/`.fusion_call.hdr.tsv`):
  the generated `sample.tsv` and `fusions.tsv` landed under
  `<output_folder>/scratch/<random>/`, the SNV/CNV VCF conversion inside
  it reproduced the same documented behavior confirmed in round 1
  (PASS-only SNV row, correctly paired CNV rows), MSI/MSI_THR and the
  fusion rows matched expectations, and `clear_scratch()` correctly
  removed both the random-named scratch subfolder and its now-empty
  `scratch/` parent afterward.
- The no-data failure path (an empty mock run folder) was confirmed to
  raise the expected `ValueError` with the "wrong vendor?" hint, and to
  still clean up its scratch folder before raising.
- `python varan.py --help` was run to confirm `--pipeline {native,
  guardant}` parses correctly alongside every pre-existing flag.

Not verified: an actual real Guardant S3 run through the real `aws` CLI
(no AWS credentials/real bucket available here, on top of the pre-existing
lack of real vendor example files), and the full downstream `varan()`
pipeline stages (VEP/vcf2maf/OncoKB/etc. — those require external tools
not installed in this environment and were out of scope for both rounds).
The dispatch layer itself (argument resolution, adapter invocation,
scratch lifecycle, error handling) was exercised directly, which is the
part that changed this round; `varan()`'s own downstream stages are
unchanged code, already relied upon by every existing native-mode run.

## Why this doc replaces `GUARDANT_INTEGRATION_NOTES.md`

The original `create_Varan_input.py` was written and debugged
specifically against Guardant360 (FPG360) output. The goal now is
multi-vendor ingestion, with Guardant as the first (reference)
implementation rather than the only one. This doc covers the
vendor-agnostic architecture; Guardant-specific behavior, bug history, and
open questions are called out explicitly in their own section below so
nothing from the original notes gets lost.

## Architecture

```
vendor_adapters/
    __init__.py    # ADAPTERS registry + get_adapter(name) lookup
    common.py      # vendor-agnostic building blocks (see below)
    guardant.py    # Guardant360 adapter (reference implementation)
create_Varan_input.py   # thin CLI wrapper: argparse + dispatch only
```

**`vendor_adapters/common.py`** holds everything that is not specific to
any one vendor's file formats:
- `run_cmd`/`list_s3_files` — generic S3 listing helpers.
- `load_oncotree_dict` — the `dict.csv` (diagnosis-name → ONCOTREE_CODE)
  lookup. This is a Varan-side concept, not a vendor one, so it's shared.
- `SampleRow` — a dataclass holding one sample's row of `sample.tsv` data,
  with an `extra: dict` field for vendor-specific passthrough columns
  (see "Golden extra fields" below).
- `write_sample_tsv` — writes a batch of `SampleRow`s to a `sample.tsv` in
  Varan's `Templates/sample.tsv` column order, unioning any `extra` keys
  across rows into extra trailing columns.
- `FUSION_TABLE_HEADER`/`append_fusions_to_table` — writes into Varan's
  plain `FUSIONS/*.tsv` ingestion path (`fill_fusion_from_temp()` in
  `walk.py`). This is Varan's contract, not a vendor's, and every adapter
  should reuse it rather than re-implementing a fusion-table writer —
  see the docstring in `common.py` for why routing fusions this way
  (rather than through a synthesized `CombinedVariantOutput.tsv` /
  `comb_path`) matters: `_walk_setup()` in `walk.py` detects "does this
  run have CombinedOutput" by folder existence, not per-sample, so a
  non-empty `comb_path` for even one sample would silently switch the
  *entire run's* MSI/TMB handling away from the sample.tsv-driven
  VALUE+THR precedence path for every sample in the batch. This hazard is
  vendor-agnostic, which is why the fix lives in `common.py`.

**Each vendor adapter module** (e.g. `vendor_adapters/guardant.py`)
implements:
- `NAME: str` — the vendor identifier, used as the `--vendor` CLI value.
- `run(folder=None, selection=None, **kwargs) -> Optional[dict]` — a
  plain function, not a script: no `argparse`/`sys.argv` and no
  module-level mutable state involved in the conversion logic itself
  (`argparse` lives only in `create_Varan_input.py`'s CLI wrapper).
  Returns `{"report_path", "fusion_table_path", "report_data"}` or `None`.
  All vendor-specific defaults (file-system/S3 paths, VCF header
  templates, XML/TSV parsing) are plain function parameters with
  vendor-specific defaults, not hardcoded globals baked into the logic —
  this is what lets the function be called directly from a test, a
  notebook, or eventually a Snakemake rule with different paths, without
  editing the module.

**Adding a second vendor**: write `vendor_adapters/<vendor>.py`
implementing the interface above, add one line to `ADAPTERS` in
`vendor_adapters/__init__.py`. No changes to `guardant.py` are required,
and none of Guardant's verified bug fixes are at risk of being touched.

## `create_Varan_input.py` remains available standalone

Superseded as the *only* option by round 2's `varan.py --pipeline`, not
replaced by it.

Round 1's reasoning below (why vendor dispatch wasn't wired into
`varan.py`'s own CLI) was superseded by round 2 above: `varan.py
--pipeline <vendor>` now runs a vendor adapter automatically, so
`create_Varan_input.py` is no longer the *only* way to invoke one. It's
still kept around, for two reasons: (1) it's the only way to drive the
cross-run `selection.tsv` batch mode (`-s`) — round 2 only wired the
single-run-folder mode (`-f`'s equivalent) into `varan.py -i`, see
"Known limitation" above; (2) it's useful on its own for inspecting or
hand-editing a generated `sample.tsv`/`fusions.tsv` before feeding it to
`varan.py -i` manually, e.g. to debug a conversion without also running
the full downstream pipeline. Both entry points call the exact same
`vendor_adapters.<name>.run()` — neither duplicates the other's logic.

Original (round 1) reasoning for keeping `varan.py`'s own CLI free of
vendor-dispatch argument-parsing logic, which is why round 2's
`--pipeline` flag was designed as a thin dispatcher rather than
expanding `varan.py`'s argument surface further: `varan.py`'s own CLI
already owns a dense set of single-letter flags (`-f`, `-s`, `-c`, ...)
and its core pipeline logic is actively evolving with **no test suite**
backing it — so round 2 deliberately added exactly one new flag
(`--pipeline`) plus one new `varan.py`-local function
(`run_vendor_adapter()`), calling straight into the unchanged
`vendor_adapters` interface, rather than growing a parallel argument
surface inside `vendor_adapters` itself.

The `run()` functions are still plain, explicit-input/output functions
(not scripts hardwired to `sys.argv` or globals) specifically so this
doesn't fight a future Snakemake-DAG migration (`stabilize-dag2`, out of
scope here) — a rule could call `guardant.run(...)` (or a future
vendor's `run(...)`) directly.

## Verification performed (and its limits)

**No real vendor example data is available on this machine** — the
example run folder the original Guardant notes referenced
(`Varan/es. output run/`) no longer exists here, and no other real vendor
files were provided. This refactor was **not** re-verified against real
Guardant files end-to-end. What was actually done instead:

- Every function moved from `create_Varan_input.py` into
  `vendor_adapters/guardant.py` was moved with **no logic changes** —
  same branching, same string handling, same field names — only:
  relocating shared/generic pieces into `common.py`, converting the
  return shape from a positional list to `SampleRow`, and turning
  hardcoded module-level path constants into function parameters with
  those same values as defaults.
- Synthetic fixtures (not real vendor data — fabricated specifically to
  exercise the documented bug scenarios, not presented as real) were
  built and run through `process_vcf`, `get_msi_data`, `load_fusions`,
  `append_fusions_to_table`, and `write_sample_tsv` directly:
  - A synthetic VCF with a `FILTER=FAIL` SNV, a `FILTER=PASS` SNV, an
    `SVTYPE=BND` row, and two `SVTYPE=CNV` rows, paired against a 2-row
    synthetic `.cnv_call.hdr.tsv` (`call=1` and `call=0`). Result:
    exactly one SNV row survives (the PASS one), the BND row appears in
    neither output, and the two CNV rows are correctly paired by
    sequential index with `call=1` → `PASS`/`<DUP>` and `call=0` →
    `FAIL`/`.` — matching the documented, verified behavior exactly.
  - A synthetic `.msi_call.hdr.tsv` (`msi_score=27`, `msi_status=MSI-H`)
    → `MSI=27`, `MSI_THR=Unstable`, matching the real example's known
    values from the original notes.
  - A synthetic `.fusion_call.hdr.tsv` with two `call=1` rows and one
    `call=0` row → only the two `call=1` fusions are kept.
  - `SampleRow`/`write_sample_tsv` round-tripped through the "golden
    extra field" mechanism (an ad-hoc `BRCA_REVERSION` column set on one
    row, blank on another) and produced the expected sample.tsv shape.

This confirms the refactor preserves the documented logic against
representative synthetic inputs, but it is **not** a substitute for
re-running against the real Guardant example files that were used to
find and fix the original bugs — those files are gone from this machine.
Treat the Guardant adapter as logically equivalent to the reviewed
original script, not as freshly re-validated against real data.

## Guardant-specific notes (carried forward from `GUARDANT_INTEGRATION_NOTES.md`)

### Verified bug fixes (unchanged, now living in `vendor_adapters/guardant.py`)

1. **Germline contamination.** SNV VCF output only carries `FILTER=PASS`
   rows. In the real example VCF this was originally letting through 225
   `FILTER=FAIL`/germline rows against 118 real somatic `FILTER=PASS`
   calls.
2. **Dead filter removed.** The `if cn_value == 2.0: continue` check in
   `load_cnv_tsv_ordered` never fired (copy_number is continuous, e.g.
   2.07/1.84/3.17) and is gone; the real filtering is the `call` column
   check in `process_vcf`.
3/4. **`SVTYPE=BND` (fusion breakend) rows are explicitly skipped** before
   the structural/SNV split in `process_vcf`, so they neither leak into
   the SNV VCF nor silently consume a `cnv_index` slot meant for a real
   CNV row (which would otherwise misattribute copy numbers to the wrong
   gene if a future file interleaves BND rows among CNV rows).
5. **Real MSI numeric score.** `MSI` gets the real `msi_score`; `MSI_THR`
   is derived from `msi_status` via case-insensitive substring match
   (`MSI-H`/`UNSTABLE` → `Unstable`, `MSS`/`MSI-L`/`STABLE` → `Stable`),
   left blank (not guessed) for anything unrecognized.
6. **Fusions routed through the plain `FUSIONS/*.tsv` path**
   (`common.append_fusions_to_table`), never through a synthesized
   `CombinedVariantOutput.tsv`/`comb_path` — `comb_path` is always blank
   for Guardant samples.

### Open questions — NOT resolved, still need a human call

These are carried forward unchanged from the original notes. Nothing in
this refactor silently resolved any of them:

- **`fill_fusion_from_temp()`'s hardcoded `min_read_count = 15`** in
  `walk.py` (currently at line ~1530) is applied independent of
  `conf.ini`'s `THRESHOLD_FUSION` (used only by the CombinedOutput fusion
  path) and independent of Guardant's own `call=1` confidence flag. In
  the real example data, `RET-NCOA4` was a real Guardant-confirmed fusion
  (`call=1`) with only 10 supporting molecules — it would be **silently
  dropped** by this hardcoded threshold even though Guardant's own
  algorithm already validated it. Whether 15 is the right cutoff for
  Guardant's molecule-count scale (not necessarily comparable to
  Illumina's read-count scale, which the CombinedOutput path was
  presumably tuned against) is unresolved. If/when this gets addressed,
  do it as a deliberate, documented decision (e.g. a vendor-aware
  threshold or reading it from `conf.ini`), not a quiet number change.
- **`ONCOTREE_CODE` mapping via `dict.csv`** (`load_oncotree_dict` in
  `common.py`, called from `get_xml_data` in `guardant.py`): unchanged
  from the original script's logic, not reviewed further this pass either.
- **"Golden" extra fields (BRCA Reversion, AutoQC metrics).** The
  mechanism now exists structurally — `SampleRow.extra` /
  `write_sample_tsv` will emit any additional column a vendor populates,
  and Varan's own `write_clinical_sample()` already carries unknown
  extra `sample.tsv` columns straight through into
  `data_clinical_sample.txt` with no new Varan-side code needed. But this
  is **not wired into `guardant.py`'s actual parsing** — the real
  board-summary/AutoQC file (BRCA Reversion, germline pathogenic count,
  etc.) has still never been provided to verify column names/formats
  against. Extend `guardant.py`'s `process_single_sample` to populate
  `SampleRow.extra` once that file is available.
- **No local test data.** The real example run folder referenced by the
  original notes (`Varan/es. output run/`) no longer exists on this
  machine. `Templates/sample_guardant_example.tsv` and
  `Templates/sample_guardant_example_fusions.tsv` remain the only
  fixtures reflecting real (redacted) example values. Do not fabricate
  fake "real vendor" data and pass it off as such — the synthetic
  fixtures used to verify this refactor (see above) were built
  specifically to probe the documented bug scenarios, not to simulate a
  full real run.

### How to run Varan on Guardant samples

**Preferred (round 2/3):** set `[Vendor] PIPELINE = guardant` in conf.ini
(no CLI flag - see "Round 3" above), then run `python varan.py -i
<s3_run_folder_or_selection.tsv> -o <out> -c <cancer>` — runs the Guardant
adapter automatically (whole-folder or samplesheet mode, auto-detected
from whether `-i`'s argument is a file or a folder) and continues straight
into the full pipeline.

**Manual two-step (round 1, always available, e.g. to inspect/hand-edit
the generated sample.tsv before running the full pipeline):** run `python
create_Varan_input.py --vendor guardant -f <s3_run_folder>` (or `-s
<selection.tsv>`) to produce a `sample.tsv`-shaped report and a companion
`{run_id}_fusions.tsv` yourself, then pass those to `varan.py -i
sample.tsv "" fusions.tsv` with `[Vendor] PIPELINE` left blank in conf.ini
(this is the unspecified/legacy path once the files already exist). See
`Templates/sample_guardant_example.tsv` and
`Templates/sample_guardant_example_fusions.tsv` for the expected shape.

Column-by-column source mapping (produced by either path — both call the
same `vendor_adapters.guardant.run()`) is unchanged from the original
notes:

| Column | Guardant source | Notes |
|---|---|---|
| `SAMPLE_ID` | `finalmetadata.xml` AccessionId / `_metadata.xml` AccessionId | |
| `PATIENT_ID` | `finalmetadata.xml` SubjectId (or AccessionId again for the RUO/QCI-only path, which has no SubjectId) | |
| `ONCOTREE_CODE` | `finalmetadata.xml` Diagnosis, mapped via `dict.csv` | unreviewed open question, see above |
| `snv_path` | `{sample}.snv.vcf` | PASS-only, no germline contamination |
| `cnv_path` | `{sample}.cnv.vcf` | correctly gene-attributed |
| `comb_path` | — | always blank |
| `MSI` | `.msi_call.hdr.tsv` `msi_score` | real number |
| `TMB` | — | blank — Guardant360 CDx has no TMB field in either metadata XML variant checked |
| `MSI_THR` | `.msi_call.hdr.tsv` `msi_status`, mapped to Stable/Unstable | both VALUE and THR populated — Varan's VALUE+THR precedence trusts THR, logs a warning by design |
| `TMB_THR` | — | blank, same reasoning as `TMB` |

Fusions: separate `{run_id}_fusions.tsv`, passed as the third `-i`
argument to `varan.py`. See the `min_read_count=15` open question above
before trusting it blindly for a real batch.
