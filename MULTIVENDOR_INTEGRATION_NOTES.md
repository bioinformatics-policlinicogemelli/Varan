# Multi-vendor input generation for Varan — architecture and status

## Round 2: `varan.py -i` runs a vendor's raw input directly

The first round (below) built `vendor_adapters/` but left it a manual,
two-step workflow: run `create_Varan_input.py --vendor guardant ...`
yourself, then feed its output to `varan.py -i sample.tsv patient.tsv
fusions.tsv`. This round wires the adapter directly into `varan.py`, so
pointing it at a vendor's raw input and running it "just works", the same
way the existing native (Illumina/CombinedOutput) path already does.

### Selection mechanism: explicit, not auto-detected

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

**Naming deviation, explained.** The user's own suggested spelling was
`--pipeline illumina/guardant`. This implementation uses `"native"`
instead of `"illumina"` for the default/non-vendor value, because the
existing default path isn't inherently Illumina-specific — it's simply
"input already in Varan's own canonical sample.tsv/folder shape",
regardless of which sequencer or vendor produced the files that shape
wraps. Calling it `"illumina"` would misleadingly vendor-label the one
option that is explicitly *not* vendor-specific, which cuts against the
whole point of this rename effort (moving away from vendor-specific
naming, guardant → multivendor). If this reasoning doesn't hold up under
real-world usage expectations, renaming `"native"` back to `"illumina"`
(or something else) is a one-string change in `varan.py`'s `--pipeline`
`choices=[...]` plus its help text and `resolve_pipeline()`'s docstring —
nothing structural depends on the exact spelling.

### How the raw-folder-to-canonical dispatch works

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
strictly before the existing `varan(...)` call. `--pipeline native` (or
omitting the flag with no `[Vendor] PIPELINE` set in conf.ini) reaches
`varan(...)` with `varan_input` completely untouched from what argparse
parsed — identical to every prior behavior, including anyone still using
the manual two-step Guardant flow from round 1 (that flow still works
unchanged: it just produces a sample.tsv/fusions.tsv pair you can pass to
`-i` yourself, exactly as before, with `--pipeline` simply never entering
the picture).

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

**Preferred (round 2):** `python varan.py -i <s3_run_folder> -o <out> -c
<cancer> --pipeline guardant` — runs the Guardant adapter automatically
against the raw S3 run folder and continues straight into the full
pipeline. See "Round 2" above for the mechanics.

**Manual two-step (round 1, still available for a single run folder or
for the cross-run `selection.tsv` batch mode):** run `python
create_Varan_input.py --vendor guardant -f <s3_run_folder>` (or `-s
<selection.tsv>`) to produce a `sample.tsv`-shaped report and a companion
`{run_id}_fusions.tsv` yourself, then pass those to `varan.py -i
sample.tsv "" fusions.tsv` (no `--pipeline` needed - this is the native
path once the files exist). See `Templates/sample_guardant_example.tsv`
and `Templates/sample_guardant_example_fusions.tsv` for the expected
shape.

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
