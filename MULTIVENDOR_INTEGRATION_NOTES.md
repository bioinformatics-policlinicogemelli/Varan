# Multi-vendor input generation for Varan — architecture and status

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

## Why this stays "extra tooling", not wired into `varan.py`

`varan.py`'s own CLI already owns a dense set of single-letter flags
(`-f`, `-s`, `-c`, ... — see its `-i`/`--varan_input`, which is exactly
what this preprocessing step's output feeds into) and its core pipeline
logic is actively evolving with **no test suite** backing it. Wiring
vendor dispatch directly into `varan.py`'s argument parser would add risk
to that already-complex, untested surface for no real benefit — this
preprocessing step runs once, before `-i`, as a separate CLI invocation,
and there's no workflow reason to merge the two. `create_Varan_input.py`
remains a standard, documented preprocessing step run ahead of
`varan.py -i sample.tsv patient.tsv [fusions.tsv]`, now organized as a
proper package instead of a single monolithic script.

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

### How to compile `sample.tsv` for Guardant samples

Run `python create_Varan_input.py --vendor guardant -f <s3_run_folder>` (or
`-s <selection.tsv>`) — it produces a `sample.tsv`-shaped report and a
companion `{run_id}_fusions.tsv`. See
`Templates/sample_guardant_example.tsv` and
`Templates/sample_guardant_example_fusions.tsv` for the expected shape.
Column-by-column source mapping is unchanged from the original notes:

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
