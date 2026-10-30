#Copyright 2025 bioinformatics-policlinicogemelli

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

"""CLI entry point exposing varan.py's create-flow stages 2-5 (FILTER /
CONCATENATE / MAKE AND POPULATE TABLES / VALIDATION - see varan.py's own
numbered comments) as independent Snakemake-invokable steps.

Unlike walk_stage.py's four WALK stages (CNV/SNV/fusion/clinical), these
four are NOT mutually independent - each genuinely needs the previous
stage's output on disk (filter needs walk's maf/ output, concatenate needs
filter's output, tables needs the concatenated data, validate needs
everything). Splitting them into separate rules here isn't about running
them in parallel with each other - it's so Snakemake gets real per-stage
resume-from-failure (a validation failure doesn't force re-running filter/
concatenate/tables from scratch) and so `filter` can start as soon as
walk_snv finishes, instead of waiting for walk_cnv/walk_fusion/
walk_clinical too.

Every stage reads its options from walk_stage.py setup's pickled
WalkContext (--ctx) instead of taking them again on the command line: the
real study folder is the versioned one _walk_setup() created (e.g.
`out_v1`, not the `out` given to setup), and with a vendor pipeline the
input is the adapter's converted sample.tsv, not the raw -i - only the
context knows either. Done markers go next to the context file, like
walk_stage.py's.

No pipeline logic lives here - every stage calls straight through to
varan.py's own filter_main/concatenate_main/meta_case_main/validate_output
with the same arguments and the same -t/--analysis_type gating varan()
applies, so this can never drift from what `python varan.py -i ...` does
end to end.

Usage (see Snakefile for the real rules):
    python pipeline_stage.py -C conf.ini filter --ctx out/.walk_ctx.pkl
    python pipeline_stage.py -C conf.ini concatenate --ctx out/.walk_ctx.pkl
    python pipeline_stage.py -C conf.ini tables --ctx out/.walk_ctx.pkl
    python pipeline_stage.py -C conf.ini validate --ctx out/.walk_ctx.pkl
"""
from __future__ import annotations

import argparse
from pathlib import Path

from loguru import logger

from walk_stage import load_ctx, prepare_process, run_reported

# varan() skips filter and concatenate for these -t values.
NO_SNV_TYPES = ["cnv", "fus", "tab"]

STAGE_LABELS = {
    "filter": "MAF filtering",
    "concatenate": "Concatenating mutation file",
    "tables": "Creating tables",
    "validate": "Validation",
}


def _cmd_filter(ctx) -> None:
    from filter_clinvar import filter_main
    if ctx.vcf_type in NO_SNV_TYPES:
        logger.info(f"Skipping MAF filtering (analysis type '{ctx.vcf_type}').")
        return
    filter_main(
        ctx.input_path, ctx.output_folder, ctx.output_folder, ctx.oncokb,
        ctx.filters, ctx.cancer, ctx.resume)


def _cmd_concatenate(ctx) -> None:
    from concatenate import concatenate_main
    if not (Path(ctx.output_folder) / "maf").exists() or ctx.vcf_type in NO_SNV_TYPES:
        logger.info("Skipping MAF concatenation (no maf/ folder or no SNV analysis).")
        return
    concatenate_main(ctx.filters, ctx.output_folder, "maf", ctx.oncokb)


def _cmd_tables(ctx) -> None:
    from Make_meta_and_cases import meta_case_main
    meta_case_main(ctx.cancer, ctx.output_folder)


def _cmd_validate(ctx) -> None:
    from ValidateFolder import validate_output
    validate_output(
        ctx.output_folder, ctx.input_path, ctx.multiple, False, ctx.cancer,
        ctx.oncokb, ctx.filters, ctx.start_time, ctx.vcf_type)
    # Only on success, matching varan.py's __main__: a failed run's vendor
    # scratch folder is left in place for inspection.
    if ctx.vendor_scratch_dir:
        from walk import clear_scratch
        clear_scratch(ctx.vendor_scratch_dir)


STAGE_FUNCS = {
    "filter": _cmd_filter,
    "concatenate": _cmd_concatenate,
    "tables": _cmd_tables,
    "validate": _cmd_validate,
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-C", "--config", default="conf.ini",
        help="Path to the conf.ini file to use for this run (default: ./conf.ini)")
    sub = parser.add_subparsers(dest="command", required=True)
    for name, label in STAGE_LABELS.items():
        stage_p = sub.add_parser(name, help=f"{label} (varan.py create flow).")
        stage_p.add_argument("--ctx", required=True)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    # Same ordering constraint as varan.py's own __main__ block: conf.ini
    # path and the reconciled sample type must be set before any Varan
    # module is imported, since most of them read conf.ini at import time -
    # hence the deferred `from X import Y` inside each _cmd_* function above.
    prepare_process(args.config)

    ctx_path = Path(args.ctx)
    ctx = load_ctx(ctx_path)
    run_reported(
        args.command, STAGE_LABELS[args.command], ctx.output_folder,
        ctx.cancer, ctx.start_time, lambda: STAGE_FUNCS[args.command](ctx))

    (ctx_path.parent / f".{args.command}.done").touch()
    logger.success(f"pipeline_stage {args.command} complete.")


if __name__ == "__main__":
    main()
