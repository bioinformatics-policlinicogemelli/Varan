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

"""CLI entry point exposing walk_folder's internal stages as independent,
Snakemake-invokable steps: setup / cnv / snv / fusion / clinical.

walk.py splits what used to be one ~400-line walk_folder() into a setup
stage plus four independent stages (CNV, SNV, fusion+splice, clinical
tables) that each depend only on setup's output, never on each other -
see WalkContext's docstring in walk.py. This script is the thin CLI shim
that lets Snakemake actually exploit that: `setup` runs once and pickles
the resulting WalkContext to disk, then cnv/snv/fusion/clinical each load
that same pickle and run their one stage, so `snakemake -j 4` can run all
four in parallel as real, independent rules instead of one opaque
shell-out to `python varan.py`.

No pipeline logic lives here - every stage call goes straight through to
the exact same _walk_process_cnv/_walk_process_snv/_walk_process_fusion/
_walk_write_clinical_tables functions the plain `python varan.py -i ...`
path uses, so this can never drift from that behavior.

Usage (paths below are illustrative - see Snakefile for the real rules):
    python walk_stage.py -C custom_conf.ini setup --ctx-out out/.walk_ctx.pkl \\
        --input in_folder --output out -c luad
    python walk_stage.py -C custom_conf.ini cnv --ctx out/.walk_ctx.pkl
    python walk_stage.py -C custom_conf.ini snv --ctx out/.walk_ctx.pkl
    python walk_stage.py -C custom_conf.ini fusion --ctx out/.walk_ctx.pkl
    python walk_stage.py -C custom_conf.ini clinical --ctx out/.walk_ctx.pkl

`-C/--config` must come before the subcommand (it's a top-level option, not
per-subcommand) and must be repeated on every one of these calls: each runs
in its own process, so each needs to be told the conf.ini path again - see
config_loader.py - and each one also re-runs varan.py's [Vendor] PIPELINE
resolution / [Sample_Type] TYPE reconciliation (see prepare_process()).
`setup` additionally runs the vendor adapter (e.g. Guardant) and varan.py's
CLI-consistency checks, exactly like `python varan.py -i ...` does.
`walk` is deliberately NOT imported at module level here:
it (like every other Varan module) reads conf.ini via config_loader at its
own import time, so importing it before `set_config_path()` has run would
silently lock this process to the default conf.ini regardless of -C,
exactly the bug that broke -C for these stages before this fix.
"""
from __future__ import annotations

import argparse
import pickle
from collections.abc import Callable
from datetime import datetime
from pathlib import Path

from loguru import logger

from config_loader import set_config_path
from log_utils import attach_file_log_sink

STAGE_NAMES = ["cnv", "snv", "fusion", "clinical"]

# Same labels varan.py's own create flow uses for write_report_failure(),
# so a failed Snakemake run's report reads the same as a plain one's.
STAGE_LABELS = {
    "setup": "Preparing study folder",
    "cnv": "Preparing study folder (CNV)",
    "snv": "Preparing study folder (SNV)",
    "fusion": "Preparing study folder (fusions)",
    "clinical": "Preparing study folder (clinical tables)",
}


def prepare_process(config_path: str) -> str:
    """Do, in this process, what varan.py's __main__ does before importing
    any conf.ini-dependent module: set the conf.ini path, then resolve
    [Vendor] PIPELINE and reconcile [Sample_Type] TYPE with it. Every stage
    process needs this, not just setup: the reconciled TYPE only lives in
    this process's in-memory conf.ini, and walk.py reads it at import time.
    Returns the resolved PIPELINE value. Shared with pipeline_stage.py."""
    set_config_path(config_path)
    from varan import resolve_and_reconcile_pipeline
    return resolve_and_reconcile_pipeline()


def load_ctx(ctx_path: Path):
    """Load setup's pickled WalkContext, restoring from it the module-level
    globals setup's process set: walk.isinputfile and
    versioning.old_version_exists. They do not survive into this process,
    so they must be restored before calling any stage function. See
    WalkContext's docstring in walk.py."""
    import versioning
    import walk

    with ctx_path.open("rb") as f:
        ctx: walk.WalkContext = pickle.load(f)
    walk.isinputfile = ctx.isinputfile
    versioning.old_version_exists = ctx.old_version_exists
    return ctx


def run_reported(
    stage_name: str, label: str, output_folder: str, cancer: str | None,
    start_time: str, func: Callable[[], None],
) -> None:
    """Run one stage with this process's log copied into the study folder,
    and write report_VARAN.html via write_report_failure() if it fails -
    the per-process equivalent of the try/except around varan.py's create
    flow. Shared with pipeline_stage.py."""
    attach_file_log_sink(
        Path(output_folder), f"Varan_{stage_name}_{{time:YYYY-MM-DD_HH-mm-ss.SS}}.log")
    try:
        func()
    except Exception as err:
        from write_report import write_report_failure
        logger.critical(f"Analysis failed during '{label}': {err}")
        write_report_failure(
            output_folder, err, cancer=cancer, start_time=start_time,
            stage=label)
        raise


def _cmd_setup(args: argparse.Namespace, pipeline: str) -> None:
    from varan import create_args_error, run_vendor_adapter
    from vendor_adapters import ADAPTERS

    args_error = create_args_error(
        args.filters, args.oncokb, args.resume, args.overwrite)
    if args_error:
        raise ValueError(args_error)

    start_time = datetime.now().astimezone().strftime("%d/%m/%Y, %H:%M:%S")

    def setup() -> None:
        import versioning
        import walk

        # Same as varan.py's __main__: a vendor pipeline's raw input is
        # converted to sample.tsv/fusion.tsv first, and its scratch folder
        # is only removed once the whole run has succeeded (see
        # pipeline_stage.py's validate stage).
        varan_input, scratch_dir = args.input, ""
        if pipeline in ADAPTERS:
            varan_input, scratch_dir = run_vendor_adapter(
                pipeline, args.input, args.output)

        ctx = walk._walk_setup(
            input_path=varan_input,
            multiple=args.multiple,
            output_folder=args.output,
            oncokb=args.oncokb,
            cancer=args.cancer,
            overwrite_output=args.overwrite,
            resume=args.resume,
            vcf_type=args.vcf_type,
            filters=args.filters or "",
            sigma=args.sigma,
        )
        ctx.start_time = start_time
        ctx.vendor_scratch_dir = scratch_dir
        ctx.old_version_exists = getattr(versioning, "old_version_exists", False)
        attach_file_log_sink(
            Path(ctx.output_folder), "Varan_setup_{time:YYYY-MM-DD_HH-mm-ss.SS}.log")

        ctx_out = Path(args.ctx_out)
        ctx_out.parent.mkdir(parents=True, exist_ok=True)
        with ctx_out.open("wb") as f:
            pickle.dump(ctx, f)
        logger.success(f"walk_stage setup complete, context written to {ctx_out}")

    # The versioned study folder doesn't exist until _walk_setup() has run,
    # so a setup failure's report goes to --output itself - the same place
    # varan.py's own create flow writes it when walk_folder() fails.
    try:
        setup()
    except Exception as err:
        from write_report import write_report_failure
        label = STAGE_LABELS["setup"]
        logger.critical(f"Analysis failed during '{label}': {err}")
        write_report_failure(
            args.output, err, cancer=args.cancer, start_time=start_time,
            stage=label)
        raise


def _cmd_stage(stage_name: str, args: argparse.Namespace) -> None:
    import walk

    stage_funcs = {
        "cnv": walk._walk_process_cnv,
        "snv": walk._walk_process_snv,
        "fusion": walk._walk_process_fusion,
        "clinical": walk._walk_write_clinical_tables,
    }

    ctx_path = Path(args.ctx)
    ctx = load_ctx(ctx_path)

    run_reported(
        stage_name, STAGE_LABELS[stage_name], ctx.output_folder, ctx.cancer,
        ctx.start_time, lambda: stage_funcs[stage_name](ctx))

    done_marker = ctx_path.parent / f".{stage_name}.done"
    done_marker.touch()
    logger.success(f"walk_stage {stage_name} complete.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-C", "--config", default="conf.ini",
        help="Path to the conf.ini file this stage should use. Must come "
        "before the subcommand.")
    sub = parser.add_subparsers(dest="command", required=True)

    setup_p = sub.add_parser("setup", help="Resolve input/output paths once.")
    setup_p.add_argument("--input", nargs="+", required=True)
    setup_p.add_argument("--output", required=True)
    setup_p.add_argument("-c", "--cancer", required=True)
    setup_p.add_argument("--multiple", action="store_true")
    setup_p.add_argument("--oncokb", action="store_true")
    setup_p.add_argument("--sigma", action="store_true")
    setup_p.add_argument("--overwrite", action="store_true")
    setup_p.add_argument("--resume", action="store_true")
    setup_p.add_argument(
        "--vcf-type", default=None,
        help="One of snv, cnv, fus, tab. Omit to auto-detect.")
    setup_p.add_argument("--filters", default="")
    setup_p.add_argument("--ctx-out", required=True)

    for stage_name in STAGE_NAMES:
        stage_p = sub.add_parser(
            stage_name, help=f"Run the {stage_name} stage from a pickled context.")
        stage_p.add_argument("--ctx", required=True)

    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    pipeline = prepare_process(args.config)

    if args.command == "setup":
        if args.vcf_type not in (None, "snv", "cnv", "fus", "tab"):
            msg = f"--vcf-type must be one of snv, cnv, fus, tab (got {args.vcf_type!r})"
            raise ValueError(msg)
        _cmd_setup(args, pipeline)
    else:
        _cmd_stage(args.command, args)


if __name__ == "__main__":
    main()
