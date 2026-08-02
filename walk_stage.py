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
    python walk_stage.py setup --ctx-out out/.walk_ctx.pkl \\
        --input in_folder --output out -c luad
    python walk_stage.py cnv --ctx out/.walk_ctx.pkl
    python walk_stage.py snv --ctx out/.walk_ctx.pkl
    python walk_stage.py fusion --ctx out/.walk_ctx.pkl
    python walk_stage.py clinical --ctx out/.walk_ctx.pkl
"""
from __future__ import annotations

import argparse
import pickle
from pathlib import Path

from loguru import logger

import walk

STAGE_FUNCS = {
    "cnv": walk._walk_process_cnv,
    "snv": walk._walk_process_snv,
    "fusion": walk._walk_process_fusion,
    "clinical": walk._walk_write_clinical_tables,
}


def _cmd_setup(args: argparse.Namespace) -> None:
    ctx = walk._walk_setup(
        input_path=args.input,
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
    ctx_out = Path(args.ctx_out)
    ctx_out.parent.mkdir(parents=True, exist_ok=True)
    with ctx_out.open("wb") as f:
        pickle.dump(ctx, f)
    logger.success(f"walk_stage setup complete, context written to {ctx_out}")


def _cmd_stage(stage_name: str, args: argparse.Namespace) -> None:
    ctx_path = Path(args.ctx)
    with ctx_path.open("rb") as f:
        ctx: walk.WalkContext = pickle.load(f)

    # isinputfile is a bare module-level global in walk.py, set once by
    # _walk_setup() in its own process - it does not survive across the
    # separate `python walk_stage.py ...` process this stage runs in, so
    # it must be restored from the pickled context before calling the
    # stage function. See WalkContext's docstring in walk.py.
    walk.isinputfile = ctx.isinputfile

    STAGE_FUNCS[stage_name](ctx)

    done_marker = ctx_path.parent / f".{stage_name}.done"
    done_marker.touch()
    logger.success(f"walk_stage {stage_name} complete.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
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

    for stage_name in STAGE_FUNCS:
        stage_p = sub.add_parser(
            stage_name, help=f"Run the {stage_name} stage from a pickled context.")
        stage_p.add_argument("--ctx", required=True)

    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    if args.command == "setup":
        if args.vcf_type not in (None, "snv", "cnv", "fus", "tab"):
            msg = f"--vcf-type must be one of snv, cnv, fus, tab (got {args.vcf_type!r})"
            raise ValueError(msg)
        _cmd_setup(args)
    else:
        _cmd_stage(args.command, args)


if __name__ == "__main__":
    main()
