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
walk_clinical too - the one real remaining parallelism gap noted in the
Snakefile's docstring after the WALK split.

No pipeline logic lives here - every stage calls straight through to
varan.py's own filter_main/concatenate_main/meta_case_main/validate_output,
so this can never drift from what `python varan.py -i ...` does end to end.

Usage (see Snakefile for the real rules):
    python pipeline_stage.py filter --input ... --output ... -c ... [-C conf.ini]
    python pipeline_stage.py concatenate --output ... [-C conf.ini]
    python pipeline_stage.py tables --output ... -c ... [-C conf.ini]
    python pipeline_stage.py validate --input ... --output ... -c ... [-C conf.ini]
"""
from __future__ import annotations

import argparse
from pathlib import Path


def _cmd_filter(args: argparse.Namespace) -> None:
    from filter_clinvar import filter_main
    filter_main(
        args.input, args.output, args.output, args.oncokb,
        args.filters, args.cancer, args.resume)


def _cmd_concatenate(args: argparse.Namespace) -> None:
    from concatenate import concatenate_main
    concatenate_main(args.filters, args.output, "maf", args.oncokb)


def _cmd_tables(args: argparse.Namespace) -> None:
    from Make_meta_and_cases import meta_case_main
    meta_case_main(args.cancer, args.output)


def _cmd_validate(args: argparse.Namespace) -> None:
    from ValidateFolder import validate_output
    validate_output(
        args.output, args.input, args.multiple, False,
        args.cancer, args.oncokb, args.filters)


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

    filter_p = sub.add_parser("filter", help="Stage 2: MAF filtering.")
    filter_p.add_argument("--input", nargs="+", required=True)
    filter_p.add_argument("--output", required=True)
    filter_p.add_argument("-c", "--cancer", required=True)
    filter_p.add_argument("--oncokb", action="store_true")
    filter_p.add_argument("--filters", default="")
    filter_p.add_argument("--resume", action="store_true")

    concat_p = sub.add_parser("concatenate", help="Stage 3: concatenate MAFs.")
    concat_p.add_argument("--output", required=True)
    concat_p.add_argument("--oncokb", action="store_true")
    concat_p.add_argument("--filters", default="")

    tables_p = sub.add_parser("tables", help="Stage 4: meta/case list files.")
    tables_p.add_argument("--output", required=True)
    tables_p.add_argument("-c", "--cancer", required=True)

    validate_p = sub.add_parser("validate", help="Stage 5: cBioPortal validation + report.")
    validate_p.add_argument("--input", nargs="+", required=True)
    validate_p.add_argument("--output", required=True)
    validate_p.add_argument("-c", "--cancer", required=True)
    validate_p.add_argument("--oncokb", action="store_true")
    validate_p.add_argument("--filters", default="")
    validate_p.add_argument("--multiple", action="store_true")

    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    # Same ordering constraint as varan.py's own __main__ block: conf.ini
    # path must be set before any Varan module is imported, since most of
    # them read it at import time - hence the deferred `from X import Y`
    # inside each _cmd_* function above instead of top-of-file imports.
    from config_loader import set_config_path
    set_config_path(args.config)

    STAGE_FUNCS[args.command](args)

    Path(args.output, f".{args.command}.done").touch()


if __name__ == "__main__":
    main()
