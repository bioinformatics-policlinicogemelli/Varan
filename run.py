#!/usr/bin/env python3

"""
Launcher for Varan's Snakemake workflow (see Snakefile).

Wraps `snakemake` with the project-scoped locking, PBS cluster submission
and run-summary conveniences used elsewhere on this cluster - adapted here
for Varan's four real operations (create/update/extract/remove, matching
varan.py's own -i/-u/-e/-r modes 1:1; see config.yaml for each one's
settings).

Examples
--------
Dry-run a new study creation (Snakemake's own dry-run - not to be confused
with varan.py's own -D/--dry-run preflight check, which validates conf.ini/
inputs/output-folder versioning and is a separate thing):
    ./run.py -w create -c config.yaml -q 4 -n

Actually build a new study (runs walk_setup, then CNV/SNV/fusion/clinical
in parallel across the 4 cores, then the rest of the pipeline):
    ./run.py -w create -c config.yaml -q 4

Merge a new batch of samples into an existing study:
    ./run.py -w update -c config.yaml -q 1

Extract/remove a subset of samples:
    ./run.py -w extract -c config.yaml -q 1
    ./run.py -w remove -c config.yaml -q 1

Print the DAG for a create run:
    ./run.py -w create -c config.yaml -q 4 --dag | dot -Tsvg > create_dag.svg

Unlock working directory after an interrupted run:
    ./run.py -c config.yaml -q 1 --unlock

Submit each job as its own PBS job instead of running everything in one
local allocation:
    ./run.py -w create -c config.yaml -q 1 -cl -qu workq -j 10

Run an internal stage directly (debugging only - see --allow-custom-target;
walk_cnv/walk_snv/walk_fusion/walk_clinical all require walk_setup's context
pickle to already exist):
    ./run.py -w walk_cnv -c config.yaml -q 1 --allow-custom-target
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import snakemake
from colorama import Fore, Style, init


try:
    import yaml
except ImportError:
    yaml = None


init(autoreset=True)

THIS_DIR = Path(__file__).resolve().parent

# PBS resource fallback for rules that don't declare their own `resources:
# mem_mb`/`runtime` -- only used when --cluster is passed; local (non-cluster)
# runs are unaffected. Starting placeholder, not a benchmarked number - VEP/
# vcf2maf (walk_snv) in particular can need considerably more than this on
# real samples; tune once real per-stage memory/runtime profiles exist (same
# "benchmark before committing to a number" caveat as SIGMA_INTEGRATION_
# FEASIBILITY.md's open question #1).
CLUSTER_DEFAULT_MEM_MB = 8000
CLUSTER_DEFAULT_RUNTIME_MIN = 120


def acquire_project_lock(configfile: Path):
    # Snakemake's own default locking locks the whole *working directory*
    # (.snakemake/locks/), not just the project a given invocation actually
    # touches. A shared cluster checkout of this repo can have several
    # people (or several projects for the same person) each with their own
    # config.yaml - two concurrent runs against DIFFERENT configs would
    # otherwise serialize every scheduling step through that one shared
    # directory-wide lock.
    #
    # Fix: take our own lock keyed on the *config file's* absolute path
    # (i.e. per-project, not per-directory) via flock, then pass lock=False
    # to bypass Snakemake's own directory-wide one entirely. flock is held
    # for this process's entire lifetime and is released automatically by
    # the kernel on any exit (clean, crashed, or killed -- no stale-lock
    # cleanup step needed, unlike Snakemake's own file-based lock which
    # requires an explicit --unlock after a killed process).
    #
    # Still refuses to run the SAME config twice concurrently -- that
    # protection is preserved, just scoped to the one project actually at
    # risk instead of the whole repo.
    lock_dir = THIS_DIR / ".smk_locks"
    lock_dir.mkdir(exist_ok=True)
    key = hashlib.sha1(str(configfile.resolve()).encode()).hexdigest()[:16]
    lock_path = lock_dir / f"{key}.lock"
    fh = open(lock_path, "w")
    try:
        fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        die(
            f"another run.py invocation is already active for config "
            f"'{configfile}' (lock: {lock_path}). Refusing to start a second "
            "one against the same project -- two concurrent invocations of "
            "the same config can race to declare/clean up the same output "
            "files."
        )
    fh.write(str(os.getpid()))
    fh.flush()
    return fh  # caller must keep this referenced for the process's lifetime


# Matches the Snakefile's own top-level rule names 1:1, which in turn match
# varan.py's own four modes (-i/-u/-e/-r). The internal per-datatype fan-out
# rules (walk_setup/walk_cnv/walk_snv/walk_fusion/walk_clinical) are
# deliberately NOT listed here - they're plumbing `create` depends on, not
# something a user picks directly (walk_cnv etc. all require walk_setup's
# context pickle to already exist). Reach them anyway, for debugging, via
# --allow-custom-target.
WORKFLOWS: Dict[str, str] = {
    "create": "Build a new study from raw input (varan.py -i): walk_setup "
              "-> {cnv, snv, fusion, clinical} in parallel -> filter/"
              "concatenate/tables/validate",
    "update": "Merge a new batch of samples into an existing study (varan.py -u)",
    "extract": "Extract a subset of samples into a new study (varan.py -e)",
    "remove": "Remove a subset of samples, producing a new study version (varan.py -r)",
}


LOGO = rf"""
{Fore.GREEN}
      __     {Fore.YELLOW}VARAN — Snakemake launcher{Fore.GREEN}
   .-\'  \'-.
  /  o    o  \\
 |  .------.  |
  \\  \\____/  /
   \'-.____.-\'
{Style.RESET_ALL}
"""


def color(text: str, colour: str) -> str:
    return f"{colour}{text}{Style.RESET_ALL}"


def die(message: str, exit_code: int = 1) -> None:
    sys.stderr.write(color(f"ERROR: {message}\n", Fore.RED))
    raise SystemExit(exit_code)


def warn(message: str) -> None:
    sys.stderr.write(color(f"WARNING: {message}\n", Fore.YELLOW))


def info(message: str) -> None:
    sys.stderr.write(color(f"{message}\n", Fore.CYAN))


def success(message: str) -> None:
    sys.stderr.write(color(f"{message}\n", Fore.GREEN))


def load_config(configfile: Path) -> Dict[str, Any]:
    if yaml is None:
        warn("PyYAML is not installed; config summary will be skipped.")
        return {}

    try:
        with configfile.open("r") as handle:
            data = yaml.safe_load(handle)
    except Exception as exc:
        warn(f"Could not parse config file for summary: {exc}")
        return {}

    return data or {}


def print_workflows() -> None:
    print(LOGO)
    print(color("Available workflow targets:", Fore.GREEN))
    print()
    for name, description in WORKFLOWS.items():
        print(f"  {color(name, Fore.YELLOW):<20} {description}")
    print()


def parse_resources(resource_args: Optional[List[str]]) -> Optional[Dict[str, int]]:
    if not resource_args:
        return None

    resources: Dict[str, int] = {}

    for item in resource_args:
        if "=" not in item:
            die(f"Invalid resource format: {item}. Use key=value, e.g. igv=1 mem_mb=64000")

        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()

        if not key:
            die(f"Invalid resource key in: {item}")

        try:
            resources[key] = int(value)
        except ValueError:
            die(f"Resource value must be an integer in: {item}")

    return resources


def build_cluster_cmd(queue: str, logdir: Path) -> str:
    # Braces doubled ({{...}}) are Snakemake's own per-job placeholders and
    # must survive this .format() call untouched; only {queue}/{logdir} are
    # substituted here. OpenPBS 20.0.1 / PBS Pro syntax: `-l
    # select=1:ncpus=X:mem=Ygb`, not the older Torque `-l nodes=1:ppn=X` form.
    template = (
        "qsub -q {queue} "
        "-l select=1:ncpus={{threads}}:mem={{resources.mem_mb}}mb "
        "-l walltime={{resources.runtime}}:00 "
        "-N smk.{{rule}}.{{jobid}} "
        "-o {logdir}/{{rule}}.{{jobid}}.out "
        "-e {logdir}/{{rule}}.{{jobid}}.err"
    )
    return template.format(queue=queue, logdir=logdir)


def build_singularity_args(args: argparse.Namespace) -> str:
    # No rule in the Snakefile declares `container:` yet, so this is inert
    # for now - kept for parity with this cluster's usual launcher shape,
    # ready to use the day a rule needs a container.
    pieces: List[str] = []

    if args.singularity_args:
        pieces.append(args.singularity_args.strip())

    for bind_path in args.bind:
        pieces.append(f"--bind {bind_path}")

    return " ".join(pieces).strip()


def _target_output_folders(targets: List[str], config: Dict[str, Any]) -> Dict[str, str]:
    """Resolve each requested target's own output_folder from config.yaml,
    for the run summary - each of create/update/extract/remove has its own
    section (see config.yaml)."""
    folders = {}
    for target in targets:
        section = config.get(target)
        if isinstance(section, dict) and "output_folder" in section:
            folders[target] = section["output_folder"]
    return folders


def print_run_summary(
    args: argparse.Namespace,
    snakefile: Path,
    configfile: Path,
    config: Dict[str, Any],
    singularity_args: str,
    resources: Optional[Dict[str, int]],
) -> None:
    print(LOGO)
    print(color("Run summary", Fore.GREEN))
    print(color("-----------", Fore.GREEN))
    print(f"{color('Snakefile:', Fore.YELLOW)}       {snakefile}")
    print(f"{color('Config:', Fore.YELLOW)}          {configfile}")
    print(f"{color('Targets:', Fore.YELLOW)}         {', '.join(args.workflow or [])}")
    print(f"{color('Cores:', Fore.YELLOW)}           {args.cores}")
    print(f"{color('Workdir:', Fore.YELLOW)}         {args.directory or os.getcwd()}")
    print(f"{color('conf.ini:', Fore.YELLOW)}        {config.get('conf_path', 'not found in config')}")
    print(f"{color('Conda env:', Fore.YELLOW)}       {config.get('conda_env', 'not found in config')}")

    for target, folder in _target_output_folders(args.workflow or [], config).items():
        print(f"{color(f'  {target} output:', Fore.YELLOW):<20} {folder}")

    print(f"{color('Dry run:', Fore.YELLOW)}         {args.dry_run}")
    print(f"{color('Use conda:', Fore.YELLOW)}       {not args.no_conda}")
    print(f"{color('Use Singularity:', Fore.YELLOW)} {not args.no_singularity}")

    if singularity_args:
        print(f"{color('Singularity args:', Fore.YELLOW)} {singularity_args}")

    if resources:
        print(f"{color('Resources:', Fore.YELLOW)}       {resources}")

    if args.forceall:
        print(color("Force all jobs:    yes", Fore.YELLOW))

    if args.rerun_incomplete:
        print(color("Rerun incomplete:  yes", Fore.YELLOW))

    if args.keep_going:
        print(color("Keep going:        yes", Fore.YELLOW))

    if args.cluster:
        print(f"{color('PBS cluster:', Fore.YELLOW)}     enabled (queue={args.queue}, max concurrent jobs={args.jobs})")
        print(f"{color('Restart times:', Fore.YELLOW)}   {args.restart_times}")

    if args.rerun_triggers:
        print(f"{color('Rerun triggers:', Fore.YELLOW)}  {', '.join(args.rerun_triggers)} (default: mtime, params, input, software-env, code)")

    print()


def validate_args(args: argparse.Namespace, snakefile: Path, configfile: Path) -> None:
    if not snakefile.exists():
        die(f"Cannot find Snakefile: {snakefile}")

    if not configfile.exists():
        die(f"Cannot find config file: {configfile}")

    if not args.unlock and not args.workflow:
        die(f"No workflow target provided. Use -w {', '.join(WORKFLOWS.keys())}.")

    if args.workflow and not args.allow_custom_target:
        invalid = [target for target in args.workflow if target not in WORKFLOWS]
        if invalid:
            valid = ", ".join(WORKFLOWS.keys())
            die(
                f"Unknown workflow target(s): {', '.join(invalid)}. "
                f"Valid targets are: {valid}. "
                "Use --allow-custom-target to run an internal stage name "
                "(walk_setup/walk_cnv/walk_snv/walk_fusion/walk_clinical) "
                "directly for debugging."
            )


def run_snakemake(args: argparse.Namespace) -> int:
    snakefile = Path(args.snakefile).resolve() if args.snakefile else THIS_DIR / "Snakefile"
    configfile = Path(args.configfile).resolve()
    workdir = str(Path(args.directory).resolve()) if args.directory else None

    validate_args(args, snakefile, configfile)

    config = load_config(configfile)
    resources = parse_resources(args.resources)
    singularity_args = build_singularity_args(args)

    if args.show_workflows:
        print_workflows()
        return 0

    # Keep a reference for the life of this process -- see
    # acquire_project_lock()'s docstring-comment for why. Bypasses Snakemake's
    # own directory-wide lock (lock=False below) in favor of this per-config
    # one.
    _project_lock = acquire_project_lock(configfile)

    cluster_cmd = None
    default_resources = None
    nodes = None
    cores = args.cores
    if args.cluster:
        cluster_logdir = THIS_DIR / "logs" / "pbs"
        cluster_logdir.mkdir(parents=True, exist_ok=True)
        cluster_cmd = build_cluster_cmd(args.queue, cluster_logdir)
        from snakemake.resources import DefaultResources

        default_resources = DefaultResources(
            [
                f"mem_mb={CLUSTER_DEFAULT_MEM_MB}",
                f"runtime={CLUSTER_DEFAULT_RUNTIME_MIN}",
            ]
        )
        nodes = args.jobs
        # Snakemake caps every rule's `threads:` to min(rule.threads, cores)
        # before rendering `{threads}` into the --cluster submit template --
        # regardless of execution mode. Use a sentinel well above any rule's
        # declared threads so nothing gets capped; -q's value is irrelevant
        # once --cluster is set.
        cores = 999

    print_run_summary(
        args=args,
        snakefile=snakefile,
        configfile=configfile,
        config=config,
        singularity_args=singularity_args,
        resources=resources,
    )

    snakemake_kwargs: Dict[str, Any] = {
        "snakefile": str(snakefile),
        "configfiles": [str(configfile)],
        "targets": args.workflow or [],
        "workdir": workdir,
        "cores": cores,
        "dryrun": args.dry_run,
        "use_conda": not args.no_conda,
        "use_singularity": not args.no_singularity,
        "singularity_args": singularity_args,
        "forceall": args.forceall,
        "force_incomplete": args.rerun_incomplete,
        "unlock": args.unlock,
        "lock": False,
        "printdag": args.dag,
        "lint": args.lint,
        "printshellcmds": args.printshellcmds,
        "keepgoing": args.keep_going,
        "latency_wait": args.latency_wait,
        "restart_times": args.restart_times,
    }

    if cluster_cmd:
        snakemake_kwargs["cluster"] = cluster_cmd
    if nodes:
        snakemake_kwargs["nodes"] = nodes
    if default_resources:
        snakemake_kwargs["default_resources"] = default_resources

    if resources:
        snakemake_kwargs["resources"] = resources

    if args.rerun_triggers:
        snakemake_kwargs["rerun_triggers"] = args.rerun_triggers

    if args.wrapper_prefix:
        snakemake_kwargs["wrapper_prefix"] = args.wrapper_prefix

    if args.conda_prefix:
        snakemake_kwargs["conda_prefix"] = args.conda_prefix

    status = snakemake.snakemake(**snakemake_kwargs)

    if status:
        success("Workflow finished successfully.")
        return 0

    die("Workflow failed.", exit_code=1)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=color("Launcher for Varan's Snakemake workflow.", Fore.GREEN),
        epilog=f"""
Examples:

  {Fore.CYAN}./run.py -w create -c config.yaml -q 4 -n{Style.RESET_ALL}
      Dry-run a new study creation (Snakemake's own dry-run).

  {Fore.CYAN}./run.py -w create -c config.yaml -q 4{Style.RESET_ALL}
      Actually build a new study (CNV/SNV/fusion/clinical run in parallel).

  {Fore.CYAN}./run.py -w update -c config.yaml -q 1{Style.RESET_ALL}
      Merge a new batch of samples into an existing study.

  {Fore.CYAN}./run.py --list-workflows{Style.RESET_ALL}
      Show available workflow targets.

  {Fore.CYAN}./run.py -w create -c config.yaml -q 4 --dag | dot -Tsvg > create_dag.svg{Style.RESET_ALL}
      Export the create DAG as SVG.

  {Fore.CYAN}./run.py -w create -c config.yaml -q 1 -cl -qu workq -j 10{Style.RESET_ALL}
      Submit each job as its own PBS job via qsub instead of running
      everything locally in one allocation (this cluster: OpenPBS 20.0.1).
""",
    )

    required = parser.add_argument_group("Required for normal execution")
    required.add_argument(
        "-w",
        "--workflow",
        nargs="+",
        help=f"Workflow target(s) to run. Main targets: {', '.join(WORKFLOWS.keys())}.",
    )
    required.add_argument(
        "-c",
        "--configfile",
        help="Path to the YAML Snakemake config file (see config.yaml).",
    )
    required.add_argument(
        "-q",
        "--cores",
        type=int,
        default=1,
        help="Number of CPU cores available to Snakemake.",
    )

    execution = parser.add_argument_group("Execution mode")
    execution.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Perform a dry-run without executing jobs (Snakemake's own "
        "dry-run - not varan.py's own -D/--dry-run preflight check).",
    )
    execution.add_argument(
        "-f",
        "--forceall",
        action="store_true",
        help="Force execution of all jobs even if outputs already exist.",
    )
    execution.add_argument(
        "-ri",
        "--rerun-incomplete",
        action="store_true",
        help="Rerun jobs with incomplete outputs.",
    )
    execution.add_argument(
        "-u",
        "--unlock",
        action="store_true",
        help="Unlock the working directory after an interrupted run.",
    )
    execution.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue independent jobs after one job fails.",
    )
    execution.add_argument(
        "--latency-wait",
        type=int,
        default=60,
        help="Seconds to wait for output files to appear on slow filesystems.",
    )

    cluster = parser.add_argument_group("Cluster (PBS)")
    cluster.add_argument(
        "-cl",
        "--cluster",
        action="store_true",
        help="Submit each rule instance as its own PBS job via qsub, instead of "
        "running everything locally in one allocation (default off; plain -q N "
        "still runs everything locally with N cores, unchanged).",
    )
    cluster.add_argument(
        "-qu",
        "--queue",
        type=str,
        default="workq",
        help="PBS queue to submit to when --cluster is set.",
    )
    cluster.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=20,
        help="Max number of concurrent PBS jobs when --cluster is set.",
    )
    cluster.add_argument(
        "-rt",
        "--restart-times",
        type=int,
        default=2,
        help="Number of times to automatically retry a job that fails (e.g. "
        "killed for exceeding its resources), each retry with progressively "
        "more mem_mb/runtime for rules whose resources: use a scaled() "
        "helper, once any declare one (i.e. up to N+1 total attempts). "
        "Rules with no declared resources fall back to "
        f"mem_mb={CLUSTER_DEFAULT_MEM_MB}/runtime={CLUSTER_DEFAULT_RUNTIME_MIN} "
        "regardless of attempt number.",
    )

    reporting = parser.add_argument_group("Reporting and debugging")
    reporting.add_argument(
        "--list-workflows",
        "--show-workflows",
        dest="show_workflows",
        action="store_true",
        help="Print available high-level workflow targets and exit.",
    )
    reporting.add_argument(
        "-d",
        "--dag",
        action="store_true",
        help='Print DAG in DOT format. Example: ./run.py -w create ... --dag | dot -Tsvg > dag.svg',
    )
    reporting.add_argument(
        "--printshellcmds",
        "-p",
        action="store_true",
        help="Print shell commands before executing them.",
    )

    reporting.add_argument(
        "-l",
        "--lint",
        default=None,
        help='Run Snakemake lint. Common values: "text" or "json".',
    )

    environment = parser.add_argument_group("Environment and paths")
    environment.add_argument(
        "-di",
        "--directory",
        type=str,
        default=None,
        help="Working directory for Snakemake execution.",
    )
    environment.add_argument(
        "--snakefile",
        type=str,
        default=None,
        help="Custom Snakefile path. Default: Snakefile next to this run.py.",
    )
    environment.add_argument(
        "--no-conda",
        action="store_true",
        help="Disable Snakemake conda integration (every rule currently "
        "shares one conda_env from config.yaml - see its own comment).",
    )
    environment.add_argument(
        "--conda-prefix",
        type=str,
        default=None,
        help="Optional Snakemake conda prefix directory.",
    )
    environment.add_argument(
        "--no-singularity",
        action="store_true",
        help="Disable Snakemake Singularity/Apptainer integration (currently "
        "inert either way - no rule declares a container: yet).",
    )
    environment.add_argument(
        "--bind",
        action="append",
        default=[],
        help="Singularity bind path. Can be used multiple times. No default "
        "assumed - pass whatever this cluster's real data mount is.",
    )
    environment.add_argument(
        "--singularity-args",
        type=str,
        default="",
        help="Additional raw Singularity arguments.",
    )
    environment.add_argument(
        "-wr",
        "--wrapper-prefix",
        type=str,
        default=None,
        help="Prefix for a local Snakemake wrapper repository.",
    )

    advanced = parser.add_argument_group("Advanced")
    advanced.add_argument(
        "--resources",
        nargs="+",
        default=None,
        help="Custom Snakemake resources as key=value pairs, e.g. --resources igv=1 mem_mb=64000.",
    )
    advanced.add_argument(
        "--allow-custom-target",
        action="store_true",
        help="Allow running an internal stage name directly (walk_setup/"
        "walk_cnv/walk_snv/walk_fusion/walk_clinical), not just the four "
        "main targets - for debugging one stage in isolation. These all "
        "require walk_setup's context pickle to already exist.",
    )
    advanced.add_argument(
        "--rerun-triggers",
        nargs="+",
        default=None,
        choices=["mtime", "params", "input", "software-env", "code"],
        help="Override Snakemake's rerun-trigger set (default when omitted: mtime "
        "params input software-env code -- Snakemake's own default, reruns a job "
        "if ANY of those changed since last time, even if its output already "
        "exists). Pass --rerun-triggers mtime to fall back to classic make-style "
        "mtime-only comparison -- useful when adding new samples to an existing "
        "project config and an unrelated historical code/params edit elsewhere in "
        "the pipeline would otherwise force recomputing already-finished samples. "
        "Only use mtime-only when you're sure existing outputs are still valid "
        "under the current rule code: it will NOT detect a genuine bug fix to a "
        "rule's shell command/script as a reason to rerun, so a real fix could "
        "silently leave old, wrong outputs in place. Left unset by default so no "
        "project's behavior changes unless this is passed explicitly.",
    )

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.show_workflows:
        print_workflows()
        return 0

    if not args.configfile:
        die("Missing config file. Use -c config.yaml")

    return run_snakemake(args)


if __name__ == "__main__":
    raise SystemExit(main())
