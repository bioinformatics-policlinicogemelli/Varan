"""Snakemake wrapper around the Varan CLI.

`create` is now a real per-stage DAG, not a single opaque shell-out:

    walk_setup -> {walk_cnv, walk_snv, walk_fusion, walk_clinical} (parallel)
               -> create (rest of the pipeline: filter/concatenate/tables/
                  validation, unchanged - see varan.py's numbered stages
                  2-5)

walk_setup/walk_cnv/walk_snv/walk_fusion/walk_clinical shell out to
walk_stage.py, a thin CLI wrapper around walk.py's internal
_walk_setup/_walk_process_cnv/_walk_process_snv/_walk_process_fusion/
_walk_write_clinical_tables functions (see walk.py's WalkContext
docstring for why those four stages are safe to run independently: none
of them depends on another's output, only on setup's). This means
`snakemake --cores 4 create` genuinely runs CNV, SNV, fusion and the
clinical tables in parallel instead of sequentially.

Known trade-off (documented rather than hidden): the final `create` rule
re-invokes `python varan.py -i ... -R ...` to run stages 2-5 (filter,
concatenate, tables, validation), which are not yet split into their own
rules. That second call re-enters walk_folder() too, but with resume=True
and the output folder's "temp" marker already present, so it skips
`create_folder` (no wipe) and skips re-running vcf2maf (the slowest step,
already done by walk_snv) - it does however harmlessly re-run the CNV/
fusion/clinical stages a second time (idempotent overwrite, just wasted
work). Splitting stages 2-5 out the same way is a natural next step, not
attempted tonight.

update/extract/remove do not go through walk_folder at all (they use their
own update_*/extract_*/delete_* functions, already one function per file
type - see Update_functions.py/ExtractSamples_functions.py/
Delete_functions.py) - they are left as single shell-out rules here; giving
them the same per-file-type parallel treatment as create is a separate,
smaller future step.

All rules currently point at the same shared conda environment
(see conda_env in config.yaml) as an explicit placeholder - swap in
per-stage environments once they exist; nothing else in this file needs to
change to pick them up.

Usage (from the repository root, with conf.ini and config.yaml filled in):
    snakemake --cores 4 --use-conda create
    snakemake --cores 1 --use-conda update
    snakemake --cores 1 --use-conda extract
    snakemake --cores 1 --use-conda remove

Or via run.py, which adds per-config locking, PBS cluster submission and a
run summary on top of the same four targets - see run.py --help:
    ./run.py -w create -c config.yaml -q 4
"""

configfile: "config.yaml"

_c = config["create"]
_ctx_pkl = f"{_c['output_folder']}/.walk_ctx.pkl"


def _create_setup_args() -> str:
    args = (
        f"--input {' '.join(_c['input'])} --output {_c['output_folder']} "
        f"-c {_c['cancer']} --ctx-out {_ctx_pkl}"
    )
    if _c.get("overwrite"):
        args += " --overwrite"
    if _c.get("resume"):
        args += " --resume"
    if _c.get("multiple"):
        args += " --multiple"
    if _c.get("oncokb"):
        args += " --oncokb"
    if _c.get("sigma"):
        args += " --sigma"
    if _c.get("vcf_type"):
        args += f" --vcf-type {_c['vcf_type']}"
    if _c.get("filters"):
        args += f" --filters {_c['filters']}"
    return args


def _create_rest_args() -> str:
    # Deliberately never passes -w here (the folder was already created, or
    # not, by walk_setup - re-wiping it now would destroy walk_cnv/walk_snv/
    # walk_fusion/walk_clinical's output) and always passes -R (the folder's
    # "temp" marker already exists by this point, so -R makes _walk_setup
    # skip create_folder entirely and skip re-running vcf2maf - see module
    # docstring above).
    args = f"-i {' '.join(_c['input'])} -o {_c['output_folder']} -c {_c['cancer']} -R"
    if _c.get("multiple"):
        args += " -m"
    if _c.get("oncokb"):
        args += " -k"
    if _c.get("sigma"):
        args += " -g"
    if _c.get("vcf_type"):
        args += f" -t {_c['vcf_type']}"
    if _c.get("filters"):
        args += f" -f {_c['filters']}"
    return args


rule walk_setup:
    """Resolve input/output paths once (create_folder, input extraction,
    resume-state checks) - see _walk_setup() in walk.py."""
    output:
        ctx=_ctx_pkl,
    log:
        f"{_c['output_folder']}/service/logs/snakemake_walk_setup.log",
    conda:
        config["conda_env"]
    params:
        args=_create_setup_args(),
        conf=config["conf_path"],
    shell:
        "python walk_stage.py -C {params.conf} setup {params.args} > {log} 2>&1"


rule walk_cnv:
    """CNV calls + BRCA exon-level CNV table - independent of walk_snv/
    walk_fusion, see WalkContext's docstring in walk.py."""
    input:
        ctx=_ctx_pkl,
    output:
        done=f"{_c['output_folder']}/.cnv.done",
    log:
        f"{_c['output_folder']}/service/logs/snakemake_walk_cnv.log",
    conda:
        config["conda_env"]
    params:
        conf=config["conf_path"],
    shell:
        "python walk_stage.py -C {params.conf} cnv --ctx {input.ctx} > {log} 2>&1"


rule walk_snv:
    """SNV calls -> vcf2maf - independent of walk_cnv/walk_fusion."""
    input:
        ctx=_ctx_pkl,
    output:
        done=f"{_c['output_folder']}/.snv.done",
    log:
        f"{_c['output_folder']}/service/logs/snakemake_walk_snv.log",
    conda:
        config["conda_env"]
    params:
        conf=config["conf_path"],
    shell:
        "python walk_stage.py -C {params.conf} snv --ctx {input.ctx} > {log} 2>&1"


rule walk_fusion:
    """RNA fusions + splice variants -> data_sv.txt - independent of
    walk_cnv/walk_snv."""
    input:
        ctx=_ctx_pkl,
    output:
        done=f"{_c['output_folder']}/.fusion.done",
    log:
        f"{_c['output_folder']}/service/logs/snakemake_walk_fusion.log",
    conda:
        config["conda_env"]
    params:
        conf=config["conf_path"],
    shell:
        "python walk_stage.py -C {params.conf} fusion --ctx {input.ctx} > {log} 2>&1"


rule walk_clinical:
    """data_clinical_patient.txt / data_clinical_sample.txt (MSI/TMB, exon,
    HRD) - independent of walk_cnv/walk_snv/walk_fusion."""
    input:
        ctx=_ctx_pkl,
    output:
        done=f"{_c['output_folder']}/.clinical.done",
    log:
        f"{_c['output_folder']}/service/logs/snakemake_walk_clinical.log",
    conda:
        config["conda_env"]
    params:
        conf=config["conf_path"],
    shell:
        "python walk_stage.py -C {params.conf} clinical --ctx {input.ctx} > {log} 2>&1"


rule create:
    """Build a new study folder from raw input (varan.py -i), fanning the
    walk stage out into walk_setup + {walk_cnv, walk_snv, walk_fusion,
    walk_clinical} above, then running stages 2-5 (filter/concatenate/
    tables/validation - see module docstring for the resume=True trade-off)."""
    input:
        cnv=f"{_c['output_folder']}/.cnv.done",
        snv=f"{_c['output_folder']}/.snv.done",
        fusion=f"{_c['output_folder']}/.fusion.done",
        clinical=f"{_c['output_folder']}/.clinical.done",
    output:
        report=f"{_c['output_folder']}/report_VARAN.html",
    log:
        f"{_c['output_folder']}/service/logs/snakemake_create.log",
    conda:
        config["conda_env"]
    params:
        args=_create_rest_args(),
        conf=config["conf_path"],
    shell:
        "python varan.py {params.args} -C {params.conf} > {log} 2>&1"


rule update:
    """Merge a new batch of samples into an existing study (varan.py -u)."""
    output:
        report=f"{config['update']['output_folder']}/report_VARAN.html",
    log:
        f"{config['update']['output_folder']}/service/logs/snakemake_update.log",
    conda:
        config["conda_env"]
    params:
        args=config["update"]["args"],
        conf=config["conf_path"],
    shell:
        "python varan.py {params.args} -C {params.conf} > {log} 2>&1"


rule extract:
    """Extract a subset of samples into a new study (varan.py -e)."""
    output:
        report=f"{config['extract']['output_folder']}/report_VARAN.html",
    log:
        f"{config['extract']['output_folder']}/service/logs/snakemake_extract.log",
    conda:
        config["conda_env"]
    params:
        args=config["extract"]["args"],
        conf=config["conf_path"],
    shell:
        "python varan.py {params.args} -C {params.conf} > {log} 2>&1"


rule remove:
    """Remove a subset of samples, producing a new study version (varan.py -r)."""
    output:
        report=f"{config['remove']['output_folder']}/report_VARAN.html",
    log:
        f"{config['remove']['output_folder']}/service/logs/snakemake_remove.log",
    conda:
        config["conda_env"]
    params:
        args=config["remove"]["args"],
        conf=config["conf_path"],
    shell:
        "python varan.py {params.args} -C {params.conf} > {log} 2>&1"


## ============================================================================
## Completion notification -- opt-in via `notify_email:` in config.yaml.
## Omit it and this silently no-ops (harmless for setups that don't want
## it). Uses the system mail/mailx command with a 15s timeout so a broken
## mail relay can't hang the pipeline's completion. Fires once per
## run.py/snakemake invocation, not once per rule/target.
## ============================================================================

def _notify_email(subject: str, body: str) -> None:
    email = config.get("notify_email")
    if not email:
        return
    import shutil
    import subprocess
    mail_bin = shutil.which("mail") or shutil.which("mailx")
    if not mail_bin:
        return
    try:
        subprocess.run(
            [mail_bin, "-s", subject, email],
            input=body.encode(), check=False, timeout=15)
    except Exception:
        pass


onsuccess:
    _notify_email(
        "[Varan/Snakemake] SUCCESS",
        "The Varan Snakemake run completed successfully.\n"
        f"Config file: {config.get('conf_path', 'conf.ini')}\n")


onerror:
    _notify_email(
        "[Varan/Snakemake] FAILED",
        "The Varan Snakemake run failed. Each target's own study folder "
        "has a report_VARAN.html explaining what happened (written even on "
        "failure - see write_report_failure()), and the per-rule Snakemake "
        "logs live under that same folder's service/logs/.\n"
        f"Config file: {config.get('conf_path', 'conf.ini')}\n")
