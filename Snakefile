"""Snakemake wrapper around the Varan CLI.

`create` is now a real per-stage DAG end to end, not a single opaque
shell-out:

    walk_setup -> {walk_cnv, walk_snv, walk_fusion, walk_clinical} (parallel)
               -> filter -> concatenate -> tables -> validate -> create

walk_setup/walk_cnv/walk_snv/walk_fusion/walk_clinical shell out to
walk_stage.py, a thin CLI wrapper around walk.py's internal
_walk_setup/_walk_process_cnv/_walk_process_snv/_walk_process_fusion/
_walk_write_clinical_tables functions (see walk.py's WalkContext
docstring for why those four stages are safe to run independently: none
of them depends on another's output, only on setup's).

filter/concatenate/tables/validate shell out to pipeline_stage.py, a thin
CLI wrapper around varan.py's own filter_main/concatenate_main/
meta_case_main/validate_output (varan.py's numbered stages 2-5) - unlike
the four WALK stages, these four are NOT mutually independent: each needs
the previous one's output on disk (filter needs the maf/ folder,
concatenate needs filter's output, tables needs the concatenated data,
validate needs everything). Splitting them into separate rules here isn't
about running them in parallel with *each other* - it's:
  1. Real per-stage resume-from-failure: a validation failure no longer
     forces re-running filter/concatenate/tables from scratch.
  2. `filter` only depends on walk_snv's output, not on walk_cnv/
     walk_fusion/walk_clinical too - so with enough cores, filter can
     start running while CNV/fusion/clinical are still in progress,
     instead of waiting for all four WALK stages to finish. This is the
     one real remaining parallelism gap from the WALK-only split.

`create` is now a thin final rule that only waits on `validate`'s output -
no more re-invoking `python varan.py -R` as a catch-all for stages 2-5.

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
    if _c.get("vcf_type"):
        args += f" --vcf-type {_c['vcf_type']}"
    if _c.get("filters"):
        args += f" --filters {_c['filters']}"
    return args


def _filter_args() -> str:
    args = (
        f"--input {' '.join(_c['input'])} --output {_c['output_folder']} "
        f"-c {_c['cancer']} --resume"
    )
    if _c.get("oncokb"):
        args += " --oncokb"
    if _c.get("filters"):
        args += f" --filters {_c['filters']}"
    return args


def _tables_args() -> str:
    return f"--output {_c['output_folder']} -c {_c['cancer']}"


def _validate_args() -> str:
    args = f"--input {' '.join(_c['input'])} --output {_c['output_folder']} -c {_c['cancer']}"
    if _c.get("multiple"):
        args += " --multiple"
    if _c.get("oncokb"):
        args += " --oncokb"
    if _c.get("filters"):
        args += f" --filters {_c['filters']}"
    return args


rule walk_setup:
    """Resolve input/output paths once (create_folder, input extraction,
    resume-state checks) - see _walk_setup() in walk.py."""
    output:
        ctx=_ctx_pkl,
    log:
        "Logs/snakemake_walk_setup.log",
    conda:
        config["conda_env"]
    params:
        args=_create_setup_args(),
        conf=config["conf_path"],
    shell:
        "python walk_stage.py setup {params.args} > {log} 2>&1"


rule walk_cnv:
    """CNV calls + BRCA exon-level CNV table - independent of walk_snv/
    walk_fusion, see WalkContext's docstring in walk.py."""
    input:
        ctx=_ctx_pkl,
    output:
        done=f"{_c['output_folder']}/.cnv.done",
    log:
        "Logs/snakemake_walk_cnv.log",
    conda:
        config["conda_env"]
    shell:
        "python walk_stage.py cnv --ctx {input.ctx} > {log} 2>&1"


rule walk_snv:
    """SNV calls -> vcf2maf - independent of walk_cnv/walk_fusion."""
    input:
        ctx=_ctx_pkl,
    output:
        done=f"{_c['output_folder']}/.snv.done",
    log:
        "Logs/snakemake_walk_snv.log",
    conda:
        config["conda_env"]
    shell:
        "python walk_stage.py snv --ctx {input.ctx} > {log} 2>&1"


rule walk_fusion:
    """RNA fusions + splice variants -> data_sv.txt - independent of
    walk_cnv/walk_snv."""
    input:
        ctx=_ctx_pkl,
    output:
        done=f"{_c['output_folder']}/.fusion.done",
    log:
        "Logs/snakemake_walk_fusion.log",
    conda:
        config["conda_env"]
    shell:
        "python walk_stage.py fusion --ctx {input.ctx} > {log} 2>&1"


rule walk_clinical:
    """data_clinical_patient.txt / data_clinical_sample.txt (MSI/TMB, exon,
    HRD) - independent of walk_cnv/walk_snv/walk_fusion."""
    input:
        ctx=_ctx_pkl,
    output:
        done=f"{_c['output_folder']}/.clinical.done",
    log:
        "Logs/snakemake_walk_clinical.log",
    conda:
        config["conda_env"]
    shell:
        "python walk_stage.py clinical --ctx {input.ctx} > {log} 2>&1"


rule filter:
    """Stage 2: MAF filtering (varan.py's filter_main). Depends only on
    walk_snv, not on walk_cnv/walk_fusion/walk_clinical - the one real
    parallelism gap this branch closes, see module docstring."""
    input:
        snv=f"{_c['output_folder']}/.snv.done",
    output:
        done=f"{_c['output_folder']}/.filter.done",
    log:
        "Logs/snakemake_filter.log",
    conda:
        config["conda_env"]
    params:
        args=_filter_args(),
        conf=config["conf_path"],
    shell:
        "python pipeline_stage.py -C {params.conf} filter {params.args} > {log} 2>&1"


rule concatenate:
    """Stage 3: concatenate per-sample MAFs (varan.py's concatenate_main)."""
    input:
        filter=f"{_c['output_folder']}/.filter.done",
    output:
        done=f"{_c['output_folder']}/.concatenate.done",
    log:
        "Logs/snakemake_concatenate.log",
    conda:
        config["conda_env"]
    params:
        oncokb="--oncokb" if _c.get("oncokb") else "",
        filters=f"--filters {_c['filters']}" if _c.get("filters") else "",
        output=_c["output_folder"],
        conf=config["conf_path"],
    shell:
        "python pipeline_stage.py -C {params.conf} concatenate "
        "--output {params.output} {params.oncokb} {params.filters} > {log} 2>&1"


rule tables:
    """Stage 4: meta/case list files (varan.py's meta_case_main). Depends
    only on concatenate - CNV/fusion/clinical are already on disk from the
    WALK stages, this stage doesn't re-touch them."""
    input:
        concat=f"{_c['output_folder']}/.concatenate.done",
        cnv=f"{_c['output_folder']}/.cnv.done",
        fusion=f"{_c['output_folder']}/.fusion.done",
        clinical=f"{_c['output_folder']}/.clinical.done",
    output:
        done=f"{_c['output_folder']}/.tables.done",
    log:
        "Logs/snakemake_tables.log",
    conda:
        config["conda_env"]
    params:
        args=_tables_args(),
        conf=config["conf_path"],
    shell:
        "python pipeline_stage.py -C {params.conf} tables {params.args} > {log} 2>&1"


rule validate:
    """Stage 5: cBioPortal validation + report (varan.py's validate_output)."""
    input:
        tables=f"{_c['output_folder']}/.tables.done",
    output:
        done=f"{_c['output_folder']}/.validate.done",
        report=f"{_c['output_folder']}/report_VARAN.html",
    log:
        "Logs/snakemake_validate.log",
    conda:
        config["conda_env"]
    params:
        args=_validate_args(),
        conf=config["conf_path"],
    shell:
        "python pipeline_stage.py -C {params.conf} validate {params.args} > {log} 2>&1"


rule create:
    """Thin final target for the whole create flow - everything real
    happens in walk_setup/walk_cnv/walk_snv/walk_fusion/walk_clinical/
    filter/concatenate/tables/validate above."""
    input:
        report=f"{_c['output_folder']}/report_VARAN.html",
    output:
        touch(f"{_c['output_folder']}/.create.done"),


rule update:
    """Merge a new batch of samples into an existing study (varan.py -u)."""
    output:
        report=f"{config['update']['output_folder']}/report_VARAN.html",
    log:
        "Logs/snakemake_update.log",
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
        "Logs/snakemake_extract.log",
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
        "Logs/snakemake_remove.log",
    conda:
        config["conda_env"]
    params:
        args=config["remove"]["args"],
        conf=config["conf_path"],
    shell:
        "python varan.py {params.args} -C {params.conf} > {log} 2>&1"
