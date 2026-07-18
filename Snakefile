"""Lightweight Snakemake wrapper around the varan.py CLI.

Wraps each of Varan's four top-level operations (create, update, extract,
remove) as its own rule, by shelling out to the existing varan.py CLI rather
than reimplementing or restructuring its internals. This keeps the actual
pipeline logic exactly as it is today, and layers Snakemake's dependency
tracking, per-run logging, and conda-environment management on top of it.

Scope note: this is a starting point, not a full internal DAG decomposition
of walk -> filter -> concatenate -> table -> validate into separate,
independently re-runnable Snakemake rules. varan.py does not currently
expose those stages as separate CLI entry points (they are only reachable
from inside the varan() function), so splitting them out is a larger,
separate change - not attempted here to avoid restructuring pipeline code
without being able to test it end to end.

All four rules currently point at the same shared conda environment
(see conda_env in config.yaml) as an explicit placeholder, per instruction -
swap in the four per-stage environments once they exist; nothing else in
this file needs to change to pick them up.

NOT YET RUN AGAINST A REAL SNAKEMAKE INSTALLATION: no Snakemake install or
real cluster conda environment was available while writing this - review
and test on the cluster before relying on it.

Usage (from the repository root, with conf.ini and config.yaml filled in):
    snakemake --cores 1 --use-conda create
    snakemake --cores 1 --use-conda update
    snakemake --cores 1 --use-conda extract
    snakemake --cores 1 --use-conda remove
"""

configfile: "config.yaml"


rule create:
    """Build a new study folder from raw input (varan.py -i)."""
    output:
        report=f"{config['create']['output_folder']}/report_VARAN.html",
    log:
        "Logs/snakemake_create.log",
    conda:
        config["conda_env"]
    params:
        args=config["create"]["args"],
        conf=config["conf_path"],
    shell:
        "python varan.py {params.args} -C {params.conf} > {log} 2>&1"


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
