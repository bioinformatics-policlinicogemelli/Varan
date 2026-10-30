# Running SigMA on an HPC cluster via conda (no Docker needed)

This is the cluster-friendly alternative to the Dockerfile's apt+R
recipe, for users who run Varan through a conda/mamba environment on an
HPC cluster rather than through the Docker image. Same end result (R +
SigMA installed, `run_sigma.R` runnable), different install path.

**Status: tested this round.** An `environment.yml` with this exact
package list was actually solved and installed with `micromamba` (conda-
forge/bioconda's real package resolver and repodata - the same channels
a cluster's own `conda`/`mamba` would use) in a clean Ubuntu 22.04
container, then `devtools::install_github("parklab/SigMA")` and
`library(SigMA)` were run for real inside that environment - not just a
documentation translation of the Dockerfile. See "What was actually
verified" at the end for the exact scope of that test.

## 1. Environment file

Save as `environment.yml` (e.g. alongside `conf.ini`):

```yaml
name: varan-sigma
channels:
  - conda-forge
  - bioconda
dependencies:
  - r-base
  - bioconductor-bsgenome
  - bioconductor-bsgenome.hsapiens.ucsc.hg19
  - bioconductor-variantannotation
  - bioconductor-genomicranges
  - bioconductor-iranges
  - r-gbm
  - r-nnls
  - r-reshape2
  - r-rmisc
  - r-dt
  - r-gridextra
  - r-ggplot2
  - r-devtools
  - r-remotes
```

Every one of these package names was checked against the real
conda-forge/bioconda package index before writing this list (not
assumed) - all of them exist as prebuilt binaries, including
`bioconductor-bsgenome.hsapiens.ucsc.hg19` and `r-rmisc`, which were
flagged in `SIGMA_INTEGRATION_FEASIBILITY.md` as "may or may not have a
conda package, check for real." Both do.

**`r-remotes` is listed explicitly, not left implicit.** Testing this
round found that conda-forge's `r-devtools` package does **not** pull in
`r-remotes` as a hard dependency, even though `devtools::install_github()`
calls straight into `remotes` internally and fails outright without it
("The package \"remotes\" is required.") - confirmed by an actual failed
run, then fixed by adding `r-remotes` to the environment. `devtools` is
still kept in the list too (not swapped out) because SigMA's own
`DESCRIPTION` lists `devtools` as a hard runtime `Imports`, not just an
install-time convenience - `library(SigMA)` itself needs it present.

**None of SigMA's declared dependencies need a CRAN/GitHub fallback
inside the conda env** - every package in the Dockerfile's
`BiocManager::install(...)` list has a bioconda/conda-forge equivalent.
The only thing that has to come from GitHub either way is **SigMA
itself** - it isn't on CRAN, Bioconductor, or bioconda, so
`devtools::install_github()` is unavoidable regardless of which install
path (Docker or conda) is used.

## 2. Create the environment

```bash
# with conda:
conda env create -f environment.yml
# or, much faster (recommended if available - this is what was
# actually used to test this round, resolved the full 287-package
# environment in ~20 seconds vs. conda's usual solver time):
mamba env create -f environment.yml
```

Then install SigMA itself into that environment (not a conda/bioconda
package - see above):

```bash
conda run -n varan-sigma Rscript -e \
  'devtools::install_github("parklab/SigMA", dependencies = FALSE, upgrade = "never")'
```

`dependencies = FALSE` is deliberate: every one of SigMA's declared
`Imports` (BSgenome, VariantAnnotation, GenomicRanges, IRanges, gbm,
nnls, reshape2, Rmisc, DT, gridExtra, ggplot2, devtools) is already
installed by the `environment.yml` step above as a conda binary -
letting `install_github` try to also resolve them from source via CRAN
would be redundant at best, and risks pulling in a mismatched/newer
source build alongside the conda-provided binary at worst.

Verify it loads:

```bash
conda run -n varan-sigma Rscript -e 'library(SigMA)'
```

## 3. Point Varan at this environment's `Rscript`

Varan's `run_sigma.R` subprocess is invoked via whatever `[Paths] RSCRIPT`
in `conf.ini` resolves to (default: plain `Rscript`, resolved via the
calling process's `$PATH` - see `sigma_runner.py`). On a cluster this
conda environment's `Rscript` is very unlikely to be the one on the
default `$PATH` (that's typically the system R, if any exists at all),
so **this is exactly what `[Paths] RSCRIPT` exists to override**:

```ini
[Paths]
RSCRIPT = /home/<user>/miniforge3/envs/varan-sigma/bin/Rscript
```

Find the exact path with:

```bash
conda run -n varan-sigma which Rscript
# or, without activating anything:
echo "$(conda info --base)/envs/varan-sigma/bin/Rscript"
```

This is a plain absolute path, not an activation command - `sigma_runner.py`
calls it directly via `subprocess.run([rscript_bin, ...])` (never
`shell=True`), so no `conda activate` step is needed or possible in that
call; pointing `RSCRIPT` straight at the env's own `bin/Rscript` is
sufficient and is the standard way to invoke a specific conda
environment's binary without activating it in the parent shell.

## 4. Whole story? Yes, with one caveat

For a single user, on a single cluster node, running Varan directly
(not through the `Snakefile`): **yes, this is the whole story** - create
the conda env once, install SigMA into it once, set `[Paths] RSCRIPT` in
`conf.ini` once, then run `python varan.py -g ...` as normal. Nothing
else in the pipeline changes; `sigma_runner.py`/`run_sigma.R` don't know
or care whether R came from Docker, conda, or a system package - they
just exec whatever `RSCRIPT` points to.

The caveat, flagged here for completeness and explicitly **out of scope
for this round** (per direction received): once there are *multiple*
per-stage conda environments (e.g. a separate env for VEP/vcf2maf vs.
this SigMA env vs. whatever else), manually tracking and setting the
right `RSCRIPT`/tool paths per stage stops scaling, and the `Snakefile`'s
own per-rule `conda:` directive (already present as a placeholder -
see its module docstring: "same environment for all four rules for
now... replace with the four per-stage environments once they exist")
is the right mechanism to manage that automatically. That's a
Snakemake-wiring change, not attempted here.

## What was actually verified this round

- The exact `environment.yml` above was solved and installed for real
  with `micromamba` (conda-forge's own minimal, spec-compatible
  installer/resolver - not a stand-in or simulation) against live
  conda-forge/bioconda repodata, in a clean `ubuntu:22.04` container:
  287 packages resolved, ~540MB of conda packages downloaded, all
  prebuilt binaries (no source compilation, no missing system headers -
  contrast with the Dockerfile's apt+CRAN-source path, where this same
  round's testing found and fixed missing `libuv1-dev`/`libharfbuzz-dev`/
  `libfribidi-dev` system headers that were silently breaking
  `devtools`'s own dependency chain). Separately, `bioconductor-
  bsgenome.hsapiens.ucsc.hg19`'s post-link script then downloads the
  actual ~677MB hg19 sequence data straight from Bioconductor's own CDN
  (this is normal bioconda behavior for large annotation-data packages -
  the conda package itself is just the R wrapper code; expect this
  same ~677MB download either way, conda or Docker, since the Dockerfile
  path fetches the identical tarball from the identical source).
- `devtools::install_github("parklab/SigMA", dependencies = FALSE)` run
  for real inside that environment (after the `r-remotes` fix above -
  the first attempt failed exactly as described, confirming this was a
  real gap and not a hypothetical one).
- `library(SigMA)` loads successfully inside that environment.
- **`run_sigma.R` itself run end-to-end against this conda-provided R**,
  not just `library(SigMA)`: `micromamba run -n varan-sigma Rscript
  run_sigma.R --maf ... --tumor-type breast --do-mva TRUE ...` against
  the same synthetic filtered MAF used in the Docker-based smoke test,
  producing a full result row (`Signature_3_mva = 0.871`, `pass_mva =
  TRUE`, `pass_mva_strict = TRUE`, `categ = "Signature_3_hc"`) -
  byte-for-byte the same numeric result as the Docker-based run of the
  identical input, as expected since it's the same SigMA version and
  the same input spectrum either way.

**Not verified this round**: an actual `conda`/`mamba` (as opposed to
`micromamba`) run on a real HPC login node - `micromamba` uses the same
solver/package format and the same conda-forge/bioconda channels, so
this should be representative, but a real cluster may have channel-
priority defaults, proxy/firewall restrictions, or a pre-existing
`.condarc` that changes resolution behavior in ways a clean container
can't reproduce.
