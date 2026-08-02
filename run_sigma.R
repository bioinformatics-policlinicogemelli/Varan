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

# run_sigma.R - thin CLI wrapper around the SigMA R package (parklab/SigMA),
# invoked once per sample from walk.py's _walk_process_snv() via
# subprocess.run() with an argument list (never shell=True - same fix this
# codebase already applies to vcf2maf_constructor's vcf-query call).
#
# Contract:
#   Rscript run_sigma.R \
#     --maf <path>            single-sample MAF (TSV), already run through
#                              sigma_filter.prepare_sigma_maf() - PASS +
#                              population-AF + VAF-exclude-bands filtered,
#                              still standard MAF columns (Chromosome,
#                              Start_Position, End_Position,
#                              Reference_Allele, Tumor_Seq_Allele1/2,
#                              Tumor_Sample_Barcode, Variant_Type).
#     --sample-id <str>        used for logging and to force
#                              Tumor_Sample_Barcode if the MAF's own value
#                              looks wrong/blank.
#     --tumor-type <str>       SigMA tumor_type, resolved Python-side by
#                              sigma_cancer_type_map.get_sigma_call_params().
#     --data-platform <str>    SigMA `data` argument, e.g. "msk"
#                              ([SigMA] DATA_PLATFORM in conf.ini).
#     --do-mva <TRUE|FALSE>    resolved Python-side (forced FALSE for
#                              medullo/ewing - see
#                              DO_MVA_UNSAFE_FOR_PANEL_DATA).
#     --check-msi <TRUE|FALSE> [SigMA] CHECK_MSI.
#     --catalog-name <str>     SigMA run()'s required catalog_name, e.g.
#                              "cosmic_v2_inhouse" (resolved Python-side
#                              from [SigMA] COSMIC_VERSION).
#     --lite-format <TRUE|FALSE>
#     --snv-cutoff <int>       SigMA's own run(snv_cutoff=...); walk.py
#                              already screens samples against this same
#                              value before ever invoking this script, this
#                              is a defensive second check, not the primary
#                              gate.
#     --ref-genome <hg19|hg38> make_matrix()'s ref_genome_name. Varan is
#                              hg19-only today (see
#                              SIGMA_INTEGRATION_FEASIBILITY.md section 3).
#     --output <path>          where the one-row result CSV is written.
#
# Exit codes: 0 on success (output file written), non-zero on any failure
# (missing package, malformed MAF, SigMA's own stop()/error, below
# snv_cutoff). walk.py treats any non-zero exit or missing/empty output
# file as "skip this sample's SigMA result, log a warning, keep going" -
# never a fatal error for the whole batch.
#
# --- Known upstream quirk this wrapper works around (NOT independently
# executed against a live SigMA install before this comment was written -
# this analysis was later verified by an actual smoke test, see
# SIGMA_INTEGRATION_FEASIBILITY.md / commit history for the smoke-test
# results this round) ---
# SigMA's own make_matrix(file_type = "maf") dispatches to
# .make_matrix_from_maf_list(). Reading that source directly
# (R/make_matrix.R on GitHub): when given exactly ONE maf file AND that
# file contains only one distinct Tumor_Sample_Barcode (i.e. exactly
# Varan's per-sample use case), the function's own branch
# (`if(length(tumors) > 1){ ... }`) never assigns `matrix_snvs` at all,
# so the function goes on to fail evaluating `rownames(matrix_snvs)`.
# SigMA's own shipped example (examples/test_maf.R) sidesteps this because
# it deliberately uses a 50-sample MAF file, not a single-sample one - so
# this bug is not something the upstream project's own tests would catch.
# The `is_list = TRUE` code path in make_matrix() dispatches into the
# *other*, working branch of .make_matrix_from_maf_list() (the multi-file
# branch, which builds one spectrum column per list entry regardless of
# how many Tumor_Sample_Barcode values are inside each individual file) as
# long as the list has 2+ entries - so this wrapper passes the same
# single-sample MAF path twice via is_list = TRUE, then keeps only the
# first (identical) column. This uses only SigMA's public, documented
# `is_list` parameter - no internal/unexported SigMA function is touched -
# but it is still a workaround for what reads as an upstream edge-case gap,
# not a documented, sanctioned usage pattern. Flagged for a human reviewer;
# re-check against a newer SigMA release if this wrapper ever starts
# failing on the "single tumor in file" case.

suppressPackageStartupMessages({
  ok <- requireNamespace("SigMA", quietly = TRUE)
})
if (!ok) {
  message("ERROR: SigMA package is not installed in this R environment.")
  quit(status = 2)
}
suppressPackageStartupMessages(library(SigMA))

parse_args <- function(argv) {
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument (expected --flag value): %s", key))
    }
    key <- sub("^--", "", key)
    key <- gsub("-", "_", key)
    if (i == length(argv)) {
      stop(sprintf("Missing value for --%s", key))
    }
    out[[key]] <- argv[[i + 1]]
    i <- i + 2
  }
  out
}

as_r_bool <- function(x) {
  toupper(trimws(x)) %in% c("TRUE", "T", "1", "YES")
}

argv <- commandArgs(trailingOnly = TRUE)
args <- tryCatch(parse_args(argv), error = function(e) {
  message(sprintf("ERROR: could not parse run_sigma.R arguments: %s", conditionMessage(e)))
  quit(status = 2)
})

required <- c("maf", "sample_id", "tumor_type", "data_platform", "do_mva",
              "check_msi", "catalog_name", "lite_format", "snv_cutoff",
              "ref_genome", "output")
missing <- setdiff(required, names(args))
if (length(missing) > 0) {
  message(sprintf("ERROR: missing required argument(s): %s",
                   paste0("--", missing, collapse = ", ")))
  quit(status = 2)
}

maf_path <- args$maf
sample_id <- args$sample_id
tumor_type <- args$tumor_type
data_platform <- args$data_platform
do_mva <- as_r_bool(args$do_mva)
check_msi <- as_r_bool(args$check_msi)
catalog_name <- args$catalog_name
lite_format <- as_r_bool(args$lite_format)
snv_cutoff <- suppressWarnings(as.integer(args$snv_cutoff))
ref_genome <- args$ref_genome
output_path <- args$output

if (is.na(snv_cutoff)) snv_cutoff <- 5L

if (!file.exists(maf_path)) {
  message(sprintf("ERROR: MAF file not found: %s", maf_path))
  quit(status = 2)
}

maf <- tryCatch(
  read.delim(maf_path, sep = "\t", header = TRUE, comment.char = "#",
             stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) {
    message(sprintf("ERROR: could not read MAF %s: %s", maf_path, conditionMessage(e)))
    NULL
  })
if (is.null(maf) || nrow(maf) == 0) {
  message(sprintf("ERROR: MAF %s is empty or unreadable after filtering - nothing for SigMA to run on for sample %s.",
                   maf_path, sample_id))
  quit(status = 3)
}

# Defensive: prepare_sigma_maf()'s output should already carry
# Tumor_Sample_Barcode from vcf2maf, but force it to the sample_id we were
# told to use if it's missing/blank/inconsistent, since make_matrix()'s
# maf-list branch (see workaround above) treats the whole file as one
# sample's spectrum regardless of this column's content - this is only
# used for the output row's own bookkeeping / sanity, not the spectrum
# computation itself.
if (!("Tumor_Sample_Barcode" %in% colnames(maf)) ||
    length(unique(maf$Tumor_Sample_Barcode)) == 0 ||
    all(trimws(maf$Tumor_Sample_Barcode) == "")) {
  maf$Tumor_Sample_Barcode <- sample_id
}

# --- Build the 96-dim trinucleotide spectrum ---
# See the workaround note at the top of this file: pass the same
# single-sample MAF path twice via is_list = TRUE to force make_matrix()
# into its working multi-file branch instead of the single-file/
# single-tumor branch that (per source reading) never assigns its result.
genomes_matrix <- tryCatch(
  SigMA::make_matrix(
    directory = c(maf_path, maf_path),
    file_type = "maf",
    is_list = TRUE,
    ref_genome_name = ref_genome),
  error = function(e) {
    message(sprintf("ERROR: SigMA::make_matrix() failed for sample %s: %s",
                     sample_id, conditionMessage(e)))
    NULL
  })
if (is.null(genomes_matrix)) quit(status = 3)

genomes <- SigMA::conv_snv_matrix_to_df(genomes_matrix)
genomes <- genomes[1, , drop = FALSE]
genomes$tumor <- sample_id

total_snvs <- sum(as.numeric(genomes[1, 1:96]))
message(sprintf("Sample %s: %d SNVs in SigMA-ready spectrum (snv_cutoff=%d).",
                 sample_id, total_snvs, snv_cutoff))

result <- tryCatch(
  SigMA::run(
    input_df = genomes,
    data = data_platform,
    tumor_type = tumor_type,
    catalog_name = catalog_name,
    do_assign = TRUE,
    do_mva = do_mva,
    check_msi = check_msi,
    lite_format = lite_format,
    snv_cutoff = snv_cutoff,
    return_df = TRUE),
  error = function(e) {
    message(sprintf("ERROR: SigMA::run() failed for sample %s (tumor_type=%s, data=%s, do_mva=%s): %s",
                     sample_id, tumor_type, data_platform, do_mva, conditionMessage(e)))
    NULL
  })

if (is.null(result) || nrow(result) == 0) {
  message(sprintf("ERROR: SigMA::run() produced no result for sample %s - likely below snv_cutoff (%d) or an unsupported tumor_type/data combination.",
                   sample_id, snv_cutoff))
  quit(status = 3)
}

# IMPORTANT, confirmed by reading R/run.R directly (and by an actual smoke
# test run this round - see SIGMA_INTEGRATION_FEASIBILITY.md / commit
# history): run()'s own `lite_format` argument is a no-op when
# `return_df = TRUE` - SigMA's run() only calls lite_df() on the internal
# "write to file" path (the `else` branch after `if(return_df) return(...)`
# in its source), never on the returned data frame. Since this wrapper
# always uses return_df = TRUE (to avoid SigMA's own output-file-naming
# logic and read it back ourselves), lite_df() must be called explicitly
# here to actually honor --lite-format - otherwise every column
# (including all the individual Signature_3_c{N}_ml / 96-context columns)
# comes back regardless of the flag, and the summary `categ` column (which
# only lite_df() computes) is silently never present at all.
if (lite_format) {
  result <- tryCatch(
    SigMA::lite_df(result),
    error = function(e) {
      message(sprintf("WARNING: SigMA::lite_df() failed for sample %s (%s) - falling back to the full (non-lite) result.",
                       sample_id, conditionMessage(e)))
      result
    })
}

result$SAMPLE_ID <- sample_id
write.csv(result, output_path, row.names = FALSE)

if (!file.exists(output_path) || file.info(output_path)$size == 0) {
  message(sprintf("ERROR: output file %s was not written or is empty.", output_path))
  quit(status = 3)
}

message(sprintf("SigMA result for sample %s written to %s", sample_id, output_path))
quit(status = 0)
