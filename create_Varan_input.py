#!/usr/bin/env python3
"""Standard preprocessing step for converting a third-party panel-
sequencing vendor's raw run output into the sample.tsv / patient.tsv /
{run_id}_fusions.tsv shape Varan's own `varan.py -i` flag consumes.

This is a thin CLI wrapper: all vendor-specific conversion logic lives
under the `vendor_adapters/` package (one module per vendor - see
`vendor_adapters/__init__.py` for the adapter interface and
`vendor_adapters/guardant.py` for the reference implementation). This
script itself only parses arguments and dispatches to the selected
vendor's `run()` function - it intentionally does not contain any
vendor-specific parsing so a second vendor never requires touching this
file's logic, only registering a new adapter module.

It remains "extra tooling" alongside Varan rather than wired into
`varan.py`'s own CLI: `varan.py` already owns a dense set of single-letter
flags (-f, -s, -c, ... - see its own -i/--varan_input, which is what this
script's output feeds into) and is actively evolving core pipeline logic
with no test suite backing it. Keeping vendor conversion as a separate,
standard preprocessing step run before `-i` avoids adding vendor-dispatch
risk to that already-complex surface. See MULTIVENDOR_INTEGRATION_NOTES.md
for the full reasoning and for vendor-specific open questions that still
need a human decision.

Usage:
    python create_Varan_input.py --vendor guardant -f <s3_run_folder>
    python create_Varan_input.py --vendor guardant -s <selection.tsv>
"""

import argparse
import sys

from vendor_adapters import ADAPTERS, get_adapter


def main():
    parser = argparse.ArgumentParser(
        description="Varan multi-vendor input generator: converts a "
                     "vendor's raw run output into sample.tsv/fusions.tsv "
                     "for `varan.py -i`.")
    parser.add_argument(
        "--vendor", default="guardant", choices=sorted(ADAPTERS),
        help="Which vendor adapter to use (default: guardant).")
    parser.add_argument("-f", "--folder", help="Path to the vendor's run folder (e.g. an S3 folder).")
    parser.add_argument("-s", "--selection", help="TSV file listing samples to process across multiple run folders.")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if not args.folder and not args.selection:
        parser.error("one of -f/--folder or -s/--selection is required")

    adapter = get_adapter(args.vendor)
    adapter.run(folder=args.folder, selection=args.selection)


if __name__ == "__main__":
    main()
