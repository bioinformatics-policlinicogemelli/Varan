"""Registry of vendor adapters for converting third-party panel-sequencing
vendor output into Varan's `sample.tsv` / `patient.tsv` /
`{run_id}_fusions.tsv` input shape (see `create_Varan_input.py` for the
CLI entry point, and `MULTIVENDOR_INTEGRATION_NOTES.md` for the design
rationale).

Each adapter module registered here is expected to expose:

  - `NAME: str` - the vendor's short identifier, used as the `--vendor`
    CLI value and as this dict's key.
  - `run(folder=None, selection=None, **kwargs) -> Optional[dict]` - a
    plain function (no argparse/sys.argv involved) that performs the
    conversion and returns `{"report_path", "fusion_table_path",
    "report_data"}`, or None if nothing was produced. `folder` is a
    single vendor run folder to process; `selection` is a TSV of
    `sample_id`/`s3_path_run` pairs spanning multiple run folders -
    both vendor-defined in their exact shape, but every adapter accepts
    at least these two ways of specifying input so the CLI wrapper can
    stay vendor-agnostic.

To add a second vendor: write `vendor_adapters/<vendor>.py` implementing
the interface above (reusing `vendor_adapters/common.py`'s S3 helpers,
`SampleRow`, `write_sample_tsv`, and `append_fusions_to_table` wherever
its file formats allow), then add one entry to `ADAPTERS` below. No
changes to `guardant.py` or to this file's existing entries are needed.
"""

from vendor_adapters import guardant

ADAPTERS = {
    guardant.NAME: guardant,
}


def get_adapter(name: str):
    """Look up a registered vendor adapter module by name.

    Raises KeyError with the list of known vendors if `name` isn't
    registered, rather than failing with an opaque AttributeError later.
    """
    try:
        return ADAPTERS[name]
    except KeyError:
        known = ", ".join(sorted(ADAPTERS)) or "(none registered)"
        raise KeyError(f"Unknown vendor {name!r}. Known vendors: {known}") from None
