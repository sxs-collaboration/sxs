# Metadata fields

This page describes the fields that may appear in the metadata of an
individual simulation.  It is generated automatically from
`sxs.metadata.metadata_docs`, which can also be inspected interactively:

```python
from sxs.metadata.metadata_docs import metadata_fields_dict
metadata_fields_dict["reference_time"]
```

Each field lists its type, the range of versions of the metadata format in
which it is present, and the kinds of simulations for which it is present.
For the versions, "2+" means the field was introduced in version 2 and is
still current, while "0–1" means it was introduced in version 0 and deprecated
in version 2.  Fields marked *computed by `sxs`* are filled in by the `sxs`
package rather than recorded directly by the simulation code.  Deprecated
fields are listed separately in each section, in a collapsed "Deprecated
fields" block.

<!--
  NOTE: Do not remove or edit the placeholder comment below.  At build time,
  the MkDocs hook in mkdocs_metadata_fields.py (registered under `hooks:` in
  mkdocs.yml) replaces it with the tables of metadata fields generated from
  sxs/metadata/metadata_docs.py, plus the Bibliography for any citations.
  Without it, this page will contain only the text above.
-->
<!-- metadata-fields-tables -->
