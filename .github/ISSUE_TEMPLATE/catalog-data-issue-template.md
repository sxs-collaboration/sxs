---
name: Catalog data issue template
about: For users to file issues about catalog data or metadata
title: "[DATA] <SXS ID>v<version>/Lev<Lev> <problem>"
labels: catalog-data
assignees: ''

---

<!--
Please set the title above to include the SXS ID of the simulation, the catalog version number, the resolution Lev, and if appropriate, the extrapolation order.
-->

**Which data has an issue?**
Include here the catalog tag (if you have not specified it, you can get it from `sxs.load("simulations").tag`).
Include here the SXS ID of the simulation, the catalog version number, the resolution Lev, and if appropriate, the extrapolation order. Is the issue with metadata, strain data, Psi_4 data, or apparent horizons data?

**Describe the issue**
A clear and concise description of what the issue is.

**To Reproduce**
Please provide a minimal example using the `sxs` package to reproduce the issue. Please include the version of the `sxs` package you're using that experience the problem (from `sxs.__version__`).

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Additional context**
Add any other context about the problem here.
