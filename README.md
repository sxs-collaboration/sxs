# SXS Simulations

This is a branch of the [`sxs` package](https://github.com/sxs-collaboration/sxs/)
repo that just contains a JSON file with all the metadata from all simulations
published by the SXS collaboration.  This is intended to be loaded by the `sxs`
package as
```python
import sxs

simulations = sxs.load("simulations")
```
See the main branch's documentation for details.
