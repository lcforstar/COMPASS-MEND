# Native output profile viewer

Run from the project root after a model run:

```matlab
addpath('PlotResults');
fig = plot_profiles('data/output/changwu');
```

The five panels show water, total NH4/NO3, model SOC and free DOC at the last
written accepted step, using the units and geometry in the native NetCDF.
PNG/PDF files are written inside that run's `figures/` directory.

This is a native-state viewer, not the observation-matching evaluation used for
the accepted annual example. Node values, layer inventories and20cm dry-mass
means have different meanings; do not relabel the carbon bulk-volume units as
mg/kg, add QOC to free DOC, or calculate observation R2 from unmatched values.
