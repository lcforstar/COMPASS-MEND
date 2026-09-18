# Global foliage clumping index

This directory stores the global foliage clumping-index dataset reserved for
future regional COMPASS simulations.

## Model-ready file

- `He_et_al_2012_1x1.nc`
- Variable: `ci`
- Grid: global 1 degree latitude-longitude grid
- Intended missing-data policy: replace missing ocean or unavailable cells with
  `1.0` only when materializing a land-model parameter field
- ClimaLand artifact tree: `b849eb95c09190095e7bf021494ddeda8858af01`
- Downloaded from the Julia package artifact mirror on 2026-07-16
- NetCDF SHA-256: `B0802D9F6FE3467BBED5DB1EB06C4E48330AD6B125FA53EB29364F8172860F73`
- Archive SHA-256: `32F168082923ECB11B9691482638396E0AEF73C1ED420E76FD181740A25D82DC`

This is the exact model-ready file referenced by ClimaLand's
`spatially_varying_parameters.jl`. It is retained as a regional-model input and
is not used by the current single-site Changwu run.

## Scientific source

The gridded product is derived from the global MODIS foliage clumping-index map:

> He, L., Chen, J. M., Pisek, J., Schaaf, C. B., and Strahler, A. H. (2012).
> Global clumping index map derived from the MODIS BRDF product. *Remote Sensing
> of Environment*, 119, 118-130. https://doi.org/10.1016/j.rse.2011.12.008

The archived 500 m source product is:

> He, L., Chen, J. M., Pisek, J., Schaaf, C., and Strahler, A. H. (2017).
> Global 500-m Foliage Clumping Index Data Derived from MODIS BRDF, 2006.
> ORNL DAAC. https://doi.org/10.3334/ORNLDAAC/1531

The ORNL product contains one global GeoTIFF, uses a scale factor of `0.01`,
uses `255` as NoData, and covers 2006. The 500 m archive is about 360 MB and
requires NASA Earthdata authentication, so it is not duplicated here.

## Selection priority

Current point simulations use:

1. a site value when `clumping_index_source = "manual"`;
2. otherwise the selected COMPASS PFT value.

Future regional simulations may add a third, higher-priority spatial source by
interpolating `ci` at each grid cell. That reader is deliberately not active yet.

The PFT constants in `vegetation.getPFTParameterSet` are a separate parameter
source: Chen et al. (2012), *Global Biogeochemical Cycles*, 26, GB1019,
Table 1, https://doi.org/10.1029/2010GB003996. They must not be described as
values sampled from this MODIS grid.

## Changwu reference value

The Changwu site is located at `107.6877 E, 35.2497 N`. Nearest-neighbour
sampling of this file selects the cell centered at `107.5 E, 35.5 N`, where
`ci = 0.7193`. The site configuration rounds this landscape-scale effective
value to `0.72`. This does not replace the Apple PFT prior of `0.80`, which
comes from the apple-orchard measurement reported by Liu et al. (2013),
Section 3.4.4 and Figure 10, https://doi.org/10.1016/j.scienta.2013.10.009.
