# Bundled species GEM assets

RegCompassR 1.8.3 distributes two transformed, compressed metabolic models so
canonical human and mouse workflows can start without network access:

| Species | Upstream model | Release | Installed file |
|---|---|---|---|
| Human | SysBioChalmers/Human-GEM | v2.0.0 | `Human2_2.0.0_regcompass.rds` |
| Mouse | SysBioChalmers/Mouse-GEM | v1.8.0 | `Mouse_1.8.0_regcompass.rds` |

The upstream YAML model, reaction annotations, and gene annotations were
converted to the RegCompass sparse GEM structure. Reaction roles, GPR tables,
external database identifiers, model source, release, checksum, and citation
metadata are retained. The files are modified/converted assets rather than
unmodified upstream archives.

Both upstream repositories distribute their model content under the Creative
Commons Attribution 4.0 International license (CC BY 4.0). Users should retain
attribution to the original model projects and cite the corresponding model
publication. The installed `manifest.tsv` records release, source, RDS checksum,
file size, citation DOI, and license.

The package default is offline loading:

```r
gem <- rc_prepare_gem("human")
```

To explicitly require the installed asset:

```r
gem <- rc_prepare_gem("human", source = "bundled")
```

To rebuild from an official pinned release or prepare an updated release:

```r
gem <- rc_prepare_gem(
  "human",
  version = "2.0.0",
  source = "download",
  force_download = TRUE
)
```

`scripts/build-bundled-gems.R` is the reproducible maintainer build script for
these installed assets.
