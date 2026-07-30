
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dsGeospatial : 

**Server-side package for privacy-preserving geospatial visualisation and analysis for health
and environmental data in DataSHIELD**

## Installation

You can install the development version of dsGeospatial from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("FederatedMethods/dsGeospatial")
```
**dsGeospatial** is currently under active development.
For a full list of development branches, checkout https://github.com/FederatedMethods/dsGeospatial/branches

## About

**dsGeospatial** is a DataSHIELD (https://www.datashield.org)) package that provides federated methods
for the visualisation and spatial analysis of health and environmental
data while preserving the privacy of individual-level records.

The package has been developed as part of the **GROVE** project within
the **DARE UK** programme, with the aim of enabling secure,
collaborative geospatial analyses across distributed datasets without
transferring sensitive household-level data.

Unlike conventional geospatial workflows, **dsGeospatial** performs all
computations within the secure DataSHIELD environment. Only
non-disclosive summary statistics and visualisations are returned to the
analyst. A key point to highlight is that the dsGeospatial package (https://github.com/FederatedMethods/dsGeospatial/) needs to be used in conjunction with the dsGeospatialClient package (https://github.com/FederatedMethods/dsGeospatialClient) - trying to use one without the other makes no sense.


## Contributing

Contributions, feature requests and bug reports are welcome.

Please see .github/:

- `CONTRIBUTING.md`
- `ISSUE_template.md`
- `DISCLOSUREissue_template.md`

for information on contributing to the project.
