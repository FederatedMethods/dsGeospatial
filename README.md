# dsGeospatial


**Privacy-preserving geospatial visualisation and analysis for health and environmental data in DataSHIELD**

---

## Overview

**dsGeospatial** is a DataSHIELD package that provides federated methods for the 
visualisation and spatial analysis of health and environmental data while preserving 
the privacy of individual-level records.

The package has been developed as part of the **GROVE** project within the **DARE UK** programme, 
with the aim of enabling secure, collaborative geospatial analyses across distributed 
datasets without transferring sensitive household-level data.

Unlike conventional geospatial workflows, **dsGeospatial** performs all 
computations within the secure DataSHIELD environment. 
Only non-disclosive summary statistics and visualisations are returned to the 
analyst, allowing organisations to collaboratively investigate spatial patterns 
and environmental-health relationships while maintaining data governance and confidentiality.

---

## Motivation

Many important public health questions require linking health outcomes with environmental exposures such as:

* Greenspace accessibility
* Air pollution
* Noise exposure
* Urban infrastructure
* Socio-economic indicators
* Household location characteristics

These datasets are often held by different organisations and cannot be combined 
because of privacy, governance or legal constraints.

**dsGeospatial** provides a federated analytical framework that enables these 
data to be explored safely without exposing household-level information.

---

## Features

Current and planned functionality includes:

### Geospatial visualisation

* Heatmaps of region-level summaries
* Choropleth maps - Geographic hotspot identification
* Spatial correlation analysis

### Federated analytics

* Distributed computation across multiple DataSHIELD servers
* Automatic aggregation of study-level results
* Privacy-preserving disclosure controls
* Compatible with Opal, Armadillo and DSLite environments

---

## Example applications

The package can be used to investigate questions such as:

* How does asthma prevalence vary with proximity to greenspace?
* Which neighbourhoods have the highest average air pollution exposure?
* Are environmental exposures distributed equally across socio-economic groups?
* Can spatial trends be reproduced across multiple organisations without sharing household-level data?

---

## Privacy and disclosure protection

All analyses are performed within the DataSHIELD framework.

No individual-level health or environmental records leave the secure data repositories. 
Results returned to the analyst consist only of disclosure-controlled aggregated outputs that comply with DataSHIELD disclosure rules.


---

## Development status

**dsGeospatial** is currently under active development.

Planned future capabilities include:

---

## Contributing

Contributions, feature requests and bug reports are welcome.

Please see:

* `CONTRIBUTING.md`

for information on contributing to the project.

---

## License

This package is released under the license specified in the `DESCRIPTION` file.

