# mfd_core
The scripts in this repository are part of the [Microflora Danica project](https://github.com/cmc-aau/mfd_wiki/wiki). 


Be advised, that for the metagenomic-derived data, the term "OTU" is only used due to format requirements by ampvis2, and they do not represent classical OTUs. 
The generated profiles can be thought of as taxonomic bins. 

## Scripts
### Metagenomic 16S data 
`scripts/core_microbes.R` identifies core genera using a mean relative abudane filter of 0.1% and a habitat-specific prevalence of 50% across all levels of the MFD Ontology. 


`scripts/core_summaries.R` performs a summaries of the identified core genera. 


`scripts/core_analysis.R` recreates the analysis of the core genera across the MFDO1 ontology level. The script also renders the figures used in the manuscript. 

## Data
The scripts rely on data files available from the MFD [github](https://github.com/cmc-aau/mfd_metadata) and the MFD Zenodo [repo](https://zenodo.org/records/12605769), from where the output files used in the manuscript are also available. 
