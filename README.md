# Identifying POLE-mutant CRCs in a routine pathologic workflow
POLE-mutated CRCs in TCGA and their co-mutations in small NGS-panel approaches

Please refer to our publication: currently under review

### Project Overview

This repository contains a short R script for analyzing POLE-mutated colorectal cancers (CRC) within TCGA-COAD and TCGA-READ as an external validation cohort. The analysis includes a simulation of a small NGS panel (Illumina Focus Panel) to assess the feasibility of targeted sequencing in this context. We hypothesized that POLE-mutated CRCs are characterized by an increased numbers of mutations even when tested with a small panel, that does not cover the POLE gene itself. As small NGS-based panel approaches are used currently in the routine work-up of CRCs, this could be helpful in screening for the small but relevant subset of POLE-mutated CRCs. POLE-mutated CRCs show an exceptional respone to immunotherapeutic approaches.

For the clinical response of POLE-mutated CRCs to immune checkpoint inhibitors, please refer to [Ambrosini et al, Annals of Oncology 2024](https://linkinghub.elsevier.com/retrieve/pii/S0923-7534(24)00104-2), doi: 10.1016/j.annonc.2024.03.009.

### Purpose of This Script

This script serves as framework for how we accessed and analyzed TCGA data in our study.
It is not identical to the exact code used to generate all figures in the publication.
Specific file paths, data filtering, and processing steps may differ, especially for our own internal cohort.
For more detailed information, please contact the authors.

### Analysis Workflow

- TCGA Data Query & Processing.
- Retrieves colonic (COAD) and rectal (READ) adenocarcinoma mutation data.
- Extracts POLE-mutated CRC cases.
- Simulation of Small NGS Panel.
- Uses Illumina Focus Panel gene list.
- Filters for predefined SNVs.
- Identifies concurrent alterations in POLE-mutated cases.
- Generates oncoplots and mutation interaction maps.

### Citation

If you use this script, please cite also the relevant [TCGA-COAD & TCGA-READ](https://www.nature.com/articles/nature11252) and [TCGAbiolinks](https://academic.oup.com/nar/article/44/8/e71/2465925?login=true) publications:

Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D, Sabedot TS, Malta TM, Pagnotta SM, Castiglioni I, Ceccarelli M, Bontempi G, Noushmehr H. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. Nucleic Acids Res. 2016 May 5;44(8):e71. doi: 10.1093/nar/gkv1507. 

Cancer Genome Atlas Network. Comprehensive molecular characterization of human colon and rectal cancer. Nature. 2012 Jul 18;487(7407):330-7. doi: 10.1038/nature11252.

