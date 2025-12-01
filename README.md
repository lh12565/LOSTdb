# 🫁 LOSTdb
<p align="center">
  <a href="http://lostdbcancer.com:8080">
    <img src="https://img.shields.io/badge/LOSTdb-1.0-blue">
  </a>
  <a href="http://lostdbcancer.com:8080">
    <img src="https://img.shields.io/badge/status-active-success">
  </a>
  <br>
  <a href="https://zenodo.org/communities/lungcancer/records">
    <img src="https://img.shields.io/badge/Zenodo-v1-blue">
  </a>
  <a href="https://groups.google.com/g/lostdb">
    <img
      src="https://img.shields.io/badge/Google%20Group-Discussions-green?logo=google&logoColor=white" 
      alt="Google Group Discussions"
    >
  </a>
  <br>
  <a href="https://cn.vuejs.org/">
    <img src="https://img.shields.io/badge/vue-3.5.13-green" alt="Vue 3.5.13">
  </a>
  <a href="https://vuetifyjs.com/">
    <img src="https://img.shields.io/badge/vuetify-3.x-blue" alt="Vuetify 3.x">
  </a>
</p>

[LOSTdb](http://lostdbcancer.com:8080/): a manually curated multi-omics database for lung cancer research

![github](https://github.com/user-attachments/assets/cf803575-c019-4519-9936-b981b6b2cfba)

## Features

- **Multi-Omics Data Integration**: Seamlessly explore genomic, transcriptomic(Bulk, Single-cell), proteomic, and epigenomic data
- **Cross-Species Analysis**: Compare data from clinical samples, mouse models, and cell lines
- **Interactive Visualizations**: Dynamic charts and tables for data exploration
- **Subtype-Specific Analysis**: Focus on LUAD, LUSC, and SCLC subtypes
- **Tools Analysis**: Integrated analysis, Metadata analysis, Gene-Metadata analysis, Target and drug analysis
- **User-Friendly Interface**: Clean, modern design for easy navigation

This repository contains partial code and tables relating to LOSTdb.
## Scripts:
- RNA.R: Script containing RNA-seq analysis functions
- Proteomics.R: Script containing proteomics analysis functions
- 🧬Mutation.R: Script containing genomics analysis functions
- Methylation.R: Script containing methylation analysis functions
- scRNA.R: Script containing scRNA-seq analysis functions
- 🏥Sigclin.R: Script for significance analysis between clinical/meta data.
- 🎯Target.R: Script for target scoring.

## 📝Tables:
- datasets.xlsx: Dataset information in LOSTdb.
- subtypes.xlsx: Detailed information of molecular subtypes, including classical subtypes and multi-omics meta-program (MP) subtypes.
- clinical-metadata.xlsx: Manually curated clinical/meta data.
