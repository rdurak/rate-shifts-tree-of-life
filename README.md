# Site-Specific Rate Shifts Across the Tree of Life
This repository contains the computational workflows used to investigate site-specific evolutionary rate shifts across orthologous protein families spanning the tree of life.

## Repository Structure

```text
scripts/

├── 01_data_preparation/
├── 02_alignments/
├── 03_raser/
├── 04_structural_mapping/
├── 05_statistics/
├── 06_structural_clustering/
├── 07_phylogenetic_effects/
└── hpc/
```

## Main Dependencies

### Python

* Python >= 3.10
* Biopython
* pandas
* numpy
* matplotlib
* requests
* ete3
* scipy
* statsmodels

### R

* base R
* readr
* stats

### External Software

* MAFFT
* MUSCLE5
* IQ-TREE2
* RASER
* SGEDTools
* Bio++

## Repository Overview

## 01_data_preparation

Preparation of orthogroup datasets and identification of experimentally resolved protein structures.

### Scripts

* `01_create_orthogroup_fastas.py`
  Generate one FASTA file per orthogroup.

* `02_check_pdb_availability.py`
  Screen orthogroups for experimentally resolved structures.

* `03_collect_pdb_orthogroups.py`
  Retrieve all experimental PDB entries associated with each orthogroup.

---

## 02_alignments

Generation and processing of multiple sequence alignments.

### Scripts

* `01_create_alignment_directories.py`
  Organize MAFFT and MUSCLE5 alignments into orthogroup directories.

* `02_create_alnscorebpp_config.py`
  Generate AlnScores.bpp configuration files.

* `03_create_seqmanbpp_config.py`
  Generate SeqMan.bpp configuration files.

* `04_run_alnscorebpp.sh`
  Run Bio++ alignment quality assessment.

* `05_run_seqmanbpp.sh`
  Run Bio++ alignment processing.

---

## 03_raser

Preparation and analysis of RASER rate-shift inference runs.

### Scripts

* `01_create_raser_inputs.py`
  Create orthogroup-specific RASER working directories.

* `02_create_parameter_files.sh`
  Generate orthogroup-specific RASER parameter files.

* `03_run_lrt.py`
  Perform likelihood-ratio tests between RASER and null models.

---

## 04_structural_mapping

Mapping inferred rate-shift sites onto experimentally resolved protein structures.

### Scripts

* `01_sged_pipeline.sh`
  Execute the SGEDTools structural mapping workflow.

---

## 05_statistics

Statistical analyses of orthogroup-level and site-level rate-shift patterns.

### Scripts

* `01_family_level_nonparametric.py`
  Perform non-parametric analyses of orthogroup-level rate-shift prevalence.

* `02_partial_correlation.R`
  Perform residual-based partial correlation analyses.

* `03_glmm_baseline.py`
  Fit a baseline site-level logistic mixed model.

* `04_glmm_interaction.py`
  Fit logistic mixed models with ancestry–structure interactions.

* `05_glmm_kingdom_level.py`
  Fit separate logistic mixed models for major evolutionary domains.

---

## 06_structural_clustering

Analysis of spatial clustering of heterotachous residues within protein structures.

### Scripts

* `01_prepare_mapped_heterotachy_observed.py`
  Prepare observed structurally mapped heterotachous residues.

* `02_prepare_mapped_heterotachy_random.py`
  Generate randomized residue datasets.

* `03_compute_contact_subgraphs_observed.py`
  Compute observed residue contact subgraphs.

* `04_compute_contact_subgraphs_random.py`
  Compute randomized residue contact subgraphs.

* `05_run_nsub_analysis.py`
  Perform permutation-based structural clustering analyses.

* `06_analyze_nsub_results.py`
  Summarize and visualize structural clustering results.

---

## 07_phylogenetic_effects

Evaluation of the phylogenetic consequences of site-specific rate shifts.

### Scripts

* `01_select_phylogenetic_effect_families.py`
  Select orthogroups for phylogenetic effect analyses.

* `02_prepare_phylogenetic_effect_dataset.py`
  Prepare datasets for phylogenetic effect analyses.

* `03_make_modified_alignments.py`
  Generate modified alignments with altered rate-shift sites.

* `04_compute_rf_bl_bs.py`
  Compute topological, branch-length, and support differences among trees.

* `05_lm_rf_bl_bs.py`
  Model factors associated with phylogenetic effects.

---

## hpc

SLURM submission scripts used for large-scale analyses on HPC systems.

### Scripts

* `01_run_mafft_alignments.slurm`
  Run MAFFT alignments.

* `02_run_muscle5_alignments.slurm`
  Run MUSCLE5 alignments.

* `03_run_iqtree_pdb_families.slurm`
  Infer phylogenetic trees with IQ-TREE2.

* `04_run_raser.slurm`
  Execute large-scale RASER analyses.

---

## Notes

Large intermediate datasets are not included in this repository. The data supporting this study are available on Edmond (Durak and Dutheil 2026): https://doi.org/10.17617/3.HCC2ME

## Citation

If you use this repository, please cite:

Durak et al. 2026. The Shifting Tempo of Evolution: Mapping Site-Specific Rate
Shifts Across the Tree of Life




