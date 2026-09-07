# *Trichoderma harzianum* metabolic-model reconstruction raw scripts

This repository contains all the scripts used for the reconstruction of the genomic scale metabolic model of *Trichoderma harzianum*.

The main organism is *T. harzianum* strain THM10 (`THM10.fasta`) obtained from:
https://mycocosm.jgi.doe.gov/TriharM10_1/TriharM10_1.home.html

## Workflow

The numbered scripts are intended to be run in ascending order. Run them

| Number | Script | Purpose | Main outputs |
| --- | --- | --- | --- |
| 0 | `s0_models_to_mat.m` | Import `.xml`/`.sbml` COBRA models and save them as MATLAB structures. | `<model>.mat` |
| 1 | `s1_get_genes_models.m` | Export the genes stored in each `.mat` model. | `<model>_genes.csv` |
| 2 | `s2_clean1.m` | Match model genes against GTF and protein FASTA records and create normalized FASTA files from standarnd nomenclature of models (locus tag). | `<model>_protein_clean_1.faa`, `<model>_protein_clean_reduced.faa` `<model>.fasta`  |
| 2a | `s2_clean2.m` through `s2_clean7.m` | Special-case FASTA cleaning for individual models whose identifiers require manual handling (not using locus tag). | `<model>_protein_clean_1.faa`, `<model>_protein_clean_reduced.faa` `<model>.fasta` |
| 3 | `s3_blast.m` | Run BLAST from `THM10.fasta` against the cleaned protein sets. Filters use e-value `<= 1e-5`, alignment length `>= 20`, and identity `>= 30`. | `blast_<X>.mat` |
| 4 | `s4_blast_export_csv.m` | Export BLAST hits tables to CSV. Reduced results are written to `hits_reduced/`. | `Hits_THM10_vs_*.csv` |
| 5 | `s5_hits_analysis.ipynb` | Interactive exploration of BLAST hits. | Notebook analyses |
| 6 | `hits_reduced/s6_rxns_subsystem.m` | Map reduced BLAST model IDs to model reactions and subsystems. Run from `hits_reduced/`. | Reaction/subsystem mapping CSV files |
| 7 | `s7_check_genes_rxns.m` | Calculate gene and reaction coverage for reduced BLAST hits. | `resultados_1.csv` |
| 8 | `s8_genes_rxns_analysis.ipynb` | Explore gene/reaction coverage results. | Notebook analyses |
| 9 | `s9_compare_mat.m` | Compare COBRA model fields and reaction identifiers against a reference model and BiGG data. | Command-window comparison report |
| 10 | `s10_biomasrxn.m` | Inspect biomass reactions and their metabolite composition in selected models. | MATLAB tables `T1`, `T2`, `T3`, ... |
| 11 | `s11_DraftModel.m` | Build draft models by homology from BLAST results and prepare a template-based reconstruction. | `workspaces/workspaceDraftModels_THM10.mat` and draft structures |
| 11a | `s11_draftmodels_analysis.ipynb` | Analyze draft-model statistics and results. | Notebook analyses |
| 12 | `s12_API.ipynb` | Produce or inspect API-derived pathway resources. | HTML files in `API/` |
| 13 | `s13_KEGG.m` | Join THM10 EC/KEGG annotations with BiGG reaction and metabolite tables. | MATLAB tables for in-model and missing KEGG entries |

## Requirements

- MATLAB with the COBRA Toolbox installed and initialized with `initCobraToolbox()`.
- COBRA Toolbox functions used by the scripts, including `readCbModel`, `buildRxnGeneMat`, `generateRules`, `getModelFromHomology`, `findExcRxns`, `findMetsFromRxns`, `optimizeCbModel`, and `createMetComp`.
- MATLAB Bioinformatics Toolbox functions such as `fastaread` and `fastawrite`.
- MATLAB parallel support for the `parfor` sections in `s3_blast.m` and `s11_DraftModel.m`.
- Jupyter support for the `.ipynb` notebooks. The notebooks may require a Python environment and packages not specified by this repository.

Before running the pipeline, add the COBRA Toolbox to the MATLAB path and verify that `initCobraToolbox()` completes successfully.

## Input data

The repository includes several classes of input and reference data:

- `THM10.fasta`: reference protein FASTA for *T. harzianum* THM10.
- `<model>.xml` and `<model>.mat`: COBRA metabolic models. XML files are ignored by Git, while representative MAT files are kept in the working tree.
- `<model>_protein.faa`: original protein FASTA files for the model organisms.
- `<model>_gtf.gtf`: optional annotations used to map `protein_id` values to `locus_tag` values.
- `<model>_genes.csv`: model gene exports.
- `bigg_models_reactions.txt` and `bigg_models_metabolites.txt`: BiGG reference tables used by `s9_compare_mat.m` and `s13_KEGG.m`.
- `THM10_EC_info.csv`, `THM10_info.csv`, `ProteinID_EC_info.csv`, and `ProteinIDs_faltantes.csv`: THM10 annotation and missing-identifier tables.
- `Biomass_analysis.xlsx`: biomass-related analysis data.

Model files are expected to use a consistent basename. For example, `iCNG99.mat`, `iCNG99_protein.faa`, and `iCNG99_gtf.gtf` belong to the same model. The cleaning scripts rely on the suffixes `_protein.faa` and `_gtf.gtf`.

## Running the core pipeline

```matlab
cd('path/to/t_harzianum_work')
initCobraToolbox()

% Run each stage in the MATLAB editor or command window:
s0_models_to_mat
s1_get_genes_models
s2_clean1
s3_blast
s4_blast_export_csv
s7_check_genes_rxns
```

For draft-model reconstruction, configure and review `s11_DraftModel.m` first. It currently assumes `THM10.fasta`, uses files discovered in the current directory, and selects a configured template model near the end of the script. Run `s6_rxns_subsystem.m` from `hits_reduced/` after the reduced hit CSV files have been generated.

Several scripts modify or regenerate files in place. Keep a copy of important model files and inspect the generated tables before using them in downstream reconstruction.

## Helper functions and utilities

- `reconstructionOpt2.m`: transfers homology-derived gene-reaction rules to a template model, rebuilds gene matrices/rules, adds ATP maintenance when needed, and prepares an optimized model.
- `addReaction2.m`: COBRA-style helper for adding reactions while checking duplicate reactions and preserving confidence metadata.
- `exchangeSingleModel.m`: reports exchange reactions, bounds, and optimized fluxes for a model.
- `KlebVarii.m`: earlier/experimental workflow for reciprocal BLAST-style comparisons and model linking.
- `s_check_mat.m`: normalizes model fields such as `metComps`, `description`, and `id` across MAT files.

## Generated files and Git behavior

The `.gitignore` intentionally excludes:

- CSV files generated during analysis.
- Large BLAST workspaces (`blast_THM10.mat` and `blast_reduced_THM10.mat`).
- XML model inputs.
- MATLAB autosave files.
- `workspaces/` and `API/` generated/working directories.

If a generated result is needed for publication or review, export it explicitly and document the script, input files, and parameter settings used to create it.

## Reproducibility notes

- Run scripts from the directory they expect; many discover inputs with `dir`, `pwd`, or relative paths.
- Review hard-coded filenames and model indices before rerunning, especially in `s9_compare_mat.m`, `s10_biomasrxn.m`, `s11_DraftModel.m`, and the `s2_clean2`–`s2_clean7` scripts.
- BLAST thresholds are embedded in the scripts. Record any changes to e-value, alignment length, identity, or homology-reconstruction parameters.
- Confirm that model gene identifiers match FASTA headers. Models that do not use locus tags may require manual identifier correction.
- The pipeline has no automated test suite; validation is currently performed through the generated summaries, model-structure checks, and notebook analyses.

## Repository structure

```text
.
├── s0_*.m ... s13_*.m       Numbered reconstruction and analysis stages
├── *.m                       Reusable helpers and experimental workflows
├── *.mat                     COBRA models and saved intermediate data
├── *.faa, *.fasta            Protein sequences
├── *.gtf, *_genes.csv        Gene annotations and model gene exports
├── hits_reduced/              Reduced-hit analysis script and results
├── *.ipynb                    Interactive analyses and API workflow
├── API/                       Generated pathway/API HTML resources
└── workspaces/                Local MATLAB intermediate workspaces
```

