# *Trichoderma harzianum M10* metabolic-model reconstruction raw scripts

This repository contains all the scripts used for the reconstruction of the genomic scale metabolic model of *Trichoderma harzianum*.

The main fasta file of *T. harzianum* strain THM10 (`THM10.fasta`) was obtained from:
https://mycocosm.jgi.doe.gov/TriharM10_1/TriharM10_1.home.html

## Workflow (Main scripts)

The numbered scripts are intended to be run in ascending order. Run them

| Number | Script | Purpose | Main outputs |
| --- | --- | --- | --- |
| 0 | `s0_models_to_mat.m` | Import `.xml`/`.sbml` COBRA models and save them as MATLAB structures. | `<model>.mat` |
| 1 | `s1_get_genes_models.m` | Export the genes stored in each `.mat` model. | `<model>_genes.csv` |
| 2 | `s2_clean1.m` | Match model genes against GTF and protein FASTA records and create normalized FASTA files from standarnd nomenclature of models (locus tag). | `<model>_protein_clean_1.faa (ALL protein anotations)`, `<model>_protein_clean_reduced.faa (ONLY protein anotations of genes inside models.mat)` |
| 2a | `s2_clean2.m` through `s2_clean7.m` | Special-case FASTA cleaning for individual models whose identifiers require manual handling (not using locus tag nomenclature). | `<model>_protein_clean_1.faa`, `<model>_protein_clean_reduced.faa` |
| 3 | `s3_blast.m` | Run BLAST from `THM10.fasta` against the cleaned protein sets. Filters use e-value `<= 1e-5`, alignment length `>= 20`, and identity `>= 30`. | `blast_<X>.mat` |
| 4 | `s4_blast_export_csv.m` | Export BLAST hits tables to CSV. Reduced results are written to `hits_reduced/`. | `Hits_THM10_vs_*.csv`, `hits_reduced/Hits_THM10_vs_*.csv`|
| 5 | `s5_hits_analysis.ipynb` | Exploration of BLAST hits. | `N/A Notebook analyses` |
| 6 | `hits_reduced/s6_rxns_subsystem.m` | Map reduced BLAST model IDs to model reactions and subsystems. Run from `hits_reduced/`. | `rxns_subsystems_*.csv (Reaction/subsystem mapping CSV files)`  |
| 7 | `s7_check_genes_rxns.m` | Calculate relation between gene -reaction coverage for the reduced BLAST hits. | `ModelsInfo_1.csv` |
| 8 | `s8_genes_rxns_analysis.ipynb` | Explore gene/reaction coverage results. | `N/A Notebook analyses` |
| 9 | `s9_compare_mat.m` | Compare models fields and reaction identifiers against a reference model (iML1515.mat)and BiGG data. | `N/A Command-window comparison report` |
| 10 | `s10_biomasrxn.m` | Inspect biomass reactions and their metabolite composition in selected models. | `N/A (Tables inside matlab has T1, T2 ...)` |
| 11 | `s11_DraftModel.m` | Build draft models by homology from BLAST results. | `workspaces/workspaceDraftModels_THM10.mat, results/THM10_draftmodel_v1.mat`  |
| 11a | `s11_draftmodels_analysis.ipynb` | Analyze draft-model statistics and results. | Notebook analyses |
| 12 | `s12_API.ipynb` | Web scraping of subsystems,Rxns,KEGG ids asociated to T. Harzianum in mycocom website: https://mycocosm.jgi.doe.gov/cgi-bin/metapathways?db=TriharM10_1 . | HTML files in `API/` |
| 13 | `s13_KEGG.m` | Analyze KEGG presence reactions of HTML files ids inside/not inside THM10_draftmodel_v1.mat | MATLAB tables for in-model and missing KEGG entries |

## Input data

The repository includes several classes of input and reference data:

- `THM10.fasta`: reference protein FASTA for *T. harzianum* THM10.
- `<model>.xml` and `<model>.mat`: COBRA metabolic models. 
- `<model>_protein.faa`: original protein FASTA files for the model organisms (.mat or .xml files).
- `<model>_gtf.gtf`: optional annotations used to map `protein_id` values to `locus_tag` values for the model organisms (.mat or .xml files).
- `bigg_models_reactions.txt` and `bigg_models_metabolites.txt`: BiGG reference tables obtained from web.


