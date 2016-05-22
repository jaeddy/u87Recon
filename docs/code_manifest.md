# Code

Organized by language. Version tags in parentheses indicate the workflow iteration for which the code was used.

## Matlab scripts/functions

### Reading/formatting reference model

+ **`convertSBMLToCobra_Recon2.m`** (`1.0`, `0.5`) This is a modified version of the `convertSBMLToCobra` function in the COBRA toolbox (which, at the time, I had trouble getting to work with the Recon 2 SBML/XML file).

+ **`removeReconGeneMets.m`** (`1.0`, `0.5`) The Recon 2 SBML file includes genes as metabolites; this script strips those metabolites out.

+ **`HR2_addAstrocyte.m`** (`0.5`) This script adds reactions and metabolites from the astrocyte model in **GLU_norm** (published by Nathan Lewis) to Recon 2 (specifically, the model in **HR2_CbModel_Dec2012.mat**).

+ **`Recon203_updater.m`:** (`1.5`) This script was used to clean up / fill in some missing metabolite annotation information in the Recon2.03 model, based on data from a separate publication.

### Compiling/formatting evidence for mCADRE

+ **`U87GEODataProcessing.m`:** (`1.0`, `0.5`) This script was used to compile and format GEO gene expression evidence for use with mCADRE; the first part of the script, when not commented out, processed and normalizes raw data from CEL files. The final output of the script is stored in **U87_GEO_evidence.mat** and represents the expression *calls* from MAS5 (0: absent, 1: marginal, 2: present) for the intersection list of genes across all datasets. Note: I have an updated version of GEO data collection code and results (see `importU87AffyData.R`, `compileGEOData.m` scripts and `geo_new` data folder).

+ **`compileGEOData.m`:** (`1.5`) This script processes gene expression data from GEO (not normalization or presence/absence calls, which happens in `importU87AffyData.R`) to create inputs for mCADRE. In addition to the steps covered by `U87GEODataProcessing.m`, this script also calculates ubiquity scores.

### Preparing inputs for mCADRE

+ **`read_biomass.m`:** (`1.0`, `0.5`, `1.5`)This script reads the biomass equation from a single line in a text file and adds the corresponding reaction to a COBRA model.

+ **`write_biomass.m`:** (`1.0`, `0.5`, `1.5`) This script extracts biomass components and coefficients from an Excel spreadsheet and writes the formatted biomass equation as a single line in a text file.

+ **`primeRecon2forU87mCADRE.m` (`1.0` ) / `primeRecon203forU87mCADRE.m` (`1.5`):** These scripts share the same goal - modifying/augmenting/testing the generic Recon 2 model prior to running mCADRE to generate a U87-specific model. I'm fairly confident that `primeRecon2forU87mCADRE.m` was used to create both **Recon201.mat** and **Recon201_U87.mat** in iteration `1.0` of the workflow; however, I seem to have continued using and making changes to the script with later versions of Recon 2, so bits of code have been added or commented out. Unfortunately, I don't have a record of exactly what commands/parameters were used to produce **Recon201_U87.mat**.
    + Note: from what I can tell, these scripts never evolved past the point of being 'interactive' (a mix of automated processing steps and some debugging/exploratory steps) - you can see at the bottom that there's code for testing and troubleshooting production of biomass and problematic biomass components/precursors.
    + In terms of logic/approach, `primeRecon2forU87mCADRE.m` seems to be focused more around modifying and testing Recon 2 capabilities based on components in DMEM media; it also includes more direct modification of Recon 2 fluxes and flux bounds, largely related to fatty acid and lipid synthesis.
    + On the other hand, `primeRecon203forU87mCADRE.m` seems to examine functionality of Recon 2 based more on what are known to be "essential" macromolecules (amino acids, fatty acids, vitamins, etc.); some of these essential metabolite lists also appear to be used in `primeRecon2`, but I suspect these were added later. The latter script, `primeRecon203` makes no mention of and does not appear to consider DMEM.
    + As best as I can tell, `primeRecon2forU87mCADRE` was used to add the 4 reactions described in the first part of the script to Recon 2 (in **HR2_CbModel_Dec2012.mat**), and this model was saved as **Recon201.mat**. The script then added the U87 biomass reaction and adjusted exchange reactions to simulate DMEM media conditions; this model was then saved as **Recon201_U87.mat**.
    + I'm pretty sure that I was still experimenting with `primeRecon203forU87mCADRE.m` when I decided to move on to other projects, so I don't believe this script was ever used to produce a saved output.

+ **`U87_mCADRE_preprocessing.m`:** (`0.5`) This script represents an earlier version of the approach used in the `primeRecon2...` scripts. In this case, the starting point for Recon 2 is the model with astrocyte-specific reactions and metabolites added. Note: while I don't seem to have used any output from this script for subsequent steps in the workflow (i.e., I didn't generate a U87 model from the Recon 2 + astrocyte model), *the steps depicted are probably closer to the state of `primeRecon2forU87mCADRE.m` when it was used to generate* ***Recon201_U87*** *than what is currently found in that script*.

+ **`createU87mCADREInputs.m`:** (`1.0`) This script puts together all evidence from HPA, GEO, and an in-house RNA-seq experiment, defines the high-confidence core list of genes, and calculates ubiquity scores for mCADRE. These inputs are combined with the augmented Recon 2 model (**Recon201_U87.mat**) from `primeRecon2forU87mCADRE.m` and saved in **U87mCADREInputs.mat**.

**Note:** I can't seem to find the script/command I used to actually run mCADRE, which is unfortunate, because it probably means the code is lost on a cluster somewhere.

### Working with mCADRE output

+ **`compileU87mCADRE.m`:** (`1.0`) This script compiles, formats, and summarizes outputs from mCADRE to yield the `U87` model in `U87_CbModel.mat`. The script also generates some statistics such as the number of core reactions that are blocked in the final model and distribution of metabolites across cellular compartments.

+ **`U87DMEMGrowth.m`:** (`1.0`) This appears to test growth of the mCADRE-produced U87 model on DMEM; it resets bounds on exchange reactions to match available inputs from DMEM (similar to the preprocessing step in `primeRecon2forU87mCADRE` with the generic Recon 2 model).

+ **`formulateDMEM.m`:** (`1.0`) This is a function that performs pretty much the same exchange reaction adjustment steps as in `U87DMEMGrowth.m`; the function is called by `U87InhibitorScreen.m`.

+ **`U87InhibitorScreen.m`:** (`1.0`) This script uses enzyme inhibitor results from the Hockenbery lab to perform a sort of knock-out screen / validation test.

### Generic helper functions

+ **`textArrayRead.m`:** utility/wrapper function I wrote years ago to make it somewhat easier to read in text files (Matlab probably has something better than this by now).

+ **`textArrayWrite.m`:** see above for `textArrayRead`.

## Python scripts

### Compiling/formatting evidence for mCADRE

+ **`getHPA.py`:** (`1.0`, `0.5`, `1.5`) This script was used to scrape protein expression data from one of three GBM cell lines (U-87, U-138, or U-251) from HTML tables on the Human Protein Atlas (HPA) website; no idea if it still works (seems doubtful).

+ **`getHPA_neg.py`:** (`1.0`, `0.5`, `1.5`) This script was a modified version of `getHPA` used to collect 'negative' staining results from HPA (for possible use in validation).

+ **`HPAensembl2gene.py`:** (`1.0`, `0.5`, `1.5`) This script uses the ID table in `Homo_sapiens.gene_info` (from NCBI) to convert Ensembl IDs from HPA to NCBI Entrez gene ID; used to produce `U-87_MG_genes.txt`.

## R scripts

### Compiling/formatting evidence for mCADRE

+ **`importU87AffyData.R`:** (`1.5`) uses bioconductor packages and biomaRt to parse GEO datasets, preprocess and save data to text files, and create mapping between probe IDs and Ensembl ID.
