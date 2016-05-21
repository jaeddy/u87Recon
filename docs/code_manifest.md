# Code

Organized by language.

## Matlab scripts/functions

+ **`textArrayRead.m`:** utility/wrapper function I wrote years ago to make it somewhat easier to read in text files (Matlab probably has something better than this by now).

+ **`textArrayWrite.m`:** see above for `textArrayRead`.

+ **`Recon203_updater`:** used to clean up / fill in some missing information in the Recon2.03 model, based on data from a separate publication.

+ **`primeRecon203forU87mCADRE.m` / `primeRecon2forU87mCADRE.m`:** These scripts all seem to share the same goal - modifying/augmenting/testing the generic Recon2 model prior to running mCADRE to generate a U87-specific model. `primeRecon2` appears to be the most recent version, but might want to refer to `primeRecon203` to transfer some comments/description/logic  
    + Note: from what I can tell, these scripts never evolved past the point of being 'interactive' (a mix of automated processing steps and some debugging/exploratory steps), so I'm not sure if they were ever used to produce a saved output.
    + Also: the scripts appear to use some functions from the Warburg modeling paper; I'll need to include these (and probably check permissions). However, it's possible I rewrote these functions myself at some point...

+ **`U87InhibitorScreen.m`:** This script uses enzyme inhibitor results from the Hockenbery lab to perform a sort of knock-out screen / validation test.

+ **`formulateDMEM.m`:** This appears to be a precursor to the `primeRecon2` scripts above. I believe the difference is that I used this script during an iteration of the project where I tried to include media and test for growth *after* running mCADRE (but I'll have to double-check); this function is called by `U87InhibitorScreen.m`.

+ **`U87DMEMGrowth.m`:** This looks like it does more or less the same thing as `formulateDMEM`, but is a script rather than a function.

+ **`compileU87mCADRE.m`:** This script compiles, formats, and summarizes outputs from mCADRE to yield the `U87` model in `U87_CbModel.mat` - appears to be based on Recon2.01.

**Note:** I can't seem to find the script/command I used to actually run mCADRE, which is unfortunate, because it probably means the code is lost on a cluster somewhere.

+ **`createU87mCADREInputs.m`:** This script takes care of the mCADRE pre-processing prior to running mCADRE with Recon1...? Or not: code seems to refer to Recon2.01. Need to sort this out!

+ **`read_biomass.m`:** This script reads the biomass equation from a single
line in a text file and adds the corresponding reaction to a COBRA model.

+ **`write_biomass.m`: This script extracts biomass components and coefficients
from an Excel spreadsheet and writes the formatted biomass equation as a single
line in a text file.

+ **`U87GEODataProcessing.m`:** another pre-processing step - compiling and formatting GEO gene expression evidence for use with mCADRE. Note: I appear to have an updated version of GEO data collection code and results (see `compileGEOData` script and `geo_new` data folder).

+ **`compileGEOData.m`:** process gene expression data from GEO to create inputs for mCADRE; uses outputs from R script `importU87AffyData`

+ **`convertSBMLToCobra_Recon2.m`**

+ **`removeReconGeneMets.m`**

+ **`HR2_addAstrocyte.m`**

## Python scripts

+ **`getHPA.py`:** This script was used to scrape protein expression data from one of three GBM cell lines (U-87, U-138, or U-251) from HTML tables on the Human Protein Atlas (HPA) website; no idea if it still works (seems doubtful).

+ **`getHPA_neg.py`:** This script was a modified version of `getHPA` used to collect 'negative' staining results from HPA (for possible use in validation).

+ **`getHPA_sub.py`:** just seems to be testing out some stuff; should discard.

+ **`HPAensembl2gene.py`:** This script uses the ID table in `Homo_sapiens.gene_info` (from NCBI) to convert Ensembl IDs from HPA to NCBI Entrez gene ID; used to produce `U-87_MG_genes.txt`.

## R scripts

+ **`importU87AffyData.R`:** uses bioconductor packages and biomaRt to parse GEO datasets, preprocess and save data to text files, and create mapping between probe IDs and Ensembl ID.
