# Workflows

Data processing steps for model generation, with inputs and outputs for each script used, organized by iteration/version.

## Iteration `1.0`

This is the most complete version of the workflow, which actually resulted in a draft U87 model and included some subsequent simulations. I believe that the problem I ultimately ran into with this iteration was that imbalances in Recon 2 allowed for "free" growth (i.e., production of biomass with *no* inputs). After putting the project on hold for several months after my defense, I made another push with version `1.5` described below. The more recent iteration, however, also ended due to frustration with the reference human model (Recon 2).

### Collect/prepare evidence

+ **Collect HPA protein expression evidence**
  + Code: `getHPA.py`
  + Input: None (HPA website) 
  + Output: **U-87_MG.txt**
+ **Format HPA data**
  + Code: `HPAensembl2gene.py`
  + Input: **U-87_MG.txt** 
  + Output: **U-87\_MG\_genes.txt**
+ **Collect GEO gene expression evidence**
  + Code: `U87GEODataProcessing.m`
  + Input: CEL files from GEO
  + Output: **all_datasets.mat**
+ **Compile/format GEO expression evidence**
  + Code: `U87GEODataProcessing.m`
  + Input: **all_datasets.mat**  
  + Output: **U87\_GEO\_evidence.mat**
+ **Format U87 biomass equation**
  + Code: `write_biomass.m`
  + Input:  **U87 biomass.xlsx**
  + Output: **U87_biomass.txt**

### Format/prepare reference model

+ **Import Recon2 model**  
  + Code: `convertSBMLToCobra_Recon2.m`
  + Input:  **recon2_model.xml**
  + Output: **HR2\_CbModel\_Dec2012\_wEnzymes.mat**
+ **Format Recon2 model**
  + Code: `removeReconGeneMets.m`
  + Input: **HR2\_CbModel\_Dec2012\_wEnzymes.mat**
  + Output: **HR2\_CbModel\_Dec2012.mat**
+ **Augment Recon2 model**
  + Code: `primeRecon2forU87mCADRE.m`
  + Input: **HR2\_CbModel\_Dec2012.mat** 
  + Output: **Recon201.mat**, **Recon201_U87.mat**

### Run mCADRE & collect results

+ **Create mCADRE inputs**
  + Code: `createU87mCADREInputs.m`
  + Input:  
     + **U-87\_MG\_genes.txt**
     + **U87\_GEO\_evidence.mat**
     + **U87_fpkm**
     + **symbols2id**
     + **Recon201_U87.mat**
  + Output: **U87mCADREInputs.mat**
+ **Run mCADRE**
+ **Compile mCADRE results**
  + Code: `compileU87mCADRE.m` 
  + Input: **U87mCADREResults.mat** 
  + Output: **U87_CbModel.mat**, **modelStats.xlsx**

### Model simulations

+ **Simulate growth on DMEM**
  + Code: `U87DMEMGrowth.m` 
  + Input:  **U87_CbModel.mat**
  + Output: NA
+ **Perform knockout simulation / inhibitor screen**
  + Code: `U87InhibitorScreen.m` 
  + Input: **U87_CbModel.mat**, **u87inhibitors.mat**
  + Output: NA

## Iteration `0.5`

I'm including this earlier version, just because it seemed like an interesting strategy at the time. The only difference between this and `1.0` was that I added (supposedly) astrocyte-specifc reactions from Nathan Lewis' model, prior to any U87-related preprocessing. I think this was a hold-over from some earlier attempts with Human Recon 1, and I might have eventually realized/decided that the astrocyte reactions/metabolites were already accounted for in Recon 2.

### Collect/prepare evidence

Same as `1.0`

### Format/prepare reference model

+ **Import Recon2 model** (same as `1.0`)  
+ **Format Recon2 model** (same as `1.0`)   
+ **Add astrocyte model info**
  + Code: `HR2_addAstrocyte.m`
  + Input: **HR2\_CbModel\_Dec2012.mat**, **GLU_norm.mat**
  + Output: **HR2plusA_CbModel.mat**, **HR2plusA_info.mat**
+ **Augment Recon2 model**
  + Code: `U87_mCADRE_preprocessing.m`
  + Input: **HR2plusA_CbModel.mat** 
  + Output: NA

### Run mCADRE & collect results

NA

### Model simulations

NA

## Iteration `1.5`

Probably the most notable change in this version was the more systematic and better documented approach to collecting gene expression evidence from GEO. I also spent a lot of time trying to clean up and resolve issues/discrepancies with a slightly newer versions of Recon 2 (v2.02); some of these changes were included in the updated version 2.03, as you can see in the acknowledgement on the Recon 2 website. Recon 2.03 was still troublesome, and at one point I even tried to merge/reconcile information between 2.03 and a version labelled 2.1 (published by Keiran Smallbone), but this proved too complicated. I won't get into the Recon 2 work here, but I can try to provide more information, if needed. These efforts should hopefully be unneeded for more recent Recon 2 versions.

### Collect/prepare evidence

+ **Collect HPA protein expression evidence** (same as `1.0`)  
+ **Format HPA data** (same as `1.0`)  
+ **Collect GEO gene expression evidence**
  + Code: `importU87AffyData.R`
  + Input: **sampleList.txt** 
  + Output: **GSE[x].txt**, **GSE[x]_map.txt**
+ **Compile/format GEO expression evidence**
  + Code: `compileGEOData.m`
  + Input: **GSE[x].txt**, **GSE[x]_map.txt** 
  + Output: **U87\_GEO\_Jan2014.mat**
+ **Format U87 biomass equation** (same as `1.0`)  

### Format/prepare reference model

+ **Clean up Recon2 model**
  + Code: `Recon203_updater.m`
  + Input: **Recon2.v03.mat**, **1317714181111037_add2 2.tsv**
  + Output: **HR2v03\_CbModel\_Jan2014.mat** 
+ **Augment Recon2 model**
  + Code: `primeRecon203forU87mCADRE.m`
  + Input: **HR2v03\_CbModel\_Jan2014.mat** 
  + Output: NA

### Run mCADRE & collect results

NA

### Model simulations

NA