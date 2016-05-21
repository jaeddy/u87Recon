# Workflows

Data processing steps for model generation, with inputs and outputs for each script used, organized by iteration/version.

## Iteration `1.0`

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
  + Output: **HR2v03_CbModel_Jan2014.mat** 
+ **Augment Recon2 model**
  + Code: `primeRecon203forU87mCADRE.m`
  + Input: **HR2v03_CbModel_Jan2014.mat** 
  + Output: NA

### Run mCADRE & collect results

NA

### Model simulations

NA