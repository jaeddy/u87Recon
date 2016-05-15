# Workflows

Data processing steps for model generation, with inputs and outputs for each script used, organized by iteration/version.

## Iteration 1.0

### Collect/prepare evidence

+ **Collect HPA protein expression evidence**
  + Code: `getHPA.py`
  + Input:  
  + Output: 
+ **Format HPA data**
  + Code: `HPAensembl2gene.py`
  + Input:  
  + Output: 
+ **Collect GEO gene expression evidence**
  + Code: `U87GEODataProcessing.m`
  + Input:  
  + Output: 
+ **Compile/format GEO expression evidence**
  + Code: `U87GEODataProcessing.m`
  + Input:  
  + Output: 
+ **Format U87 biomass equation**
  + Code: `write_biomass.m`
  + Input:  
  + Output: 

### Format/prepare reference model

+ **Import Recon2 model**  
  + Code: `convertSBMLToCobra_Recon2.m`
  + Input:  
  + Output: 
+ **Format Recon2 model**
  + Code: `removeReconGeneMets.m`
  + Input:
  + Output:
+ **Augment Recon2 model**
  + Code: `primeRecon2forU87mCADRE.m`
  + Input:  
  + Output: 

### Run mCADRE & collect results

+ **Create mCADRE inputs**
  + Code: `createU87mCADREInputs.m`
  + Input:  
  + Output: 
+ **Run mCADRE**
+ **Compile mCADRE results**
  + Code: `compileU87mCADRE.m` 
  + Input:  
  + Output: 

### Model simulations

+ **Simulate growth on DMEM**
  + Code: `U87DMEMGrowth.m` 
  + Input:  
  + Output: 
+ **Perform knockout simulation / inhibitor screen**
  + Code: `U87InhibitorScreen.m` 
  + Input:  
  + Output: 

## Iteration 0.5

### Collect/prepare evidence

Same as `1.0`

### Format/prepare reference model

+ **Import Recon2 model** (same as `1.0`)  
+ **Format Recon2 model** (same as `1.0`)   
+ **Add astrocyte model info**
  + Code: `HR2_addAstrocyte.m`
  + Input:
  + Output:
+ **Augment Recon2 model**
  + Code: `U87_mCADRE_preprocessing.m`
  + Input:  
  + Output: 

### Run mCADRE & collect results

NA

### Model simulations

NA

## Iteration 1.5

### Collect/prepare evidence

+ **Collect HPA protein expression evidence** (same as `1.0`)  
+ **Format HPA data** (same as `1.0`)  
+ **Collect GEO gene expression evidence**
  + Code: `importU87AffyData.R`
  + Input:  
  + Output: 
+ **Compile/format GEO expression evidence**
  + Code: `compileGEOData.m`
  + Input:  
  + Output: 
+ **Format U87 biomass equation** (same as `1.0`)  

### Format/prepare reference model

+ **Import Recon2 model**  
  + Code: `convertSBMLToCobra_Recon2.m`
  + Input:  
  + Output: 
+ **Format Recon2 model**
  + Code: `removeReconGeneMets.m`
  + Input:
  + Output:
+ **Clean up Recon2 model**
  + Code: `Recon203_updater.m`
  + Input:
  + Output: 
+ **Augment Recon2 model**
  + Code: `primeRecon203forU87mCADRE.m`
  + Input:  
  + Output: 

### Run mCADRE & collect results

NA

### Model simulations

NA