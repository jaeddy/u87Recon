# Data

Organized by subfolder.

## `evidence`

Nominally, these files are ultimately used to produce inputs to mCADRE (or MBA, in previous iterations). I've stored them in a separate folder, as they might be useful in their own right. Files are further grouped into subfolders based on source:

### `hpa`

### `geo`

### `lab`

### `geo_new`


## `simulation`

This folder includes relevant information for defining parameters and constraints in model simulations (whether for testing, pre-processing, model-building, or validation). Files are further grouped into subfolders based on type of information:

### `biomass`

### `media`

### `experiments`


## `models`

This folder contains different versions of both the generic (Human Recon) and U87-specific COBRA models, along with some related metadata. Files are grouped into subfolders based on generic or cell-type specific model:

### `recon2`

### `u87`


## `mcadre`

This folder contains the final inputs for and outputs from mCADRE.


## `reference`

This folder contains some reference data from NCBI, mostly for tasks that could be done much more easily with `biomaRt`...


## other

+ **`Workbook3.xlsx`:** I have no idea what this file is, but I'm going to hang onto it until I'm confident it's not important.
