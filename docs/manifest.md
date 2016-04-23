
## Code

+ **`textArrayRead.m`:** utility/wrapper function I wrote years ago to make it somewhat easier to read in text files (Matlab probably has something better than this by now).

+ **`textArrayWrite.m`:** see above for `textArrayRead`.

+ **`Recon203_updater`:** used to clean up / fill in some missing information in the Recon2.03 model, based on data from a separate publication.

+ **`primeRecon202forU87mCADRE.m` / `primeRecon203forU87mCADRE.m` / `primeRecon2forU87mCADRE.m`:** These scripts all seem to share the same goal - modifying/augmenting/testing the generic Recon2 model prior to running mCADRE to generate a U87-specific model. `primeRecon2` appears to be the most recent version, but might want to refer to `primeRecon203` to transfer some comments/description/logic (`primeRecon202` can probably be discarded).
    + Note: from what I can tell, these scripts never evolved past the point of being 'interactive' (a mix of automated processing steps and some debugging/exploratory steps), so I'm not sure if they were ever used to produce a saved output.
    + Also: the scripts appear to use some functions from the Warburg modeling paper; I'll need to include these (and probably check permissions). However, it's possible I rewrote these functions myself at some point...

+ **`U87InhibitorScreen.m`:** This script uses enzyme inhibitor results from the Hockenbery lab to perform a sort of knock-out screen / validation test.

+ **`formulateDMEM.m`:** This appears to be a precursor to the `primeRecon2` scripts above. I believe the difference is that I used this script during an iteration of the project where I tried to include media and test for growth *after* running mCADRE (but I'll have to double-check).

+ **`U87DMEMGrowth.m`:** This looks like it does more or less the same thing as `formulateDMEM`, but is a script rather than a function.

+ **`U87_growth_OLD.m`:** another script for testing growth and other capabilities of the U87-specific metabolic model. I think this uses a model produced from Recon2.01, but it mentions Recon1 in some places.

+ **`compileU87mCADRE.m`:** This script compiles, formats, and summarizes outputs from mCADRE to yield the `U87` model in `U87_CbModel.mat` - appears to be based on Recon2.01.

**Note:** I can't seem to find the script/command I used to actually run mCADRE, which is unfortunate, because it probably means the code is lost on a cluster somewhere.

+ **`createU87mCADREInputs.m`:** This script takes care of the mCADRE pre-processing prior to running mCADRE with Recon1...? Or not: code seems to refer to Recon2.01. Need to sort this out!

+ **`U87GEODataProcessing.m`:** another pre-processing step - compiling and formatting GEO gene expression evidence for use with mCADRE.

+ **`processU87genes.m`:** I can't quite tell what this one is doing. I'm pretty sure it's working with a earlier version of the U87 model (probably based on Recon1..?). This might have been used to format inputs for DIRAC.
