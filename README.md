# u87Recon

This repository contains code, data, and results from efforts to build a COBRA model for the **U-87 MG glioblastoma (GBM)** cell line.

## Background

I worked on this project off and on during my time as a PhD student and post-doc in the lab of Nathan Price. Some of my earliest efforts on the project date as far back as 2009 and came to an end in early 2014. Over this time, a variety of strategies, data, tools, and reference models were used - more so than I could reasonably hope to share or describe in any sensible way. None of these efforts resulted in a final, validated, and published model of U87; however, some pieces of the workflow could potentially be of use to others.

## Workflow

The most complete iteration of the model building workflow started with Human Recon 2, used **mCADRE** (an extension of **MBA** that I co-developed with Yuliang Wang in the Price lab) to produce a U87-specific draft model, and even included some downstream simulation for testing and validation. This workflow is summarized [**here**](docs/workflow.md) and is labeled as iteration `1.0`. Two other closely related iterations of the workflow (labeled as `0.5` and `1.5`) are also described, though neither were carried as far towards completion.

For more details/text related to the workflow used in iteration `1.0`, you can check out **Chapter 6** of my [**PhD thesis**](docs/Eddy_James_Thesis.pdf), included in the **`docs`** folder.

For an overview of my more recent progress and line of thinking towards the end of the project, I've also included a [**poster**](docs/COBRAPoster2014.pdf) that I presented at the 2014 COBRA Conference in Virginia.

## Code

The bulk of the model building workflow was done in **Matlab**, using either COBRA toolbox functions or custom scripts. Some **Python** and **R** scripts were also used to collect and format external data, which was used as evidence for mCADRE (or **MBA**, in earlier versions). All code in the repository is described [**here**](docs/code_manifest.md) and is grouped into folders by language.

I have not attempted to directly clean up, test, or update any of the code here. While I would prefer to share a complete, well-documented pipeline, where one could run (and understand) all code in order to reproduce the data and results, this probably isn't going to happen for a few reasons:  

1. I haven't touched any of this code in over 2 years, making it hard even for me to always figure out exactly what's going on; not to mention that past me wasn't as organized with code and data as I might have liked. Altogether, the idea of updating and cleaning up all the code is incredibly daunting and not very appealing.
2. I haven't regularly coded in Matlab for almost 2 years, exacerbating the effects of (1).
3. I don't really think anyone should try to use my code or reproduce my results at this point. Nothing (besides a section of my thesis) was published from this work, and the project wasn't really successful in terms of producing a validated U87-specific model. I also think many aspects could be done differently or better (e.g., there are newer versions of Human Recon and likely better alternatives to mCADRE at this point). 

Instead, I'm more hopeful that some of the data, evidence, and resources I identified and collected over the course of the project might prove useful to someone working on their own U87/GBM metabolic model reconstruction. Some of this information is contained directly in the files of this repository, while in other cases, the logic or even just the targets of certain scripts might be more useful at this point.

## Data

As with the code, I've done my best to describe all of the data files included in this repo [**here**](docs/data_manifest.md). The data can be broken into the following categories:

+ **Evidence:** U87-specifc gene and protein measurements used as inputs for model building with mCADRE - e.g., microarray gene expression data from GEO, protein expression data from HPA, and RNA-seq gene expression data from in-house experiments in the Price lab.
+ **Simulation:** data used to constrain, parameterize, or validate properties and simulation results of the U87 model - e.g., biomass composition, growth media definition, growth curves and results an of enzyme inhibitor screen.
+ **mCADRE/Models:** COBRA models and related data for Human Recon 2 and U87.
+ **Reference:** look-up tables and variables used for converting between different gene/ID nomenclatures (probably less useful these days with the availability of tools like `biomaRt`).

The relationship between data files and code is summarized in the documentation for [**workflow iterations**](docs/workflow.md).

## Conclusion

Putting together the material in this repo is my attempt to effectively wrap up and hand off my work on this project. By doing so, I hope that I can at least provide some useful information or ideas to help with their own project. I've moved on from academia, the Price lab, and metabolic modeling research, so I don't currently have the time (or motivation) to contribute further to the code or analysis.

I am happy to answer questions (or do my best) about anything contained here.

If you are interested in using any of the code, data, or ideas here, I only ask that you give me a heads up (and maybe an acknowledgement somewhere). I'm not concerned about authorship, but if this work is used as part of a publication, you should contact [**Nathan Price**](https://price.systemsbiology.org/bio/nathan-price/).