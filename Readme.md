# Boundary stacking interactions enable cross-TAD enhancer-promoter communication during limb development

[![DOI](https://zenodo.org/badge/666579572.svg)](https://zenodo.org/badge/latestdoi/666579572)

Tzu-Chiao Hung<sup>a</sup>, David M. Kingsley<sup>a, b</sup>, Alistair N. Boettiger<sup>a</sup>, 

<sup>a</sup> Department of Developmental Biology, Stanford University School of Medicine, Stanford, CA 94305, USA

<sup>b</sup> Howard Hughes Medical Institute, Stanford University School of Medicine, Stanford, CA 94305, USA

Correspondence: boettiger@stanford.edu 


This repository contains original code and model parameters to used in the analyses and simulations described in our work investigating chromatin folding at the Pitx1 locus in mammalian limb development using Optical Reconstruction of Chromatin Architecture.

The intention of this repository is to make the transparent the workflow behind the data analysis and behind the simulations.  

The chromatin tracing image data is deposited at the [4DN data portal](https://data.4dnucleome.org/) and can be retrieved there, accession number: 4DNES4TC13IL.  

The "Simulations" folder contains python scripts to run both the Loop Extrusion based polymer simulations and the micro-phase-separation based simulations.  To run the polymer simulations, this code requires the open-access [polychrom repository](https://github.com/open2c/polychrom) developed by the open2c team.  We added several additional functions to the repository to allow non-uniform loading of Cohesin, though this ability was not ultimately used in the simulations described in the current manuscript, the enclosed code does require the expanded repository to be installed to run correctly.  These additions can be found in our open-access fork of the project, [https://github.com/BoettigerLab/polychrom](https://github.com/BoettigerLab/polychrom).  Simulation data used in our study can be found at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8148723.svg)](https://doi.org/10.5281/zenodo.8148723). Running the code described here will produce similar results, non-identical due to the stochastic nature of the simulations.  

Also included are analysis scripts to reproduce the figures presented in our work.  These scripts depend on several functions, included here and in the public release of our analysis toolbox: [https://github.com/BoettigerLab/ORCA-public](https://github.com/BoettigerLab/ORCA-public)
A few supplemental functions not included in ORCA-public are provided in the "Functions" folder within the simulation folder. The scripts are labeled by figure.  Most of the figure scripts also produce associated images that are present in the Extended Data related to the figure.  A subset of Extended Data figures are not included in the Main Figures and are presented instead as stand-alone scripts, labeled by their ED figure number. 

## Project Abstract

While long-range enhancers and their target promoters are frequently contained within a TAD, many developmentally important genes have their promoter and enhancers within different TADs. Hypotheses about molecular mechanisms enabling such cross-TAD interactions remain to be assessed. To test these hypotheses, we used Optical Reconstruction of Chromatin Architecture (ORCA) to characterize the conformations of the Pitx1 locus on thousands of single chromosomes in developing mouse limbs. Our data supports a model in which neighboring boundaries are stacked with each other as a result of loop-extrusion, bringing boundary-proximal cis-elements into contact. This stacking interaction also contributes to the appearance of architectural stripes in the population average maps (e.g. Hi-C data). Through molecular dynamics simulations, we further propose that increasing boundary strengths facilitates the formation of the stacked boundary conformation, counter-intuitively facilitating border bypass. This work provides a revised view of the TAD borders’ function, both facilitating as well as preventing cis-regulatory interactions, and introduces a framework to distinguish border-crossing from border-respecting enhancer-promoter pairs. 