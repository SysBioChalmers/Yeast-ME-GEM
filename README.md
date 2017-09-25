# Yeast-ME

**We updated the metabolic model by including protein biosynthesis and degradation reactions and coupled metabolic fluxes with protein abundances.**

- Abstract:

Genome-scale metabolic models (GEMs) elucidate how a cell consumes nutrients to build its biomass. The main challenge of the current GEMs is to predict the metabolic shifts or  metabolism overflow when the cell reshapes its metabolism from the respiration phase to the respiro-fermentative phase. We are reconstructing the next-generation yeast GEM (the first model for a eukaryotic cell) to overcome this limitation by using further information on how energy is consumed in the cell, such as protein translation and degradation reactions that take up approximately of 80% of total cell energy consumption. Additional intracellular constraints have been also included, such as cell volume, that can be coupled to protein abundances in the model. At high growth rates, the cell needs more amount of ribosomes and metabolic enzymes, so the cell proteome can reach to one of space capacity constraints, such as the total cytosolic protein, mitochondrial membrane protein, or any other protein constraint. Therefore, the model can predict a metabolic shift, because the protein cost in respiration pathway is more expensive. Model predictions under different conditions such as nitrogen limitation (reducing in protein content) and high temperatures (increasing in misfolded proteins) can decide the best constraint. Integration model predictions with absolute quantitative proteomics (~2000 protein abundances) and experimental data can subsequently guide the model towards predicting when a metabolic shift may occur and validate this constraint causing the metabolic shift. This model can be used to interpret the yeast adaptive responses as well as to determine the optimal metabolic engineering strategies for metabolite or protein productions. 

- Model KeyWords:Yeast metabolic models, genome-scale metabolic models, metabolic shift, metabolism overflow.

**GEM Category:** Species; **Utilisation:** maximising product growth, minimising product growth, predictive simulation, experimental data reconstruction; **Field:** metabolic engineering, thermodynamic modeling, metabolic-network reconstruction; **Type of Model:** Init, curated, reconstruction, dynamic-kinetic simulation; **Model Source:** HPA, HMR, iHepatocytes2322, YeastMetabolicNetwork; **Omic Source:** Transcriptomics, Proteomics, Metabolomics; **Enzymatically Constrained** ; **Taxonomy:** Saccharomyces cerevisiae; **Metabolic System:** Mitochondrial metabolism, Fatty Acid Metabolism, Glucose Metabolism, General Metabolism; **Strain:** str.xxxx; **Condition:** Glucose-limit, Rich Media;




## Installation

### Required Software:

* *_PROGRAMMING LANGUAGE/Version_*  (e.g.):
  *  You need a functional Matlab installation of **Matlab_R_2015_b**  (MATLAB 7.3 and higher)
  * The [RAVEN](https://github.com/SysBioChalmers/RAVEN) toolbox for MATLAB. An up-to-date version from COBRA GitHub repository is strongly recommended . Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)

### Dependencies - Recommended Software:
* SoPlex (http://soplex.zib.de/).


### Installation Instructions
* Clone [model](https://github.com/SysBioChalmers/) branch from [SysBioChalmers GitHub](https://github.com/SysBioChalmers)
* Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)


## Contributors
- Ibrahim Elsemman (http://sysbio.se/people.html), Technical University of Denmark, Kgs.Lyngby, Denamrk

## License
The MIT License (MIT)

> Copyright (c) 2017 Systems and Synthetic Biology
>
> Chalmers University of Technology Gothenburg, Sweden
>
>Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
>
>The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
>
>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
