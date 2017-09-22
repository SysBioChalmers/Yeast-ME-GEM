# Yeast-ME

**We updated the metabolic model by including protein biosynthesis and degradation reactions and coupled metabolic fluxes with protein abundances.**

- Abstract:

Flux Balance Analysis diverts the glucose uptake rate into many pathways to illustrate how a cell builds its biomass. FBA does not consider the enzymatic cost in each pathway, so all predicted flux distributions are only in respiration state. In this computational framework (Yeast-ME), we estimated the cost of protein biosynthesis and degradation. We included also in silico ribosome and proteasome. Finally, for each peptide added three reactions for translation, degradation, and dilutions. Using the new reactions, the model was to estimate a protein abundance, and the metabolic fluxes were coupled with the protein abundances. By constrained protein abundance per cell and other 22 constraints, the model becomes able to predict the metabolic shifts. 

- Model KeyWords:Yeast metabolic models, genome-scale metabolic models, metabolic shiftm metabolism overflow.

**GEM Category:** Species/Community/Collection; **Utilisation:** maximising product growth, minimising product growth, predictive simulation, experimental data reconstruction; **Field:** metabolic engineering, bacterial community, thermodynamic modeling, metabolic-network reconstruction; **Type of Model:** Init, curated, reconstruction, dynamic-kinetic simulation; **Model Source:** HPA, HMR, iHepatocytes2322, YeastMetabolicNetwork, [kbase](https://kbase.us/); **Omic Source:** Transcriptomics, Proteomics, Metabolomics; **Enzymatically Constrained** ; **Taxonomy:** Homo sapiens, Saccharomyces cerevisiae, Bacteria; **Metabolic System:** Mitochondrial metabolism, Fatty Acid Metabolism, Glucose Metabolism, General Metabolism; **Tissue:** NameOfTissue; **Bioreactor** ; **Cell Type:** CellType; **Cell Line:** Cancer_NumberOfLine; **Strain:** str.xxxx; **Condition:** Cancer, Malnourishment, Starvation, Oxidative Stress, Rich Media;




## Installation

### Required Software:

* *_PROGRAMMING LANGUAGE/Version_*  (e.g.):
  *  You need a functional Matlab installation of **Matlab_R_2015_b**  (MATLAB 7.3 and higher)
  * The [RAVEN](https://github.com/SysBioChalmers/RAVEN) toolbox for MATLAB. An up-to-date version from COBRA GitHub repository is strongly recommended . Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)

### Dependencies - Recommended Software:
* libSBML MATLAB API (version [5.13.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/)  is recommended).


### Installation Instructions
* Clone [model](https://github.com/SysBioChalmers/) branch from [SysBioChalmers GitHub](https://github.com/SysBioChalmers)
* Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)


## Contributors
- [Contributor1](https://www.chalmers.se/en/Staff/Pages/.aspx), Chalmers University of Technology, Gothenburg Sweden
- [Contributor2](https://www.chalmers.se/en/staff/Pages/.aspx), Chalmers University of Technology, Gothenburg Sweden

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
