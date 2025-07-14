# SAS_MoCa_Proteo
Multiscale scattering model and analysis tool for large unilamellar proteoliposomes (pLUVs)

## Small-angle scattering analysis for proteoliposomes

Current supported model:
* _Lipids_: POPC and DLPC
* _Protein_: Outer membrane phospholipase A (OmpLA)
 
### Minimization algorithm: Adaptive Thermodynamic Simulated Annealing $^1$

## Installation

The package works without installation. Just download the entire folder and run the scripts *Fit_TSA.py* or *Fit_TSA_mp.py* as described in the documentation.<br>

It is necessary to have Python installed on your computer in order to run the scripts, along with the modules:
- yaml
- numpy
- scipy
- pandas
- matplotlib
- tqdm
- multiprocess

#### Recommended: use conda environment

As an alternative it is recommended to run SAS-MoCa_proteo package within a **conda environment**.
This will ensure to have all the right Python installation and module dependencies in a separated environemnt.

0. Download and extract **SAS-MoCa_proteo** package in a folder of preference (ideally in the same drive where Python will be installed).
1. Install Miniconda (or Anaconda) following the instruction at https://conda.io/projects/conda/en/latest/user-guide/install/index.html.
2. Go to a command-line prompt (using "Anaconda Prompt (miniconda3)” or “Anaconda Powershell Prompt (miniconda3)" on Windows) and create the, e.g., **SAS-MoCa** environment by typing:

```
conda env create -f <path-to>/environment.yml
```
3. Before running *Fit_TSA.py* or *Fit_TSA_mp.py* script activate the conda environment:
```
conda activate SAS-MoCa
```

For more information about managing conda environments go to https://conda.io/projects/conda/en/latest/user-guide/tasks/index.html.

## Documentation

### Setting up the analysis: Input file (*.yml file)

The input file contains all the information needed to configure and initialize data analysis, including the data file path, minimization options, scattering model, and initial parameters.
This is a *.yml file, so any changes must comply with YAML syntax (e.g., indentation, etc.). You can include/remove comments and descriptions (text preceded by #).<br>
This is how it looks like: 

```
config:
   datafile: "/filepath/filename.dat" 
   save-folder: "MySample"
   qrange: [0.005,0.6] 
   model: IPLUV_Stk_SDP_POPCbase_RecBuf 
   temperature-init: 400
   temperature-gain: 50 
   target-X2: 12.0
   data-binning: 10
   error-scale: 1
   iterations: 42
   processes: 24
   plot-only: off 
   state: monomer

parameters:
    parameter_1: [4.2, off, Null, 1.0, 10.]
    parameter_3: [2.4, on, Null, 0.8, 4.8] 
    parameter_2: [420, on, 0.042, Null, Null]
    .
    .
    .
    parameter_N: [420, on, 0.042, 0, 10000]
```

#### Configuration

In the _config_ block there are all the configuration options.

**datafile** <br> Path of scattering data. The datafile should not contain text in header and footer sections, but just numbers is 3 columns: an array for $q_i$, one for the scattering intensity $I_i$, and the last for the associated error $\sigma_i$.

**save-folder** <br> String of characters used to create the folder for saving the output (do not uses empty spaces): *RES_"save-folder"*; the position of the folder is the current working directory.

**qrange** <br> Min. and max. _q_-values to define the range of data to fit.

**model** <br> Scattering model to use (chosen among the available ones in Model_list.py).

**temperature-init** <br> Initial temperature of the Simulated Annealing algorithm; ideally it should be about 100 times bigger than the expected best $\chi^2$.

**temperature-gain** <br> Factor that influences the variation of the simulated annealing temperature, the higher the gain the faster the minimum temperature will be reached. To avoid *crystallization* the minimum temperature is dynamically set to $[1 + \log_{10}(\chi^2_{min})]$.

**data-binning:** <br> Reduce the dataset size by a factor *n*, keeping data point with a "frequency" of *n* and discarding the others.

**error-scale(0..1):** <br> Scale factor for $\sigma_i$; in the case of scaling values in [0:1), $\sigma_i$ values do not become smaller that 2\% of $I_i$.

**iterations** <br> Number of iterations to build up statistics; each iteration consists of one simulated annealing run, but alone does not provide error estimations.

**processes** <br> Number of processes to handle different iterations with parallel computing; if this value is set to 1 serial computation will be used. 

**plot-only** <br> on/off entry: only plotting (on) or data fitting (off).

**state** <br> Assume either 'monomer' or 'dimer' state of OmpLA.

#### Initialize model parameters

In the second block there is the list of the parameters required by a given scattering model.
```
    parameter_1: [4.2, off, Null, 1.0, 10.]
    parameter_3: [2.4, on, Null, 0.8, 4.8] 
    parameter_2: [420, on, 0.042, Null, Null]
    .
    .
    .
    parameter_N: [420, on, 0.042, 0, 10000]
```
* **label:** <br> Name of the parameters, it is just a label for visualization and should be short and without spaces; use the description are after *#* to include a longer description if needed. The order of parameters matters: check the model to verify the order of the parameters. 

* **List of values** <br> It contains, in the order, 1) initialization value, 2) fix/free parameter option, 3) prior information, 4) low hard-boundary, 5) high hard-boundary.

    * **initialization** <br> Starting values to initialize the fitting routine. 

    * **fix/free parameter option** <br> set off to have a fixed value, set on to adjust the parameter during the minimization run.

    * **prior information** <br> Set the relative $\sigma_{prioir}$ value of a Gaussian prior pdf centered at the initialization value. If the value is _Null_, no informative prior is set for the given parameter.

    * **low and high hard-boundaries** <br> Set the lower and higher boundaries accessible to the adjustable parameters. The _Null_ option is only valid if the prior information is not _Null_, as the lower and higher boundaries are automatically set to $\pm5\times\sigma_{prioir}$.

> [!NOTE]
> To simplify the creation of the fisrt input file, the repository contains a template file *Template_parameter-file.yml*.

### Run the fitting routine

To fit data o preview the outcome of the chosen scattering model, open a command-line terminal (idealyl on Linux Ubuntu), create and/or move to the folder where you want to save the results and type:

```
>python <path to SAS-MoCa_proteo>/Fit_TSA_mp.py ./input_parameter-file.yml
```

In this example *input_parameter-file.yml* is located in the same working folder.
<br>

> [!WARNING]
> The script *Fit_TSA_mp.py* requires the python modules **multiprocessing** and **tqdm**. If multiple processing causes problems, serial computation can be used by either setting *processes=1* or running the script *Fit_TSA.py*.<br>
> The current version of *Fit_TSA_mp.py* was not tested on MAC OS.

### Output

The fitting results (or preview) are saved in the configured folder. Here the list of saved files:

* **Plot.png** <br> Plot of data and fitted scattering curve, it includes a comparison between relative experimental error and relative deviation from the model.

* **Plot_histograms.png** <br> Histograms showing the distribution obtained for each single adjustable parameter (only saved if the number of iterations is higher than 10). The plots include a kernel-density visualization of distributions, marks pointing the mean ands standard deviations, and the Gaussian prior profiles when present.

* **Plot_histogram_X2.png** <br> Histograms showing the distribution obtained for the best $\chi^2$ values. The plot includes a kernel-density visualization of the distribution.

* **Results_collection.dat** <br> List of the set of resulting parameter for each iterations. 

* **Results_intensity.dat** <br> File containing the data points used for fitting and the resulting scattering intensity; column 1: $q_i$, column 2:  $I_i^{exp}$, column 3: $\sigma_i$, column 4: $I^{model}(q_i)$.

* **Results_metadata.dat** <br> Summary of the configuration used to fit the data.

* **Results_pearsonr.dat** <br> Symmetric matrix containing the Pearson correlation coefficient for the adjustable parameters  (only saved if the number of iterations is higher than 10).

* **Results_recap.dat** <br> Summary of the input parameter table (from name to boiundaries) that includes two new columns of results: mean and standard deviation for each adjustable parameter.

* **Results_X2_mean.dat** <br> Single $\chi^2$-value retrieved from the calculated mean results.

### Scattering models

The available scattering models are listed in *Model_list.py* and the relative scripts are in the subfolder *models*.

<br>

The curretly tested modules are:
* **PLUV_POPC.py** classes:  
    * *IPLUV_Stk_SDP_POPCbase_RecBuf* : scattering model for LUVs including homogeneously distributed OmpLA monomers or dimers, suspended in buffer consisiting of 20 mM TRIS and 2 mM EDTA. The fuction makes use of a SDP $^{2,3}$ modelling for 95:5 mol/mol POPC/POPG $^{4,5}$ bilayer combined with separated form-factor model $^6$ including polydispersity.
* **PLUV_DLPC.py** classes:  
    * *IPLUV_Stk_SDP_DLPCbase_RecBuf* : scattering model for LUVs including homogeneously distributed OmpLA monomers or dimers, suspended in buffer consisiting of 20 mM TRIS and 2 mM EDTA. The fuction makes use of a SDP $^{2,3}$ modelling for 95:5 mol/mol DLPC/DLPG $^{4,5}$ bilayer combined with separated form-factor model $^6$ including polydispersity.    


## For details about the scattering models and the minimization algorithm check:

Ref


## References


1. de Vicente, J., Lanchares, J., & Hermida, R. (2003). Placement by thermodynamic simulated annealing. Physics Letters A, 317(5–6), 415–423. https://doi.org/10.1016/j.physleta.2003.08.070
2. Kučerka, N., Nagle, J. F., Sachs, J. N., Feller, S. E., Pencer, J., Jackson, A., & Katsaras, J. (2008). Lipid Bilayer Structure Determined by the Simultaneous Analysis of Neutron and X-Ray Scattering Data. Biophysical Journal, 95(5), 2356–2367. https://doi.org/10.1529/biophysj.108.132662
3. Frewein, M. P. K., Doktorova, M., Heberle, F. A., Scott, H. L., Semeraro, E. F., Porcar, L., & Pabst, G. (2021). Structure and Interdigitation of Chain-Asymmetric Phosphatidylcholines and Milk Sphingomyelin in the Fluid Phase. Symmetry, 13(8), 1441. https://doi.org/10.3390/sym13081441
4. Kučerka, N., Nieh, M.-P., & Katsaras, J. (2011). Fluid phase lipid areas and bilayer thicknesses of commonly used phosphatidylcholines as a function of temperature. Biochimica et Biophysica Acta (BBA) - Biomembranes, 1808(11), 2761–2771. https://doi.org/10.1016/j.bbamem.2011.07.022
5. Pan, J., Heberle, F. A., Tristram-Nagle, S., Szymanski, M., Koepfinger, M., Katsaras, J., & Kučerka, N. (2012). Molecular structures of fluid phase phosphatidylglycerol bilayers as determined by small angle neutron and X-ray scattering. Biochimica et Biophysica Acta (BBA) - Biomembranes, 1818(9), 2135–2148. https://doi.org/10.1016/j.bbamem.2012.05.007
6. Pencer, J., Krueger, S., Adams, C. P., & Katsaras, J. (2006). Method of separated form factors for polydisperse vesicles. Journal of Applied Crystallography, 39(3), 293–303. https://doi.org/10.1107/S0021889806005255


**For details about the scattering models and the minimization algorithm check:**

Ref

## If you use _SAS-MoCa_proteo_ repository please cite:

Ref 
