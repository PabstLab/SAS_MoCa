# _SAS_MoCa_

## A stochastic analysis tool for studying the structure of lipid membrane systems using small-angle scattering (SAS).

### Currently supported models
Currently supported models of Large unilamellar vesicles (LUVs) and protein-containinn LUVs (pLUVs, or proteoliposomes) are based on the combination of scattering density profile (SDP) $^1$ and separated form factor (SFF) $^2$ scattering models.
Recent updates for bound-water molecules and thickness fluctuaions are included. $^3$ </br>

* **Large unilamellar vesicles (LUVs)**: 
     * **LUV_DMPC** -> _LUVs_: DMPC/DMPG 95:5 mol/mol mixture; suspension in pure water.
     * **LUV_POPC** -> _LUVs_: POPC/POPG 95:5 mol/mol mixture; suspension in pure water.
     * **LUV_POPE** -> _LUVs_: POPE/POPG 90:10 mol/mol mixture; suspension in pure water.
* **Proteoliposomes (pLUVs)**: 
    * **pLUV_DLPC_OmpLA_RecBuf** -> _Hosting LUV_: DLPC/DLPG 95:5 mol/mol mixture; _protein_: Outer membrane phospholipase A (OmpLA) monomer/dimer; suspension in 20 mM TRIS, 2mM EDTA.
    * **pLUV_POPC_OmpLA_RecBuf** -> _Hosting LUVs_: POPC/POPG 95:5 mol/mol mixture; _protein_: Outer membrane phospholipase A (OmpLA) monomer/dimer; suspension in 20 mM TRIS, 2mM EDTA;

 
**Abbreviations**:
* *DLPC* -> 1,2-dilauroyl-sn-glycero-3-phosphatidylcholine
* *DMPC* -> 1,2-dimyristoyl-sn-glycero-3-phosphocholine
* *POPC* ->	1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine 
* *POPE* ->	1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine 
* *DLPG* -> 1,2-dilauroyl-sn-glycero-3-phosphoglycerol
* *DMPG* -> 1,2-dimyristoyl-sn-glycero-3-phosphoglycerol
* *POPG* ->	1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoglycerol 

### Minimization algorithm: 
* Adaptive **Thermodynamic Simulated Annealing** (TSA) $^{4,5}$

## Installation

The package works without installation. We recommend to run _sasmoca_ on a dedicated conda environment. Installation instructions are in the _./conda-env_ folder.

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
   model: pLUV_POPC_OmpLA_RecBuf 
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
> To simplify the the initial use of _SAS_MoCa_ see the _examples_ folder for working examples and templates.

---

### Run the fitting routine

To fit data o preview the outcome of the chosen scattering model, open a command-line terminal (idealyl on Linux Ubuntu), create and/or move to the folder where you want to save the results and type:

```
>python <path-to-SAS_MoCa>/sasmoca/sasmoca.py ./input_parameter-file.yml
```

---

### Output
The fitting results (or preview) are saved in the configured folder. Here the list of saved files:

#### Text files
Each of the following files includes a header containing a timestamp, the version used, and a summary of the configuration used to fit the data.

* **correlations.dat** <br> It contains a symmetric matrix of Pearson's correlation coefficients for the adjustable parameters. _This file is only created if the number of iterations is greater than 10!_.

* **iterations_collection.dat** <br> Collection of the result parameters (means only) for each iteration. Fixed parameters are also included.

* **plot_intensity.dat** <br> It contains the data points used for fitting and the resulting scattering intensity; column 1: $q_i$, column 2:  $I_i^{exp}$, column 3: $\sigma_i$, column 4: $I^{model}(q_i)$.

* **plot_SDP-SLD.dat**  <br> It contains the calculate probability density for each quasi-molecular group of the (p)LUV model (SDP profile), as well as the corresponding SLD contrast profile. The $z$-axis units are $\AA$, and the SLD units are $\AA^{-2}$.

* **results_recap.dat** <br> Summary of the input parameter table (from name to boundaries) that includes two new columns for the results: **mean** and **standard deviation** for each adjustable parameter. Fixed parameters are also included. Additionally, there is the list of calculated values (_marked with *_) such as area per lipid $A_L$, Luzzati thickness $D_B$, phosphate-to-phosphate distance $D_{pp}$ (in the case of SAXS), and estimate of number of water molecules per headgroup $n_W$.

#### Plots

* **Plot.png** <br> Plot of data and fitted scattering curve, it includes a comparison between relative experimental error and relative deviation from the model.

* **Plot_histograms.png** <br> Histograms showing the distribution obtained for each single adjustable parameter (only saved if the number of iterations is higher than 10). The plots include a kernel-density visualization of distributions, marks pointing the mean ands standard deviations, and the Gaussian prior profiles when present.

* **Plot_histogram_X2.png** <br> Histograms showing the distribution obtained for the best $\chi^2$ values. The plot includes a kernel-density visualization of the distribution.

* **plot_SDP-SLD.png**  <br> Calculated SDP and SLD contrast profiles. The plots include visualizations of acyl-chain thickness ($D_C$) and Luzzati length ($D_B$).

## If you use _SAS_MoCa_ repository please cite:

* Semeraro, E. F., & Pabst, G. SAS_MoCa (Version 1.3.0) [Computer software]. https://github.com/PabstLab/SAS_MoCa
* Semeraro E. F., et al., in preparation 
* Semeraro E. F., et al., in preparation 

## License
**SAS_MoCa** is free and open source software, distributed under the [BSD 3‑Clause “New” or “Revised” License](https://opensource.org/license/BSD-3-Clause).<br> 
For the full legal text, see the `LICENSE` file in this repository.

## References

1. Kučerka, N., Nagle, J. F., Sachs, J. N., Feller, S. E., Pencer, J., Jackson, A., & Katsaras, J. (2008). Lipid Bilayer Structure Determined by the Simultaneous Analysis of Neutron and X-Ray Scattering Data. Biophysical Journal, 95(5), 2356–2367. https://doi.org/10.1529/biophysj.108.132662
2. Pencer, J., Krueger, S., Adams, C. P., & Katsaras, J. (2006). Method of separated form factors for polydisperse vesicles. Journal of Applied Crystallography, 39(3), 293–303. https://doi.org/10.1107/S0021889806005255
3. Frewein, M. P. K., Doktorova, M., Heberle, F. A., Scott, H. L., Semeraro, E. F., Porcar, L., & Pabst, G. (2021). Structure and Interdigitation of Chain-Asymmetric Phosphatidylcholines and Milk Sphingomyelin in the Fluid Phase. Symmetry, 13(8), 1441. https://doi.org/10.3390/sym13081441
4. de Vicente, J., Lanchares, J., & Hermida, R. (2003). Placement by thermodynamic simulated annealing. Physics Letters A, 317(5–6), 415–423. https://doi.org/10.1016/j.physleta.2003.08.070
5. Semeraro et al., in submission 
