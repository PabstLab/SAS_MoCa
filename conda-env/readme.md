# Use and installation

The package works without installation. Just download the package folder and run the script *sasmoca.py* as described in the documentation.<br>

To run the script it is necessary to have Python installed (ideally 3.12 version), along with the modules:
  - pyyaml
  - scipy
  - numpy
  - tqdm
  - pandas
  - python
  - matplotlib
  - multiprocess (tested on version 0.70)

#### Recommended: use conda environment!

As an alternative it is recommended to run **SAS_MoCa** package within a **conda environment**.
This will ensure to have all the right Python installation and module dependencies in a separated environemnt.

0. Download and extract **SAS-MoCa** package in a folder of preference (ideally in the same drive where Python will be installed).
1. Install Miniconda (or Anaconda) correctly by following the instructions at https://conda.io/projects/conda/en/latest/user-guide/install/index.html.
2. Go to a command-line prompt (using "Anaconda Prompt (miniconda3)” or “Anaconda Powershell Prompt (miniconda3)" on Windows) and create the, e.g., **SAS-MoCa** environment by typing:

```
conda env create -f <path-to>/environment.yml
```
3. Before running *sasmoca.py* script activate the conda environment:
```
conda activate SAS-MoCa
```

For more information about managing conda environments go to https://conda.io/projects/conda/en/latest/user-guide/tasks/index.html.
