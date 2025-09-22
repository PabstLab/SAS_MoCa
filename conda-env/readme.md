# Use and installation

The package works without installation. Just download the package folder and run the script *sasmoca.py* as described in the documentation.<br>

In order to run the script, you will need to have Python installed (ideally version 3.12), along with the following modules:
  - pyyaml
  - scipy
  - numpy
  - tqdm
  - pandas
  - python
  - matplotlib
  - multiprocess (tested on version 0.70)

#### Recommended: Use conda environment!

As an alternative, it is recommended to run the  **SAS_MoCa** package within a **Conda environment**.
This ensures that the correct Python installation and module dependencies are in place in a separate environment.

0. Download and extract the **SAS-MoCa** package into a folder of your choice. Ideally, this should be on the same drive where Python is installed.
1. Install Miniconda (or Anaconda) correctly by _strictly_ following the instructions at https://conda.io/projects/conda/en/latest/user-guide/install/index.html.
2. Open a command-line prompt. On Windows, open either "Anaconda Prompt (Miniconda3)” or “Anaconda Powershell Prompt (Miniconda3)" and create the, e.g., **SAS-MoCa** environment by typing:

```
conda env create -f <path-to>/environment.yml
```
3. Before running *sasmoca.py* script activate the conda environment:
```
conda activate SAS-MoCa
```

For more information about managing conda environments go to https://conda.io/projects/conda/en/latest/user-guide/tasks/index.html.
