---
layout: default
title: Getting Started
nav_order: 1
---

# Getting Started

🚧 __*Under construction*__

Content coming soon.

---

## Installation and prerequisites

**SAS MoCa** is a Python package that works without installation, that is it runs directly from source (no ``` pip install``` required). It is sufficient to download the [last release](https://github.com/PabstLab/SAS_MoCa/releases/latest) and be sure to have the right Python installation along with the following modules:

- python (ideally version >3.12)
- pyyaml
- scipy
- numpy
- tqdm
- pandas
- matplotlib
- multiprocess (tested on version 0.70)

### Recommended: Use conda environment!

We strongly recommend running **SAS_MoCa** within a [**Conda environment**](https://docs.conda.io/en/latest/) to ensure dependency management and reproducibility.
Alternatives like [_Mamba_](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) also work.
In the following, we give a guide in the case of _Conda_.

#### Step 1: Download the Package

Download the latest release from [GitHub Releases](https://github.com/PabstLab/SAS_MoCa/releases/latest) and extract it to your desired folder.

#### Step 2: Set Up the Conda Environment (_recommended_)

First, install Miniconda (or Anaconda) by _**strictly**_ following the instructions at [https://conda.io/projects/conda/en/latest/user-guide/install/index.html](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
Then, open a command-line terminal:
- **Linux/macOS**: Terminal with conda initialized
- **Windows**: "Anaconda Prompt (Miniconda3)" or "Anaconda Powershell Prompt (Miniconda3)"

Navigate to the _SAS\_MoCa_ folder and create the environment:
```
cd <path-to-SAS_MoCa>
conda env create -f conda-env/environment.yml 
```
💡 If the environment creation gets stuck at "solving environment", try instead :
```
conda env create -f conda-env/environment.yml --solver libmamba
```

#### Step 3: Activate the Environment

```
conda activate SAS-MoCa
```

You're now ready to run SAS_MoCa!

## Quick start and verification
Show users how to run something immediately after install

## Input file basics
Link template folder

## Next steps
Point to User Guide, Examples folder




