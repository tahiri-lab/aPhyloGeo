﻿﻿﻿﻿﻿﻿﻿<h1  align="center"> aPhyloGeo <p align='center'>
        [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
        [![Contributions](https://img.shields.io/badge/contributions-welcome-blue.svg)](https://pysd.readthedocs.io/en/latest/development/development_index.html)
        [![Py version](https://img.shields.io/pypi/pyversions/pysd.svg)](https://pypi.python.org/pypi/pysd/)
        [![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Ftahiri-lab%2FaPhylogeo&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
        [![GitHub release](https://img.shields.io/github/v/release/tahiri-lab/aPhylogeo.svg?maxAge=3600)](https://github.com/tahiri-lab/aPhylogeo/releases/)
        </p>


<h2  align="center"> 🌳 Multi-platform application for analyze phylogenetic trees with climatic parameters</h2>

<details open>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#-about-the-project">About the project</a>
    </li>
    <li>
      <a href="#%EF%B8%8F-installation">Installation</a>
      <ul>
        <li><a href="#linux-unix-mac-os--windows-versions">Linux/UNIX, Mac OS and Windows versions</a></li>
      </ul>
    </li>
     <li>
      <a href="#-settings">Settings</a>
    </li>
    <li>
      <a href=#%EF%B8%8F-references">References</a>
    </li>
    <li>
      <a href="#-contact">Contact</a>
    </li>
  </ol>
</details>


# 📝 About the project

`aPhyloGeo` is a powerful and versatile bioinformatics pipeline specifically designed for the analysis of phylogeography. Developed by the dedicated team led by Professor [Nadia Tahiri](https://tahirinadia.github.io/) at the University of Sherbrooke in Quebec, Canada, this open-source multi-platform application is implemented in Python. It serves as a valuable tool for researchers and scientists interested in unraveling the intricate relationships between genetic sequences and geographic locations.

💡 If you are using our algorithm in your research, please cite our recent papers:

1️⃣  Li, W. & Tahiri, N. (2023). aPhyloGeo-Covid: A Web Interface for Reproducible Phylogeographic Analysis of SARS-CoV-2 Variation using Neo4j and Snakemake.
[Proceeding in SciPy 2023, Auxtin, TX, USA](https://conference.scipy.org/proceedings/scipy2023/pdfs/nadia_tahiri.pdf)

2️⃣ Koshkarov, A., Li, W., Luu, M. L., & Tahiri, N. (2022). Phylogeography: Analysis of genetic and climatic data of SARS-CoV-2.
[Proceeding in SciPy 2022, Auxtin, TX, USA](https://conference.scipy.org/proceedings/scipy2022/pdfs/nadia_tahiri.pdf)

# ⚒️ Installation

## Linux UNIX & Windows versions 
`aPhyloGeo` is available as a Python script.

### Prerequisites
This package use ```Poetry``` dependency management and packaging tool for Python. Poetry installation guide can be found here: [Poetry Install](https://python-poetry.org/docs/#installation)
⚠️ For windows installation it's recommended to launch powershell in **Administrator mode**.

Once Poetry is installed, you can clone the repository and install the package using the following commands:

```
poetry install
```

### Usage
Poetry will handle the virtual environment automatically. if you want to use the virtual environment manually, you can use the following command:

```
poetry shell
```

⚠️ Assuming Python 3.8 or higher is installed on the machine, these scripts should run well with the libraries installed.

You can also launch the package using the `make` command from your terminal when you are in the `root`. This command will use the `Makefile` to run the script. If you use the command `make clean`, it will erase the `output.csv` file previously created with the first command.

# 🚀 Settings
The `aPhyloGeo` software can be encapsulated in other applications and applied to other data by providing a YAML file. This file will include a set of parameters for easy handling (see [Wiki documentation](https://github.com/tahiri-lab/aPhyloGeo/wiki/User-Preference)).


# ✔️ References

1️⃣ Calculation of distance between phylogenetic tree: `Least Square metric`
+ [Cavalli-Sforza, L. L., & Edwards, A. W. (1967). Phylogenetic analysis. Models and estimation procedures. American journal of human genetics, 19(3 Pt 1), 233.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1706274/)
+ [Felsenstein, J. (1997). An alternating least squares approach to inferring phylogenies from pairwise distances. Systematic biology, 46(1), 101-111.](https://pubmed.ncbi.nlm.nih.gov/11975348/)
+ [Makarenkov, V., & Lapointe, F. J. (2004). A weighted least-squares approach for inferring phylogenies from incomplete distance matrices. Bioinformatics, 20(13), 2113-2121.](https://pubmed.ncbi.nlm.nih.gov/15059836/)

2️⃣ Calculation of distance between phylogenetic tree: `Robinson-Foulds metric`
+ [Robinson, D.F. and Foulds, L.R., 1981. Comparison of phylogenetic trees. Mathematical biosciences, 53(1-2), pp.131-147.](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub)

3️⃣ Dataset full description: `Analysis of genetic and climatic data of SARS-CoV-2`
+ [Koshkarov, A., Li, W., Luu, M. L., & Tahiri, N. (2022). Phylogeography: Analysis of genetic and climatic data of SARS-CoV-2.](https://conference.scipy.org/proceedings/scipy2022/nadia_tahiri.html)
+ [Li, W., & Tahiri, N. (2023). aPhyloGeo-Covid: A Web Interface for Reproducible Phylogeographic Analysis of SARS-CoV-2 Variation using Neo4j and Snakemake.](https://conference.scipy.org/proceedings/scipy2023/nadia_tahiri.html)

# 📧 Contact
Please email us at: <Nadia.Tahiri@USherbrooke.ca> for any questions or feedback.
