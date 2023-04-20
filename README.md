# adata_and_scanpy_tools
Mostly scanpy wrappers and tools that use adata objects to stream line single cell analysis 

# Getting started

## Clone this repo

clone this repo into into to the directory where you keep your cloned python packages

```
cd home/ubuntu/projects/github_repos
git clone https://github.com/megadesk/adata_and_scanpy_tools.git
```

## Code to add to the start of a jupyter notebook in order to run the package
``` 
# option 1 use an absolute path to direcotry containing python packages from github
repo_parent_dir='home/ubuntu/projects/github_repos'
# option 2 use a relative path that will work as long as a the direcotry containing python packages is two levels above the directory containing the jupyter notebook
repo_parent_dir='../../'

import sys
if repo_parent_dir not in sys.path:
    sys.path.append(repo_parent_dir)

# import the package
import adata_and_scanpy_tools as adsctl

``` 
## Directory structure

I recommend placing all your github repos (packages and analysis_project_repos) in the same parent directory that way you can use relative paths to access packages from your analysis_project_repos. This way the same notebook will work on your aws enviroment as as well as your local envirmoment ( without changing path variables) 

>github_repos
    >adata_and_scanpy_tools
    >another_package_repo_directory
    >analysis_project_repo_zero
        >analysis_project_A
          >00_notebook
          >01_notebook
          >02_notebook



# adata_and_scanpy_tools package functions

## Getting start will the package functions
### Check the doc strings of the functions for more info on each

`print(functionOrClass.__doc__) #(double underscores)`

`print(adsctl.pp.basic_filitering.__doc__)`

`print(adsctl.pp.annotate_n_view_adata_raw_counts.__doc__)`

`print(adsctl.pp.filter_cells_by_anotated_QC_gene.__doc__)`

`print(adsctl.pp.remove_genes.__doc__)`

`print(adsctl.pp.process2scaledPCA.__doc__)`

`print(adsctl.pp.leiden_clustering.__doc__)`

`print(adsctl.pl.silhouette_score_n_plot.__doc__)`

# Example notebook using the PBMC3k data set (most usefull starting point)

try out adsctl_PBMC3k_test.ipynb in PBMC3k_example_nb folder

# I'll add more to read me later
