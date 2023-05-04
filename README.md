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

I recommend placing all your github repos (packages and analysis_project_repos) in the same parent directory that way you can use relative paths to access packages from your analysis_project_repos.\
This way the same notebook will work on your aws enviroment as as well as your local envirmoment ( without changing path variables).\
For example adding, adding '../../' to your sys.path from any of the notebooks in the example direcotry tree below will give you access to both the adata_and_scanpy_tools and the second_favorite_package_repo_directory regardless of where the "github_repos" is located.
> github_repos
>> adata_and_scanpy_tools
>> 
>> second_favorite_package_repo_directory
>> 
>> singlecell_analysis_project_repo
>>> singlecell_analysis_project_A
>>>> 00_A_notebook\
>>>> 01_A_notebook\
>>>> 02_A_notebook
>>>
>>> singlecell_analysis_project_B
>>>> 00_B_notebook\
>>>> 01_B_notebook\
>>>> 02_B_notebook



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


### latest update is a modified version of scanpy's injest funciton
`print(adsctl.tl.ingest_verbose.__doc__)`



# most usefull  (most usefull starting point)
## check out the eample notebook folders
* preprocessing_PBMC3k_example_nbs
   * adsctl_gex_class_preprocessing_PBMC3k.ipynb # single class object for preprocessing
   * adsctl_functions_preprocessing_PBMC3k_.ipynb # multiple functions preprocess an adata object 
      * more memory efficient



# I'll add more to read me later


## Subpackage - adata_and_scanpy_tools - preprocessing - ANsc.pp.

`adsctl.pp.basic_filitering()`
```
basic_filitering(adata,
                     filter_cells_min_counts=0,
                      filter_cells_min_genes=200,
                     filter_genes_min_cells=3,
                     filter_genes_min_counts=0,
                     **parameters):
```


`adsctl.pp.annotate_QC_genes(adata)`

```
annotate the group of QC genes
### double HB gene annoation works.... maybe just give it a list of genes
#### code 
adata.var['mt'] = adata.var_names.str.startswith("MT-")  # mitochondrial genes as 'mt'
adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL")) # ribosomal genes genes as 'ribo'
adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)(S)]")) & ~adata.var_names.str.contains(("HBEGF")) 
# "^HB[^(P)" changed to "^HB[^(P)(S)" and  & ~adata_test.var_names.str.contains(("HBEGF")) added to remove HBS1L and HBEGF which are NOT memoglobin genes
adata.var['malat1'] = adata.var_names.str.contains(("MALAT1"))  # MALAT1 genes as 'malat1'
return adata
```


`adsctl.pp.calculate_qc_metrics(adata)`
```
calculate_qc_metrics
# add code to check if genes already annotated 
#### code 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) # mitocohndrial  genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True) # ribosomal genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['hb'], percent_top=None, log1p=False, inplace=True) # hemoglobin genes.
sc.pp.calculate_qc_metrics(adata, qc_vars=['malat1'], percent_top=None, log1p=False, inplace=True) # MALAT1 gene.
return adata
```

`adsctl.pp.plot_QC_metrics_scatter(adata)`
```
#### code 
figQC, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(1 ,5,figsize=(20,4), gridspec_kw={'wspace':0.9})
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',ax=ax1, show=False) # plot number of dected genes vs total counts 
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',ax=ax2, show=False) #percent mt counts vs total counts
sc.pl.scatter(adata, x='total_counts', y='pct_counts_ribo',ax=ax3, show=False) #percent ribo counts vs total counts
sc.pl.scatter(adata, x='total_counts', y='pct_counts_malat1',ax=ax4, show=False) #percent HB counts vs total count
sc.pl.scatter(adata, x='total_counts', y='pct_counts_hb',ax=ax5, show=False) #percent HB counts vs total counts
```

`adsctl.pp.plot_QC_metrics_violin(adata)`
```
#### code 
fig1, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(1 ,6,figsize=(20,4), gridspec_kw={'wspace':0.9})
sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4,ax=ax1, show=False)
sc.pl.violin(adata, ['total_counts'], jitter=0.4 ,ax=ax2, show=False)
sc.pl.violin(adata, [ 'pct_counts_mt'], jitter=0.4,ax=ax3, show=False) # mitocohndrial  genes
sc.pl.violin(adata, [ 'pct_counts_ribo'], jitter=0.4,ax=ax4, show=False) # ribosomal genes
sc.pl.violin(adata, [ 'pct_counts_malat1'], jitter=0.4,ax=ax5, show=False) # hemoglobin genes.
sc.pl.violin(adata, [ 'pct_counts_hb'], jitter=0.4,ax=ax6, show=False) # hemoglobin genes.  
```

`adsctl.pp.plot_qc_metrics(adata)`
```
plot_qc_metrics of Annotated technical gene groups  and top 20 highly expressed
#### code 
plot_QC_metrics_violin(adata)  
plot_QC_metrics_scatter(adata) 
sc.pl.highest_expr_genes(adata, n_top=20, )
```

`adsctl.pp.annotate_n_view_adata_raw_counts(adata)`
```
Annotate technical gene groups  and calculate qc metrics
#### code 
annotate_QC_genes(adata)
calculate_qc_metrics(adata)
plot_qc_metrics(adata)
```


`adsctl.pp.filter_cells_by_anotated_QC_gene()`
```
filter_cells_by_anotated_QC_gene(adata,
                                    filter_ncount=True,
                                    n_genes_bycounts=10000,
                                    filter_pct_mt=True,
                                    percent_mt=20,
                                    over_percent_mt=0,
                                    filter_pct_ribo=False,
                                    percent_ribo=100,
                                    over_percent_ribo=0,
                                    filter_pct_hb=False,
                                    percent_hb=100,
                                    over_percent_hb=0,
                                    filter_pct_malat1=False,
                                    percent_malat1=100,
                                    over_percent_malat1=0,
                                      **parameters
                                    )
```

`adsctl.pp.remove_genes()`
```
remove_genes(adata,
                remove_MALAT1=False,
                remove_MT=False,
                remove_HB=False,
                remove_RP_SL=False,
                remove_MRP_SL=False,
                 **parameters
                )
```

`adsctl.pp.process2scaledPCA()`
```
 process2scaledPCA(adata,
                   normalize_total_target_sum=1e4,logarithmize=True,
                   filter_HVG=False,HVG_min_mean=0.0125, HVG_max_mean=3, HVG_min_disp=3,
                   regress_mt=False,regress_ribo=False,regress_malat1=False,regress_hb=False,
                   scale=True,PCA=True,
                      cell_cycle_score=True,
                   regress_cell_cycle_score=False,HVG_flavor='seurat',HVG_n_top_genes=1500,
                      **parameters)
```

`adsctl.pp.leiden_clustering()`
```
#### code
leiden_clustering(adata, number_of_neighbors=10,number_of_PC=40, leiden_res=1,key_added='leiden', **parameters)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=leiden_res,key_added=key_added)
sc.pl.umap(adata, color=[key_added] )
```




## Subpackage - adata_and_scanpy_tools - tools - ANsc.tl.

### latest update is a modified version of scanpy's injest funciton

`adsctl.tl.ingest_verbose(adata, adata_ref, obs=['adata_ref_obs1', 'adata_ref_obs2'])`
```
hello i modified this ingest

* more than just a single obs column (e.g. adata_ref_obs1) is added to the adata_query.obs for each key specified in obs argument
* standard ingests normally adds a cell labels column to the adata_query.obs for each key specified in obs argument
* in this case there was only a single obs key passed so standard ingest would have just added adata.obs['louvain']
* these are the addtional columns that ingest_verbose adds
* 'KNN_distances', #### this is all of the distances from the query cell to the KNN network in the adata_ref reference cells
* 'KNN_avg_distance', #### this is the average distance from the query cell to the KNN network in the adata_ref reference cells
* 'KNN_std_distance',  #### this is the standard deviation of the distances from the query cell to the KNN network in the adata_ref reference cells
* 'adata_ref_obs1_all_KNN', ### this the cell labels for all of the KNN group
* adata_ref_obs1_KNN_pct_most_common' #### frequency of the most common category label in each cells  KNN group is 
* 'adata_ref_obs1_KNN_pct', ##### the pct of each category label in each cell's KNN group is 
* the category corresponding to each pct is stored adata_query.uns['louvain_KNN_pct_cat_labels']
   these are in the same order as the adata.obs['adata_ref_obs1_KNN_pct']  values

* to run this modified ingest 
repo_parent_dir='../../'

import sys
if repo_parent_dir not in sys.path:
   sys.path.append(repo_parent_dir)

import adata_and_scanpy_tools as adsctl
from adata_and_scanpy_tools.adsctl_gex_class import *
# then this is the function call to ingest_verbose all arguements are the same as standard ingest 
# ingest verbose just has additional obs columns outputs

adsctl.tl.ingest_verbose(adata, adata_ref, obs=['adata_ref_obs1', 'adata_ref_obs2'])
################# code that was modified
def to_adata(self, inplace=False):
   """\
   Returns `adata_new` with mapped embeddings and labels.

   If `inplace=False` returns a copy of `adata_new`
   with mapped embeddings and labels in `obsm` and `obs` correspondingly.
   If `inplace=True` returns nothing and updates `adata_new.obsm`
   and `adata_new.obs` with mapped embeddings and labels.
   """
   adata = self._adata_new if inplace else self._adata_new.copy()

   adata.obsm.update(self._obsm)


   for key in self._obs:
      adata.obs[key] = self._obs[key]

   ########################### added code starts here
   ##### add a list of distances values for each ingested obs label
   adata.obs['KNN_distances']=self._distances.tolist()
   ##### add an average distance value for each ingested obs label
   adata.obs['KNN_avg_distance']=np.mean(self._distances, axis=1).tolist()
   ##### add the standard deviation of the distances for each ingested obs label
   adata.obs['KNN_std_distance']=np.std(self._distances, axis=1).tolist()
   
   
   for key in self._obs:
      ##### add a list of label values for each ingested obs label
      cat_array = self._adata_ref.obs[key].astype('category' )  # ensure it's categorical
      values = [cat_array[inds].tolist() for inds in self._indices]
      adata.obs[key+'_all_KNN'] = values

      #### addd percentage of each label in the KNN
      values_KNN_percent = [cat_array[inds].value_counts(sort=False,normalize=True).tolist() for inds in self._indices]
      adata.obs[key+'_KNN_pct'] = values_KNN_percent

      #### add ordered list of cat labels for each ingested obs label to adata.uns to identify adata.obs[key+'_KNN_percent'] values
      adata.uns[key+'_KNN_pct_cat_labels'] = cat_array.cat.categories.tolist()

      #### add the pct of the most common label in the KNN
      values_KNN_pct_most_common = [cat_array[inds].value_counts(normalize=True).tolist()[0] for inds in self._indices]
      adata.obs[key+'_KNN_pct_most_common'] = values_KNN_pct_most_common

   ########################### added code ends here
   if not inplace:
      return adata
```




## Subpackage - adata_and_scanpy_tools - plots - adsctl.pl.

`adsctl.pl.silhouette_score_n_plot(adata,leiden_res='unk',**parameters)`
```
adsctl.pl.silhouette_score_n_plot(adata,parameters,leiden_res='unk'):
> assumes ledien clusteirng to subset cells
> uses X_pca for silhoutte_scores
samples_silhoutte_scores=silhouette_samples(adata.obsm['X_pca'], adata.obs['leiden']

```

`adsctl.pl.plot_batch_obs_key_of_obs_key2(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4),flavor="pct_count")`
```
adsctl.pl.plot_batch_obs_key_of_obs_key2(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4),flavor="pct_count")
makes two side by side bar charts, each bar is a batch_obs_key='batch category and each bar is stacked and colored by obs_key2="leiden"
left bar chart is fraction on y -axis 
right  bar chart is obs/cell count on y -axis 
flavor="pct_count"  >>> both charts
flavor="pct"  >>> only pct chart
flavor="count"  >>> only count chart

```

`adsctl.pl.plot_percent_obs_key2_per_batch_obs_key(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4))`
```
adsctl.pl.plot_percent_obs_key2_per_batch_obs_key(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4)
This produce one column of individual bar charts (one chart for each catagory in batch_obs_key='batch'") # batch_obs_key="sample_ID is good to use
each bar chart show percentage of cells in "batch" assigned to obs_key2="leiden"
```




# Stand alone module - adsctl_gex class- adata_and_scanpy_tools - adsctl_gex - adsctl.adsctl_gex

`adsctl_gex(adata,**parameters)`
* see the notebook 
   * adsctl_gex_class_preprocessing_PBMC3k.ipynb # single class object for preprocessing


# more coming soon