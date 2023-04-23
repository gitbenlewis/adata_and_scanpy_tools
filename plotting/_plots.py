# module level import libraries
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples

sc.settings.set_figure_params(dpi=80, facecolor='white')



def silhouette_score_n_plot(adata,leiden_res='unk',**parameters):
    '''
    silhouette_score_n_plot(adata,parameters,leiden_res='unk'):
    > assumes ledien clusteirng to subset cells
    > uses X_pca for silhoutte_scores
    samples_silhoutte_scores=silhouette_samples(adata.obsm['X_pca'], adata.obs['leiden']
    '''
    '''
    if parameters!=None:
        leiden_res=parameters["leiden_res"]
    '''
    ##################### sillhouette scoreing
    samples_silhoutte_scores=silhouette_samples(adata.obsm['X_pca'], adata.obs['leiden'])
    adata.obs['silhoutte']=samples_silhoutte_scores.tolist()
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_res}')
    #print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of add this to adata.uns later')

    ##################### sillhouette scoreing  #####END

    ###################### umap and sillhouette scoreing graph results of final leiden resolution setting 
    fig_PP2C_cluster_scores, (UMAP_final,ax_final,UMAP_sil,pca_leiden) = plt.subplots(nrows=1, ncols=4, figsize=(20,5), gridspec_kw={'wspace':0.4})

    UMAP_final=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata), ax=UMAP_final,#palette=sc.pl.palettes.vega_20_scanpy,
                          show=False)

    cluster_silhouette_score_list=[]
    for i in range(0, cluster_number):
        cluster_silhouette_score=adata.obs['silhoutte'].loc[adata.obs['leiden']==str(i)].mean()
        cluster_silhouette_score_list.append(cluster_silhouette_score) 

    #  cluster scores
    pre_scores=cluster_silhouette_score_list
    pre_y_pos = np.arange(len(cluster_silhouette_score_list))
    ax_final.barh(pre_y_pos,pre_scores)
    #ax_final.set_yticks(pre_y_pos, labels=pre_y_pos)
    ax_final.set_yticks(pre_y_pos)
    ax_final.set_yticklabels(pre_y_pos) #new
    ax_final.invert_yaxis()  # labels read top-to-bottom
    ax_final.set_title('Cluster Scores')

    UMAP_sil=sc.pl.umap(adata, color=['silhoutte'], ax=UMAP_sil, show=False,#palette=sc.pl.palettes.vega_20_scanpy
                       )

    pca_leiden=sc.pl.pca(adata, color='leiden', ax=pca_leiden, show=False)

    #fig_PP2C_cluster_scores.savefig(dataset_figures_output_directory+'silscore.pdf')
    ###################### umap and sillhouette scoreing graph results of final leiden resolution setting  #####END
    #return adata
    return 



####################################################################################################################

def plot_batch_obs_key_of_obs_key2(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4),flavor="pct_count"):
    '''
    MD_plot_batch_obs_key_of_obs_key2(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4),flavor="pct_count")
    makes two side by side bar charts, each bar is a batch_obs_key='batch category and each bar is stacked and colored by obs_key2="leiden"
    left bar chart is fraction on y -axis 
    right  bar chart is obs/cell count on y -axis 
    flavor="pct_count"  >>> both charts
    flavor="pct"  >>> only pct chart
    flavor="count"  >>> only count chart
    '''
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    df_norm=pd.crosstab(adata.obs[obs_key2],adata.obs[batch_obs_key], normalize='index')
    print(df_norm)
    df=pd.crosstab(adata.obs[obs_key2],adata.obs[batch_obs_key] )
    print(df)
    if flavor=="pct_count":
        fig1, axes = plt.subplots(nrows=1, ncols=2,figsize=figsize)
        ax1=df_norm.plot.bar(stacked=True,ax=axes[0]).legend().set_visible(False)
        ax2=df.plot.bar(stacked=True,ax=axes[1]).legend(loc='upper right')
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
    if flavor=="pct":
        fig1, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize)
        ax1=df_norm.plot.bar(stacked=True,ax=axes).legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
    if flavor=="count":
        fig1, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize)
        ax2=df.plot.bar(stacked=True,ax=axes).legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
    return df

####################################################################################################################


####################################################################################################################

def plot_percent_obs_key2_per_batch_obs_key(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4)):
    '''
     MD_plot_percent_obs_key2_per_batch_obs_key(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4)
     This produce one column of individual bar charts (one chart for each catagory in batch_obs_key='batch'") # batch_obs_key="sample_ID is good to use
     each bar chart show percentage of cells in "batch" assigned to obs_key2="leiden"
    '''
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    obs_key2_Xtab_batch_obs_key_df=pd.crosstab(adata.obs[obs_key2],adata.obs[batch_obs_key] )

    df=obs_key2_Xtab_batch_obs_key_df

    df_col_norm=pd.DataFrame()
    for i in df.columns:
        df_col_norm[i]=list(map(lambda x:x/df[i].sum(axis=0),df[i]))
    print(df_col_norm)
    fig1, axes = plt.subplots(nrows=df_col_norm.shape[1], ncols=1,figsize=figsize)
    ax_n=0
    #obs_key2_groups=np.arange(len(df_col_norm.index))
    obs_key2_groups=adata.obs[obs_key2].cat.categories.tolist()
    for i in df_col_norm.columns:
        #ax=df_col_norm[df_col_norm.columns[ax_n]].plot.bar(ax=axes[ax_n]).legend().set_visible(True)
        axes[ax_n].barh(obs_key2_groups,df_col_norm[df_col_norm.columns[ax_n]].tolist())
        axes[ax_n].set_title(f' {df_col_norm.columns[ax_n]}')
        axes[ax_n].set_yticks(obs_key2_groups, labels=obs_key2_groups)
        axes[ax_n].invert_yaxis() 

        for bars in axes[ax_n].containers:
            axes[ax_n].bar_label(bars, label_type='center',  fmt='%.2g',padding=30,)
        ax_n=ax_n+1
    if savefig==True:
        fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
    return

####################################################################################################################




