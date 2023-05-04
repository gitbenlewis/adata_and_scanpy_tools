# ANgex.py
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

class adsctl_gex:
    def __init__(self, adata, **parameters):
        
        self.adata=adata.copy()
        self.adata_raw_counts = self.adata.copy()
        self.adata.layers["counts"] = self.adata.X.copy()  # preserve counts


        if parameters:
            for key, value in parameters.items():
                setattr(self, key, value)
            self.set_output_directories()
            print('parameters is not empty')
        else:
            print('parameters is empty')
            self.set_default_parameters(**parameters)

        ########## set scanpy settings 
        sc.settings.verbosity = 1             # verbosity: errors (0), warnings (1), info (2), hints (3)
        sc.logging.print_header()
        sc.settings.set_figure_params(dpi=80, facecolor='white')
        sc.settings.n_jobs = self.n_jobs  

        ##########
    
################################## define methods

    def set_output_directories(self):
        ##################### ######################## make output directories  
        os.makedirs(self.output_dir+self.output_prefix, exist_ok=True) ####  make output directory and dataset output director
        
        os.makedirs(self.output_dir+self.output_prefix+'/tables/', exist_ok=True) ####  make output directory for tables
        self.dataset_tables_output_directory=self.output_dir+self.output_prefix+'/tables/'

        os.makedirs(self.output_dir+self.output_prefix+'/figures/', exist_ok=True)  ####  make output directory for figures
        self.dataset_figures_output_directory=self.output_dir+self.output_prefix+'/figures/'
        
        sc.settings.figdir=self.dataset_figures_output_directory


    def set_default_parameters(self):
        ##################### set_default_parameters variables start
        self.n_jobs=24
        #### dataset specfic parameters
        self.output_prefix="Dataset_"
        self.output_dir="ProjectName_"
        
        ###Basic filters
        self.filter_genes_min_cells=3 # min  of cells a gene is detected in else gene is tossed out default 3
        self.filter_genes_min_counts=0 # min  of counts a gene must have to pass basic filter default 0
        self.filter_cells_min_genes=200 # min  of genes detected or else cell/observation is tossed out default 200
        self.filter_cells_min_counts=0 # min  of counts detected or else cell/observation is tossed out default 0
        
        ####Filter  on off switches
        self.filter_ncount=True
        self.filter_pct_mt=True
        self.filter_pct_ribo=False
        self.filter_pct_hb=False
        self.filter_pct_malat1=False
        self.filter_HVG=False
        ### less than filter percent 
        self.n_genes_bycounts=7000 # less than filter
        self.percent_mt=10 # less than filter
        self.percent_ribo=100 # less than filter
        self.percent_malat1 =100 # less than filter
        self.percent_hb=100 # less than filter
        ### Greater than filter percent
        #self.over_n_genes_bycounts=200 # greater than filter  # now filter_cells_min_genes
        self.over_percent_mt=0 # greater than filter
        self.over_percent_ribo=0 # greater than filter
        self.over_percent_malat1 =0 # greater than filter
        self.over_percent_hb=0 # greater than filter
        ### Remove gene sets  on off switches
        self.remove_MALAT1=False
        self.remove_MT=False
        self.remove_HB=False
        self.remove_RP_SL=False
        self.remove_MRP_SL=False
        #### processing parameters and options
        self.normalize_total_target_sum=1e4  # scanpy  default 1e4
        self.filter_genes_min_counts_normed=0 # min counts to keep a gene after library size normalization
        self.HVG_min_mean=0.0125  # scanpy  default 0.0125
        self.HVG_max_mean=3  # scanpy  default 3
        self.HVG_min_disp=0.5  # scanpy  default 0.5
        self.logarithmize=True  # scanpy  default True
        self.scale=True # scanpy  default True
        #### regression on off switches
        self.regress_mt=False
        self.regress_ribo=False
        self.regress_malat1=False
        self.regress_hb=False
        self.regress_cell_cycle_score=False
        ### clustering parameters for clusters
        self.number_of_PC=30 ### dataset demensionality 
        self.number_of_neighbors=10
        self.leiden_res=1 #leiden clustering resolution
        # UMAP graph parameters
        self.umap_marker_gene=False
        self.umap_marker_gene_list=[#'IL7R',
            'CD14','LYZ', 'MS4A1','CD8A','GNLY','NKG7','FCGR3A','MS4A7','FCER1A','CST3','PPBP'] # from PBMC 3k
        # cluster naming parameters
        self.rename_cluster=False
        self.new_cluster_names=False
        ##################### set_default_parameters variables END
        self.set_output_directories()


    def basic_filitering(self):
        """ Basic Filtering """
        print(f'number of Cells BEFORE Basic Filtering : {self.adata.n_obs}')
        sc.pp.filter_cells(self.adata, min_genes=self.filter_cells_min_genes,)  #min_genes=self.over_n_genes_bycounts
        print(f'Filtering cells pp.filter_cells(adata, min_cells=filter_cells_min_genes)  Cells remaining : {self.adata.n_obs}')
        print(f'min_cells=filter_cells_min_genes = ', self.filter_cells_min_genes )
        sc.pp.filter_cells(self.adata,min_counts=self.filter_cells_min_counts,)  # cells / observations must have min # of coutns
        print(f'Filtering cells pp.filter_cells(adata, min_cells=filter_cells_min_counts)  Cells remaining : {self.adata.n_obs}')
        print(f'min_counts=filter_cells_min_counts = ',self.filter_cells_min_counts )
        print(f'number of GENES BEFORE Basic Filtering : {self.adata.n_vars}')
        sc.pp.filter_genes(self.adata, min_cells=self.filter_genes_min_cells ) #genes must be present in min # of cells / observations
        print(f'Filtering genes pp.filter_genes(adata, min_cells=filter_genes_min_cells)  Genes remaining : {self.adata.n_vars}')
        print(f'min_cells=filter_genes_min_cells = ',self.filter_genes_min_cells )
        sc.pp.filter_genes(self.adata, min_counts=self.filter_genes_min_counts ) #genes must have min # of counts for gene to be kept
        print(f'Filtering genes pp.filter_genes(adata, min_cells=filter_genes_min_counts)  Genes remaining :  {self.adata.n_vars}')
        print(f'min_counts=filter_genes_min_counts = ',self.filter_genes_min_counts )
        return self.adata

    def annotate_QC_genes(self):
        """ annotate the group of QC genes """
        self.adata.var['mt'] = self.adata.var_names.str.startswith("MT-")  # mitochondrial genes as 'mt'
        self.adata.var['ribo'] = self.adata.var_names.str.startswith(("RPS","RPL")) # ribosomal genes genes as 'ribo'
        self.adata.var['hb'] = self.adata.var_names.str.contains(("^HB[^(P)(S)]")) & ~self.adata.var_names.str.contains(("HBEGF")) 
        # "^HB[^(P)" changed to "^HB[^(P)(S)" and  & ~adata_test.var_names.str.contains(("HBEGF")) added to remove HBS1L and HBEGF which are NOT memoglobin genes
        self.adata.var['malat1'] = self.adata.var_names.str.contains(("MALAT1"))  # MALAT1 genes as 'malat1'
        return self.adata

    def calculate_qc_metrics(self):
        """ calculate_qc_metrics"""
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) # mitocohndrial  genes
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True) # ribosomal genes
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['hb'], percent_top=None, log1p=False, inplace=True) # hemoglobin genes.
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['malat1'], percent_top=None, log1p=False, inplace=True) # MALAT1 gene.
        return self.adata

    
    def annotate_n_view_adata_raw_counts(self):
        """  Annotate technical gene groups  and calculate qc metrics"""
        self.adata=self.annotate_QC_genes()
        self.adata=self.calculate_qc_metrics()
        self.plot_qc_metrics()
        return
    
    
    def filter_cells_by_anotated_QC_gene(self):
        """  Remove cells that have too many mitochondrial genes expressed or too many total counts:""" 

        print(f' {self.filter_ncount} keep cells with less than {self.n_genes_bycounts} (n_genes_bycounts) dected genes ')

        print(f' {self.filter_pct_mt} keep cells with less than {self.percent_mt} (percent_mt) mitochondiral gene counts ')
        print(f' {self.filter_pct_mt} keep cells with greater than {self.over_percent_mt} (percent_mt) mitochondiral gene counts ')

        print(f' {self.filter_pct_ribo} keep cells with less than {self.percent_ribo} (percent_ribo) ribosomal protein gene counts ')
        print(f' {self.filter_pct_ribo} keep cells with greater than {self.over_percent_ribo} (percent_ribo) ribosomal protein gene counts ')


        print(f' {self.filter_pct_hb} keep cells with less than {self.percent_hb} (percent_hb) hemoglobin protein gene counts ')
        print(f' {self.filter_pct_hb} keep cells with greater than {self.over_percent_hb} (percent_ribo) ribosomal protein gene counts ')

        # Actually do the filtering by slicing the `AnnData` object.
        print(f'number of Cells BEFORE pct Filtering : {self.adata.n_obs}')
        if self.filter_ncount ==True:
            self.adata = self.adata[self.adata.obs.n_genes_by_counts <= self.n_genes_bycounts, :]  # by n_genes_bycounts
            print(f'number of Cells AFTER n_genes_bycounts Filtering : {self.adata.n_obs}')
        if self.filter_pct_mt ==True:
            self.adata = self.adata[self.adata.obs.pct_counts_mt <= self.percent_mt, :]   # by percent_mt
            print(f'number of Cells AFTER percent_mt Filtering : {self.adata.n_obs}')
            self.adata = self.adata[self.adata.obs.pct_counts_mt >= self.over_percent_mt, :]   # by percent_mt
            print(f'number of Cells AFTER over_percent_mt Filtering : {self.adata.n_obs}')
        if self.filter_pct_ribo ==True:
            self.adata = self.adata[self.adata.obs.pct_counts_ribo <= self.percent_ribo, :]   # by percent_ribo
            print(f'number of Cells AFTER percent_ribo Filtering : {self.adata.n_obs}')
            self.adata = self.adata[self.adata.obs.pct_counts_ribo >= self.over_percent_ribo, :]   # by percent_ribo
            print(f'number of Cells AFTER over_percent_ribo Filtering : {self.adata.n_obs}')
        if self.filter_pct_malat1 ==True:
            self.adata = self.adata[self.adata.obs.pct_counts_malat1 <= self.percent_malat1, :]    # by percent_hb
            print(f'number of Cells AFTER percent_malat1 Filtering : {self.adata.n_obs}')
            self.adata = self.adata[self.adata.obs.pct_counts_malat1 >= self.over_percent_malat1, :]    # by percent_hb
            print(f'number of Cells AFTER over_percent_malat1 Filtering : {self.adata.n_obs}')        
        if self.filter_pct_hb ==True:
            self.adata = self.adata[self.adata.obs.pct_counts_hb <= self.percent_hb, :]    # by percent_hb
            print(f'number of Cells AFTER percent_hb Filtering : {self.adata.n_obs}')
            self.adata = self.adata[self.adata.obs.pct_counts_hb >= self.over_percent_hb, :]    # by percent_hb
            print(f'number of Cells AFTER over_percent_hb Filtering : {self.adata.n_obs}')
        return self.adata

    def remove_genes(self):
        """ ################################# Remove Filter out genes with ""techincal bias""
        ### Remove gene sets  on off switches""" 
        print(f'remove_MALAT1 {self.remove_MALAT1}')
        print(f'remove_MT {self.remove_MT}')
        print(f'remove_HB {self.remove_HB}')
        print(f'remove_RP_SL {self.remove_RP_SL}')
        print(f'remove_MRP_SL {self.remove_MRP_SL} ')


        print(f' BEFORE filtering for specific gene : number of Cells {self.adata.n_obs}, number of genes {self.adata.n_vars}')

        nothing = self.adata.var_names.str.startswith('NO_GENES_HAVE_THIS_NAME')
        remove = np.add(nothing, nothing)
        print(len((nothing)))

        if self.remove_MALAT1==True:
            malat1 = self.adata.var_names.str.startswith('MALAT1')
            remove = np.add(remove, malat1)
        # we need to redefine the mito_genes since they were first 
        # calculated on the full object before removing low expressed genes.

        if self.remove_MT==True:
            mito_genes = self.adata.var_names.str.startswith('MT-')
            remove = np.add(remove,mito_genes)
        if self.remove_HB==True:
            hb_genes = (self.adata.var_names.str.startswith('HB')& ~self.adata.var_names.str.contains(("HBEGF"))  & ~self.adata.var_names.str.contains(("HBS1L"))  & ~self.adata.var_names.str.contains(("HBP1"))) # HBEGF,HBS1L, HBP1 not a hemeoglobin genes 
            remove = np.add(remove,hb_genes)
        if self.remove_RP_SL==True:
            RP_SL_genes = self.adata.var_names.str.startswith(("RPS","RPL"))
            remove = np.add(remove,RP_SL_genes )
        if self.remove_MRP_SL==True:
            MRP_SL_genes = self.adata.var_names.str.startswith(("MRPS","MRPL"))
            remove = np.add(remove,MRP_SL_genes )    

        keep = np.invert(remove)
        self.adata = self.adata[:,keep]

        print(f' AFTER filtering for specific gene : number of Cells {self.adata.n_obs}, number of genes {self.adata.n_vars}')
        return self.adata

    def regress_out_anotated_QC_genes(self):
        ################################# and Regression 
        print(f'we are regressing out  total_counts {self.regress_mt}')
        print(f'we are regressing out  pct_counts_mt {self.regress_mt}')
        print(f'we are regressing out  pct_counts_ribo {self.regress_ribo}')
        print(f'we are regressing out  pct_counts_malat1 {self.regress_malat1}')
        print(f'we are regressing out  pct_counts_hb {self.regress_hb}')

        ################ Do the regression 
        if self.regress_mt ==True:
            # by total_counts
            sc.pp.regress_out(self.adata, ['total_counts' ])
        if self.regress_mt ==True:
            # by percent_mt
            sc.pp.regress_out(self.adata, ['pct_counts_mt' ])
        if self.regress_ribo ==True:
            # by percent_ribo
            sc.pp.regress_out(self.adata, ['pct_counts_ribo' ])
        if self.regress_malat1 ==True:
            # by percent_hb
            sc.pp.regress_out(self.adata, ['pct_counts_malat1'])
        if self.regress_hb ==True:
            # by percent_hb
            sc.pp.regress_out(self.adata, ['pct_counts_hb' ])
        return self.adata

    def cell_cycle_score_and_regress(self):
        """ ################# cell cycle score  and (True/False) Regress out effects of cell cylce score"""
        # Import cell cycle list and split into s_genes and g2m_genes
        s_genes=['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
        g2m_genes=['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1',         'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']
        cell_cycle_genes=g2m_genes+s_genes
        print(f' there are {len(s_genes)} s_genes   {len(g2m_genes)} g2m_genes  {len(cell_cycle_genes)} cell_cycle_genes')
        cell_cycle_genes = [x for x in cell_cycle_genes if x in self.adata.var_names]
        print(f' there are {len(cell_cycle_genes)} cell_cycle_genes in the dataset')    
         ## do scoring 
        sc.tl.score_genes_cell_cycle(self.adata, s_genes=s_genes, g2m_genes=g2m_genes)
       # plot the cell cycle scores 
        sc.pl.violin(self.adata, ['S_score', 'G2M_score'],jitter=0.4,rotation=45)
     # Regress out effects of cell cylce score   # cell cycle regress
        print(f'we are regressing out  cell_cycle_score {self.regress_cell_cycle_score}')
        if self.regress_cell_cycle_score ==True:
            ## do scoring 
            #do regression 
            sc.pp.regress_out(self.adata, ['S_score', 'G2M_score'])
        return self.adata


### plot functions
    def plot_QC_metrics_scatter(self):
        figQC, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(1 ,5,figsize=(20,4), gridspec_kw={'wspace':0.9})
        sc.pl.scatter(self.adata, x='total_counts', y='n_genes_by_counts',ax=ax1, show=False) # plot number of dected genes vs total counts 
        sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_mt',ax=ax2, show=False) #percent mt counts vs total counts
        sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_ribo',ax=ax3, show=False) #percent ribo counts vs total counts
        sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_malat1',ax=ax4, show=False) #percent HB counts vs total count
        sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_hb',ax=ax5, show=False) #percent HB counts vs total counts 
        return

    def plot_QC_metrics_violin(self):
        fig1, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(1 ,6,figsize=(20,4), gridspec_kw={'wspace':0.9})
        sc.pl.violin(self.adata, ['n_genes_by_counts'], jitter=0.4,ax=ax1, show=False)
        sc.pl.violin(self.adata, ['total_counts'], jitter=0.4 ,ax=ax2, show=False)
        sc.pl.violin(self.adata, [ 'pct_counts_mt'], jitter=0.4,ax=ax3, show=False) # mitocohndrial  genes
        sc.pl.violin(self.adata, [ 'pct_counts_ribo'], jitter=0.4,ax=ax4, show=False) # ribosomal genes
        sc.pl.violin(self.adata, [ 'pct_counts_malat1'], jitter=0.4,ax=ax5, show=False) # hemoglobin genes.
        sc.pl.violin(self.adata, [ 'pct_counts_hb'], jitter=0.4,ax=ax6, show=False) # hemoglobin genes.
        return self.adata

    def silhouette_score_n_plot(self):
        ##################### sillhouette scoreing
        samples_silhoutte_scores=silhouette_samples(self.adata.obsm['X_pca'], self.adata.obs['leiden'])
        self.adata.obs['silhoutte']=samples_silhoutte_scores.tolist()
        silhouette_score_adata=silhouette_score(self.adata.obsm['X_pca'], self.adata.obs['leiden'],)
        cluster_number=len(set(self.adata.obs['leiden'].tolist()))
        print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {self.leiden_res}')
        ##################### sillhouette scoreing  #####END


        ###################### umap and sillhouette scoreing graph results of final leiden resolution setting 
        fig_PP2C_cluster_scores, (UMAP_final,ax_final,UMAP_sil,pca_leiden) = plt.subplots(nrows=1, ncols=4, figsize=(20,5), gridspec_kw={'wspace':0.4})

        UMAP_final=sc.pl.umap(self.adata, color='leiden',title='Leiden.res='+str(self.leiden_res)+' Avg.sil.='+str(silhouette_score_adata), ax=UMAP_final,#palette=sc.pl.palettes.vega_20_scanpy,
                              show=False)

        cluster_silhouette_score_list=[]
        for i in range(0, cluster_number):
            cluster_silhouette_score=self.adata.obs['silhoutte'].loc[self.adata.obs['leiden']==str(i)].mean()
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

        UMAP_sil=sc.pl.umap(self.adata, color=['silhoutte'], ax=UMAP_sil, show=False,#palette=sc.pl.palettes.vega_20_scanpy
                           )

        pca_leiden=sc.pl.pca(self.adata, color='leiden', ax=pca_leiden, show=False)

        fig_PP2C_cluster_scores.savefig(self.dataset_figures_output_directory+'silscore.pdf')

        ###################### umap and sillhouette scoreing graph results of final leiden resolution setting  #####END
        return self.adata

    def plot_qc_metrics(self):
        """ plot_qc_metrics of Annotated technical gene groups  and top 20 highly expressed"""
        self.plot_QC_metrics_violin()  
        self.plot_QC_metrics_scatter() 
        sc.pl.highest_expr_genes(self.adata, n_top=20, )
    

    
    
    def PP(self):
        """ adsctl_gex.PP() doc string to add later"""
        sc.logging.print_header()
        
        print(f'############################################################################# # reset self.adata to raw counts')
        self.adata=self.adata_raw_counts.copy() # reset self.adata to raw counts
        print(f'################################################ # preserve counts ########## # self.adata.layers["counts"] = self.adata.X.copy() ')
        self.adata.layers["counts"] = self.adata.X.copy()  # preserve counts

        print(f'############################################################################# preprocessing start')
        print(f'#################################################### Basic Filtering')
        #####  Basic Filtering
        self.adata=self.basic_filitering()
        print(f'#################################################### Annotate technical gene groups  and calculate qc metrics')
        print(f'####################################################  plot_qc_metrics before QC metric filters and gene removal ')
        ####### Annotate technical gene groups  and calculate qc metrics
        self.annotate_n_view_adata_raw_counts()
        print(f'####################################################  filter by QC metrics')
        ################# filter by QC metrics 
        self.adata=self.filter_cells_by_anotated_QC_gene()
        print(f'####################################################  remove_genes')
        ################# remove specfic gene groups 
        self.adata=self.remove_genes()
        print(f'####################################################  plot_qc_metrics after filters and gene removal ')
        # #  re annoate / re calc  and  after fitlering and removing gene groups  plot_qc_metrics
        self.plot_qc_metrics()
        print(f'####################################################  library-size correct  the data')
        ################################## library-size correct  the data:
        sc.pp.normalize_total(self.adata, target_sum=self.normalize_total_target_sum)  # normalize to depth normalize_total_target_sum usually 1e4 ...1e6 equals CPM
        print(f'library-size correct the observations sc.pp.normalize_total(adata, target_sum=self.normalize_total_target_sum)   depth normalize_total_target_sum : {self.normalize_total_target_sum}')
        if self.filter_genes_min_counts_normed!=0:
            sc.pp.filter_genes(self.adata, min_counts=self.filter_genes_min_counts_normed ) # cells / observations must have min # of coutns
            print(f'Filtering genes pp.filter_genes(adata, min_counts=filter_genes_min_counts_normed)  Genes remaining : {self.adata.n_vars}')
        else:
             print(f'filter_genes_min_counts_normed = 0 ... skipping > Filtering genes pp.filter_genes(adata, min_counts=filter_genes_min_counts_normed')
        if self.logarithmize==True:
            print(f'####################################################  Logarithmize  the data')
            ################################## Logarithmize
            sc.pp.log1p(self.adata)  # logaritmize
        print(f'####################################################  Identify highly-variable genes and plot')
        ########################  Identify highly-variable genes.
        sc.pp.highly_variable_genes(self.adata, min_mean=self.HVG_min_mean, max_mean=self.HVG_max_mean, min_disp=self.HVG_min_disp)
        print(f'############################# the number of highly varriable gens are = ',sum(self.adata.var.highly_variable))
        sc.pl.highly_variable_genes(self.adata) #### plot HVGs
        
        ################ save "raw" data before regression
        # for later use in differential testing and visualizations of gene expression. 
        # This simply freezes the state of the AnnData object.
        if self.logarithmize==True:
            print(f'############################# to adata.raw save filtered, normalized and logarithmized gene expression and plot')
        else:
            print(f'############################# to adata.raw save filtered and normalized gene expression and plot')
        self.adata.raw = self.adata
        ################################# Filtering for HVG
        if self.filter_HVG ==True:
            print(f' Before  filtering for highly_variable genes : number of Cells {self.adata.n_obs}, number of genes {self.adata.n_vars}')
            print(f' filter_HVG = True ... only highly_variable gene will be kept ')
            self.adata = self.adata[:, self.adata.var.highly_variable] # Keep only highly variable genes
            print(f' AFTER  filtering for highly_variable genes: number of Cells {self.adata.n_obs}, number of genes {self.adata.n_vars}')
        else:
            print(f' filter_HVG = False ... all genes will be kept ')
        print(f'####################################################  regress_out_anotated_QC_genes ')
        self.adata=self.regress_out_anotated_QC_genes()

        if self.scale==True:
            print(f'####################################################  Scale the data (each gene to unit variance)')
        ############################### Scale each gene to unit variance. (do before cell cylce regression)  
            sc.pp.scale(self.adata, max_value=10)  # Scale and Clip values exceeding standard deviation 10.    
        print(f'####################################################  regress out cell cycle score and (True/False) Regress out score')
        self.adata=self.cell_cycle_score_and_regress()
        print(f'####################################################  Principal component analysis ')
        ###################################### Principal component analysis
        sc.tl.pca(self.adata, svd_solver='arpack')
        sc.pl.pca_variance_ratio(self.adata, log=True)
        print(f'####################################################  Graph based clustering ')
        print(f'## Computing the neighborhood graph ')
        print(f'number_of_neighbors parameter uset to {self.number_of_neighbors} number_of_PC parameter uset to {self.number_of_PC}')
        sc.pp.neighbors(self.adata, n_neighbors=self.number_of_neighbors, n_pcs=self.number_of_PC)
        print(f'## Embedding the neighborhood graph ')
        sc.tl.umap(self.adata)
        ################################################################ Leiden based clustering
        print(f'## Leiden based clustering ')
        print(f'leiden_res parameter set to {self.leiden_res}')   # Print the leiden resolution parameter used 
        sc.tl.leiden(self.adata,resolution=self.leiden_res)  # perofrom the leiden clustering
        
        #########################################  preprocessing complete
        print(f'#################################################### preprocessing complete')
        
        
        print(f'#################################################### ploting sillhouette scoreing ')
        ##################### sillhouette scoreing
        self.adata=self.silhouette_score_n_plot()
        

        ##################### cluster remnameing 
        if self.rename_cluster==True:
            cluster_numbers_len=len(set(self.adata.obs['leiden'].tolist()))
            if cluster_numbers_len==len(self.new_cluster_names):
                self.adata.obs['leiden_renamed']=self.adata.obs['leiden']
                self.adata.rename_categories('leiden_renamed', self.new_cluster_names)
        ##################### cluster remnameing #####END
        
        print(f'#################################################### ploting umap of marker genes ')
        ######################## plotting umap 
        if self.umap_marker_gene==True:
            cluster_numbers_len=len(set(self.adata.obs['leiden'].tolist()))
            if ((self.rename_cluster==True) & (cluster_numbers_len==len(self.new_cluster_names))):
                print(f'## leiden clusters_renamed ')
                sc.pl.umap(self.adata, color=self.umap_marker_gene_list+['leiden','leiden_renamed'],save=self.output_prefix+'umap_markergenes.pdf')
            else:
                sc.pl.umap(self.adata, color=self.umap_marker_gene_list+['leiden'],
                           save=self.output_prefix+'umap_markergenes.pdf')
            
        else:
            sc.pl.umap(self.adata ,color='leiden',
                       save=self.output_prefix+'umap_markergenes.pdf')

        ######################## ANgex.PP(self)  #####END


    
    def df_loadings_ordered_byPC(self,ascending=False,save_table=False):
        """
        df_loadings_ordered_byPC(self,ascending=False,save_table=False)
        ascending=False gives ...
        ########### idea from https://github.com/scverse/scanpy/issues/836
        """
        os.makedirs(self.output_dir+self.output_prefix, exist_ok=True)
        os.makedirs(self.output_dir+self.output_prefix+'/tables/', exist_ok=True)

        dataset_tables_output_directory=self.output_dir+self.output_prefix+'/tables/'

        df_loadings = pd.DataFrame(self.adata.varm['PCs'], index=self.adata.var_names)
        df_loadings_ordered_byPC=pd.DataFrame()
        
        for i in df_loadings.columns:
            df_loadings_ordered_byPC['PC_'+str(i+1)+'_n']=df_loadings[df_loadings.columns[i]].sort_values(ascending=ascending).index.tolist()
            df_loadings_ordered_byPC['PC_'+str(i+1)+'_val']=df_loadings[df_loadings.columns[i]].sort_values(ascending=ascending).tolist()
        if save_table==True:
            if ascending==False:
                df_loadings_ordered_byPC.to_csv(self.dataset_tables_output_directory+self.output_prefix+"PC_embedings_POS.csv")
            if ascending==True:
                df_loadings_ordered_byPC.to_csv(self.dataset_tables_output_directory+self.output_prefix+"PC_embedings_NEG.csv")
        return df_loadings_ordered_byPC

    
    
    
    def rank_genes(self,wilcox=True,logreg=True,t_test=True,rank_use_raw=True,obs_key="leiden"):
        """
        adsctl_gex.rank_genes(self,
        adata=self.adata
        output_dir=self.output_dir, # use same output_dir as in the parameters["output_dir"] 
        output_prefix=self.output_prefix, # use same output_prefix as in as in the parameters["output_prefix"] 
        wilcox=True,logreg=True,t_test=True, ####  which test to run 
        rank_use_raw=True, # if set to false only uses the highly varrible genes 
        obs_key="leiden", adata.obs key to use to find differentially expressed genes
        n_jobs=self.n_jobs # number of threads
        )
        dataframes saved to variables 
        self.rank_genes_groups_wilcox
        self.rank_genes_groups_logreg
        self.rank_genes_groups_t_test
        dataframes saved to  .csv files y
        self.output_dir+self.output_prefix+'/tables/'+self.output_prefix+"rank_genes_groups_wilcox.csv"
        self.output_dir+self.output_prefix+'/tables/'+self.output_prefix+"rank_genes_groups_logreg.csv"
        self.output_dir+self.output_prefix+'/tables/'+self.output_prefix+"rank_genes_groups_t_test.csv"
        
        """
        import os
        import numpy as np
        import pandas as pd
        import scanpy as sc
        sc.settings.verbosity = 1             # verbosity: errors (0), warnings (1), info (2), hints (3)
        sc.logging.print_header()
        sc.settings.set_figure_params(dpi=80, facecolor='white')
        sc.settings.n_jobs = self.n_jobs  

        os.makedirs(self.output_dir+self.output_prefix, exist_ok=True)


        os.makedirs(self.output_dir+self.output_prefix+'/tables/', exist_ok=True)
        dataset_tables_output_directory=self.output_dir+self.output_prefix+'/tables/'

        os.makedirs(self.output_dir+self.output_prefix+'/figures/', exist_ok=True)
        dataset_figures_output_directory=self.output_dir+self.output_prefix+'/figures/'

        sc.settings.figdir=dataset_figures_output_directory

        
        rank_genes_groups_wilcox=pd.DataFrame()
        rank_genes_groups_logreg=pd.DataFrame()
        rank_genes_groups_t_test=pd.DataFrame()
        
        #bug work around found on github
        #needed for next cell to run
        self.adata.uns['log1p']["base"] = None
        if wilcox==True:
            #########################  Wilcox
            sc.tl.rank_genes_groups(self.adata, obs_key, method='wilcoxon', use_raw=rank_use_raw)
            sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False,
                                    save=self.output_prefix+'wilcoxon_topgenes.pdf')

            result = self.adata.uns['rank_genes_groups']
            groups = result['names'].dtype.names
            rank_genes_groups_wilcox=pd.DataFrame(
                {group + '_' + key[:15]: result[key][group]
                for group in groups for key in ['names', 'scores','pvals', 'pvals_adj','logfoldchanges']})
            self.rank_genes_groups_wilcox=rank_genes_groups_wilcox
            rank_genes_groups_wilcox.to_csv(dataset_tables_output_directory+self.output_prefix+"rank_genes_groups_wilcox.csv")

            #########################  Wilcox
        if logreg==True:
            #########################  logical reggression
            sc.tl.rank_genes_groups(self.adata, obs_key, method='logreg',use_raw=rank_use_raw)
            sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False,
                                    save=self.output_prefix+'logreg_topgenes.pdf')

            rank_genes_groups_logreg=pd.DataFrame(self.adata.uns['rank_genes_groups']['names'])
            self.rank_genes_groups_logreg=rank_genes_groups_logreg

            rank_genes_groups_logreg.to_csv(dataset_tables_output_directory+self.output_prefix+"rank_genes_groups_logreg.csv")

            #########################  logical reggression

        if t_test==True:                                
            ######################### t-test
            sc.tl.rank_genes_groups(self.adata, obs_key, method='t-test',use_raw=rank_use_raw)
            sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False,
                                   save=self.output_prefix+'t_test_topgenes.pdf')
            result = self.adata.uns['rank_genes_groups']
            groups = result['names'].dtype.names
            rank_genes_groups_t_test=pd.DataFrame(
                {group + '_' + key[:15]: result[key][group]
                for group in groups for key in ['names', 'scores','pvals','pvals_adj','logfoldchanges']})
            self.rank_genes_groups_t_test=rank_genes_groups_t_test
            rank_genes_groups_t_test.to_csv(dataset_tables_output_directory+self.output_prefix+"rank_genes_groups_t_test.csv")
             ######################### t-test
        #return (rank_genes_groups_wilcox,rank_genes_groups_logreg,rank_genes_groups_t_test)



    def rank_genes_obscat1_vs_obscat2(self,wilcox=True,logreg=True,t_test=True,rank_use_raw=True,
                                      obs_key="leiden",obscat1='0',obscat2='1'):
        """
        adsctl_gex.rank_genes_obscat1_vs_obscat2(
        adata=self.adata,
        output_dir=self.output_dir  # use same output_dir as in the parameters["output_dir"] used in MD_PP2C(adata,parameters)
        output_prefix=self.output_prefix #  use same output_prefix as in as in the parameters["output_prefix"] used in MD_PP2C(adata,parameters)
        wilcox=True,logreg=True,t_test=True, ####  which test to run 
        rank_use_raw=True, # if set to false only uses the highly varrible genes 
        n_jobs=self.n_jobs # number of threads
        obs_key="leiden", adata.obs key to use to find differentially expressed genes
        obscat1='0' # diffenretioally expressed genes in adata[obs_key]=obscat1 vs adata[obs_key]=obscat2
        obscat2='1'
        )
        """
        import os
        import numpy as np
        import pandas as pd
        import scanpy as sc
        sc.settings.verbosity = 1             # verbosity: errors (0), warnings (1), info (2), hints (3)
        sc.logging.print_header()
        sc.settings.set_figure_params(dpi=80, facecolor='white')
        sc.settings.n_jobs = self.n_jobs  

        os.makedirs(self.output_dir+self.output_prefix, exist_ok=True)


        os.makedirs(self.output_dir+self.output_prefix+'/tables/', exist_ok=True)
        dataset_tables_output_directory=self.output_dir+self.output_prefix+'/tables/'

        os.makedirs(self.output_dir+self.output_prefix+'/figures/', exist_ok=True)
        dataset_figures_output_directory=self.output_dir+self.output_prefix+'/figures/'

        sc.settings.figdir=dataset_figures_output_directory

        
        rank_genes_groups_wilcox=pd.DataFrame()
        rank_genes_groups_logreg=pd.DataFrame()
        rank_genes_groups_t_test=pd.DataFrame()
        
        #bug work around found on github
        #needed for next cell to run
        self.adata.uns['log1p']["base"] = None
        if wilcox==True:
            #########################  Wilcox
            sc.tl.rank_genes_groups(self.adata,obs_key, groups=[obscat1], reference=obscat2, method='wilcoxon', use_raw=rank_use_raw)
            #sc.tl.rank_genes_groups(adata, obs_key, method='wilcoxon', use_raw=rank_use_raw)
            sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False,
                                    save=self.output_prefix+obs_key+'_'+obscat1+'_VS_'+obscat2+'_'+'wilcoxon_topgenes.pdf')

            result = self.adata.uns['rank_genes_groups']
            groups = result['names'].dtype.names
            rank_genes_groups_wilcox=pd.DataFrame(
                {group + '_' + key[:15]: result[key][group]
                for group in groups for key in ['names', 'scores','pvals','pvals_adj','logfoldchanges']})
            rank_genes_groups_wilcox.to_csv(dataset_tables_output_directory+self.output_prefix+obs_key+'_'+obscat1+'_VS_'+obscat2+'_'+"rank_genes_groups_wilcox.csv")

            #########################  Wilcox
        if logreg==True:
            #########################  logical reggression
            sc.tl.rank_genes_groups(self.adata,obs_key, groups=[obscat1], reference=obscat2, method='logreg', use_raw=rank_use_raw)
           # sc.tl.rank_genes_groups(adata, obs_key, method='logreg',use_raw=rank_use_raw)
            sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False,
                                    save=self.output_prefix+obs_key+'_'+obscat1+'_VS_'+obscat2+'_'+'logreg_topgenes.pdf')

            rank_genes_groups_logreg=pd.DataFrame(self.adata.uns['rank_genes_groups']['names'])

            rank_genes_groups_logreg.to_csv(dataset_tables_output_directory+self.output_prefix+obs_key+'_'+obscat1+'_VS_'+obscat2+'_'+"rank_genes_groups_logreg.csv")

            #########################  logical reggression

        if t_test==True:                                
            ######################### t-test
            sc.tl.rank_genes_groups(self.adata,obs_key, groups=[obscat1], reference=obscat2, method='t-test', use_raw=rank_use_raw)
            #sc.tl.rank_genes_groups(adata, obs_key, method='t-test',use_raw=rank_use_raw)
            sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False,
                                   save=self.output_prefix+obs_key+'_'+obscat1+'_VS_'+obscat2+'_'+'t_test_topgenes.pdf')
            result = self.adata.uns['rank_genes_groups']
            groups = result['names'].dtype.names
            rank_genes_groups_t_test=pd.DataFrame(
                {group + '_' + key[:15]: result[key][group]
                for group in groups for key in ['names', 'scores','pvals','pvals_adj','logfoldchanges']})
            rank_genes_groups_t_test.to_csv(dataset_tables_output_directory+self.output_prefix+obs_key+'_'+obscat1+'_VS_'+obscat2+'_'+"rank_genes_groups_t_test.csv")    
        return (rank_genes_groups_wilcox,rank_genes_groups_logreg,rank_genes_groups_t_test)
    
    
    
    
    
    
    def GSEA_enrichr_all_clusters(self,test_library_names=['GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021'], top_nth=10):
        """This is the doc string
        This functions take the tables produced by the MD_rank_genes(adata,output_dir,output_prefix) function and perfroms GSEA analysis using the gseapy enrichr package 

        default arguements:
        adsctl_gex.GSEA_enrichr_all_clusters(
        output_dir:self.output_dir # set this to same output_dir used for MD_rank_genes(adata,output_dir,output_prefix)
        output_prefix:self.output_prefix # set this to same output_prefix directory used for MD_rank_genes(adata,output_dir,output_prefix)
        test_library_names=['GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021'], # pick from list below
        top_nth=10, # set the top_nth percentile of the backgorund list to be used as the foregorund list 
        ##top_nth=10 (default) means the foregourd list is the top 10% of the background list
        )


        List of avalible gene sets to look for enrichement. These can go into  test_library_names=[] list .
          ['ARCHS4_Cell-lines',
         'ARCHS4_IDG_Coexp',
         'ARCHS4_Kinases_Coexp',
         'ARCHS4_TFs_Coexp',
         'ARCHS4_Tissues',
         'Achilles_fitness_decrease',
         'Achilles_fitness_increase',
         'Aging_Perturbations_from_GEO_down',
         'Aging_Perturbations_from_GEO_up',
         'Allen_Brain_Atlas_10x_scRNA_2021',
         'Allen_Brain_Atlas_down',
         'Allen_Brain_Atlas_up',
         'Azimuth_Cell_Types_2021',
         'BioCarta_2013',
         'BioCarta_2015',
         'BioCarta_2016',
         'BioPlanet_2019',
         'BioPlex_2017',
         'CCLE_Proteomics_2020',
         'CORUM',
         'COVID-19_Related_Gene_Sets',
         'COVID-19_Related_Gene_Sets_2021',
         'Cancer_Cell_Line_Encyclopedia',
         'CellMarker_Augmented_2021',
         'ChEA_2013',
         'ChEA_2015',
         'ChEA_2016',
         'Chromosome_Location',
         'Chromosome_Location_hg19',
         'ClinVar_2019',
         'DSigDB',
         'Data_Acquisition_Method_Most_Popular_Genes',
         'DepMap_WG_CRISPR_Screens_Broad_CellLines_2019',
         'DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019',
         'Descartes_Cell_Types_and_Tissue_2021',
         'DisGeNET',
         'Disease_Perturbations_from_GEO_down',
         'Disease_Perturbations_from_GEO_up',
         'Disease_Signatures_from_GEO_down_2014',
         'Disease_Signatures_from_GEO_up_2014',
         'DrugMatrix',
         'Drug_Perturbations_from_GEO_2014',
         'Drug_Perturbations_from_GEO_down',
         'Drug_Perturbations_from_GEO_up',
         'ENCODE_Histone_Modifications_2013',
         'ENCODE_Histone_Modifications_2015',
         'ENCODE_TF_ChIP-seq_2014',
         'ENCODE_TF_ChIP-seq_2015',
         'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
         'ESCAPE',
         'Elsevier_Pathway_Collection',
         'Enrichr_Libraries_Most_Popular_Genes',
         'Enrichr_Submissions_TF-Gene_Coocurrence',
         'Enrichr_Users_Contributed_Lists_2020',
         'Epigenomics_Roadmap_HM_ChIP-seq',
         'FANTOM6_lncRNA_KD_DEGs',
         'GO_Biological_Process_2013',
         'GO_Biological_Process_2015',
         'GO_Biological_Process_2017',
         'GO_Biological_Process_2017b',
         'GO_Biological_Process_2018',
         'GO_Biological_Process_2021',
         'GO_Cellular_Component_2013',
         'GO_Cellular_Component_2015',
         'GO_Cellular_Component_2017',
         'GO_Cellular_Component_2017b',
         'GO_Cellular_Component_2018',
         'GO_Cellular_Component_2021',
         'GO_Molecular_Function_2013',
         'GO_Molecular_Function_2015',
         'GO_Molecular_Function_2017',
         'GO_Molecular_Function_2017b',
         'GO_Molecular_Function_2018',
         'GO_Molecular_Function_2021',
         'GTEx_Aging_Signatures_2021',
         'GTEx_Tissue_Expression_Down',
         'GTEx_Tissue_Expression_Up',
         'GWAS_Catalog_2019',
         'GeneSigDB',
         'Gene_Perturbations_from_GEO_down',
         'Gene_Perturbations_from_GEO_up',
         'Genes_Associated_with_NIH_Grants',
         'Genome_Browser_PWMs',
         'HDSigDB_Human_2021',
         'HDSigDB_Mouse_2021',
         'HMDB_Metabolites',
         'HMS_LINCS_KinomeScan',
         'HomoloGene',
         'HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression',
         'HuBMAP_ASCTplusB_augmented_2022',
         'HumanCyc_2015',
         'HumanCyc_2016',
         'Human_Gene_Atlas',
         'Human_Phenotype_Ontology',
         'InterPro_Domains_2019',
         'Jensen_COMPARTMENTS',
         'Jensen_DISEASES',
         'Jensen_TISSUES',
         'KEA_2013',
         'KEA_2015',
         'KEGG_2013',
         'KEGG_2015',
         'KEGG_2016',
         'KEGG_2019_Human',
         'KEGG_2019_Mouse',
         'KEGG_2021_Human',
         'Kinase_Perturbations_from_GEO_down',
         'Kinase_Perturbations_from_GEO_up',
         'L1000_Kinase_and_GPCR_Perturbations_down',
         'L1000_Kinase_and_GPCR_Perturbations_up',
         'LINCS_L1000_Chem_Pert_down',
         'LINCS_L1000_Chem_Pert_up',
         'LINCS_L1000_Ligand_Perturbations_down',
         'LINCS_L1000_Ligand_Perturbations_up',
         'Ligand_Perturbations_from_GEO_down',
         'Ligand_Perturbations_from_GEO_up',
         'MCF7_Perturbations_from_GEO_down',
         'MCF7_Perturbations_from_GEO_up',
         'MGI_Mammalian_Phenotype_2013',
         'MGI_Mammalian_Phenotype_2017',
         'MGI_Mammalian_Phenotype_Level_3',
         'MGI_Mammalian_Phenotype_Level_4',
         'MGI_Mammalian_Phenotype_Level_4_2019',
         'MGI_Mammalian_Phenotype_Level_4_2021',
         'MSigDB_Computational',
         'MSigDB_Hallmark_2020',
         'MSigDB_Oncogenic_Signatures',
         'Microbe_Perturbations_from_GEO_down',
         'Microbe_Perturbations_from_GEO_up',
         'Mouse_Gene_Atlas',
         'NCI-60_Cancer_Cell_Lines',
         'NCI-Nature_2015',
         'NCI-Nature_2016',
         'NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions',
         'NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions',
         'NIH_Funded_PIs_2017_Human_AutoRIF',
         'NIH_Funded_PIs_2017_Human_GeneRIF',
         'NURSA_Human_Endogenous_Complexome',
         'OMIM_Disease',
         'OMIM_Expanded',
         'Old_CMAP_down',
         'Old_CMAP_up',
         'Orphanet_Augmented_2021',
         'PPI_Hub_Proteins',
         'PanglaoDB_Augmented_2021',
         'Panther_2015',
         'Panther_2016',
         'Pfam_Domains_2019',
         'Pfam_InterPro_Domains',
         'PheWeb_2019',
         'PhenGenI_Association_2021',
         'Phosphatase_Substrates_from_DEPOD',
         'ProteomicsDB_2020',
         'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO',
         'RNAseq_Automatic_GEO_Signatures_Human_Down',
         'RNAseq_Automatic_GEO_Signatures_Human_Up',
         'RNAseq_Automatic_GEO_Signatures_Mouse_Down',
         'RNAseq_Automatic_GEO_Signatures_Mouse_Up',
         'Rare_Diseases_AutoRIF_ARCHS4_Predictions',
         'Rare_Diseases_AutoRIF_Gene_Lists',
         'Rare_Diseases_GeneRIF_ARCHS4_Predictions',
         'Rare_Diseases_GeneRIF_Gene_Lists',
         'Reactome_2013',
         'Reactome_2015',
         'Reactome_2016',
         'SILAC_Phosphoproteomics',
         'SubCell_BarCode',
         'SysMyo_Muscle_Gene_Sets',
         'TF-LOF_Expression_from_GEO',
         'TF_Perturbations_Followed_by_Expression',
         'TG_GATES_2020',
         'TRANSFAC_and_JASPAR_PWMs',
         'TRRUST_Transcription_Factors_2019',
         'Table_Mining_of_CRISPR_Studies',
         'TargetScan_microRNA',
         'TargetScan_microRNA_2017',
         'Tissue_Protein_Expression_from_Human_Proteome_Map',
         'Tissue_Protein_Expression_from_ProteomicsDB',
         'Transcription_Factor_PPIs',
         'UK_Biobank_GWAS_v1',
         'Virus-Host_PPI_P-HIPSTer_2020',
         'VirusMINT',
         'Virus_Perturbations_from_GEO_down',
         'Virus_Perturbations_from_GEO_up',
         'WikiPathway_2021_Human',
         'WikiPathways_2013',
         'WikiPathways_2015',
         'WikiPathways_2016',
         'WikiPathways_2019_Human',
         'WikiPathways_2019_Mouse',
         'dbGaP',
         'huMAP',
         'lncHUB_lncRNA_Co-Expression',
         'miRTarBase_2017']
        """
        import scanpy as sc
        import pandas as pd
        import gseapy as gp
        import numpy as np
        import matplotlib.pyplot as plt
        import os
        ################## supress FutureWarning
        import warnings
        warnings.simplefilter(action='ignore', category=FutureWarning)
        ##################

        os.makedirs(self.output_dir+self.output_prefix, exist_ok=True)


        os.makedirs(self.output_dir+self.output_prefix+'/tables/', exist_ok=True)
        dataset_tables_output_directory=self.output_dir+self.output_prefix+'/tables/'

        os.makedirs(self.output_dir+self.output_prefix+'/figures/', exist_ok=True)
        dataset_figures_output_directory=self.output_dir+self.output_prefix+'/figures/'

        sc.settings.figdir=dataset_figures_output_directory


        os.makedirs(self.output_dir+self.output_prefix+"/GSEA_out/", exist_ok=True)
        dataset_GESA_output_directory=self.output_dir+self.output_prefix+"/GSEA_out/"

        #total_cluster_number=len(set(adata.obs['leiden'].tolist()))

        top_percentile=(top_nth)/100


        ############################## logical regression test GSEA
        test="logreg"
        full_table = pd.read_csv(dataset_tables_output_directory+self.output_prefix+"rank_genes_groups_"+test+".csv",header=0,index_col=0)
        total_cluster_number=len(full_table.columns) # set total cluster number to column # of loreg rank table
        background_list_len=full_table.shape[0]
        print(f'logreg: the full_table  is {full_table.shape[0]} genes long by {full_table.shape[1]} columns for {total_cluster_number} clusters')
        foreground_list_len=len((full_table[ :int(background_list_len * top_percentile)]))
        print(f'logreg: the foreground list is {foreground_list_len} genes long')

        for i in range(0, total_cluster_number):
            test_cluster_number=i

            os.makedirs(dataset_GESA_output_directory+test+"_top_"+str(top_nth)+"pct_"+"cluster_"+str(test_cluster_number), exist_ok=True)
            cluster_gsea_output_dir=dataset_GESA_output_directory+test+"_top_"+str(top_nth)+"pct_"+"cluster_"+str(test_cluster_number)

            background_list=full_table[full_table.columns[test_cluster_number]].squeeze().str.strip().tolist()
            foreground_list=(background_list[ :int(background_list_len * top_percentile)])
            print(f"<CLUSTER {test_cluster_number}> for {test} gene rank top 3 background genes {background_list[:3]}, bottom 3 background genes {background_list[-3:]} ")
            print(f"<CLUSTER {test_cluster_number}> for {test} gene rank top 3 foreground genes {foreground_list[:3]}, bottom 3 foreground genes {foreground_list[-3:]} ")
            # run enrichr
            # list, dataframe, series inputs are supported
            try:
                enr = gp.enrichr(gene_list=foreground_list,
                                 background=background_list,
                                 gene_sets=test_library_names,
                                 organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
                                 #description=test+"_top_"+str(top_nth)+"pct_"+"cluster_"+str(test_cluster_number),
                                 outdir=cluster_gsea_output_dir,
                                 # no_plot=True,
                                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                                )
            except:
                print("Something  went wrong")
        ############################## logical regression test GSEA END

        ############################## wilcox regression test GSEA
        test="wilcox"

        full_table = pd.read_csv(dataset_tables_output_directory+self.output_prefix+"rank_genes_groups_"+test+".csv",header=0,index_col=0)
        total_cluster_number=int(len(full_table.columns)/5) # set total cluster number to column # / 2 of wilcox or t test rank table
        background_list_len=full_table.shape[0]
        print(f'wilcox: the full_table  is {full_table.shape[0]} genes long by {full_table.shape[1]} columns for {total_cluster_number} clusters')
        foreground_list_len=len((full_table[ :int(background_list_len * top_percentile)]))
        print(f'wilcox: the foreground list is {foreground_list_len} genes long')

        for i in range(0, total_cluster_number):
            test_cluster_number=i

            os.makedirs(dataset_GESA_output_directory+test+"_top_"+str(top_nth)+"pct_"+"cluster_"+str(test_cluster_number), exist_ok=True)
            cluster_gsea_output_dir=dataset_GESA_output_directory+test+"_top_"+str(top_nth)+"pct_"+"cluster_"+str(test_cluster_number)

            background_list=full_table[full_table.columns[(test_cluster_number*5)]].squeeze().str.strip().tolist()
            foreground_list=(background_list[ :int(background_list_len * top_percentile)])
            print(f"<CLUSTER {test_cluster_number}> for {test} gene rank top 3 background genes {background_list[:3]}, bottom 3 background genes {background_list[-3:]} ")
            print(f"<CLUSTER {test_cluster_number}> for {test} gene rank top 3 foreground genes {foreground_list[:3]}, bottom 3 foreground genes {foreground_list[-3:]} ")
            # run enrichr
            # list, dataframe, series inputs are supported
            try:
                enr = gp.enrichr(gene_list=foreground_list,
                                 background=background_list,
                                 gene_sets=test_library_names,
                                 organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
                                 #description=test+"_top_"+str(top_nth)+"pct_"+"cluster_"+str(test_cluster_number),
                                 outdir=cluster_gsea_output_dir,
                                 # no_plot=True,
                                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                                )
            except:
                print("Something went wrong")
        ############################## wilcox regression test GSEA END


        ############################## t_test regression test GSEA
        test="t_test"

        full_table = pd.read_csv(dataset_tables_output_directory+self.output_prefix+"rank_genes_groups_"+test+".csv",header=0,index_col=0)
        total_cluster_number=int(len(full_table.columns)/5) # set total cluster number to column # / 2 of wilcox or t test rank table
        background_list_len=full_table.shape[0]
        print(f't_test: the full_table  is {full_table.shape[0]} genes long by {full_table.shape[1]} columns for {total_cluster_number} clusters')
        foreground_list_len=len((full_table[ :int(background_list_len * top_percentile)]))
        print(f't_test: the foreground list is {foreground_list_len} genes long')

        for i in range(0, total_cluster_number):
            test_cluster_number=i

            os.makedirs(dataset_GESA_output_directory+test+"_top_"+str(top_nth)+"pct_"+"cluster_"+str(test_cluster_number), exist_ok=True)
            cluster_gsea_output_dir=dataset_GESA_output_directory+test+"_top_"+str(top_nth)+"pct_"+"cluster_"+str(test_cluster_number)

            background_list=full_table[full_table.columns[(test_cluster_number*5)]].squeeze().str.strip().tolist()
            foreground_list=(background_list[ :int(background_list_len * top_percentile)])
            print(f"<CLUSTER {test_cluster_number}> for {test} gene rank top 3 background genes {background_list[:3]}, bottom 3 background genes {background_list[-3:]} ")
            print(f"<CLUSTER {test_cluster_number}> for {test} gene rank top 3 foreground genes {foreground_list[:3]}, bottom 3 foreground genes {foreground_list[-3:]} ")
            # run enrichr
            # list, dataframe, series inputs are supported
            try:
                enr = gp.enrichr(gene_list=foreground_list,
                                 background=background_list,
                                 gene_sets=test_library_names,
                                 organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
                                 #description=test+"_top_"+str(top_nth)+"pct_"+"cluster_"+str(test_cluster_number),
                                 outdir=cluster_gsea_output_dir,
                                 # no_plot=True,
                                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                                )
            except:
                print("Something  went wrong")
        ############################## t_test regression test GSEA END
        return
    
    