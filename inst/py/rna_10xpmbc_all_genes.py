#!/usr/bin/env python
# coding: utf-8

# # scRNA-seq analysis

# Here we will use 10X PBMCs scRNA-seq dataset as an example to illustrate how SIMBA performs scRNA-seq analysis

# In[1]:


import os
import simba as si
si.__version__


# In[2]:


workdir = 'result_simba_rnaseq'
si.settings.set_workdir(workdir)


# In[3]:


si.settings.set_figure_params(dpi=80,
                              style='white',
                              fig_size=[5,5],
                              rc={'image.cmap': 'viridis'})


# In[4]:


# make plots prettier
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina')


# In[ ]:





# ### load example data

# In[5]:


adata_CG = si.datasets.rna_10xpmbc3k()


# In[6]:


adata_CG


# In[7]:


adata_CG.var


# ### preprocessing

# In[8]:


#si.pp.filter_genes(adata_CG)#n_cells=3)


# In[9]:


si.pp.cal_qc_rna(adata_CG)


# In[10]:


si.pl.violin(adata_CG,list_obs=['n_counts','n_genes','pct_mt'])


# Filter out cells if needed:
# 
# ```python
# si.pp.filter_cells_rna(adata,min_n_genes=100)
# ```

# In[11]:


si.pp.normalize(adata_CG,method='lib_size')


# In[12]:


si.pp.log_transform(adata_CG)


# In[ ]:





# Optionally, variable gene selection step can be also performed. 
# 
# ```python
# si.pp.select_variable_genes(adata_CG, n_top_genes=2000)
# si.pl.variable_genes(adata_CG,show_texts=True)
# ```
# 
# This will speed up the training procedure as only variable genes are encoded into the graph. But we won't obtain the embeddings of non-variable genes.

# In[ ]:





# ### discretize RNA expression

# In[13]:


si.tl.discretize(adata_CG,n_bins=5)


# In[14]:


si.pl.discretize(adata_CG,kde=False)


# In[ ]:





# ### generate graph

# In[15]:


si.tl.gen_graph(list_CG=[adata_CG],
                copy=False,
                use_highly_variable=False,
                dirname='graph0')


# In[ ]:





# ### PBG training

# Before PBG training, let's take a look at the parameters:

# In[ ]:





# In[16]:


si.settings.pbg_params


# In[ ]:





# If no parameters need to be adjusted, the training can be simply done with:
# 
# ```python
# si.tl.pbg_train(auto_wd=True, save_wd=True, output='model')
# ```

# In[ ]:





# Here we show how to adjust training-related parameters if needed. In general, weight decay `wd` is the only parameter that might need to be adjusted based on the following pbg metric plots. However, in almost all the cases, the automatically decided `wd` (enabling it by setting `auto_wd=True`) works well.
# 
# E.g. we want to change `wd_interval`:

# In[17]:


# modify parameters
dict_config = si.settings.pbg_params.copy()
# dict_config['wd'] = 0.015521 
dict_config['wd_interval'] = 10 # we usually set `wd_interval` to 10 for scRNA-seq datasets for a slower but finer training
dict_config['workers'] = 12 #The number of CPUs.

## start training
si.tl.pbg_train(pbg_params = dict_config, auto_wd=True, save_wd=True, output='model')


# In[ ]:





# > If `wd` is specified by users instead of being automatically decided, then make sure to update it in simba setting:
# ```python
# si.settings.pbg_params = dict_config.copy()
# ```

# In[ ]:





# The trained result can be loaded in with the following steps:
# 
# By default, it's using the current training result stored in `.setting.pbg_params`
# ```python
# # load in graph ('graph0') info
# si.load_graph_stats()
# # load in model info for ('graph0')
# si.load_pbg_config()
# ```
# Users can also specify different pathss
# ```python
# # load in graph ('graph0') info
# si.load_graph_stats(path='./result_simba_rnaseq/pbg/graph0/')
# # load in model info for ('graph0')
# si.load_pbg_config(path='./result_simba_rnaseq/pbg/graph0/model/')
# ```

# In[ ]:





# plotting training metrics to make sure the model is not overfitting

# In[18]:


si.pl.pbg_metrics(fig_ncol=1)


# In[ ]:





# ### post-training analysis

# In[19]:


palette_celltype={'B':'#1f77b4',
                  'CD4 T':'#ff7f0e', 
                  'CD8 T':'#279e68',
                  'Dendritic':"#aa40fc",
                  'CD14 Monocytes':'#d62728',
                  'FCGR3A Monocytes':'#b5bd61',
                  'Megakaryocytes':'#e377c2',
                  'NK':'#8c564b'}


# In[20]:


dict_adata = si.read_embedding()


# In[21]:


dict_adata


# In[22]:


adata_C = dict_adata['C']  # embeddings for cells
adata_G = dict_adata['G']  # embeddings for genes


# In[23]:


adata_C


# In[24]:


adata_G


# In[ ]:





# visualize embeddings of cells

# In[25]:


## Add annotation of celltypes (optional)
adata_C.obs['celltype'] = adata_CG[adata_C.obs_names,:].obs['celltype'].copy()
adata_C


# In[26]:


si.tl.umap(adata_C,n_neighbors=15,n_components=2)


# In[27]:


si.pl.umap(adata_C,color=['celltype'],
           dict_palette={'celltype': palette_celltype},
           fig_size=(6,4),
           drawing_order='random')


# In[ ]:





# visualize embeddings of cells and genes

# In[28]:


# embed cells and genes into the same space
adata_all = si.tl.embed(adata_ref=adata_C,list_adata_query=[adata_G])


# In[29]:


adata_all.obs.head()


# In[30]:


## add annotations of cells and genes
adata_all.obs['entity_anno'] = ""
adata_all.obs.loc[adata_C.obs_names, 'entity_anno'] = adata_all.obs.loc[adata_C.obs_names, 'celltype']
adata_all.obs.loc[adata_G.obs_names, 'entity_anno'] = 'gene'


# In[31]:


adata_all.obs.head()


# In[32]:


si.tl.umap(adata_all,n_neighbors=15,n_components=2)


# In[33]:


palette_entity_anno = palette_celltype.copy()
palette_entity_anno['gene'] = "#607e95"


# In[34]:


si.pl.umap(adata_all,color=['id_dataset','entity_anno'],
           dict_palette={'entity_anno': palette_entity_anno},
           drawing_order='original',
           fig_size=(6,5))


# In[ ]:





# In[35]:


marker_genes_ini = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']


# In[36]:


marker_genes = ['IL7R']#, 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
    #            'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
       #         'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']


# In[37]:


si.pl.umap(adata_all[::-1,],color=['entity_anno'],dict_palette={'entity_anno': palette_entity_anno},
           drawing_order='original',
           texts=marker_genes, #+ ['GAPDH', 'B2M'],
           show_texts=True,
           fig_size=(8,6))


# In[38]:


marker_genes


# SIMBA metrics

# In[39]:


adata_cmp = si.tl.compare_entities(adata_ref=adata_C, adata_query=adata_G)


# In[40]:


adata_cmp


# In[41]:


si.pl.entity_metrics(adata_cmp,
                     x='max',
                     y='gini',
                     show_contour=False,
                     texts=marker_genes,# + ['GAPDH', 'B2M'],
                     show_texts=True,
                     show_cutoff=True,
                     size=5,
                     text_expand=(1.3,1.5),
                     cutoff_x=1.,
                     cutoff_y=0.3,
                     save_fig=False)


# In[ ]:





# SIMBA barcode plots

# In[42]:


# add annoations of cells
adata_cmp.obs['celltype'] = adata_CG.obs.loc[adata_cmp.obs_names,'celltype']


# In[43]:


list_genes = ['CST3', 'NKG7', 'MS4A1', 'GAPDH']


# In[44]:


si.pl.entity_barcode(adata_cmp, 
                     layer='softmax',
                     entities=list_genes, 
                     anno_ref='celltype',
                     show_cutoff=True,
                     cutoff=0.001,
                     palette=palette_celltype,
                     fig_size=(6, 2.5),
                     save_fig=False)


# In[ ]:





# visualize genes on umap of cells

# In[ ]:


adata_CG.obsm['X_umap'] = adata_C[adata_CG.obs_names,].obsm['X_umap'].copy()
si.pl.umap(adata_CG,
           color=['CST3', 'NKG7', 'MS4A1', 'GAPDH'],
           drawing_order='sorted',
           size=5,
           alpha=0.9,
           fig_ncol=2,
           save_fig=False)


# In[ ]:





# SIMBA queries

# In[ ]:


si.pl.umap(adata_all,
           color=['entity_anno'],dict_palette={'entity_anno': palette_entity_anno},
           drawing_order='random',
           show_texts=False,
           fig_size=(7,5))


# In[ ]:





# In[ ]:


# find neighbor genes around the location [6, 16] on UMAP
query_result = si.tl.query(adata_all,
                           pin=[6,16],
                           obsm='X_umap',
                           use_radius=True,r=2,
                           anno_filter='entity_anno',
                           filters=['gene'])
print(query_result.shape)
query_result.head()


# In[ ]:


# show locations of pin point and its neighbor genes 
si.pl.query(adata_all,
            show_texts=False, 
            color=['entity_anno'], dict_palette={'entity_anno': palette_entity_anno},
            alpha=0.9,
            alpha_bg=0.1,
            fig_size=(7,5))


# In[ ]:





# In[ ]:


# find top 50 neighbor genes around cell "ACTCAGGATTCGTT-1" (CD14 Monocytes) in SIMBA space
query_result = si.tl.query(adata_all,
                           entity=['ACTCAGGATTCGTT-1'],
                           obsm=None,
                           use_radius=False,
                           k=50,
                           anno_filter='entity_anno',
                           filters=['gene'])
print(query_result.shape)
query_result.head()


# In[ ]:


# show locations of entity and its neighbor genes 
si.pl.query(adata_all,
            obsm='X_umap',
            color=['entity_anno'], dict_palette={'entity_anno': palette_entity_anno},
            alpha=0.9,
            alpha_bg=0.1,
            fig_size=(7,5))


# In[ ]:





# In[ ]:


# find top 50 neighbor genes for multiples cells in SIMBA space
query_result = si.tl.query(adata_all,entity=['GATGCCCTCTCATT-1', 'CTGAAGTGGCTATG-1'],
                           obsm=None,
                           use_radius=False,
                           k=50,
                           anno_filter='entity_anno',
                           filters=['gene'],
                           )
print(query_result.shape)
query_result.head()


# In[ ]:


# show locations of entities and their neighbor genes 
si.pl.query(adata_all,
            obsm='X_umap',
            alpha=0.9,
            alpha_bg=0.1,            
            fig_size=(5,5))


# In[ ]:





# In[ ]:


# find neighbor entities (both cells and genes) of a given gene on UMAP
query_result = si.tl.query(adata_all,
                           entity=['CD79A'],
                           obsm='X_umap',
                           use_radius=False,
                           k=50
                           )
print(query_result.shape)
query_result.iloc[:10,]


# In[ ]:


si.pl.query(adata_all,
            obsm='X_umap',
            color=['entity_anno'], dict_palette={'entity_anno': palette_entity_anno},
            show_texts=False,
            alpha=0.9,
            alpha_bg=0.1,
            fig_size=(7,5))


# In[ ]:


print(workdir)


# ### save results

# In[ ]:


adata_CG.write(os.path.join(workdir, 'adata_CG.h5ad'))
adata_C.write(os.path.join(workdir, 'adata_C.h5ad'))
adata_G.write(os.path.join(workdir, 'adata_G.h5ad'))
adata_all.write(os.path.join(workdir, 'adata_all.h5ad'))
adata_cmp.write(os.path.join(workdir, 'adata_cmp.h5ad'))


# In[ ]:





# Read back anndata objects
# ```python
# adata_CG = si.read_h5ad(os.path.join(workdir, 'adata_CG.h5ad'))
# adata_C = si.read_h5ad(os.path.join(workdir, 'adata_C.h5ad'))
# adata_G = si.read_h5ad(os.path.join(workdir, 'adata_G.h5ad'))
# adata_all = si.read_h5ad(os.path.join(workdir, 'adata_all.h5ad'))
# adata_cmp = si.read_h5ad(os.path.join(workdir, 'adata_cmp.h5ad'))
# ```

# In[ ]:


get_ipython().system('pip install sinfo')


# In[ ]:


from sinfo import sinfo


# In[ ]:


sinfo()


# In[ ]:


help(sinfo)


# In[ ]:


sinfo(dependencies=True)


# In[ ]:




