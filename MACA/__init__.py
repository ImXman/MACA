# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:31:17 2020

@author: Yang Xu
"""
import umap
import scipy
#import random
import anndata
import collections
import cosg as cosg
import numpy as np
import pandas as pd
import scanpy as sc
import multiprocessing

from sklearn.decomposition import PCA
from sklearn.metrics import confusion_matrix
#from sklearn.preprocessing import Normalizer
from sklearn.feature_extraction.text import TfidfTransformer

import warnings
warnings.filterwarnings("ignore")

##-----------------------------------------------------------------------------
def ensemble_labels(multi_labels=None):
    ensemble = []
    for i in range(multi_labels.shape[0]):
        ks=[]
        vs=[]
        for k,v in collections.Counter(multi_labels[i,:].tolist()).items():
            ks.append(k)
            vs.append(v)
        ensemble.append(ks[vs.index(max(vs))])
    return ensemble

def singleMACA(ad=None, cell_markers=None,n_pcs=None,res=[1,2,3],n_neis = [5,10],
               freq=0.5,use_weight=False):
    ##TF-IDF transformation
    X = ad.X.copy()
    if scipy.sparse.issparse(X):
        X = X.todense()
        
    tf_transformer = TfidfTransformer(use_idf=True).fit(X)
    X= tf_transformer.transform(X).todense()
    
    labels = pd.DataFrame(np.zeros((X.shape[0],len(cell_markers))))
    labels.columns = cell_markers.keys()
    exprsed = pd.DataFrame(np.zeros((X.shape[0],len(cell_markers))))
    exprsed.columns = cell_markers.keys()
    celltype_size = {}
    
    ##create artifical labels for each cell
    if use_weight == True:
        for k, v in cell_markers.items():
            celltype_size[k]=0
            sums=0
            n = np.zeros((X.shape[0]))
            marker_index = -1
            for i in v:
                marker_index += 1
                if i in ad.var.index:
                    expr95 = np.percentile(X[:,ad.var.index == i],95)
                    thresh = 0.25 * expr95
                    l = np.array(X[:,ad.var.index == i])
                    l[X[:,ad.var.index == i]<=thresh]=0
                    
                    ##consider marker weight
                    l = l*(1-marker_index/(len(v)*2))##default 2
                    
                    n[np.array(l>0).reshape(X.shape[0])] += 1
                    sums += 1
                    labels[k] += l.reshape(X.shape[0])
            n = n/sums
            celltype_size[k]=sums
            exprsed[k] = n.reshape(X.shape[0]) 

    else:
        for k, v in cell_markers.items():
            celltype_size[k]=0
            sums=0
            n = np.zeros((X.shape[0]))
            for i in v:
                if i in ad.var.index:
                    expr95 = np.percentile(X[:,ad.var.index == i],95)
                    thresh = 0.25 * expr95
                    l = np.array(X[:,ad.var.index == i])
                    l[X[:,ad.var.index == i]<=thresh]=0
                    n[np.array(l>0).reshape(X.shape[0])] += 1
                    sums += 1
                    labels[k] += l.reshape(X.shape[0])
            n = n/sums
            celltype_size[k]=sums
            exprsed[k] = n.reshape(X.shape[0])   
    
    assess1 = np.argmax((labels*exprsed).values,axis=1)
    vals1 = 0
    for k,v in collections.Counter(assess1).items():
        if v >= 5:
            vals1 += 1
                        
    assess1 = vals1
    
    assess2 = np.argmax((labels).values,axis=1)
    vals2 = 0
    for k,v in collections.Counter(assess2).items():
        if v >= 5:
            vals2 += 1
                        
    assess2 = vals2
    
    ass = [assess1,assess2]
    
    new_labels = [labels*exprsed,labels][ass.index(max(ass))]
    #new_labels = labels*exprsed##consider the number of expressed marker of each cell-type for each cell
    #new_labels = labels#*exprsed
    
    celltype_size = pd.DataFrame.from_dict(celltype_size,orient='index')

    ##remove cell types with > 300 and < 5 marker genes
    labels = labels.iloc[:,np.logical_and(celltype_size.values<=300,celltype_size.values>=3)]
    new_labels = new_labels.iloc[:,np.logical_and(celltype_size.values<=300,celltype_size.values>=3)]
    celltype_size = celltype_size.iloc[np.logical_and(celltype_size.values<=300,celltype_size.values>=3),:]
    print(labels.shape)
    
    ##l2 normalize
    #transformer = Normalizer().fit(labels)
    #l2 = transformer.transform(labels)
    
    ##create Label1
    labels1 = np.argmax(new_labels.values,axis=1)
    ad.obs['Label1'] = labels1
    
    ##UMAP visualization
    embedding = umap.UMAP(n_neighbors=15,min_dist=0.2,n_components=2,
                          metric='cosine').fit_transform(labels.values)#(l2)#
    embedding = pd.DataFrame(embedding)
    embedding.columns=['UMAP1','UMAP2']
    ad.obsm['X_umap'] = embedding.iloc[:,:2].values
    
    ##create Label2
    if n_pcs is not None:
        pca = PCA(n_components=n_pcs)
        scores = pca.fit_transform(labels.values)
        ad.obsm['Score']=scores
    else:
        ad.obsm['Score']=labels.values#l2#
    
    label_list = np.zeros((ad.X.shape[0],len(res)*len(n_neis))).astype('str')
    indexs = 0
    for r in res:
        for nei in n_neis:
            
            sc.pp.neighbors(ad, use_rep="Score", n_neighbors=nei,metric='cosine')
            sc.tl.louvain(ad, resolution=r, key_added = 'louvain')
    
            ##Map Label2 to Label1; remove all non-candiate cell types
            labels1_2 = list(collections.Counter(labels1.tolist()).keys())
            labels_2 = labels.copy().iloc[:,labels1_2]
            new_labels_2 = new_labels.copy().iloc[:,labels1_2]
            #celltype_size_2 = celltype_size.iloc[labels1_2,:].copy()
            celltype_names = new_labels_2.columns
            
            cm = confusion_matrix(ad.obs['louvain'].values.astype(int),
                                  np.argmax(new_labels_2.values,axis=1).astype(int))
            cm = cm[:np.max(ad.obs['louvain'].values.astype(int)),:]
            
            normed_cm = cm.copy().T
            normed_cm = normed_cm/np.sum(normed_cm,axis=0)
            normed_cm = np.nan_to_num(normed_cm)
            normed_cm = normed_cm.T
            mapping={}
            
            mapping = np.argmax(normed_cm,axis=1)
            mapmax = np.max(normed_cm,axis=1)
            mapmax = np.nan_to_num(mapmax)
            
            clustering = ad.obs['louvain'].values.astype(int)
            new_cluster = np.zeros((len(clustering)))
            for i in range(len(mapping)):
                tof = clustering==i
                if mapmax[i]>=freq:
                    new_cluster[tof]=mapping[i]
                else:
                    sub_labels = new_labels_2.values[tof,:]
                    if sub_labels.shape[0]>0:
                        freqs = []
                        for j in range(sub_labels.shape[0]):
                            zscore=scipy.stats.zscore(sub_labels[j,:])
                            orderi = np.array([g for g in range(sub_labels.shape[1])])[zscore>3].tolist()
                            orderj = np.argsort(sub_labels[j,:])[::-1][:3].tolist()#default 3
                            a = [len(orderi),len(orderj)]
                            freqs += [orderi,orderj][a.index(max(a))]
                        vals = 0
                        ks = 0
                        for k,v in collections.Counter(freqs).items():
                            if v >= vals:
                                ks = k
                                vals = v
                        
                        if vals/sub_labels.shape[0]>=freq:
                            new_cluster[tof]=ks
                        else:
                            new_cluster[tof]=-1
    
            mapped =new_cluster.astype('int')
            mapped=mapped.astype('str')
            ad.obs['Mapped'] = mapped
            
            cell_dict = {}
            for k,v in collections.Counter(new_cluster.tolist()).items():
                if int(k)>=0:
                    cell_dict[int(k)]=labels_2.columns[int(k)]
                else:
                    cell_dict[int(k)]="unassigned"
            label_list[:,indexs]=ad.obs['Mapped'].values
            indexs+=1
    
    ensemble = ensemble_labels(label_list)
    ensemble = np.array(ensemble)
    
    annotations=[]
    for e in ensemble:
        if int(e)>=0:
            annotations.append(celltype_names[int(e)])
        else:
            annotations.append("unassigned")

    return ad, np.array(annotations)
    
def gene2cell(ad=None, cell_markers=None,use_weight=False):
    ##TF-IDF transformation
    X = ad.X.copy()
    if scipy.sparse.issparse(X):
        X = X.todense()
        
    tf_transformer = TfidfTransformer(use_idf=True).fit(X)
    X= tf_transformer.transform(X).todense()
    
    labels = pd.DataFrame(np.zeros((X.shape[0],len(cell_markers))))
    labels.columns = cell_markers.keys()
    exprsed = pd.DataFrame(np.zeros((X.shape[0],len(cell_markers))))
    exprsed.columns = cell_markers.keys()
    celltype_size = {}
    
    ##create artifical labels for each cell
    if use_weight == True:
        for k, v in cell_markers.items():
            celltype_size[k]=0
            sums=0
            n = np.zeros((X.shape[0]))
            marker_index = -1
            for i in v:
                marker_index += 1
                if i in ad.var.index:
                    expr95 = np.percentile(X[:,ad.var.index == i],95)
                    thresh = 0.25 * expr95
                    l = np.array(X[:,ad.var.index == i])
                    l[X[:,ad.var.index == i]<=thresh]=0
                    
                    ##consider marker weight
                    l = l*(1-marker_index/(len(v)*2))##default 2
                    
                    n[np.array(l>0).reshape(X.shape[0])] += 1
                    sums += 1
                    labels[k] += l.reshape(X.shape[0])
            n = n/sums
            celltype_size[k]=sums
            exprsed[k] = n.reshape(X.shape[0]) 

    else:
        for k, v in cell_markers.items():
            celltype_size[k]=0
            sums=0
            n = np.zeros((X.shape[0]))
            for i in v:
                if i in ad.var.index:
                    expr95 = np.percentile(X[:,ad.var.index == i],95)
                    thresh = 0.25 * expr95
                    l = np.array(X[:,ad.var.index == i])
                    l[X[:,ad.var.index == i]<=thresh]=0
                    n[np.array(l>0).reshape(X.shape[0])] += 1
                    sums += 1
                    labels[k] += l.reshape(X.shape[0])
            n = n/sums
            celltype_size[k]=sums
            exprsed[k] = n.reshape(X.shape[0])        
    
    assess1 = np.argmax((labels*exprsed).values,axis=1)
    vals1 = 0
    for k,v in collections.Counter(assess1).items():
        if v >= 5:
            vals1 += 1
                        
    assess1 = vals1
    
    assess2 = np.argmax((labels).values,axis=1)
    vals2 = 0
    for k,v in collections.Counter(assess2).items():
        if v >= 5:
            vals2 += 1
                        
    assess2 = vals2
    
    ass = [assess1,assess2]
    
    new_labels = [labels*exprsed,labels][ass.index(max(ass))]
    #new_labels = labels*exprsed##consider the number of expressed marker of each cell-type for each cell
    #new_labels = labels#*exprsed
    
    celltype_size = pd.DataFrame.from_dict(celltype_size,orient='index')

    ##remove cell types with > 300 and < 5 marker genes
    labels = labels.iloc[:,np.logical_and(celltype_size.values<=300,celltype_size.values>=3)]
    new_labels = new_labels.iloc[:,np.logical_and(celltype_size.values<=300,celltype_size.values>=3)]
    celltype_size = celltype_size.iloc[np.logical_and(celltype_size.values<=300,celltype_size.values>=3),:]
    print(labels.shape)
    
    return labels, new_labels

def multiMACA(labels=None,new_labels=None):
    
    ad = anndata.AnnData(X=labels)
    ##create Label1
    labels1 = np.argmax(new_labels.values,axis=1)
    ad.obs['Label1'] = labels1
    
    ##create Label2
    #if n_pcs is not None:
    #    pca = PCA(n_components=n_pcs)
    #    scores = pca.fit_transform(labels.values)
    #    ad.obsm['Score']=scores
    #else:
    #    ad.obsm['Score']=labels.values
    
    ##l2 normalize
    #transformer = Normalizer().fit(new_labels)
    #l2 = transformer.transform(new_labels)
    
    ad.obsm['Score']=labels.values#l2#
    sc.pp.neighbors(ad, use_rep="Score", n_neighbors=5,metric='cosine')
    sc.tl.louvain(ad, resolution=2, key_added = 'louvain')##default 1

    ##Map Label2 to Label1; remove all non-candiate cell types
    labels1_2 = list(collections.Counter(labels1.tolist()).keys())
    #labels_2 = labels.copy().iloc[:,labels1_2]
    new_labels_2 = new_labels.copy().iloc[:,labels1_2]
    #celltype_size_2 = celltype_size.iloc[labels1_2,:].copy()
    celltype_names = new_labels_2.columns
    
    cm = confusion_matrix(ad.obs['louvain'].values.astype(int),
                          np.argmax(new_labels_2.values,axis=1).astype(int))
    cm = cm[:np.max(ad.obs['louvain'].values.astype(int)),:]
    
    normed_cm = cm.copy().T
    normed_cm = normed_cm/np.sum(normed_cm,axis=0)
    normed_cm = np.nan_to_num(normed_cm)
    normed_cm = normed_cm.T
    mapping={}
    
    mapping = np.argmax(normed_cm,axis=1)
    mapmax = np.max(normed_cm,axis=1)
    mapmax = np.nan_to_num(mapmax)
    
    clustering = ad.obs['louvain'].values.astype(int)
    new_cluster = np.zeros((len(clustering)))
    for i in range(len(mapping)):
        tof = clustering==i
        if mapmax[i]>=0.5:
            new_cluster[tof]=mapping[i]
        else:
            sub_labels = new_labels_2.values[tof,:]
            if sub_labels.shape[0]>0:
                freqs = []
                for j in range(sub_labels.shape[0]):
                    zscore=scipy.stats.zscore(sub_labels[j,:])
                    orderi = np.array([g for g in range(sub_labels.shape[1])])[zscore>3].tolist()
                    orderj = np.argsort(sub_labels[j,:])[::-1][:3].tolist()##default 3
                    a = [len(orderi),len(orderj)]
                    freqs += [orderi,orderj][a.index(max(a))]
                    vals = 0
                    ks = 0
                    for k,v in collections.Counter(freqs).items():
                        if v >= vals:
                            ks = k
                            vals = v
                    
                    if vals/sub_labels.shape[0]>=0.5:
                        new_cluster[tof]=ks
                    else:
                        new_cluster[tof]=-1
        
        mapped=[]
        for e in new_cluster:
            if int(e)>=0:
                mapped.append(celltype_names[int(e)])
            else:
                mapped.append("unassigned")
        
    return np.array(mapped)

def parallel(scores=None,labels=None,batch_size=20000,repeats=9,n_core=10):
    index = np.array([i for i in range(scores.shape[0])])
    label_list = np.zeros((scores.shape[0],repeats)).astype('str')
    
    for i in range(repeats):
        r = np.random.permutation(scores.shape[0])
        r_index = index[r]
        r_scores = scores.iloc[r,:]
        r_label1 = labels.iloc[r,:]
        scores_list= []
        for j in range(scores.shape[0]//batch_size+1):
            scores_list.append((r_scores.iloc[j*batch_size:(j+1)*batch_size,:],
                                r_label1.iloc[j*batch_size:(j+1)*batch_size,:]))
            
        pool = multiprocessing.Pool(processes=n_core)
        mapped = pool.starmap(multiMACA, scores_list)
        merged = []
        for m in mapped:
            merged += m.tolist()
        merged = np.array(merged)
        merged = merged[r_index.argsort()]
        label_list[:,i]=merged
        
    return label_list

##-----------------------------------------------------------------------------
##new updates
def downsample_to_smallest_category(adata,column="cell_type",random_state=None,
                                    min_cells=15,keep_small_categories=False):
    
    counts = adata.obs[column].value_counts(sort=False)
    min_size = min(counts[counts >= min_cells])
    sample_selection = None
    for sample, num_cells in counts.items():
        if num_cells <= min_cells:
            if keep_small_categories:
                sel = adata.obs.index.isin(
                    adata.obs[adata.obs[column] == sample].index)
            else:
                continue
        else:
            sel = adata.obs.index.isin(
                adata.obs[adata.obs[column] == sample]
                .sample(min_size, random_state=random_state)
                .index
            )
        if sample_selection is None:
            sample_selection = sel
        else:
            sample_selection |= sel
    return adata[sample_selection].copy()

def marker_identification(source_data=None,diff_method="wilcoxon",
                          repeats=10,num_sub_cells=500):
    
    cell_markers_freq ={}
    cell_markers_rank ={}
    all_cell_types = list(set(source_data.obs['cell_type']))
    for celltype in all_cell_types:
        cell_markers_freq[celltype] = {}
        cell_markers_rank[celltype] = {}
    
    cell_markers = []
    for i in range(repeats):#default 50
        subsample = source_data.copy()
        subsample = downsample_to_smallest_category(subsample, 'cell_type', min_cells=num_sub_cells, 
                                                    keep_small_categories=True)#default 500
        if diff_method == 'cosg':
            cosg.cosg(subsample,key_added='cosg',mu=1,n_genes_user=50,groupby='cell_type')
            cellmarkers = pd.DataFrame(subsample.uns['cosg']['names']).iloc[:50,:]#default 50
        else:
            sc.tl.rank_genes_groups(subsample, 'cell_type', method=diff_method)#,tie_correct=True)##Seurat wilcoxon
            cellmarkers = pd.DataFrame(subsample.uns['rank_genes_groups']['names']).iloc[:50,:]#default 50
        cell_markers.append(cellmarkers)
        
    for i in range(repeats):
        for celltype in cell_markers[i].columns:
            markers = cell_markers[i][celltype].values.tolist()
            for m in markers:
                if m not in cell_markers_freq[celltype].keys():
                    cell_markers_freq[celltype][m]=1
                    cell_markers_rank[celltype][m]=markers.index(m)
                else:
                    cell_markers_freq[celltype][m]+=1
                    cell_markers_rank[celltype][m]+=markers.index(m)
    
    cell_markers = {}
    for k,v in cell_markers_freq.items():
        freq = pd.DataFrame.from_dict(v,orient='index')
        rank = pd.DataFrame.from_dict(cell_markers_rank[k],orient='index')
        fr = pd.concat((freq,rank),1)
        fr.columns =["freq","rank"]
        fr = fr[fr["freq"]==repeats]
        fr = fr.sort_values(by=['rank'])
        cell_markers[k]=fr.index.tolist()
        
    num_markers=[]
    for k,v in cell_markers.items():
        num_markers.append(len(v))
    min_num = min(num_markers)
    for k,v in cell_markers.items():
        cell_markers[k]=v[:min_num]
    
    return cell_markers
