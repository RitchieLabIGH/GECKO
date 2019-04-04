'''
Mainfold generation from sample data (type input text file matrix) generate figure and json data
@agr1 : comon start path for AG folders, or list of them separte by a coma

'''

from time import time

import json,os
import numpy as np
import pandas as pd
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
else:
    try:
        import tkinter
    except ImportError:
        mpl.use('Agg')
import sys
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter

from sklearn import manifold
from sklearn.preprocessing import normalize
from sklearn import metrics
from sklearn import decomposition

'''
mainfold_organism  : generate mainfold result for each individuals
Arguments:
    - file_path : path to input text matrix
    - n_neighbors : n_neighbors parameter for concern mainfold technics ( recommend value 30)
    - maxkmer : number max of kmer to analyze ( recommend : number of km by individual * number of individuals)
    - outerfilepath : path to output analyze
    - nbkm : number of kmer by individual
    
'''
def mainfold_organism(file_path,n_neighbors=30,maxkmer=1000,outerfilepath="",nbkm=10):
    # Variables for manifold learning.

    if os.stat(file_path).st_size == 0:
        print(file_path+" is empty")
        return
    Y = np.genfromtxt(file_path,delimiter="\t")
    labels = np.genfromtxt(file_path, delimiter="\t", dtype='unicode')
    labels = labels[ 1,1:].T
    startkmer=2
    i=1
    while (startkmer+nbkm-1)<np.min([Y.shape[0],maxkmer+2]):
        X=Y[startkmer:startkmer+nbkm,1:].T
        mainfold(X,labels,file_path+"_orga_"+str(i), n_neighbors, outerfilepath,i)
        startkmer+=nbkm
        i+=1


'''
mainfold_organism  : generate mainfold result for each individuals
Arguments:
    - file_path : path to input text matrix
    - n_neighbors : n_neighbors parameter for concern mainfold technics ( recommend value 30)
    - maxkmer : number max of kmer to analyze
    - outerfilepath : path to output analyze
    - startkmer : kmer indice to start (base 1)

'''
def mainfoldtest(file_path,n_neighbors,maxkmer=1000,outerfilepath="",startkmer=1):
    # Variables for manifold learning.

    if os.stat(file_path).st_size == 0:
        print(file_path+" is empty")
        return
    X = np.genfromtxt(file_path,delimiter="\t")
    if len(X.shape)==2:
        labels = np.genfromtxt(file_path, delimiter="\t", dtype='unicode')
        labels = labels[1,1:].T
        X=X[startkmer+1:np.min([X.shape[1],maxkmer+1]),1:].T

        if X.shape[1]>=2:
            mainfold(X, labels,file_path,n_neighbors, outerfilepath)
        else:
            print("not enought kmer dimension < 2 "+file_path)
    else:
        print("not enought kmer, no kmer " + file_path)




'''
old version of function than mainfoldtest before input file translation
'''
def mainfoldtest_beforetranslation(file_path,n_neighbors,maxkmer=1000,outerfilepath="",startkmer=1):
    # Variables for manifold learning.

    if os.stat(file_path).st_size == 0:
        print(file_path+" is empty")
        return
    X = np.genfromtxt(file_path,delimiter=",")
    if len(X.shape)==2:
        labels = np.genfromtxt(file_path, delimiter=",", dtype='unicode')
        labels = labels[1:, 0]
        X=X[1:,startkmer:np.min([X.shape[1],maxkmer+1])]
        if X.shape[1]>=2:
            mainfold(X, labels,file_path,n_neighbors, outerfilepath)
        else:
            print("not enought kmer dimension < 2 "+file_path)
    else:
        print("not enought kmer, no kmer " + file_path)

'''
old version of function than mainfold_organism before input file translation
'''

def mainfold_organism_beforetranslation(file_path,n_neighbors,maxkmer=1000,outerfilepath="",nbkm=10):
    # Variables for manifold learning.

    if os.stat(file_path).st_size == 0:
        print(file_path+" is empty")
        return
    Y = np.genfromtxt(file_path,delimiter=",")
    labels = np.genfromtxt(file_path, delimiter=",", dtype='unicode')
    labels = labels[1:, 0]
    startkmer=1
    ######
    ###wrong nbkm remove nborga
    i=1
    while (startkmer+nbkm-1)<Y.shape[1]:
        X=Y[1:,startkmer:startkmer+nbkm]
        mainfold(X,labels,file_path+"_orga_"+str(i), n_neighbors, outerfilepath,i)
        startkmer+=nbkm
        i+=1


'''
mainfold  : mainfold analyze 
Arguments:
    - X : count matrix
    - file_path : path to input text matrix
    - n_neighbors : n_neighbors parameter for concern mainfold technics ( recommend value 30)
    - outerfilepath : path to output analyze
    - organ : individuals number ( for output label)

'''
def mainfold(X,labels,file_path, n_neighbors,  outerfilepath="",organ=0):
    restxt = open(file_path + "MainfoldNeig" + str(n_neighbors) + ".txt", 'w')
    print("compute mainfolds : "+file_path + "MainfoldNeig" + str(n_neighbors) + ".txt")
    X=normalize(X, axis=1, copy=False)

    initial_data_y = labels
    num_tmp_data_y = np.zeros(len(initial_data_y))
    cat = np.unique(initial_data_y)
    i = 0
    for grptmp in initial_data_y:
        num_tmp_data_y[i] = np.where(cat == grptmp)[0][0] + 1
        i += 1
    colors = num_tmp_data_y.astype(int)
    if (outerfilepath != "" )& os.path.exists(outerfilepath)& os.path.isfile(outerfilepath + "/_outerSelected.csv"):
        outter = np.genfromtxt(outerfilepath + "/_outerSelected.csv", delimiter=",")*-1
        outter = outter==1
    else:
        outter=np.zeros(colors.shape,dtype=bool)
        #colors = colors * ((outter+1)*-1)

    fig = plt.figure(figsize=(15, 8))
    plt.suptitle("Manifold Learning with %i points, %i neighbors\n%s"
                 % (X.shape[1], n_neighbors,file_path), fontsize=14)


    # Perform Locally Linear Embedding Manifold learning
    methods = ['standard', 'ltsa']#, 'modified', 'hessian'
    labels = ['LLE', 'LTSA']#, 'Modified LLE', 'Hessian LLE'
    
    for i, method in enumerate(methods):

        t0 = time()
        try:
            trans_data = manifold.LocallyLinearEmbedding(n_neighbors, 2,method=method).fit_transform(X).T
       

            t1 = time()
            # print("%s: %.2g sec" % (methods[i], t1 - t0))
    
            ax = fig.add_subplot(241 + i)
            if np.isnan(trans_data[0,0])==False:
                plt.scatter(trans_data[0,np.logical_not(outter)], trans_data[1,np.logical_not(outter)], c=colors[np.logical_not(outter)], cmap=plt.cm.rainbow, alpha=0.7)
                plt.scatter(trans_data[0,outter], trans_data[1,outter],marker='^', c=colors[outter], cmap=plt.cm.rainbow, alpha=0.7)
    
                score = metrics.silhouette_score(trans_data.T, initial_data_y, metric='euclidean')
                score2 = metrics.calinski_harabaz_score(trans_data.T, initial_data_y)
                # plt.title("t-SNE (%.2g sec), score= %.2g" % (t1 - t0,score))
                plt.title("%s, score=%.2g/%.2g" % (labels[i], score,score2))
    
                print("%ssilhouette,%g" % (labels[i],score),file=restxt)
                print("%scalinski,%g" % (labels[i],score2), file=restxt)
            ax.xaxis.set_major_formatter(NullFormatter())
            ax.yaxis.set_major_formatter(NullFormatter())
            plt.axis('tight')
        except OSError as err:
            print("OS error: {0}".format(err))
        except ValueError:
            print("Could not perform"+method+".")
        except:
            print("Unexpected error:", sys.exc_info()[0])



    # Perform PCA.
    try:
        t0 = time()
        pca = decomposition.PCA(n_components=2)
        trans_data = pca.fit_transform(X).T
        explained_variance_ratio_ = pca.explained_variance_ratio_;
        t1 = time()
        # print("%s: %.2g sec" % ('ISO', t1 - t0))

        ax = fig.add_subplot(243)
        plt.scatter(trans_data[0, np.logical_not(outter)], trans_data[1, np.logical_not(outter)],
                    c=colors[np.logical_not(outter)], cmap=plt.cm.rainbow, alpha=0.7)
        plt.scatter(trans_data[0, outter], trans_data[1, outter], marker='^', c=colors[outter], cmap=plt.cm.rainbow,
                    alpha=0.7)
        score = metrics.silhouette_score(trans_data.T, initial_data_y, metric='euclidean')
        score2 = metrics.calinski_harabaz_score(trans_data.T, initial_data_y)
        plt.title("%s, score= %.2g/%.2g" % ('PCA', score, score2))
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        plt.axis('tight')
    except OSError as err:
        print("OS error: {0}".format(err))

    except ValueError:
        print("Could not perform PCA.")
    except:
        print("Unexpected error:", sys.exc_info()[0])

    ax = fig.add_subplot(244)
   # plt.bar(range(1,len(explained_variance_ratio_)+1) , explained_variance_ratio_)
    plt.bar(1, explained_variance_ratio_.sum(),label="Y")
    plt.bar(1, explained_variance_ratio_[0],label="X")
    plt.ylim([0,1])
    plt.title("%s,sum=%.2g" % ('PCA explained_variance_ratio_',explained_variance_ratio_.sum()))
    ax.xaxis.set_major_formatter(NullFormatter())
    print("PCAsilhouette,%g" % ( score), file=restxt)
    print("PCAcalinski,%g" % (score2), file=restxt)
    print("PCAvarratiocum,%g" % (explained_variance_ratio_.sum()), file=restxt)

    # Perform Isomap Manifold learning.
    t0 = time()
    iso = manifold.Isomap(n_neighbors, n_components=2)
    trans_data=    iso.fit_transform(X).T
    reconstruction_error=iso.reconstruction_error();
    t1 = time()
    # print("%s: %.2g sec" % ('ISO', t1 - t0))
    
    ax = fig.add_subplot(245)
    plt.scatter(trans_data[0, np.logical_not(outter)], trans_data[1, np.logical_not(outter)],
                c=colors[np.logical_not(outter)], cmap=plt.cm.rainbow, alpha=0.7)
    plt.scatter(trans_data[0, outter], trans_data[1, outter], marker='^', c=colors[outter], cmap=plt.cm.rainbow,
                alpha=0.7)

    score = metrics.silhouette_score(trans_data.T, initial_data_y, metric='euclidean')
    score2 = metrics.calinski_harabaz_score(trans_data.T, initial_data_y)
    plt.title("%s, score= %.2g/%.2g\nreconstruction_error=%.2g" % ('Isomap',score,score2,reconstruction_error))
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    plt.axis('tight')
    print("ISOsilhouette,%g" % (score), file=restxt)
    print("ISOcalinski,%g" % (score2), file=restxt)
    print("ISOreconstruction_error,%g" % (reconstruction_error), file=restxt)
    
    # Perform Multi-dimensional scaling.
    t0 = time()
    mds = manifold.MDS(2, max_iter=100, n_init=1)
    try:
        trans_data = mds.fit_transform(X).T
    
        t1 = time()
        # print("MDS: %.2g sec" % (t1 - t0))
        
        ax = fig.add_subplot(246)
        plt.scatter(trans_data[0, np.logical_not(outter)], trans_data[1, np.logical_not(outter)],
                    c=colors[np.logical_not(outter)], cmap=plt.cm.rainbow, alpha=0.7)
        plt.scatter(trans_data[0, outter], trans_data[1, outter], marker='^', c=colors[outter], cmap=plt.cm.rainbow,
                    alpha=0.7)
    
        score = metrics.silhouette_score(trans_data.T, initial_data_y, metric='euclidean')
        score2 = metrics.calinski_harabaz_score(trans_data.T, initial_data_y)
        plt.title("MDS, score= %.2g/%.2g\nStress = %.2g" % (score,score2,mds.stress_))
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        plt.axis('tight')
    except OSError as err:
        print("OS error: {0}".format(err))
    except ValueError:
        print("Could not Perform Multi-dimensional scaling.")
    except:
        print("Unexpected error:", sys.exc_info()[0])
    
    
    # Perform Spectral Embedding.
    try:
        t0 = time()
        se = manifold.SpectralEmbedding(n_components=2,
                                        n_neighbors=n_neighbors)
        trans_data = se.fit_transform(X).T
        t1 = time()
        # print("Spectral Embedding: %.2g sec" % (t1 - t0))
        
        ax = fig.add_subplot(247)
        plt.scatter(trans_data[0, np.logical_not(outter)], trans_data[1, np.logical_not(outter)],
                    c=colors[np.logical_not(outter)], cmap=plt.cm.rainbow, alpha=0.7)
        plt.scatter(trans_data[0, outter], trans_data[1, outter], marker='^', c=colors[outter], cmap=plt.cm.rainbow,
                    alpha=0.7)
    
        score = metrics.silhouette_score(trans_data.T, initial_data_y, metric='euclidean')
        score2 = metrics.calinski_harabaz_score(trans_data.T, initial_data_y)
        plt.title("Spectral Embedding, score =%.2g / %.2g" % (score,score2))
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        plt.axis('tight')
        print("SpectralEmbeddingsilhouette,%g" % (score), file=restxt)
        print("SpectralEmbeddingcalinski,%g" % (score2), file=restxt)
    except OSError as err:
        print("OS error: {0}".format(err))
    except ValueError:
        print("Could not Perform Spectral Embedding.")
    except:
        print("Unexpected error:", sys.exc_info()[0])


    # Perform t-distributed stochastic neighbor embedding.
    t0 = time()
    tsne = manifold.TSNE(n_components=2, init='pca', random_state=0,perplexity=30.0,
                         early_exaggeration=12.0, learning_rate=200.0, n_iter=1000,  n_iter_without_progress=300)
    trans_data = tsne.fit_transform(X).T
    t1 = time()
    # print("t-SNE: %.2g sec\nkl divergence=%.2g" % (t1 - t0,tsne.kl_divergence_))
    
    ax = fig.add_subplot(2, 4, 8)
    plt.scatter(trans_data[0, np.logical_not(outter)], trans_data[1, np.logical_not(outter)],
                c=colors[np.logical_not(outter)], cmap=plt.cm.rainbow, alpha=0.7)
    plt.scatter(trans_data[0, outter], trans_data[1, outter], marker='^', c=colors[outter], cmap=plt.cm.rainbow,
                alpha=0.7)

    score = metrics.silhouette_score(trans_data.T, initial_data_y, metric='euclidean')
    score2 = metrics.calinski_harabaz_score(trans_data.T, initial_data_y)
    plt.title("t-SNE, score= %.2g / %.2g\nkl divergence=%.2g" % (score,score2,tsne.kl_divergence_))
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    plt.axis('tight')

    print("tSNEsilhouette,%g" % (score), file=restxt)
    print("tSNEcalinski,%g" % (score2), file=restxt)

    print("tSNEkl_divergence_,%g" % (tsne.kl_divergence_), file=restxt)


    restxt.close()
    #ax.legend(np.unique(colors),cat)
    #plt.legend(cat)
    #plt.savefig(file_path+ "MainfoldNeig" +str( n_neighbors)+"_km"+str(X.shape[1])+".png", dpi=300)
    plt.savefig(file_path + "MainfoldNeig" + str(n_neighbors) + ".png", dpi=300)



    plt.close()

    '''
    Generate json data for web interface , interpreted by Highcharts-6.0.4
    '''
    yinner = initial_data_y[np.logical_not(outter)]
    tdatainner = trans_data[:,np.logical_not( outter)]
    youtter = initial_data_y[outter]
    tdataoutter = trans_data[:, outter]
    i=0
    icat=1
    jsond=""
    jsond = 'mainfoldorga'+str(organ)+'={'
    if organ>0:
        jsond += 'title:"t-SNE of winner organism <b>rank'+str(organ)+'</b>, '+str(X.shape[1])+'Kmers",'
    else:
        jsond += 'title:"t-SNE ' + str(X.shape[1]) + 'k-mers",'
    jsond += 'file_path:"'+file_path+'",'
    jsond += 'subtitle:"score silhouette=%.2f, score calinski=%.2f / kl divergence=%.2f",'%(score,score2,tsne.kl_divergence_)
    # jsond += 'subtitle:"t-SNE, score: silhouette=' + str(score) + ', calinski=' + str( score2) + ', kl divergence=' + str(tsne.kl_divergence_) + '",'
    jsond +="data:["
    data=[]
    name=[]
    color=[]
    shape=[]
    for c in cat:
        aa = np.array([tdatainner[0, yinner == c], tdatainner[1, yinner == c]]).T
        # data.append( aa.tolist())
        # # data.append( [tdatainner[0, yinner==c].tolist(), tdatainner[1, yinner==c].tolist()])
        # name.append( c)
        # color.append(icat)
        # shape.append("circle")
        if i!=0:
            jsond += ","
        jsond += "{data:" + json.dumps(aa.tolist()) + ",\n"
        jsond += "colorIndex:" + json.dumps(icat) + ",\n"
        jsond += "name:" + json.dumps(c) + ",\n"
        jsond += 'marker: {  symbol:"circle" } }'
        i+=1
        # data.append( [tdataoutter[0, youtter==c], tdataoutter[1, youtter==c]].tolist())
        aa = np.array([tdataoutter[0, youtter == c], tdataoutter[1, youtter == c]]).T
        # data.append( aa.tolist())
        # name.append( c+" outter")
        # color.append( icat)
        # shape.append("triangle")
        if len(aa)>0:
            jsond += ",{data:" + json.dumps(aa.tolist()) + ",\n"
            jsond += "colorIndex:" + json.dumps(icat) + ",\n"
            jsond += "name:" + json.dumps(c+" outter") + ",\n"
            jsond += 'marker: {  symbol:"triangle" } }'
            i+=1
        icat+=1

    # jsond ="datamainfold="+json.dumps(data)+"\n"
    # jsond +="colors="+json.dumps(color)+"\n"
    # jsond += "name=" + json.dumps(name) + "\n"
    jsond += "]}"
    with open(file_path +  ".json", 'w') as f:
        f.write( jsond)






'''
main function for test purpose only
@arg1 : matrix input
@arg2 : output path
'''
if __name__ == "__main__":
    file_path = "0_16/fig/countkmerwin0.8.txt_SampleMat.csv"
    file_path = "BEAUTYtrinegRes/logoutersample/0_17/fig/countkmerwin0.8.txt_SampleMat.csv"
    output_path = 'BEAUTYtrinegRes/logoutersample'
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
    if len(sys.argv) > 2:
        output_path = sys.argv[2]
    n_neighbors = 30
    #
    for i in [10,20,30]:
        mainfoldtest(file_path,i,500,output_path)

    mainfold_organism(file_path, n_neighbors, maxkmer=1000, outerfilepath="", nbkm=10)