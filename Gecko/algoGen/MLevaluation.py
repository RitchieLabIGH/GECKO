from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

#if (sys.version_info > (3, 0)):
import pickle
#else:
#import _pickle as cPickle


import math
import numpy as np
import pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import  AdaBoostClassifier #RandomForestClassifier,
from sklearn.preprocessing import normalize
from sklearn.metrics import accuracy_score
from sklearn.utils import class_weight
from sklearn import  svm, pipeline
from sklearn import linear_model
from sklearn.neighbors import KNeighborsClassifier
from sklearn.kernel_approximation import (RBFSampler,
                                          Nystroem)
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
import itertools

import sys
import os.path
import time,gc
import struct

from sklearn import model_selection
#from sklearn import cross_validation

import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
else:
    try:
        import tkinter
    except ImportError:
        mpl.use('Agg')

import matplotlib.pyplot as plt

import seaborn as sns

# from kerascallback import EarlyStoppingTestData
# from keras.models import Sequential
# from keras.layers import Dense,GaussianNoise, Activation, normalization,Dropout
# #from keras.utils.visualize_util import plot
# from keras.utils.np_utils import to_categorical
# from keras.callbacks import EarlyStopping
# from keras import optimizers
# from keras.callbacks import LearningRateScheduler

#import matplotlib.pyplot as plt

#nbcallin=0
################################################
# MLevaluation
# class to evaluate/analyse classifier from subdataset using crossvalidation able to generate full report
# We recommand to use the static method generateNlzReportClassifierQuality for complete report on GA solutions
# advanced treatement like removing samples, or export/load train model are also available
class MLevaluation:
    samplename =[]
    classname =[]
    kmerName=[]
    countMat =[]
    hideNeurons = 5
    classfierType =["linSVC"]
    kmByIndiv=0
    crossvalNlyzData = pd.DataFrame([])
    def setSavePath(self,path):
        self.savepath = path
        if os.path.isdir(self.savepath) == False:
            os.makedirs(self.savepath)

    ############################################
    # loadData
    # @arg1 : datafile : data input file name
    # @arg2 : hideNeurons : number of hiddenneurons for mpl ( neural network) classifier (average betwen number of input and otput could be a good default value to start with)
    # @arg3 : classfierType : main classifier type ( default : ["linSVC"])
    # @arg4 : kmByIndiv : number of kmer by individuals
    # @arg5 : outerdirpath : folder path to the genetic algorithme root result where to find outer file selection (_outerSelected.csv)
    # @arg6 : extendsavepath : folder or subfolder name to save image and html report of the analysis

    def loadData(self,datafile,hideNeurons=5,classfierType=["linSVC"],kmByIndiv=0,outerdirpath="",extendsavepath=""):
        self.hideNeurons =hideNeurons
        self.classfierType = classfierType
        # import io
        # ...
        # s = open(datafile).read().replace(',', '\t')
        # textinfo = np.loadtxt(io.StringIO(s),  dtype='S')
        textinfo = np.loadtxt(datafile, delimiter='\t', comments='#', dtype='S')
        if len(textinfo.shape)==1 :
            textinfo = np.loadtxt(datafile, delimiter=',', comments='#', dtype='S')

        self.setSavePath( (os.path.splitext(datafile)[0])+extendsavepath)


        self.datafilename=datafile
        self.samplename = textinfo[0,1:]#.astype('U64')
        self.classname = textinfo[1,1:].astype('U64')
        self.kmerName = textinfo[2:,0].T.astype('U64')
        self.countMat = textinfo[2:,1:].astype(np.float)
        self.outerdirpath=outerdirpath

        uniqSamplName,count=np.unique(self.samplename.astype('U64'), return_counts = True)

        for isamplename in np.where(count>1)[0] :
            print("Warning you got "+str(count[isamplename])+" samples with the name : "+uniqSamplName[isamplename])
            formultid=np.where(self.samplename.astype('U64')==uniqSamplName[isamplename])[0]
            for j,idsamp in  enumerate(formultid):
                #self.samplename[idsamp]= self.samplename[idsamp]+"-"+str(j).dtype('U64')
                self.samplename[idsamp] = self.samplename[idsamp].decode("utf-8") + "-{}".format(j)

        self.samplename=self.samplename.astype('U64')

        if kmByIndiv==0:
            self.kmByIndiv = self.kmerName.size
        else:
            self.kmByIndiv =kmByIndiv

        self.infoDatasetLoad()



    ##############################
    # infoDatasetLoad
    # display infos size and classes for the data loaded
    def infoDatasetLoad(self):
        #resume
        nbsamplebyclass=[]
        print("DATA file : ",self.datafilename)
        print("# of samples = ", self.samplename.size)
        print("# of features = ", self.kmerName.size)
        print("# of class = ", np.unique(self.classname).size)
        print("kmByIndiv = ", self.kmByIndiv)
        uniq_classes = np.unique(self.classname)
        for classnamei in uniq_classes:
            ans = np.in1d(self.classname, classnamei).sum()
            print("Number of sample in class :",classnamei ,"\t = ",ans)
            nbsamplebyclass.append(ans)

        return nbsamplebyclass
    ########################################################
    # pairPlotDisplayModel
    # display kmcount pair plot and distribution
    # @arg1 : idindiv : position of the individuals in the input file (base 0 )
    # @arg2 : maxkmfeature : number max of feature to display ( default : 15)
    def pairPlotDisplayModel(self,idindiv=0,maxkmfeature=15):

        data = self.countMat[ idindiv * self.kmByIndiv:(idindiv + 1) * self.kmByIndiv,:].T
        if data.shape[1]>maxkmfeature:
            data=data[:,:maxkmfeature]
        data = normalize(data, axis=1, copy=False)
        dataf=pd.DataFrame( data)

        classname = pd.DataFrame({'classname':self.classname})
        dataf=dataf.join(classname)
        try:
            sns.pairplot(dataf, hue="classname",diag_kind="kde")
        except np.linalg.linalg.LinAlgError as err:
            print("Singular matrix error catch : pair plot af the feature are not completely generate")


        plt.savefig(os.path.join(self.savepath, "featuredistrib.png"))
        plt.close()
        return os.path.join(self.savepath, "featuredistrib.png")
    #######################################
    # sklearn build in cross validation
    # @arg1 : clf : model object sklearn
    # @arg2 : idindiv : position of the individuals in the input file (base 0 )
    # @arg3 : nrot : number of train/test rotation for cross evaluate the model
    def sklearcrossval(self,clf='',idindiv=0,nrot=50):
        if clf == '':
            clf= self.initClassifier(self.classfierType[0])
        data = self.countMat[ idindiv * self.kmByIndiv:(idindiv + 1) * self.kmByIndiv,:].T
        data = normalize(data, axis=1, copy=False)
        data = pd.DataFrame(data)
        res = model_selection.cross_validate(clf, data, y=self.classname, cv=nrot, groups=None, scoring="accuracy") #scoring="hamming"
        print("cross validation sklearn pre-build")
        print("Train score mean = {}".format(np.mean(res["train_score"])))
        print("Test score mean = {}".format(np.mean(res["test_score"])))
        return [np.mean(res["train_score"]), np.std(res["train_score"]), np.mean(res["test_score"]),
                              np.std(res["test_score"])]

        #return pd.DataFrame([[np.mean(res["train_score"]),np.std(res["train_score"]),np.mean(res["test_score"]),np.mean(res["test_score"])]],columns=["train_score","std_train_score","test_score","std_test_score"])
    #########################
    # evaluateModel
    # cross validation analyse of the model with several technics
    # @arg1 : clf : model object sklearn
    # @arg2 : idindiv : position of the individuals in the input file (base 0 )
    # @arg3 : nrot : number of train/test rotation for cross evaluate the model
    # @arg4 : ratio test sample size
    # @arg5 : eliminate : boolean to keep one sample out at everey try in order to mesure the impact of each sample on the model (long computational time, desactivate by default=False)

    def evaluateModel(self,clf='',idindiv=0,nbrot=50 , test_size = 0.33, eliminate = False):
        self.crossvalNlyzData=pd.DataFrame(0.0,columns=np.append(np.unique(self.classname),["scoreInInner","scoreInTest","cptInInner","cptInTest","scoreOut"]),
                                           index=self.samplename.astype('U64'))
        self.crossvalNlyzData.loc[:, "group" ] = pd.Series(self.classname.astype('U64'), index=self.crossvalNlyzData.index)
        self.crossvalNlyzData.loc[:, "falseNeg"] = pd.Series(0.0,
                                                          index=self.crossvalNlyzData.index)
        if clf=='':
            clf= self.initClassifier(self.classfierType[0])
        data = self.countMat[ idindiv * self.kmByIndiv:(idindiv + 1) * self.kmByIndiv,:].T
        data = normalize(data, axis=1, copy=False)
        data = pd.DataFrame(data)
        seriesclassname = pd.Series(self.classname)

        seed = 7
        result=np.array([])
        metaKillLoop = 1
        if eliminate:
            oriData=data
            metaKillLoop = self.samplename.size


        for ikillLoop in range(metaKillLoop) :

            if eliminate :
                print("sample out {}/{}".format(ikillLoop + 1, metaKillLoop))
                seriesclassname = pd.Series(np.delete(self.classname, ikillLoop))
                data=oriData.drop(oriData.index[ikillLoop])
                samplename = np.delete(self.samplename,ikillLoop)
                removedSample=self.samplename[ikillLoop]
            else:
                samplename = self.samplename
            for i in range(nbrot):

                X_train, X_test, Y_train, Y_test = model_selection.train_test_split(data,seriesclassname, test_size=test_size)#,random_state=seed)
                clf.fit(X_train, Y_train)
                scoretest=clf.score(X_test, Y_test)
                result=np.append(result, scoretest)
                print(".", end='',flush=True)
                #print("rotation {}/{} test result = {}".format(i+1,nbrot,scoretest))
                y_pred = clf.predict(X_test)
                #self.generate_confusion_matrix(Y_test, y_pred)
                self.crossval_getnlyzdata( samplename,Y_train, Y_test , y_pred,scoretest)
                if eliminate:
                    self.crossvalNlyzData.loc[removedSample, "scoreOut"] += scoretest
                #print("nb sample in A and B train ",(Y_train=="groupA".encode("utf8")).sum(),(Y_train=="groupB".encode("utf8")).sum())
                #print("nb sample in A and B test ", (Y_test == "groupA".encode("utf8")).sum(), ( Y_test == "groupB".encode("utf8")).sum())
            print("Done")
            if eliminate:
                self.crossvalNlyzData.loc[removedSample,"scoreOut"] /= nbrot
            gc.collect()


        #average the data
        for clasid in np.unique(self.classname).astype('U64'):
            self.crossvalNlyzData.loc[:, clasid]     = self.crossvalNlyzData.loc[:,clasid] / self.crossvalNlyzData.cptInTest
        self.crossvalNlyzData.loc[:, "scoreInTest"]  = self.crossvalNlyzData.loc[:, "scoreInTest"] / self.crossvalNlyzData.cptInTest
        self.crossvalNlyzData.loc[:, "scoreInInner"] = self.crossvalNlyzData.loc[:, "scoreInInner"] / self.crossvalNlyzData.cptInInner
        self.crossvalNlyzData.loc[:, "falseNeg"]     = self.crossvalNlyzData.loc[:, "falseNeg"] / self.crossvalNlyzData.cptInTest
        self.crossvalNlyzData["diffscore"]           = self.crossvalNlyzData.loc[:, "scoreInInner"]-self.crossvalNlyzData.loc[:, "scoreInTest"]
        self.crossvalNlyzData["diffscoreOut"]           = self.crossvalNlyzData.loc[:, "scoreInInner"]-self.crossvalNlyzData.loc[:, "scoreOut"]

        if  os.path.exists(self.outerdirpath) & os.path.isfile(os.path.join( self.outerdirpath ,"_outerSelected.csv")):
            outer = np.genfromtxt(os.path.join( self.outerdirpath ,"_outerSelected.csv"), delimiter=",") * -1
            #outer = outer == 1
            self.crossvalNlyzData["outer"]=outer

        self.crossvalNlyzData                        = self.crossvalNlyzData.sort_values(by=["group","falseNeg","diffscore"], ascending=[False,False,False])
        print("evaluateModel crossvalidation rand  , score :")
        #print(result)
        print("mean = ",np.mean(result),", std = ",np.std(result))
        self.crossvalNlyzDataShow(eliminate)
        self.saveCrossvalMat()
        #permutation
        cv = model_selection.StratifiedKFold(5)

        score, permutation_scores, pvalue = model_selection.permutation_test_score(
            clf, data, seriesclassname, scoring="accuracy", cv=cv, n_permutations=100, n_jobs=1)

        print("Classification score %s (pvalue : %s)" % (score, pvalue))
        # #############################################################################
        ## cross validation sklearn compare
        if eliminate:
            data=  oriData
        res=model_selection.cross_validate(clf, data, y=self.classname, cv=nbrot,groups=None, scoring=None)
        print("cross validation sklearn pre-build")
        print(res)
        res = model_selection.cross_val_score(clf, data, y=self.classname, cv=nbrot, groups=None, scoring=None)
        print("cross validation sklearn pre-build")
        print(res)
        # #############################################################################
        # View histogram of permutation scores
        n_classes = np.unique(self.classname).size

        plt.hist(permutation_scores, 20, label='Permutation scores',
                 edgecolor='black')
        ylim = plt.ylim()
        # BUG: vlines(..., linestyle='--') fails on older versions of matplotlib
        # plt.vlines(score, ylim[0], ylim[1], linestyle='--',
        #          color='g', linewidth=3, label='Classification Score'
        #          ' (pvalue %s)' % pvalue)
        # plt.vlines(1.0 / n_classes, ylim[0], ylim[1], linestyle='--',
        #          color='k', linewidth=3, label='Luck')
        plt.plot(2 * [score], ylim, '--g', linewidth=3,
                 label='Classification Score'
                       ' (pvalue %s)' % pvalue)
        plt.plot(2 * [1. / n_classes], ylim, '--k', linewidth=3, label='Luck')

        plt.ylim(ylim)
        plt.xlim([0,1])
        plt.legend()
        plt.xlabel('Score')

        plt.savefig(os.path.join(self.savepath, "permutation_score.png"))

        plt.close()
        return np.mean(result),np.std(result)
    #############################################
    # generateConfusionMat4evaluateModel
    # generate confusion matrix from from cross validation result
    def generateConfusionMat4evaluateModel(self):
        ugroups=np.unique(self.crossvalNlyzData["group"])
        confmat=pd.DataFrame([],columns=ugroups,dtype="float")
        for grpX in ugroups:
            grpres = self.crossvalNlyzData.loc[self.crossvalNlyzData['group'] == grpX, :]
            for grpY in ugroups:
                confmat.loc[grpX, grpY] = np.mean(grpres.loc[:, grpY]) * 100

        sns.heatmap(confmat, vmin=0, vmax=100, cmap="jet", annot=True, fmt=".1f")  # YlGnBu
        plt.savefig(os.path.join(self.savepath, "evalModelConfusionMat.png"))

        plt.close()

    #######################################################
    # crossvalNlyzDataShow
    # generate and save analyse figures
    # diff_score_pairplot.png,sort_diffout_heatmap.png,full_Heatmap.png,diff_score_plot.png
    def crossvalNlyzDataShow(self,eliminate = False):
        plt.figure()
        if eliminate==False:
            self.crossvalNlyzData=self.crossvalNlyzData.drop(["diffscoreOut"],axis=1)

        try:
            sns.pairplot(self.crossvalNlyzData, vars=["diffscore", "falseNeg"], hue="group",diag_kind="kde", kind="reg")
        #plt.show()
        except np.linalg.linalg.LinAlgError as err:
            print("Singular matrix error catch : pair plot af the feature are not completely generate")

        except ValueError:
            print("Could not convert data to an integer.")
        except:
            print("Unexpected error:", sys.exc_info()[0])

        plt.savefig(os.path.join(self.savepath,"diff_score_pairplot.png"))

        plt.close()

        #sns.pairplot(self.crossvalNlyzData, hue="group")
        if eliminate:
            aaa = self.crossvalNlyzData.sort_values(by=["diffscoreOut", "group", "falseNeg", "diffscore"],
                                                    ascending=[False, False, False, False])
            sns.heatmap(aaa.drop(["group", "cptInTest", "cptInInner"], axis=1),cmap="YlGnBu",vmin=0,vmax=1)
            plt.savefig(os.path.join(self.savepath,"sort_diffout_heatmap.png"))
            plt.close()
        sns.heatmap(self.crossvalNlyzData.drop(["group", "cptInTest", "cptInInner"], axis=1),cmap="YlGnBu",center=0.5)
        plt.savefig(os.path.join(self.savepath,"full_Heatmap.png"))

        plt.close()
        # sns.heatmap(self.crossvalNlyzData.loc[:,"diffscore"], axis=1)
        plt.plot(range(self.classname.size), self.crossvalNlyzData.loc[:, "diffscore"])
        if eliminate:
            plt.plot(range(self.classname.size), self.crossvalNlyzData.loc[:, "diffscoreOut"], color="green")
        plt.savefig(os.path.join(self.savepath,"diff_score_plot.png"))

        plt.close()
    ########################
    #saveCrossvalMat
    # save analyse count and scores values in crossvalidationresult.csv
    def saveCrossvalMat(self):
        self.crossvalNlyzData.to_csv(os.path.join(self.savepath,"crossvalidationresult.csv"), sep="\t")

    ########################
    # crossval_getnlyzdata
    # store analyse count and scores values in crossvalNlyzData datafram
    def crossval_getnlyzdata(self,samplename, Y_train, Y_test , y_pred,scoretest):
        for itrain in Y_train.index:
            #print(itrain)
            self.crossvalNlyzData.loc[samplename[itrain] ,"cptInInner"] += 1
            self.crossvalNlyzData.loc[samplename[itrain] ,"scoreInInner"] += scoretest
        for i,itest in enumerate(Y_test.index):
            self.crossvalNlyzData.loc[samplename[itest] ,"cptInTest"] += 1
            self.crossvalNlyzData.loc[samplename[itest] ,"scoreInTest"] += scoretest
            self.crossvalNlyzData.loc[samplename[itest] ,y_pred[i]] += 1
            if y_pred[i]!=Y_test.iloc[i] :
                self.crossvalNlyzData.loc[samplename[itest] ,"falseNeg"] += 1


    ##################################################################################"
    #oneVsAll
    # Analysyse type one vs all
    # arg1 : clf : object model (sklearn)
    # arg2 : idindiv : position of the individuals in the input file (base 0 )
    # arg3 : number of train test rotations
    # arg4: ratio test sample size
    def oneVsAll(self, clf, idindiv=0, nbrot=5, test_size=0.33):
        data = self.countMat[idindiv * self.kmByIndiv:(idindiv + 1) * self.kmByIndiv, :].T
        data = normalize(data, axis=1, copy=False)

        from sklearn.preprocessing import label_binarize
        Y = label_binarize(self.classname, classes=np.unique(self.classname))
        uniqClasname=np.unique(self.classname)
        n_classes = Y.shape[1]

        result = np.array([])
        for i in range(nbrot):
            X_train, X_test, Y_train, Y_test = model_selection.train_test_split(data, Y,
                                                                                test_size=test_size)  # ,random_state=seed)
            classifier = OneVsRestClassifier(clf)
            classifier.fit(X_train, Y_train)
            y_score = classifier.decision_function(X_test)

            # For each class
            precision = dict()
            recall = dict()
            average_precision = dict()

            for i in range(n_classes):
                if n_classes > 1:
                    tmpscore=y_score[:, i]
                else:
                    tmpscore=y_score
                precision[i], recall[i], _ = precision_recall_curve(Y_test[:, i],tmpscore)
                average_precision[i] = average_precision_score(Y_test[:, i], tmpscore)


            # A "micro-average": quantifying score on all classes jointly
            precision["micro"], recall["micro"], _ = precision_recall_curve(Y_test.ravel(),
                                                                            y_score.ravel())
            average_precision["micro"] = average_precision_score(Y_test, y_score,
                                                                 average="micro")
            print('Average precision score, micro-averaged over all classes: {0:0.2f}'
                  .format(average_precision["micro"]))
            result = np.append(result,average_precision["micro"])

            plt.figure()
            plt.step(recall['micro'], precision['micro'], color='b', alpha=0.2,
                     where='post')
            #plt.fill_between(recall["micro"], precision["micro"], alpha=0.2, color='b', **step_kwargs)

            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.ylim([0.0, 1.05])
            plt.xlim([0.0, 1.0])
            plt.title(
                'Average precision score, micro-averaged over all classes: AP={0:0.2f}'
                    .format(average_precision["micro"]))
            plt.savefig(os.path.join(self.savepath, "average_precision_score.png"))
            plt.close()

            #Plot            Precision - Recall            curve            for each class and iso-f1 curves¶
            from itertools import cycle
            # setup plot details
            colors = cycle(['navy', 'turquoise', 'darkorange', 'cornflowerblue', 'teal'])

            plt.figure(figsize=(7, 8))
            f_scores = np.linspace(0.2, 0.8, num=4)
            lines = []
            labels = []
            for f_score in f_scores:
                x = np.linspace(0.01, 1)
                y = f_score * x / (2 * x - f_score)
                l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
                plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))

            lines.append(l)
            labels.append('iso-f1 curves')
            l, = plt.plot(recall["micro"], precision["micro"], color='gold', lw=2)
            lines.append(l)
            labels.append('micro-average Precision-recall (area = {0:0.2f})'
                          ''.format(average_precision["micro"]))

            for i, color in zip(range(n_classes), colors):
                l, = plt.plot(recall[i], precision[i], color=color, lw=2)
                lines.append(l)
                labels.append('Precision-recall for class {0} (area = {1:0.2f})'
                              ''.format(uniqClasname[i] , average_precision[i]))

            fig = plt.gcf()
            fig.subplots_adjust(bottom=0.25)
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.title('Extension of Precision-Recall curve to multi-class')
            plt.legend(lines, labels, loc=(0, -.38), prop=dict(size=14))
            plt.savefig(os.path.join(self.savepath, "Precision-Recall_curve.png"))
            plt.close()


        print("oneVsAll crossvalidation Average precision score, micro-averaged over all classes:")
        # print(result)
        print("mean = ", np.mean(result), ", std = ", np.std(result))


    ###########################################################################################
    # loadEvaluateModel( idindiv=0,filename="modeldump.sklsav")
    #  train on all sample and save model already train from file
    #  @arg1 : idindiv : position of the individuals in the data file (base 0 )
    #  @arg2 : filename : filename to the save the model
    #  return : the model object sklearn
    def saveFullyTrainModel(self, idindiv=0,filename="modeldump.sklsav", classfierType="",noisefactor=0):
        if classfierType == "":
            classfierType = self.classfierType[0]
        # copy data for the selected individual and translate it for ML processing so each row becom a sample and columns are the "kmer"
        data = self.countMat[ idindiv * self.kmByIndiv:(idindiv + 1) * self.kmByIndiv,:].T
        data = normalize(data, axis=1, copy=False)
        if noisefactor!= 0:
            data= self.addnoise(data, noisefactor)

        clf = self.initClassifier(classfierType)
        clf.fit(data, self.classname)


        # save the model to disk

        pickle.dump(clf, open(filename, 'wb'))

        # some time later...

        # load the model from disk
        loaded_model = pickle.load(open(filename, 'rb'))
        result = loaded_model.score(data, self.classname)
        print("saveFullyTrainModel model save as "+filename+" , score on training data =")
        print(result)
        return clf
        1
    ###########################################################################################
    # loadEvaluateModel( idindiv=0,filename="modeldump.sklsav")
    #   load and evalute model already train from file
    #  @arg1 : idindiv : position of the individuals in the data file (base 0 )
    #  @arg2 : filename : filename to the saved model

    def loadEvaluateModel(self, idindiv=0,filename="modeldump.sklsav"):

        # copy data for the selected individual and translate it for ML processing so each row becom a sample and columns are the "kmer"
        data = self.countMat[ idindiv * self.kmByIndiv:(idindiv + 1) * self.kmByIndiv,:].T
        data = normalize(data, axis=1, copy=False)



        # load the model from disk
        loaded_model = pickle.load(open(filename, 'rb'))
        result = loaded_model.score(data, self.classname)

        print("loadEvaluateModel model  "+filename+" loaded , score =")
        print(result)
        y_pred=loaded_model.predict(data)
        self.generate_confusion_matrix( self.classname, y_pred)



        return loaded_model

    def addnoise(self,table, noisefactor):

        return table + (np.random.normal(0, 1, table.shape) * np.std(table, axis=0) * noisefactor)
    ################################################################################################
    # initialyse classifier
    # @arg1 : classiType possible value : ['MLP','ADA','SDG','SVC','linSVC', 'KNEIG','SVMfourier','SVMnystro']
    # @arg2 : hideNeurons : number of hiddenneurons for mpl ( neural network) classifier (average betwen number of input and otput could be a good default value to start with)
    # return : model object sklearn
    def initClassifier(self,classiType,hideNeurons=10):
        if classiType == 'MLP':
            clf = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=hideNeurons, random_state=1,
                                learning_rate='adaptive', activation='logistic')
        else:
            if classiType == 'ADA':
                clf = AdaBoostClassifier()
            if classiType == 'SDG':
                clf = linear_model.SGDClassifier(class_weight='balanced')
            if classiType == 'SVC':
                clf = svm.SVC(gamma=.4, class_weight='balanced')
            if classiType == 'linSVC':
                clf = svm.LinearSVC(class_weight='balanced')
            if classiType == 'KNEIG':
                clf = KNeighborsClassifier(n_neighbors=5, algorithm='auto')

            # create pipeline from kernel approximation
            # and linear svm
            if classiType == 'SVMfourier':
                feature_map_fourier = RBFSampler(gamma=.2, n_components=hideNeurons)  # , random_state=1

                clf = pipeline.Pipeline([("feature_map", feature_map_fourier),
                                         ("svm", svm.LinearSVC(class_weight='balanced'))])
            if classiType == 'SVMnystro':
                feature_map_nystroem = Nystroem(gamma=.2, n_components=hideNeurons)  # , random_state=1
                clf = pipeline.Pipeline([("feature_map", feature_map_nystroem),
                                         ("svm", svm.LinearSVC(class_weight='balanced'))])
        return clf
    def generate_confusion_matrix(self,y_test, y_pred):
        from sklearn.metrics import confusion_matrix
        class_names=np.unique(self.classname)
        # Compute confusion matrix
        cnf_matrix = confusion_matrix(y_test, y_pred,labels=class_names)
        np.set_printoptions(precision=2)

        # Plot non-normalized confusion matrix
        plt.figure()
        self.plot_confusion_matrix(cnf_matrix, classes=class_names,
                              title='Confusion matrix')
        plt.savefig(os.path.join(self.savepath, "Confusion_matrix.png"))
        plt.close()




        plt.show()
    def geterrorbysample(self,y_test, y_pred):
        from sklearn.metrics import confusion_matrix
        class_names=np.unique(self.classname)
        # Compute confusion matrix
        cnf_matrix = confusion_matrix(y_test, y_pred,labels=class_names)
        np.set_printoptions(precision=2)

        # Plot non-normalized confusion matrix
        plt.figure()
        self.plot_confusion_matrix(cnf_matrix, classes=class_names,
                              title='Confusion matrix')



        plt.show()
    def plot_confusion_matrix(self, cm, classes,
                              normalize=False,
                              title='Confusion matrix',
                              cmap=plt.cm.Blues):
        classes = [x.decode('UTF8') for x in classes]
        """
        This function prints and plots the confusion matrix.
        Normalization can be applied by setting `normalize=True`.
        """
        if normalize:
            cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
            print("Normalized confusion matrix")
        else:
            print('Confusion matrix')

        print(cm)

        plt.imshow(cm, interpolation='nearest', cmap=cmap)
        plt.title(title)
        plt.colorbar()
        tick_marks = np.arange(len(classes))
        plt.xticks(tick_marks, classes, rotation=45)
        plt.yticks(tick_marks, classes)

        fmt = '.2f' if normalize else 'd'
        thresh = cm.max() / 2.
        for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            plt.text(j, i, format(cm[i, j], fmt),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")

        plt.ylabel('True label')
        plt.xlabel('Predicted label')
        plt.tight_layout()
    #####################################################################################################################
    # applyDataFilterBySamplesName(self,filelist,filecorrespondinglist="")
    # remove sample from sample list
    # @arg1 : filelist : file with sample list
    # @arg2 : filecorrespondinglist : file to convert name if needed (optional)
    #     tabulated text file 1st column : sample data name; 2nd column : filelist sample name

    def applyDataFilterBySamplesName(self,filelist,filecorrespondinglist=""):
        nametorevmove = np.loadtxt(filelist, delimiter='\t', comments='#', dtype='S')

        if filecorrespondinglist!="":
            textinfo = np.loadtxt(filecorrespondinglist, delimiter='\t', comments='#', dtype='S')
            convertingname = textinfo[:, 0:2]
            duplicatelistorign=nametorevmove
            for i,nameconv in enumerate(duplicatelistorign):
                arr_index = np.where(convertingname[:,1] == nameconv)
                if arr_index[0].size !=0:
                    nametorevmove[i]=convertingname[arr_index[0],0][0]
                else:
                    print("not found ",nameconv)
        masktoremove = np.ones(self.samplename.shape, dtype=bool)
        nametorevmove=nametorevmove.astype('U64')
        print("{} samples name to be remove".format(len(nametorevmove) ))
        for nameconv in nametorevmove:
            masktoremove[self.samplename==nameconv]=False
        self.countMat =self.countMat[:,masktoremove]
        self.samplename = self.samplename[masktoremove]
        self.classname = self.classname[masktoremove]
        #print("lumA cat : ",self.samplename[self.classname =="LumA".encode('UTF8')])
        print("{} samples removed".format(len(masktoremove)-masktoremove.sum()))
        print("Removed sample from : ",filelist)
        self.infoDatasetLoad()
    #####################################################################################################################
    # Static method
    #generateNlzReportClassifierQuality (fi, extendsavepath="", nrot=50, hideNeurons=5, kmByIndiv=10, idindiv=0, eliminate=False,classfierType="linSVC")
    # Generate full html report in new folder name like you input file (+ extendsavepath)
    # exemple : fi="GA/0_1/fig/BestIndiv10.csvforextractkm.count_SampleMat.csv"
    #           extendsavepath="Nlyz/indiv1/"
    #           generated report will be save in "GA/0_1/fig/BestIndiv10.csvforextractkm.count_SampleMatNlyz/indiv1/"
    # input:
    # @arg1 : fi : input file
    # @arg2 : extendsavepath : folder or subfolder name to save image and html report of the analysis
    # @arg3 : nrot : number of train/test rotation too evaluate de the model
    # @arg4 : hideNeurons : number of hiddenneurons for mpl ( neural network) classifier (average betwen number of input and otput could be a good default value to start with)
    # @arg5 : kmByIndiv : number of kmer by individuals
    # @arg6 : idindiv : position of the individuals in the input file (base 0 )
    # @arg7 : eliminate : boolean to keep one sample out at everey try in order to mesure the impact of each sample on the model (long computational time, desactivate by default=False)
    # @arg8 : classfierType : main classifier type ( default : "linSVC")
    #
    @staticmethod
    def generateNlzReportClassifierQuality(fi, extendsavepath="", nrot=0, hideNeurons=5, kmByIndiv=10, idindiv=0,
                                          eliminate=False,classfierType="linSVC"):
        slcrossres = pd.DataFrame([], columns=["train_score", "std_train_score", "test_score", "std_test_score"])


        MLevaluator = MLevaluation()
        htmlreport = '<html>\n\t<head>\n\t\t<script src="sorttable.js"></script>\n\t\t<style>table {    border-collapse: collapse;   }\n'
        htmlreport += "th, td {    text-align: left;    padding: 3px;}\n"
        htmlreport += "tr:nth-child(even) {background-color: #f2f2f2;}\n.halfwidth {width:50%; display: inline-block;}"
        htmlreport += ".title{font-size:150%; font-style:bold} \n .infos{font-style:italic; font-size:80%}\nimg{width:100%}\n .tableres * {font-size:90%}</style>\n</head>\n\t<body>\n"
        htmlreport += "\t\t\t<div class='title' style='text-align: center;'>Classifier analysis number {}, {} inputs kmers</div>\n".format(idindiv+1,kmByIndiv)
        htmlreport += "\t\t<div class='sectionres'>\n\t\t\t<div class=datafilename>" + fi + "\n\t\t\t</div>\n"
        htmlreport += "\t\t\t<div class='title'>Compare machine learning algorithms</div>\n<div class='infos' >(All table are sortable by click on the columns titles)"
        print("## classificator power ##" + fi)

        outerdirpath=os.path.join(os.path.split(fi)[0],"../../")

        MLevaluator.loadData(fi, hideNeurons, classfierType, kmByIndiv=kmByIndiv,outerdirpath=outerdirpath,
                             extendsavepath=extendsavepath)
        # MLevaluator.applyDataFilterBySamplesName("/home/sylvain/genologin/TCGA_GA/list1kmertoremoveAubin.csv")
        if nrot==0:
            nbsamplebyclass=MLevaluator.infoDatasetLoad()
            nrot=int(2*np.min(nbsamplebyclass)/3)
            print("adaptative nrot = "+str(nrot))

        commandclient = "cp sorttable.js " + os.path.join( MLevaluator.savepath,"sorttable.js")
        os.system(commandclient)
        print(" Compare classifier")
        for clstype in ["linSVC", "ADA", "SDG", "SVC", "KNEIG", "SVMfourier", "SVMnystro", "MLP"]:
            print("classifier type : " + clstype)
            model = MLevaluator.initClassifier(clstype, hideNeurons=hideNeurons)

            # slcrossres[clstype]=MLevaluator.sklearcrossval(clf=model, idindiv=0)
            slcrossres.loc[clstype, :] = MLevaluator.sklearcrossval(clf=model, idindiv=idindiv,nrot=nrot)

        # print(" Compare classifier")
        # for clstype in ["linSVC", "ADA", "SDG", "SVC", "KNEIG", "SVMfourier", "SVMnystro", "MLP"]:
        #     print("classifier type : " + clstype)
        #     slcrossres.append(MLevaluator.sklearcrossval(clf=model, idindiv=0))
        from IPython.display import display
        slcrossres.to_csv(os.path.join( MLevaluator.savepath,"multiModelCrossValResult.csv"),sep="\t")
        display(slcrossres)


        htmlreport += "\t\t<div class='tableres' id='MLevalres'>" + slcrossres.to_html(classes="sortable") + "</div>\n\t\t\t</div>"

        model = MLevaluator.initClassifier(classfierType)
        meanres, stdres = MLevaluator.evaluateModel(model, nbrot=nrot,idindiv=idindiv, eliminate=eliminate)
        htmlreport += "\t\t<div class='res' id='rexcrossrand'>Random crossvalidation (rotation {}) : mean test ={}, std test={}".format(
            50, meanres, stdres) + "</div>"

        display(MLevaluator.crossvalNlyzData)
        # htmlreport +="toto\n\n"
        if MLevaluator.crossvalNlyzData["scoreOut"].sum() == 0:
            col2drop = ["cptInInner", "cptInTest", "scoreOut"]
        else:
            col2drop = ["cptInInner", "cptInTest"]
        tableressamples = "\t\t<div class='title'>Cross validation score by samples ("+classfierType+")</div> <div class='tableres' id='evaluateModelres'>" + MLevaluator.crossvalNlyzData.drop(col2drop,
                                                                                                             axis=1).to_html(classes="sortable") + "</div>\n"
        if 'outer' in MLevaluator.crossvalNlyzData.columns:
            outerfalse = MLevaluator.crossvalNlyzData[MLevaluator.crossvalNlyzData["outer"] == 1]["falseNeg"]
            innerfalse = MLevaluator.crossvalNlyzData[MLevaluator.crossvalNlyzData["outer"] == 0]["falseNeg"]
            htmlreport += "\t\t<div class='' id='resouter'> False negative average for GA outer samples : {:.3}</br>False negative average for  non outer samples : {:.3}".format(outerfalse.mean(),innerfalse.mean()) +  "</div>\n"



        MLevaluator.generateConfusionMat4evaluateModel()
        htmlreport += "\t\t<div class='confmat halfwidth'><div style='text-align:center'> Confusion matix</div> <img src='evalModelConfusionMat.png' /></div>"



        MLevaluator.pairPlotDisplayModel(idindiv)
        htmlreport += "<div class='paiplotfeature  halfwidth'> <div style='text-align:center'> Pair plot false negative / diff score ( scoreInInner-scoreInTest)</div><img src='diff_score_pairplot.png' />\n</div>\n"
        htmlreport += "<div><div class='title'>Kmers distributions</div> <img src='featuredistrib.png'   style='width:100%'/></div>"
        htmlreport += tableressamples

        MLevaluator.oneVsAll(model, idindiv=idindiv, nbrot=1, test_size=0.20)
        htmlreport += "<div class='containeroneVsAll'><div class='title'> One vs all "+classfierType+" </div>\n"
        #<div class='oneVsAll'> <img src='" + os.path.join(
        #    MLevaluator.savepath, "Confusion_matrix.png") + "' /></div>"
        htmlreport += "<div class='oneVsAll halfwidth'> <img src='average_precision_score.png' /></div>"
        htmlreport += "<div class='oneVsAll halfwidth'> <img src='Precision-Recall_curve.png' /><div class='oneVsAll'>"
        htmlreport += "\t</body>\n</html>"
        htmlreportfile = open(os.path.join(MLevaluator.savepath, "testfile.html"), "w")
        htmlreportfile.write(htmlreport)
        htmlreportfile.close()

        return slcrossres


'''
hideNeurons = 0

classfierType = "linSVC"
detailResFile =""
conffile ="BEAUTYtrinegRes/configbase.conf"
datafile = "/home/sylvain/Documents/cleangecko/GECKO_source/algoGen/PAM50/Beauty/concatfile_emptysampleremoved_noLUMA.csv"
datafile = "/home/sylvain/Documents/cleangecko/GECKO_source/algoGen/PAM50/TCGA/concatfile.csv"
#datafile = "/home/sylvain/shenron_zfs/home/sylvain/TCGA2Beauty/extractbestindiv.matrix"
datafile = "/home/sylvain/shenron_zfs/home/sylvain/TCGA2Beauty/FromTCGAtest1_AlIAGADir_0_5_BestIndiv10.csvforextractkm.count_SampleMat.csv"
datafile = "/home/sylvain/genologin/CLL_mistery/testlength/test1_Re2Km5no2Dir/0_4/fig/BestIndiv10.csvforextractkm.count_SampleMat.csv"
datafile = "/home/sylvain/genologin/CLL_mistery/testlength/test1_Re2Km5no0.5Dir/0_5/fig/BestIndiv10.csvforextractkm.count_SampleMat.csv"

datafile = "/home/sylvain/genologin/CLL_mistery/testlength/test1_Re1Km5no0.5Dir/0_4/fig/BestIndiv10.csvforextractkm.count_SampleMat.csv"
datafile = "/home/sylvain/genologin/CLL_mistery/testlength/test1_Re2Km5no2Dir/0_4/fig/BestIndiv10.csvforextractkm.count_SampleMat.csv"
#datafile = "/home/sylvain/genologin/CLL_mistery/testlength/test1_Re2Km5no0.5Dir/0_5/fig/BestIndiv10.csvforextractkm.count_SampleMat.csv"
datafile = "/home/sylvain/genologin/CLL_mistery/testlength/test1_Re/mergedWinnerkmlis0_nbkm5000.fastq"
# parameters : idapp  hideNeurons  shufTrainTestDataType  classfierType  noisefactor  detailResFile
if len(sys.argv) > 1:
    datafile = sys.argv[1] # id application for communication
if len(sys.argv) > 2:
    conffile= sys.argv[2]
MLevaluator=MLevaluation()


MLevaluator.loadData(datafile,hideNeurons,classfierType,detailResFile,kmByIndiv=5)
MLevaluator.pairPlotDisplayModel()
'''
'''
model= MLevaluator.initClassifier(classfierType)
MLevaluator.evaluateModel(model)
model = MLevaluator.saveFullyTrainModel(noisefactor=0)

#datafile = "/home/sylvain/Documents/cleangecko/GECKO_source/algoGen/PAM50/Beauty/concatfile_emptysampleremoved_noLUMA.csv"
datafile = "/home/sylvain/shenron_zfs/home/sylvain/TCGA2Beauty/extractbestindiv.TCGAcatCompatible_cleaned.matrix"
MLevaluator.loadData(datafile,hideNeurons,classfierType,detailResFile,kmByIndiv=10)
MLevaluator.applyDataFilterBySamplesName("/home/sylvain/shenron_zfs/home/sylvain/TCGA2Beauty/_sampletoremovefromBeauty.txt","/home/sylvain/shenron_zfs/home/sylvain/TCGA2Beauty/_Beautynamecorrespondance.csv")

model = MLevaluator.loadEvaluateModel()
#MLevaluator.evaluateModel(model)
#MLevaluator.evaluateModel(model,idindiv=1)
#MLevaluator.oneVsAll(model)



kmByIndiv=0
datafile = "/home/sylvain/Documents/cleangecko/GECKO_source/algoGen/PAM50/Beauty/concatfile_emptysampleremoved_noLUMA.csv"
datafile = "/home/sylvain/Documents/cleangecko/GECKO_source/algoGen/PAM50/TCGA/concatfile.csv"




MLevaluator.loadData(datafile,hideNeurons,classfierType,detailResFile,kmByIndiv=kmByIndiv)
model= MLevaluator.initClassifier(classfierType)
MLevaluator.evaluateModel(model)
model = MLevaluator.saveFullyTrainModel(noisefactor=0)

datafile = "/home/sylvain/Documents/cleangecko/GECKO_source/algoGen/PAM50/Beauty/concatfile_emptysampleremoved_noLUMA.csv"
datafile = "/home/sylvain/Documents/cleangecko/GECKO_source/algoGen/PAM50/Beauty/concatfile_emptysampleremoved.csv"

MLevaluator.loadData(datafile,hideNeurons,classfierType,detailResFile,kmByIndiv=kmByIndiv)
model = MLevaluator.loadEvaluateModel()

'''