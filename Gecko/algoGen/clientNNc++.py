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

import sys
import os.path
import time,gc
import struct
# from kerascallback import EarlyStoppingTestData
# from keras.models import Sequential
# from keras.layers import Dense,GaussianNoise, Activation, normalization,Dropout
# #from keras.utils.visualize_util import plot
# from keras.utils.np_utils import to_categorical
# from keras.callbacks import EarlyStopping
# from keras import optimizers
# from keras.callbacks import LearningRateScheduler

#import matplotlib.pyplot as plt
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'
os.environ['TMPDIR'] = "/tmp"
from mpi4py import MPI
comm = MPI.COMM_WORLD
#nbcallin=0
def main(idapp,hideNeurons,shufTrainTestDataType,classfierType,noisefactor,detailResFile):
    #global nbcallin
    if detailResFile!="":
        detailResFile=detailResFile+str(idapp)
    countit=0
    comm_size = comm.Get_size()
    comm_rank = comm.Get_rank()
    print(comm_rank, "/", comm_size)
    unique_indiv=[]

    commandkey = ""

    file_output = str(idapp)+"_response"
    file_input = str(idapp)+"_request"
    file_stop = str(idapp) + "_stop"
    if comm_rank == 0:
        # remove file log score
        if os.path.isfile(pathLog + "AllScoreByGeneration.csv"):
            os.remove(pathLog + "AllScoreByGeneration.csv")
        if os.path.isfile(pathLog + "AllOUTTERScoreByGeneration.csv"):
            os.remove(pathLog + "AllOUTTERScoreByGeneration.csv")

    # #   if os.path.isfile(file_input):
    #  #      os.remove(file_input)

        if os.path.isfile(file_output):
            os.remove(file_output)

    while True:

        if comm_rank==0:
            countit=countit+1

            #print("client : waiting for file"+str(countit))
            while (not os.path.exists(file_input)) & (not os.path.exists(file_stop)):
                time.sleep(0.01)

            if os.path.isfile(file_stop):

                commandkey="stop"
                print("End of client : " + str(idapp))
                os.remove(file_stop)
                #sys.exit()

            if os.path.isfile(file_input):
                #print("read file client<")
                #initial_data = pd.read_csv(file_input, sep=",", header=None)
                initial_data =readBinRequest(file_input)

                #os.rename(file_input, file_input + str(countit) + ".csv")
                i_group = initial_data.shape[1] - 3
                i_traintest = initial_data.shape[1] - 2
                i_indiv = initial_data.shape[1] - 1
                init_group=initial_data.loc[:,i_group].astype(int)

                init_traintest = initial_data.loc[:, i_traintest].astype(int)
                init_indiv = initial_data.loc[:, i_indiv].astype(int)
                init_indiv = init_indiv - min(init_indiv)

                unique_indiv = np.unique(init_indiv)
                if min(init_traintest) < 0:
                    #print("hdin received",nbcallin)
                    commandkey = "hdin"
                    if min(init_traintest) == -2:
                        commandkey = "hdinlogall"
                        init_traintest=init_traintest//2

                    '''if nbcallin==1:#comment for prod
                        print("debug exit......")
                        commandkey = "stop"
                    '''
                    init_traintest=np.abs(init_traintest)
                else:
                    commandkey = "hdout"

                        #os.remove(file_input)#comment for prod

                if (idapp != "UID"): #&(commandkey != "stop")) : #& (commandkey != "hdin"):
                    os.remove(file_input)  # ucomment for prod
                    #commandkey = "hdout"

                del initial_data[i_group]
                del initial_data[i_traintest]
                del initial_data[i_indiv]

             # commandkey=comm.bcast(commandkey, root=0)
             # print(" 1 fin read file client" + commandkey)
            nb_idiv_cpu=int(np.ceil(len(unique_indiv)/comm_size))
            #print("nb_idiv_ by cpu "+str(nb_idiv_cpu))

            trainsgroups=bin(max(init_traintest))[2:]
            for i in range(0, comm_size):
                req = comm.send(trainsgroups, dest=i, tag=9)


        #trainsgroups = comm.bcast(trainsgroups, root=0)
        #if comm_rank != 0:
        while not comm.Iprobe(0, 9):
              time.sleep(0.001)
        trainsgroups=comm.recv(source=0, tag=9)
        commandkey = comm.bcast(commandkey, root=0)
        if commandkey=="stop":
            sys.exit()
        #read if multiple traintest group coded in base 2
        nbInstances = len(unique_indiv)
        resRandTree = np.zeros(nbInstances)
        for ntrainpow in range(len(trainsgroups)): #max(init_traintest)>
            if comm_rank == 0 :
                idtrainrot=np.power(2,ntrainpow)
                for i in range(1,comm_size):
                    selindiv=range(i*nb_idiv_cpu,i*nb_idiv_cpu+nb_idiv_cpu)
                    seltrain=(init_traintest & idtrainrot != idtrainrot) & (init_indiv.isin(selindiv))
                    seltest=(init_traintest & idtrainrot == idtrainrot) & (init_indiv.isin(selindiv))

                    req=comm.send(initial_data.loc[seltrain], dest=i, tag=1) #idx
                    req =comm.send(initial_data.loc[seltest], dest=i, tag=2)  #tidx
                    req =comm.send(init_group.loc[seltrain], dest=i, tag=3)  # idy
                    req =comm.send(init_group.loc[seltest], dest=i, tag=4)  # tidy
                    req =comm.send(init_indiv.loc[seltrain], dest=i, tag=5) #individu train
                    req =comm.send(init_indiv.loc[seltest], dest=i, tag=6)  # individu test

                    ## for rank0

                selindiv = range(0,  nb_idiv_cpu)
                seltrain = (init_traintest & idtrainrot != idtrainrot) & (init_indiv.isin(selindiv))
                seltest = (init_traintest & idtrainrot == idtrainrot) & (init_indiv.isin(selindiv))
                idx = initial_data.loc[seltrain]
                tidx = initial_data.loc[seltest]
                idy = init_group.loc[seltrain]
                tidy = init_group.loc[seltest]
                indiv_train=init_indiv.loc[seltrain]
                indiv_test=init_indiv.loc[seltest]


            else :
                idx = comm.recv(source=0, tag=1)
                tidx = comm.recv(source=0, tag=2)
                idy = comm.recv(source=0, tag=3)
                tidy = comm.recv(source=0, tag=4)
                indiv_train = comm.recv(source=0, tag=5)
                indiv_test = comm.recv(source=0, tag=6)


            idy=np.array(idy)
            tidy = np.array(tidy)
            indiv_train = np.array(indiv_train)
            indiv_test = np.array(indiv_test)

            unique_indiv= comm.bcast(unique_indiv,root=0)
            #commandkey = comm.bcast(commandkey, root=0)
#            if commandkey=="stop":
#                sys.exit()
            nbInstances=len(unique_indiv)
           # chunksize = int(np.floor(nbInstances / comm_size))

            resRandTreesingle = np.zeros(nbInstances)
            resRandTreeTmp=np.array([])

            #print(" fin read file client"+commandkey)


            if len(idy)!=0 :
                if "hdin"in commandkey :
                    #nbcallin+=1
                    if "KERAS" in classfierType:
                        resRandTreeTmp = KERAS_computeAccuracyScore_HDIn(idx, idy,indiv_train, tidx, tidy,indiv_test, hideNeurons, classfierType)
                    else :
                        resRandTreeTmp = computeAccuracyScore_HDIn(idx, idy,indiv_train, tidx, tidy,indiv_test, hideNeurons, classfierType)
                else:
                    if "KERAS" in classfierType:
                        resRandTreeTmp = KERAS_computeAccuracyScore_HDOut(idx, idy,indiv_train, tidx, tidy,indiv_test, hideNeurons,shufTrainTestDataType,classfierType,noisefactor,detailResFile)
                    else:
                        resRandTreeTmp = computeAccuracyScore_HDOut(idx, idy,indiv_train,tidx, tidy,indiv_test, hideNeurons,shufTrainTestDataType,classfierType,noisefactor,detailResFile )

            resRandTreeTmp=np.append(resRandTreeTmp,np.zeros(nbInstances-len(resRandTreeTmp)))
            comm.Reduce(resRandTreeTmp, resRandTreesingle, op=MPI.SUM, root=0)

            if comm_rank == 0:
                resRandTree=resRandTree+resRandTreesingle


        if comm_rank == 0:
            resRandTree = resRandTree / len(trainsgroups)
            #os.remove(file_input)
            #if commandkey == "hdin":
            #print(resRandTree)
            #os.rename(file_input, file_input+str(countit)+".csv")

            #write file response to productor
            if commandkey== "hdinlogall":  #test for outter all population history only

                directory = os.path.dirname(pathLog)
                if not os.path.exists(directory):
                    os.makedirs(directory)

                with open(pathLog + "AllOUTTERScoreByGeneration.csv", "a") as file:
                    resRandTree.tofile(file, sep=",", format="%1.3f")  # , format="%s"
                    print("", file=file)
                print(", outters = {:1.3}".format(resRandTree[-1]), end='', flush=True)
            else:
                resRandTree.tofile(file_output+".tmp", sep="\n", format="%s")
                os.rename(file_output + ".tmp", file_output)

                #log score sorted of each indidu
                if (commandkey == "hdout")& (allResultHistoryBOOL==1):
                    directory = os.path.dirname(pathLog)
                    if not os.path.exists(directory):
                        os.makedirs(directory)
                    with open(pathLog + "AllScoreByGeneration.csv", "a") as file:
                    #with open(pathLog + "AllScoreByGenerationUnsort.csv", "a") as fileu:
                        np.sort(resRandTree).tofile(file, sep=",", format="%1.3f")#, format="%s"
                        #resRandTree.tofile(fileu, sep=",", format="%1.3f")
                        print("",file=file)
                        #print("", file=fileu)


        gc.collect()

#_______________________________________________________________________________
# Compute scores for data with Highly Duplicates Out
def computeAccuracyScore_HDOut( idx, idy,indiv_train, tidx, tidy,indiv_test,  hideNeurons,noise,classiType,stdfact ,log_details):

    #distrib_group = np.sort(np.unique(idg))
    # define number of test samples :  1/3 of the smallest group

    indivuniq = np.unique( indiv_train)
    #print("nbindiv train :"+str(indivuniq.size))
    resChunk = np.zeros(np.max(indivuniq)+1)
    groupuniq=np.unique(idy[indiv_train==indiv_train[0]])


    groupcount = np.bincount(idy[indiv_train==indiv_train[0]]).astype(float)
    groupcount[0] = np.nan
    mincatechant = (int)(np.nanmin(groupcount)) #nombre d'occurence minimum parmi les groupes
    nbtest = (int)(np.round(mincatechant / 3))
    nbtrain=(int)(mincatechant-nbtest)
    #subsetid=np.zeros(np.unique(idg)*mincatechant)
    if classiType == 'MLP':
        clf = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=hideNeurons, random_state=1,
                            learning_rate='adaptive', activation='logistic')
    else:
        if classiType == 'ADA':
            clf = AdaBoostClassifier()
        if classiType == 'SDG':
            clf = linear_model.SGDClassifier(class_weight='balanced')
        if classiType == 'SVC':
            clf = svm.SVC(gamma=.4,class_weight='balanced')
        if classiType == 'linSVC':
            clf = svm.LinearSVC(class_weight='balanced')
        if classiType == 'KNEIG':
            clf = KNeighborsClassifier(n_neighbors=5, algorithm='auto')


        # create pipeline from kernel approximation
        # and linear svm
        if classiType == 'SVMfourier':
            feature_map_fourier = RBFSampler(gamma=.2,n_components=hideNeurons)#, random_state=1

            clf = pipeline.Pipeline([("feature_map", feature_map_fourier),
                                                ("svm", svm.LinearSVC(class_weight='balanced'))])
        if classiType == 'SVMnystro':
            feature_map_nystroem = Nystroem(gamma=.2,n_components=hideNeurons)#, random_state=1
            clf = pipeline.Pipeline([("feature_map", feature_map_nystroem),
                                                 ("svm", svm.LinearSVC(class_weight='balanced'))])



            #a=idx+WhiteKernel(1e-1)
     #noise = np.random.normal(0, 1, idx.shape) / 100
    reschdetail = []
    for i in indivuniq:
        reschdetail_indiv = []
        resch = []
        j=1
        if noise==3:
            j=4
        #while((j<distrib_group.size) and (j<4)):

        #training with partial train data randomly choose 2/3train, 1/3test, 3 time
        #not use anymore
        while(False): #j<=3):
            idtrain = np.zeros((int)(nbtrain*groupuniq.size),dtype=np.int)
            idtest = np.zeros((int)(nbtest*groupuniq.size),dtype=np.int)
            jgrp=0
            for idgrp in groupuniq:
                idgrpsel = np.random.choice(np.where(idy[indiv_train==i] == idgrp)[0], mincatechant, replace=False)
                idgrpselidtrain = np.random.choice(range(idgrpsel.size), nbtrain, replace=False)
                idtrain[jgrp*nbtrain:(jgrp+1)*nbtrain ] = idgrpsel[idgrpselidtrain]
                idtest[ jgrp*nbtest:(jgrp+1)*nbtest ] = np.delete(idgrpsel, idgrpselidtrain)
                jgrp+=1

            #normtrain = normalize(idx.iloc[idtrain, table[i]], axis=1, copy=False)+np.random.normal(0, 1,[len(idtrain),len(table[i])])/50

            normtrain = normalize(idx.iloc[indiv_train==i, :].iloc[idtrain, :], axis=1, copy=False)
            normtest  = normalize(idx.iloc[indiv_train==i, :].iloc[idtest, :], axis=1, copy=False)
            if noise==1 :
                sepcount=0
                normtrain_org = normtrain
                normtest_org = normtest
            else:
                sepcount = 2
            while sepcount<3:
                if noise==2 :
                    # repeat and noise training
                    normtrain= addnoise(np.repeat(normtrain,3,axis=0),stdfact)
                        #,+ np.random.normal(0, 1, [len(idtrain)*3, len(table[i])]) *np.std(normtrain_org,axis=0)*stdfact
                    normtest = addnoise(np.repeat(normtest, 3, axis=0),stdfact)
                        # + np.random.normal(0, 1, [len(idtest) * 3, len(table[i])]) *np.std(normtrain_org,axis=0)*ftdfact

                    noise_trainidy = np.repeat(idy[idtrain],3)
                    noise_testidy = np.repeat(idy[idtest], 3)

                    clf.fit(normtrain, noise_trainidy)
                else:
                    if noise==1:
                        normtrain = addnoise(normtrain_org,stdfact)
                        #+ np.random.normal(0, 1, [len(idtrain), len(table[i])])*np.std(normtrain_org,axis=0)*ftdfact
                        normtest = addnoise(normtest_org,stdfact)
                        # + np.random.normal(0, 1, [len(idtest), len(table[i])]) *np.std(normtrain_org,axis=0)*ftdfact

                    clf.fit(normtrain, idy[idtrain])
                #normtrain = normalize(idx[idg!=distrib_group[j]].iloc[:,table[i]], axis=1, copy=False)
                #normtest = normalize(idx[idg==distrib_group[j]].iloc[:,table[i]], axis=1, copy=False)
                #clf.fit(normtrain, idy[idg!=distrib_group[j]])
                #data_y_predict = clf.predict(normtest)
                #resch.append(accuracy_score(idy[idg==distrib_group[j]], data_y_predict, normalize=True))

                #Calculate accuraccy score

                #data_y_predict = clf.predict(np.concatenate((normtest,normtrain), axis=0))
                data_y_predict = clf.predict(normtest)
                if noise == 2:

                    #resch.append(np.prod(accuracy_score(np.concatenate((noise_testidy, noise_trainidy), axis=0), data_y_predict,normalize=True) * (np.mean(np.amax(clf.predict_proba(np.concatenate((normtest, normtrain), axis=0)), axis=1)))))
                    #resch.append(np.prod(accuracy_score(noise_testidy, data_y_predict,normalize=True) * (
                     #   np.mean(np.amax(clf.predict_proba(normtest), axis=1)))))
                    resch.append(accuracy_score(noise_testidy, data_y_predict, normalize=True))

                    if log_details != "":
                        reschdetail_indiv.append(accuracy_score(noise_trainidy,clf.predict(normtrain)))
                        reschdetail_indiv.append(accuracy_score(noise_testidy,clf.predict(normtest)))

                else:# noise type = 1 and 0
                    #resch.append(np.prod(accuracy_score(np.concatenate((idy[idtest], idy[idtrain]), axis=0), data_y_predict,normalize=True) * ( np.mean(np.amax(clf.predict_proba(np.concatenate((normtest, normtrain), axis=0)), axis=1)))))
                    #resch.append(np.prod(accuracy_score(idy[idtest], data_y_predict,normalize=True) * (
                    #    np.mean(np.amax(clf.predict_proba(normtest), axis=1)))))
                    resch.append(accuracy_score(idy[idtest], data_y_predict, normalize=True))

                    if log_details != "":
                        reschdetail_indiv.append(accuracy_score(idy[idtrain],clf.predict(normtrain)))
                        reschdetail_indiv.append(accuracy_score(idy[idtest],clf.predict(normtest)))

                sepcount += 1
            j+=1
        '''
        del clf

        if classiType == 'MLP':
            clf = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=hideNeurons, random_state=1,
                                learning_rate='adaptive', activation='logistic')
        else:
            if classiType == 'ADA':
                clf = AdaBoostClassifier()
        '''

        #test with all datas available for training

        normtrain = normalize(idx.iloc[indiv_train==i, :], axis=1, copy=False)
        normtest = normalize(tidx.iloc[indiv_test==i, :], axis=1, copy=False)
        if noise == 1 :
            sepcount = 0
            normtrain_org= normtrain
            normtest_org = normtest
        else:
            sepcount = 2
        while sepcount < 3:

            if noise == 2:
                normtrain =  addnoise(np.repeat(normtrain, 3, axis=0),stdfact)
                normtest =  addnoise(np.repeat(normtest, 3, axis=0),stdfact)

                noise_trainidy = np.repeat(idy[indiv_train==i], 3)
                noise_testidy = np.repeat(tidy[indiv_test==i], 3)

                clf.fit(normtrain, noise_trainidy)
            else:
                if noise == 1:
                    normtrain =  addnoise(normtrain_org,stdfact)
                    normtest = addnoise( normtest_org ,stdfact)
                clf.fit(normtrain, idy[indiv_train==i])

            data_y_predict = clf.predict(normtest)
            if noise == 2:
                #resch.append(np.prod(accuracy_score(np.concatenate((noise_testidy,noise_trainidy), axis=0), data_y_predict, normalize=True) * (np.mean(np.amax(clf.predict_proba(np.concatenate((normtest, normtrain), axis=0)), axis=1)))))
                #resch.append(np.prod( accuracy_score(noise_testidy, data_y_predict,normalize=True) * (np.mean(np.amax(clf.predict_proba(normtest), axis=1)))))
                resch.append(accuracy_score(noise_testidy, data_y_predict, normalize=True))

                if log_details != "":
                    reschdetail_indiv.append(accuracy_score(noise_trainidy, clf.predict(normtrain)))
                    reschdetail_indiv.append(accuracy_score(noise_testidy, clf.predict(normtest)))
            else:
                #resch.append(np.prod(
                #    accuracy_score(np.concatenate((tidy[indiv_test==i], idy[indiv_train==i]), axis=0), data_y_predict, normalize=True) * (
                #    np.mean(np.amax(clf.predict_proba(np.concatenate((normtest, normtrain), axis=0)), axis=1)))))
                #resch.append(np.prod(accuracy_score(tidy[indiv_test == i], data_y_predict )* (
                #    np.mean(np.amax(clf.predict_proba(normtest), axis=1)))))
                resch.append(accuracy_score(tidy[indiv_test == i], data_y_predict))

                if log_details != "":
                    reschdetail_indiv.append(accuracy_score(idy[indiv_train==i], clf.predict(normtrain)))
                    reschdetail_indiv.append(accuracy_score(tidy[indiv_test==i], clf.predict(normtest)))


            sepcount += 1



        resChunk[i] = np.mean(resch)
        if log_details != "":
            reschdetail.append(reschdetail_indiv)

    if log_details != "":
        # reschdetail.tocvs(log_details+'_'+str(hideNeurons)+'_'+str(noise)+'_'+classiType+'_'+str(stdfact))
        # f_handle = np.file(log_details+'_'+str(hideNeurons)+'_'+str(noise)+'_'+classiType+'_'+str(stdfact), 'ab')
        with open(log_details + '_' + str(hideNeurons) + '_' + str(noise) + '_' + classiType + '_' + str(stdfact),
                  'ab') as f:
            np.savetxt(f, np.asarray(reschdetail, dtype='float'), delimiter=',', newline='\n', fmt='%.10f')

    return resChunk

def KERAS_computeAccuracyScore_HDOut( idx, idy,indiv_train, tidx, tidy,indiv_test,  hideNeurons,noise,classiType,stdfact ,log_details):

    #distrib_group = np.sort(np.unique(idg))
    # define number of test samples :  1/3 of the smallest group

    #t0 = time.time() #debug


    idy = idy - np.min(idy)
   # print("tidy")
   # print(tidy)
    tidy = tidy - np.min(tidy)
    indivuniq = np.unique( indiv_train)
    resChunk = np.zeros(max(indivuniq)+1)
    groupuniq=np.unique(idy[indiv_train==indiv_train[0]])


    groupcount = np.bincount(idy[indiv_train==indiv_train[0]]).astype(float)
    groupcount[0] = np.nan
    mincatechant = (int)(np.nanmin(groupcount)) #nombre d'occurence minimum parmi les groupes
    nbtest = (int)(np.round(mincatechant / 3))
    nbtrain = (int)(mincatechant-nbtest)

    #model declaration
    clf = Sequential()
    if (noise==1) | (noise==2):
        clf.add(GaussianNoise(stdfact,input_shape=(idx.shape[1],)))

    clf.add(Dense(hideNeurons, input_dim=idx.shape[1], activation='relu'))
    if "_BN" in classfierType:
        clf.add(normalization.BatchNormalization())
    if "_DO" in classfierType:
        if "_DOs" in classfierType:
            clf.add(Dropout(0.9))
        else:
            clf.add(Dropout(0.75))
    if "_2L" in classfierType:
        clf.add(Dense(hideNeurons, input_dim=idx.shape[1], activation='relu'))
        if "_BN" in classfierType:
            clf.add(normalization.BatchNormalization())
        if "_DO" in classfierType:
            if "_DOs" in classfierType:
                clf.add(Dropout(0.9))
            else:
                clf.add(Dropout(0.75))

    # model.add(Activation('relu'))
    # model.add(Dense(1,activation='sigmoid'))
    clf.add(Dense(len(groupuniq), activation='softmax'))
    learning_rate = 0.1
    sgd = optimizers.SGD(lr=learning_rate, momentum=0.5, decay=0.0, nesterov=False)
    if len(groupuniq > 2):
        clf.compile(optimizer=sgd, loss='categorical_crossentropy', metrics=['accuracy'])
    else:
        clf.compile(optimizer=sgd, loss='categorical_crossentropy', metrics=['accuracy'])
        #clf.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
    Wsave = clf.get_weights()


    earlystop = EarlyStopping(monitor='loss', min_delta=0.001, patience=10,
                              verbose=0, mode='auto')
    #EarlyStop = EarlyStoppingTestData(testdata=[],
     #            testlabel=tidy,idmetric=1,monitor='loss', min_delta=0.01, patience=10,
      #                        verbose=0, mode='auto')

    initial_lratecopy=0.1
    if "_LR2" in classfierType:
        initial_lratecopy = 0.2
    if "_LR3" in classfierType:
        initial_lratecopy = 0.3
        print(initial_lratecopy)
    def step_decay(epoch):
        initial_lrate = initial_lratecopy
        ##print(initial_lrate)
        drop = 0.5
        epochs_drop = 10.0
        lrate = initial_lrate * math.pow(drop, math.floor((epoch)/epochs_drop))
        return lrate
    

    lrate = LearningRateScheduler(step_decay)

    callbacks_list = [earlystop]
    reschdetail=[]

    for i in indivuniq:
        resch = []
        reschdetail_indiv=[]
        j=1
        #while((j<distrib_group.size) and (j<4)):

        #training with partial train data randomly choose 2/3train, 1/3test, 3 time
        while(j<=3):
            idtrain = np.zeros((int)(nbtrain*groupuniq.size),dtype=np.int)
            idtest = np.zeros((int)(nbtest*groupuniq.size),dtype=np.int)
            jgrp=0
            for idgrp in groupuniq:
                idgrpsel = np.random.choice(np.where(idy[indiv_train==i] == idgrp)[0], mincatechant, replace=False)
                idgrpselidtrain = np.random.choice(range(idgrpsel.size), nbtrain, replace=False)
                idtrain[jgrp*nbtrain:(jgrp+1)*nbtrain ] = idgrpsel[idgrpselidtrain]
                idtest[ jgrp*nbtest:(jgrp+1)*nbtest ] = np.delete(idgrpsel, idgrpselidtrain)
                jgrp+=1

            #normtrain = normalize(idx.iloc[idtrain, table[i]], axis=1, copy=False)+np.random.normal(0, 1,[len(idtrain),len(table[i])])/50


            normtrain = normalize(idx.iloc[indiv_train == i, :].iloc[idtrain, :], axis=1, copy=False)
            normtest = normalize(idx.iloc[indiv_train == i, :].iloc[idtest, :], axis=1, copy=False)

            if noise==1 :
                sepcount=0
            else:
                sepcount = 2
            while sepcount<3:
                if noise==2 :
                    # repeat and noise training
                    normtrain= np.repeat(normtrain,3,axis=0)
                    normtest = np.repeat(normtest, 3, axis=0)

                    trainidy = np.repeat(idy[idtrain],3)
                    testidy = np.repeat(idy[idtest], 3)
                else:
                    trainidy=idy[idtrain]
                    testidy = idy[idtest]
                clf.set_weights(Wsave)
                clf.fit(normtrain, to_categorical(trainidy), epochs=450,shuffle=True, batch_size=np.min([500, len(trainidy)]),
                        callbacks=callbacks_list, verbose=0)


                #model_info =clf.fit(.....
                #plot_model_history(model_info)
                #Calculate accuraccy score


                #resch.append(clf.evaluate(np.concatenate((normtest,normtrain), axis=0),  to_categorical(np.concatenate((testidy, trainidy), axis=0)),verbose=0)[1])
                resch.append(clf.evaluate(normtest, to_categorical(testidy),verbose=0)[1])

                if log_details!= "":
                    reschdetail_indiv.append(clf.evaluate(normtrain,  to_categorical(trainidy),verbose=0)[1])
                    reschdetail_indiv.append(clf.evaluate(normtest , to_categorical(testidy),verbose=0)[1])


                sepcount += 1
            j+=1

        #####################
        #test with all datas available for training
        normtrain = normalize(idx.iloc[indiv_train==i, :], axis=1, copy=False)
        normtest = normalize(tidx.iloc[indiv_test==i, :], axis=1, copy=False)
        if noise == 1:
            sepcount = 0

        else:
            sepcount = 2
        while sepcount < 3:

            if noise == 2:
                normtrain =  np.repeat(normtrain, 3, axis=0)
                normtest =  np.repeat(normtest, 3, axis=0)

                trainidy = np.repeat(idy[indiv_train==i], 3)
                testidy = np.repeat(tidy[indiv_test==i], 3)
            else:
                trainidy = idy[indiv_train == i]
                testidy = tidy[indiv_test == i]
            clf.set_weights(Wsave)
            classweight = class_weight.compute_class_weight('balanced', np.unique(trainidy), trainidy)
            '''
            earlystop2 = EarlyStoppingTestData(testdata=normtest,
                                              testlabel=testidy, idmetric=1, monitor='acc', min_delta=0.01, patience=11100,
                                              verbose=1, mode='auto')
            callbacks_list = [earlystop2]
            '''
            callbacks_list = [earlystop]
            clf.fit(normtrain, to_categorical(trainidy), epochs=450,shuffle=True, batch_size=np.min([2000,len(trainidy)]), callbacks=callbacks_list,verbose=0,class_weight=classweight)
            #history =
            #plot_model_history_callback(history,earlystop2.history)
            #resch.append(clf.evaluate(np.concatenate((normtest, normtrain), axis=0),  to_categorical(np.concatenate((testidy, trainidy), axis=0)),verbose=0)[1])
            resch.append(clf.evaluate(normtest, to_categorical(testidy), verbose=0)[1])
            if log_details!= "":
                reschdetail_indiv.append(clf.evaluate(normtrain, to_categorical(trainidy),verbose=0)[1])
                reschdetail_indiv.append(clf.evaluate(normtest, to_categorical(testidy),verbose=0)[1])

            sepcount += 1

        resChunk[i] = np.mean(resch)
        if log_details != "":
            reschdetail.append(reschdetail_indiv)

    if log_details != "":
        #reschdetail.tocvs(log_details+'_'+str(hideNeurons)+'_'+str(noise)+'_'+classiType+'_'+str(stdfact))
        #f_handle = np.file(log_details+'_'+str(hideNeurons)+'_'+str(noise)+'_'+classiType+'_'+str(stdfact), 'ab')
        with open(log_details+'_hn'+str(hideNeurons)+'_n'+str(noise)+'_'+classiType+'_s'+str(stdfact), 'ab') as f:
            np.savetxt(f, np.asarray(reschdetail,dtype='float'),delimiter=',',newline='\n',fmt='%.10f')


    return resChunk

#_______________________________________________________________________________
# Compute scores for Highly Duplicates that are In main data
def computeAccuracyScore_HDIn( idx, idy,indiv_train, tidx, tidy,indiv_test,  hideNeurons,classiType):
    #distrib_group = np.sort(np.unique(idg))
    #global nbcallin
    indivuniq = np.unique(indiv_train)

    resChunk = np.zeros(max(indivuniq) + 1)

    if classiType == 'MLP':
        clf = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=hideNeurons, random_state=1,
                            learning_rate='adaptive', activation='logistic')
    else:
        if classiType == 'ADA':
            clf = AdaBoostClassifier()
        if classiType == 'SDG':
            clf = linear_model.SGDClassifier(class_weight='balanced')
        if classiType == 'SVC':
            clf = svm.SVC(gamma=.4,class_weight='balanced')#,tol=hideNeurons)
        if classiType == 'linSVC':
            clf = svm.LinearSVC(class_weight='balanced')
        if classiType == 'KNEIG':
            clf = KNeighborsClassifier(n_neighbors=5,algorithm='auto')
        if classiType == 'SVMfourier':
            feature_map_fourier = RBFSampler(gamma=.2,n_components=hideNeurons)#, random_state=1

            clf = pipeline.Pipeline([("feature_map", feature_map_fourier),
                                     ("svm", svm.LinearSVC(class_weight='balanced'))])
        if classiType == 'SVMnystro':
            feature_map_nystroem = Nystroem(gamma=.2,n_components=hideNeurons)#, random_state=1
            clf = pipeline.Pipeline([("feature_map", feature_map_nystroem),
                                     ("svm", svm.LinearSVC(class_weight='balanced'))])


    for i in indivuniq:
        normtrain = normalize(idx.iloc[indiv_train==i, :], axis=1, copy=False)
        normtest = normalize(tidx.iloc[indiv_test==i,:], axis=1, copy=False)


        clf.fit(normtrain, idy[indiv_train==i])
        data_y_predict = clf.predict(normtest)

        resChunk[i] =  accuracy_score(tidy[indiv_test==i], data_y_predict, normalize=True)
        #print("hdin",nbcallin," i ",i)
        #resChunk[i] =i+(100* nbcallin)
    #print(resChunk)
    return resChunk

#-------------------------------------------------------------------------------
def readBinRequest(filename):


    f = open(filename, "rb")
    sldib = np.fromfile(f, dtype=np.uint32,count=1)
    nblines = np.fromfile(f, dtype=np.uint32,count=1)
    nbcolumns = np.fromfile(f, dtype=np.uint32,count=1)
    #data = np.fromfile(f, dtype=sldib,count=-1)
    #data.reshape([nblines,nbcolumns])
    #dt = np.dtype("("+str(nblines)+","+str(nbcolumns)+")f"+str(sldib))
    #dt=np.dtype(str(sldib),shape=(nblines,nbcolumns))

   # dt = np.dtype('f' + str(int(np.sqrt(sldib[0]))), (nblines, nbcolumns))

    #dt = np.dtype('f' + str(sldib[0])+','+str( nbcolumns[0]))
    #dt = np.dtype('f16')#' + str(sldib[0]))

    #print("nblines : ",nblines,", nbcolumns :",nbcolumns,",sldib : ",sldib)
    dt = np.dtype('(' + str(nblines[0]) + ',' + str(nbcolumns[0]) + ')' + 'f' + str(sldib[0]))
    #dt = np.dtype( 'f' + str(sldib[0]))

    data = pd.DataFrame(np.fromfile(f, dtype=dt)[0]) #.reshape((nblines[0],nbcolumns[0])))
    #data = data.astype(int)
    f.close()
    return data

    '''
    with open(filename, "rb") as binary_file:
        # Read the whole file at once
        #data = binary_file.read()
        #print(data)

        # Seek position and read N bytes
        #binary_file.seek(0)  # Go to beginning
        sldib = binary_file.read(2)
        nblines = binary_file.read(2)
        nbcolumns = binary_file.read(2)

        data = binary_file.read()
        print(couple_bytes)
    '''

#_______________________________________________________________________________
# Compute scores for Highly Duplicates that are In main data
def KERAS_computeAccuracyScore_HDIn(idx, idy,indiv_train, tidx, tidy,indiv_test, hideNeurons, classiType):
    #print("idy")
    #print(idy)
    #print("tidy")
    #print(tidy)
    idy = idy - np.min(idy)
    tidy = tidy - np.min(tidy)
    indivuniq = np.unique(indiv_train)
    resChunk = np.zeros(max(indivuniq) + 1)
    groupuniq = np.unique(idy)
    clf = Sequential()
    #if (noise == 1) | (noise == 2):
    #    clf.add(GaussianNoise(stdfact, input_shape=(idx.shape[1],)))

    clf.add(Dense(hideNeurons, input_dim=idx.shape[1], activation='relu'))
    if "_BN" in classfierType:
        clf.add(normalization.BatchNormalization())
    if "_DO" in classfierType:
        if "_DOs" in classfierType:
            clf.add(Dropout(0.9))
        else:
            clf.add(Dropout(0.75))
    if "_2L" in classfierType:
        clf.add(Dense(hideNeurons, input_dim=idx.shape[1], activation='relu'))
        if "_BN" in classfierType:
            clf.add(normalization.BatchNormalization())
        if "_DO" in classfierType:
            if "_DOs" in classfierType:
                clf.add(Dropout(0.9))
            else:
                clf.add(Dropout(0.75))

    # model.add(Activation('relu'))
    # model.add(Dense(1,activation='sigmoid'))
    clf.add(Dense(len(groupuniq), activation='softmax'))
    learning_rate = 0.1
    sgd = optimizers.SGD(lr=learning_rate, momentum=0.5, decay=0.0, nesterov=False)
    if len(groupuniq > 2):
        clf.compile(optimizer=sgd, loss='categorical_crossentropy', metrics=['accuracy'])
    else:
        clf.compile(optimizer=sgd, loss='categorical_crossentropy', metrics=['accuracy'])
        # clf.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
    Wsave = clf.get_weights()

    earlystop = EarlyStopping(monitor='loss', min_delta=0.001, patience=10,
                              verbose=0, mode='auto')
    # EarlyStop = EarlyStoppingTestData(testdata=[],
    #            testlabel=tidy,idmetric=1,monitor='loss', min_delta=0.01, patience=10,
    #                        verbose=0, mode='auto')

    initial_lratecopy = 0.1
    if "_LR2" in classfierType:
        initial_lratecopy = 0.2
    if "_LR3" in classfierType:
        initial_lratecopy = 0.3
        print(initial_lratecopy)

    def step_decay(epoch):
        initial_lrate = initial_lratecopy
        ##print(initial_lrate)
        drop = 0.5
        epochs_drop = 10.0
        lrate = initial_lrate * math.pow(drop, math.floor((epoch) / epochs_drop))
        return lrate

    lrate = LearningRateScheduler(step_decay)

    callbacks_list = [earlystop]

    for i in indivuniq:
        normtrain = normalize(idx.iloc[indiv_train == i, :], axis=1, copy=False)
        normtest = normalize(tidx.iloc[indiv_test == i, :], axis=1, copy=False)
        classweight = class_weight.compute_class_weight('balanced', np.unique(idy[indiv_train == i]),
                                                        idy[indiv_train == i])
        clf.set_weights(Wsave)

        clf.fit(normtrain, to_categorical(idy[indiv_train==i]),shuffle=True, epochs=450, batch_size=np.min([2000, len(normtrain)]),
                callbacks=callbacks_list, verbose=0, class_weight=classweight)
        # model.fit(idx.iloc[:, table[i]], idy ,batch_size=10,epochs=100)#epochs=100, batch_size=
        #data_y_predict
        resChunk[i]  = clf.evaluate(normtest,to_categorical(tidy[indiv_test==i]),verbose=0)[1]

         #accuracy_score(tidy[indiv_test == i], data_y_predict, normalize=True)

    return resChunk

# _______________________________________________________________________________
# add noise to the data
def addnoise(table, noisefactor):

    return table + (np.random.normal(0, 1, table.shape) * np.std(table, axis=0) * noisefactor)
# _______________________________________________________________________________
def plot_model_history(model_history):
    fig, axs = plt.subplots(1,2,figsize=(15,5))
    # summarize history for accuracy
    axs[0].plot(range(1,len(model_history.history['acc'])+1),model_history.history['acc'])
    if 'val_acc' in model_history.history:
        axs[0].plot(range(1,len(model_history.history['val_acc'])+1),model_history.history['val_acc'])
    axs[0].set_title('Model Accuracy')
    axs[0].set_ylabel('Accuracy')
    axs[0].set_xlabel('Epoch')
    axs[0].set_xticks(np.arange(1,len(model_history.history['acc'])+1),len(model_history.history['acc'])/10)
    axs[0].legend(['train', 'val'], loc='best')
    # summarize history for loss
    axs[1].plot(range(1,len(model_history.history['loss'])+1),model_history.history['loss'])
    if 'val_loss' in model_history.history:
        axs[1].plot(range(1,len(model_history.history['val_loss'])+1),model_history.history['val_loss'])
    axs[1].set_title('Model Loss')
    axs[1].set_ylabel('Loss')
    axs[1].set_xlabel('Epoch')
    axs[1].set_xticks(np.arange(1,len(model_history.history['loss'])+1),len(model_history.history['loss'])/10)
    axs[1].legend(['train', 'val'], loc='best')
    #plt.savefig(figdir + "/plot_model_history", dpi=300)
    plt.show()

def plot_model_history_callback(model_history,model_historytest):
    fig, axs = plt.subplots(1,3,figsize=(15,5))
    # summarize history for accuracy
    axs[0].plot(range(1,len(model_history.history['acc'])+1),model_history.history['acc'])
    if 'val_acc' in model_history.history:
        axs[0].plot(range(1,len(model_history.history['val_acc'])+1),model_history.history['val_acc'])
    axs[0].set_title('Model Accuracy')
    axs[0].set_ylabel('Accuracy')
    axs[0].set_xlabel('Epoch')
    axs[0].set_xticks(np.arange(1,len(model_history.history['acc'])+1),len(model_history.history['acc'])/10)
    axs[0].legend(['train', 'val'], loc='best')
    # summarize history for loss
    axs[1].plot(range(1,len(model_history.history['loss'])+1),model_history.history['loss'])
    if 'val_loss' in model_history.history:
        axs[1].plot(range(1,len(model_history.history['val_loss'])+1),model_history.history['val_loss'])
    axs[1].set_title('Model Loss')
    axs[1].set_ylabel('Loss')
    axs[1].set_xlabel('Epoch')
    axs[1].set_xticks(np.arange(1,len(model_history.history['loss'])+1),len(model_history.history['loss'])/10)
    axs[1].legend(['train', 'val'], loc='best')

    axs[2].plot(range(1, len(model_historytest) + 1), model_historytest)

    axs[2].set_title('Model Accuracy test')
    axs[2].set_ylabel('Accuracy')
    axs[2].set_xlabel('Epoch')
    axs[2].set_xticks(np.arange(1, len(model_historytest) + 1), len(model_historytest) / 10)
    axs[2].legend(['train', 'val'], loc='best')

    #plt.savefig(figdir + "/plot_model_history", dpi=300)
    plt.show()
######## snippet compare scikitlearn and KERAS
def comparelib(clf,idx, idy,indiv_train, tidx, tidy,indiv_test, hideNeurons,normtrain,normtest,trainidy,testidy):
    earlystop = EarlyStopping(monitor='acc', min_delta=0, patience=2, verbose=0, mode='auto')
    callbacks_list = [earlystop]
    clfsci = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=hideNeurons, random_state=1,
                           learning_rate='adaptive', activation='logistic')
    t0 = time.time()
    clfsci.fit(normtrain, trainidy)
    t1 = time.time()
    data_y_predict = clfsci.predict(np.concatenate((normtest, normtrain), axis=0))
    print("scikit " + str(t1 - t0) + "s")

    print(accuracy_score(np.concatenate((testidy, trainidy), axis=0), data_y_predict, normalize=True))
    t2 = time.time()
    model_info = clf.fit(normtrain, to_categorical(trainidy), epochs=150, batch_size=len(trainidy), shuffle=True,
                         callbacks=callbacks_list,
                         verbose=0)
    t3 = time.time()

    plot_model_history(model_info)
    print("KERAS " + str(t3 - t2) + "s")
    print(clf.evaluate(np.concatenate((normtest, normtrain), axis=0),
                       to_categorical(np.concatenate((testidy, trainidy), axis=0)), verbose=0)[1])

#######################################################################################
idapp = 0
hideNeurons = 0
shufTrainTestDataType = 0
classfierType = 0
noisefactor = 0
detailResFile =""
conffile="BEAUTYtrinegRes/configbase.conf"
allResultHistoryBOOL=0
pathLog=""
# parameters : idapp  hideNeurons  shufTrainTestDataType  classfierType  noisefactor  detailResFile
if len(sys.argv) > 1:
    idapp = sys.argv[1] # id application for communication
if len(sys.argv) > 2:
    conffile= sys.argv[2]

#params=np.genfromtxt(conffile,delimiter='=')
txtconf=np.loadtxt(conffile,delimiter='=',comments='#',dtype='S')
for i in txtconf:
    key=i[0].decode("utf-8").lower()
    if key =="hiddenneuron":
        hideNeurons=int(i[1])
    if key =="shufTrainTestDataType".lower():
        shufTrainTestDataType=int(i[1])
    if key =="noisefactor":
        noisefactor=float(i[1])
    if key == "method":
        classfierType = i[1].decode("utf-8")
    if key == "detailResFile".lower():
        detailResFile = i[1].decode("utf-8")
    if key == "pathLog".lower():
        pathLog = i[1].decode("utf-8")+"Dir/"

    if key == "allResultHistoryBOOL".lower():
        allResultHistoryBOOL = int(i[1])

print("idapp:",idapp,",conffile:",conffile,",hiddenNeurons:",hideNeurons, ",shufTrainTestDataType:",shufTrainTestDataType," ,classfierType:",classfierType," ,noisefactor:",noisefactor," ,(detailResFile):",detailResFile)


main(idapp,hideNeurons,shufTrainTestDataType,classfierType,noisefactor,detailResFile)




