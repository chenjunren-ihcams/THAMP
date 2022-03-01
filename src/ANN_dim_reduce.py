from sklearn.metrics import roc_auc_score
from sklearn import metrics
from sklearn.preprocessing import label_binarize
from keras.models import Sequential
from keras.layers import Dense
from keras.models import Model
from keras.models import load_model
from sklearn.model_selection import StratifiedKFold
import tensorflow as tf
from sklearn.utils import resample
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.decomposition import PCA
import joblib

#input_data
input_file="../data-2022-02-28.xlsx"
#patient_total
patient_total=1192

default_method="main analysis"
#default_method="dim_select"

if default_method == "main analysis":
    #if boot_stapping
    if_bootstrap=1
    bs_start=0
    total_bs_rounds=100
    # if cross_validation_test
    if_cv=0
    # if cross_validation_test and output auc for each diagnosis
    detailed_auc_trigger=0
    #neurons to iterate
    neu_start=65
    neu_end=neu_start+1

if default_method == "dim_select":
    #if boot_stapping
    if_bootstrap=0
    bs_start=0
    total_bs_rounds=1
    # if cross_validation_test
    if_cv=0
    # if cross_validation_test and output auc for each diagnosis
    detailed_auc_trigger=0
    #neurons to iterate
    neu_start=30
    neu_end=100

#stablization purpose
seed=4
#epochs
epo=100
#pre_defined_mapping of diag
n_classes = 8
labels = [0, 1, 2, 3, 4, 5, 6, 7]
id_to_disease ={0:"MPN",
    1:"MDS",
    2:"MDS/MPN",
    3:"CLL",
    4:"AML",
    5:"ALL",
    6:"MM",
    7:"CML"}
cv_out_path="../output/cv_stats/"

np.random.seed(seed)
tf.random.set_seed(seed)
auc_for_each_diag = {}

for diag_id in range(n_classes):
    auc_for_each_diag[id_to_disease[diag_id]] = []

def mkdir(pa):
    if not os.path.exists(pa):
        os.makedirs(pa)
mkdir("../output/cv_stats")
mkdir("../output/bs_data")
mkdir("../output/hidden_layer_out")
mkdir("../output/2d_full_data")
mkdir("../output/big_graph")
mkdir("../output/dim_sel")
mkdir("../output/model")


bs_cv_result_file = open(cv_out_path+"cv_result.csv", "w")
bs_cv_result_file.write("hidden_unit,boostrapping_round,overall_auc,auc1,auc2,auc3,auc4,auc5\n")

#generate PCA projections and draw the patient graph(raw version)
def pca_and_pic(filename,neu_units,sampled_bs_file,bs_round):
    hidden_layer_out_file=filename
    full_name=filename.split("/")[-1]
    output_file="../output/2d_full_data/"+full_name.split(".xlsx")[0]+"_2d.xlsx"

    #read hidden layer
    df = pd.read_excel(hidden_layer_out_file)
    x=df
    x=np.array(x)
    #PCA projection
    pca_model= PCA(n_components=2,random_state=2)
    X_reduced = pca_model.fit_transform(x)
    #save pca parameters
    joblib.dump(pca_model, "../output/model/{}_bs{}_pca.m".format(neu_units,bs_round))

    data=pd.DataFrame(X_reduced)
    data.columns=["0","1"]
    #add to the end of the input data for later analysis
    all_data=pd.read_excel(sampled_bs_file).copy(deep=True)
    all_data["X"]=data["0"]
    all_data["Y"]=data["1"]
    all_data.to_excel(output_file)
    #generate pictures
    gen_pic(output_file,neu_units,bs_round)


#draw the patient graph(raw version)
def gen_pic(filename,neu_units,bs_round):
    df = pd.read_excel(filename)
    id_to_disease={
    0:"MPN",
    1:"MDS",
    2:"MDS/MPN",
    3:"CLL",
    4:"AML",
    5:"ALL",
    6:"MM",
    7:"CML"}
    plt.clf()
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.set_facecolor('white')
    plt.rcParams['font.sans-serif'] = ['Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False
    df2=df
    xx=0.3
    d=df2[df2["disease_id"]==0]
    plt.scatter(d["X"],d["Y"],c='r',s=xx,label = id_to_disease[0])
    d=df2[df2["disease_id"]==1]
    plt.scatter(d["X"],d["Y"],c='g',s=xx,label = id_to_disease[1])
    d=df2[df2["disease_id"]==2]
    plt.scatter(d["X"],d["Y"],c='b',s=xx,label = id_to_disease[2])
    d=df2[df2["disease_id"]==3]
    plt.scatter(d["X"],d["Y"],c='lime',s=xx,label = id_to_disease[3])
    d=df2[df2["disease_id"]==4]
    plt.scatter(d["X"],d["Y"],c='#fac205',s=xx,label = id_to_disease[4])
    d=df2[df2["disease_id"]==5]
    plt.scatter(d["X"],d["Y"],c='#a00498',s=xx,label = id_to_disease[5])
    d=df2[df2["disease_id"]==6]
    plt.scatter(d["X"],d["Y"],c='#3af1fe',s=xx,label = id_to_disease[6])
    d=df2[df2["disease_id"]==7]
    plt.scatter(d["X"],d["Y"],c='hotpink',s=xx,label = id_to_disease[7])
    plt.legend(markerscale=3,loc="upper right", prop={'size': 7})
    plt.title("overview -{}".format(neu_units))
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.savefig('../output/big_graph/overview_100_{}_bs_{}.jpg'.format(neu_units,bs_round), dpi=500)
    plt.close()
#main process
for neu_units in range(neu_start,neu_end):
    #bootstrapping
    for current_bs_round in range(bs_start,total_bs_rounds):
        input_data = pd.read_excel(input_file)
        input_data["disease_id"]=input_data["Diagnosis"].replace({value:key for key, value in id_to_disease.items()})
        sampled_bootstrapping_file = "../output/bs_data/{}.xlsx".format(current_bs_round)
        if if_bootstrap and current_bs_round!=0:
            bootstrapSamples = resample(input_data, n_samples=patient_total, replace=1)
            bootstrapSamples.to_excel(sampled_bootstrapping_file,index=False)
            all_data=pd.read_excel(sampled_bootstrapping_file)
        else:
            all_data=input_data
            input_data.to_excel(sampled_bootstrapping_file, index=False)
        #generate training data
        X_train = all_data
        X_train = X_train.drop(["ID"], axis=1)
        X_train = X_train.drop(["Diagnosis"], axis=1)
        X_train = X_train.drop(["Sex"], axis=1)
        X_train = X_train.drop(["disease_id"], axis=1)
        X_train = X_train.drop(["MRD (for AML)"], axis=1)
        X_train = X_train.drop(["2017 ELN risk score (for AML)"], axis=1)
        X_train = X_train.drop(["M3 (for AML)"], axis=1)
        X_train = X_train.drop(["2016 WHO category (for MDS)"], axis=1)
        X_train = X_train.drop(["IPSS-R score (for MDS)"], axis=1)
        X_train = X_train.drop(["CR after 6 cycles of chemotherapy (for MDS)"], axis=1)
        X_train = X_train.drop(["AML type"], axis=1)
        X_train = X_train.drop(["treatment"], axis=1)
        
        y_train=all_data["disease_id"]
        #test data is the same since we use all patients to generate the following projections
        X_test = X_train.copy(deep=True)
        y_test = y_train.copy(deep=True)
        #data for cross validation purpose
        cv_train=X_train
        cv_test=y_train
        #re-format y values to match softmax output
        y_train = pd.get_dummies(y_train)
        # train and save the model
        model = Sequential()
        model.add(Dense(100, activation='relu', input_dim=339))
        model.add(Dense(neu_units, activation='relu', input_dim=100))
        model.add(Dense(8, activation='softmax'))
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
        model.fit(X_train, y_train, epochs=epo, batch_size=10,shuffle=False,verbose=0)
        model.save("../output/model/{}_bs{}_mod.h5".format(neu_units,current_bs_round))
        #whether to perform cross-validation on the model
        res = []
        if if_cv or (current_bs_round==0 and default_method=="main analysis"):
            kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=7)
            #for cross-validation, the original data will be split into 5 folds.
            for train, test in kfold.split(cv_train, cv_test):
                # create model
                modelx = Sequential()
                cv_train = np.array(cv_train)
                cv_test = np.array(cv_test)
                modelx.add(Dense(100, activation='relu', input_dim=339))
                modelx.add(Dense(neu_units , activation='relu', input_dim=100))
                modelx.add(Dense(8, activation='softmax'))
                modelx.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
                yt = pd.get_dummies(cv_test[train])
                modelx.fit(cv_train[train], yt, epochs=epo, batch_size=10,verbose=0)

                #to match the format
                predx = modelx.predict(cv_train[test])
                predx = [np.argmax(i) for i in predx]
                ytt = np.array(list(cv_test[test]))
                preds = label_binarize(predx, classes=labels)
                ypreds = label_binarize(ytt, classes=labels)
                #calculate AUC
                res2 = roc_auc_score(ypreds, preds, average='weighted', multi_class='ovo')
                res.append(res2)
                #whether do cross validation for each diagnosis
                if detailed_auc_trigger or current_bs_round==0:
                    for p in range(len(labels)):
                        tt = np.array([elem[p] for elem in list(ypreds)])
                        pred_y_2 = np.array([elem[p] for elem in list(preds)])
                        fpr, tpr, thresholds = metrics.roc_curve(tt,pred_y_2, pos_label=1)
                        auroc = round(metrics.auc(fpr, tpr), 2)
                        auc_for_each_diag[id_to_disease[p]].append(auroc)
            #general AUC for cross validation
            bs_cv_result_file.write(str(neu_units ) + "," +str(current_bs_round)+"," + str(sum(res) / 5) + "," + ",".join(str(vall) for vall in res) + "\n")
            #AUCrecord for cross validation on each diagnosis
            if detailed_auc_trigger or current_bs_round==0:
                cv_out=open(cv_out_path+"units_{}_bs_{}.csv".format(neu_units,current_bs_round),"w")
                cv_out.write("diag,round1,round2,round3,round4,round5,average\n")
                cv_out.write("Weighted AUC,"+",".join(str(vall) for vall in res)+","+str(sum(res) / 5)+"\n")
                for each_diag in auc_for_each_diag:
                    cv_out.write(each_diag+","+",".join(str(vall) for vall in auc_for_each_diag[each_diag])+","+str(sum(auc_for_each_diag[each_diag]) / 5) +"\n")
                cv_out.close()

        #load the model to get hidden units for each patient
        model = load_model("../output/model/{}_bs{}_mod.h5".format(neu_units,current_bs_round))
        m1 = Model(inputs=model.input, outputs=model.layers[1].output)
        p = m1.predict(X_test)
        filename=("../output/hidden_layer_out/encode_100_{}_bs_{}.xlsx").format(neu_units,current_bs_round)
        #make sure to have index=False
        pd.DataFrame(p).to_excel(filename, index=False)
        #formatting the data to generate AUC for current round. The score is meaningless since we use all the patients to train and test
        pred = model.predict(X_test)
        pred = [np.argmax(i) for i in pred]
        y_test = np.array(list(y_test))
        pred = np.array([int(i) for i in list(pred)])
        preds = label_binarize(pred, classes=labels)
        ypreds = label_binarize(y_test, classes=labels)
        res = roc_auc_score(ypreds, preds, average='weighted', multi_class='ovo')
        print("hidden unit:{}, round: {}".format(neu_units,current_bs_round))
        #get the PCA output and the patient graph.
        pca_and_pic(filename, neu_units, sampled_bootstrapping_file, current_bs_round)

bs_cv_result_file.close()





