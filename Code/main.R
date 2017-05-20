#######################################################################################
########################################################################################
###
### Script to select, make and train a model for the classification of different
### cancer types. First the data is read in to perform a data inspection.... TODO FINISH  
### .....
### 
### Dependecies: misc.R, feature_selection.R and model_tuning.R
### input: Train_call.txt, Train_clinical.txt andunlabelled_sample.txt
### output: .... models
###
###
########################################################################################
########################################################################################

#Pulling from R

# clear environment and load source files
source("code/misc.R")

load_depend()
# get the train data and write to WEKA file
combine = Load_labeled_Data()

########################################################################################


# Feature selection by importance ()
features = filter_Var_selection(combine,30)

# use rfe featureselection (not returningfeatures still)
#features_rfe = feature_selection_rfe(combine)

# Advance feature selction using WEKA. Algorithm used describe in the report
CorrelationAttribute <- c("V2185","V2214","V2215","V2211","V1657","V1679",
                          "V675","V1668","V674","V1680","V2207","V2224","V2224","V1673",
                          "V673","V1669","V1665","V2220","V693","V696","V855","V1678",
                          "V2212","V1664","V1670","V1671","V2225","V1674","V853","V1675",
                          "V850","V1663","V672","V680","V312","V673","V696","V773","V1570",
                          "V1657","V2185","V2593","V2663","V2733","V2185")
J48_previous <- c("V2733", "V1647","V1647","V1564","V624","V624","V1679","V1679","V2548",
                  "V148","V2548","V227","V2743")
filtering <- c("V312","V673","V696","V773","V1570","V1657","V2185","V2593","V2663","V2733")
Huge_combination <- c(J48_previous, CorrelationAttribute, filtering, features)
Huge_Unic <- unique(Huge_combination)

##### Combining and arranging classses ####
### HR~TRIPLE
HRTR_f <- c("V59","V68","V313","V386","V624","V625","V692","V696", "V736","V773","V842","V843", "V1567","V1568","V1640","V1647","V1657","V1658","V1679","V1688","V1839","V1974","V2018","V2026","V2548",
            "V2026","V2548", "V2593","V2733",  "V2734", "V2735","V2740","V2749", "V2752","V2808")

# J48 in selected features. 10 fold cross-validation. 0.76 ac. 
HRTRJ48_f <- c("V2733","V1647","V1567","V842","V2548","V68")
HRTRPEAR_f <- c("V1679", "V1657", "V1673", "V1680", "V1669", "V1668", "V1678", "V1665", "V675")
HRTR_features = unique(c(HRTR_f,HRTRJ48_f, HRTRPEAR_f))

### HER+~TRIPLE-    1 acuarracy
HERTR_f <- c("V17","V33","V55","V485","V731","V905","V1001","V1671","V2026","V2105","V2114","V2126","V2167","V2185","V2220","V2248","V2609","V2654","V2774")
HERTRJ48_f <- c("V2185")
HERTR_features = unique(c(HERTR_f, HERTRJ48_f))

###  HER+~HR       1 acuarracy
HERHR_f <- c("V118","V475","V673","V853","V1092","V1574","V2169","V2184","V2185","V2214","V2277","V2380","V2663")
HERHRJ48_f <- c("V2185")
HERHRJ48_features = c(HERHR_f, HERHRJ48_f)

### HER+HR~TRIPLE-
# Best accuaracy with all the feature than with the feature selection. 0.88
HERHR_TRIPLE_J48 <- c("V2185","V2026","V1657","V736","V1688","V306","V736","V59","V487","V2762")


########################################################################################




################################## FEATURE SELECTION ####################################


# load: features, features_rfe, CorrelationAttribute, J48_previous. Huge_combination and
# Huge_Unic


############################### GENERATING TEST/TRAIN SETS ###############################


# create simple cross validation sets (works only with 3 classes)
crossval = 10
simple_cross_sets = crossval_sets(combine, crossval) 

# create double cross validation sets
dcrosssets = double_crossval_datasets(combine, crossval=5, loops = 5)

# make from 3 classes sets 3 sets of two classes
two_class_data = twoclasses(combine)

# make the 2 classes data equal in classes sample numbers
HER2 = make_equal_sets(two_class_data[[1]])
HR = make_equal_sets(two_class_data[[2]])
Triple = make_equal_sets(two_class_data[[3]])


############################## TRAINING MODELS ON TWO CLASSES (deleting the other) #########


# datasets
HER2HR = combine[combine$Subgroup != "TripleNeg", ]
HER2HR["Subgroup"] = droplevels(HER2HR$Subgroup)
HERTRIPLE = combine[combine$Subgroup != "HRPlus",]
HERTRIPLE["Subgroup"] = droplevels(HERTRIPLE$Subgroup)
HRTRIPLE = combine[combine$Subgroup != "HER2Plus",]
HRTRIPLE["Subgroup"] = droplevels(HRTRIPLE$Subgroup)

HER2HR_sets = crossval_sets2(HER2HR,crossval)
HERTRIPLE_sets = crossval_sets2(HERTRIPLE,crossval)
HRTRIPLE_sets = crossval_sets2(HRTRIPLE,crossval)

#write.arff(HER2HR, file = "HER2HR.arff")
#write.arff(HERTRIPLE, file = "HERTRIPLE.arff")
#write.arff(HRTRIPLE, file = "HRTRIPLE.arff")

accur = empty_accur()
pred = empty_pred()
models = list()

accur2 = empty_accur()
pred2 = empty_pred()
models2 = list()


# feature selection before training
for (j in 1:crossval) {
  
  # data and feature selection
  train_set = simple_cross_sets[[j]][[1]]   # 90%
  test_set = simple_cross_sets[[j]][[2]]    # 10%

  train_set2 = train_set[train_set$Subgroup != "HER2Plus",]
  train_set2["Subgroup"] = droplevels(train_set2$Subgroup)  # 60%
  test_set2 = test_set

  #features = filter_Var_selection(train_set,10)
  features = Huge_Unic
  features2 = HRTR_features
  
  # train the models and store on the trainset (feature selection before on whole set)
  models_trained = train_all_models_with_feature_selection(train_set,features)
  models_trained2 = train_all_models_with_feature_selection(train_set2,features2)
  
  
  # saving all the models
  models = save_models(models, models_trained)
  models2 = save_models(models2, models_trained2)
  
  # saving all the accuracies of the models
  accur = save_accur(accur, models_trained[[6]])
  accur2 = save_accur(accur2, models_trained2[[6]])
  
  # save predictions
  pred = save_pred(pred,models_trained, test_set)
  pred2 = save_pred(pred2,models_trained2, test_set2)
  
}


############################## TRAINING MODELS ON TWO CLASSES (renaming another class) #####



# train models for two class sets (HER2 and NOHER), (HR and NOHR) and (Triple and NOTriple)
features_HER2 = filter_Var_selection(HER2,number_of_features=30)
result_HER2 = train_all_models(HER2, features_HER2)

features_HR = filter_Var_selection(HR,number_of_features=30)
result_HR = train_all_models(HR, features_HR)

features_Triple = filter_Var_selection(Triple,number_of_features=30)
results_Triple = train_all_models(Triple, features_Triple)

summary(result_HER2)[3]$statistics$Accuracy
summary(result_HR)[3]$statistics$Accuracy
summary(results_Triple)[3]$statistics$Accuracy


##################### visualize the preformance of the models ##############################
models[[1]][[5]]
models[[1]][[5]]$coefnames
models[[1]][[4]]
models[[1]][[4]]$coefnames
models[[1]][[3]]
models[[1]][[3]]$coefnames
models[[1]][[2]]
models[[1]][[2]]$coefnames
models[[1]][[1]]
models[[1]][[1]]$coefnames
for (i in 1:10) {
  for (j in 1:5) {
    name = paste(i,"_",j,".pkl", sep="")
    model = models[[i]][[j]]
    save(model, file= name)
  }
}

# one model (3 classes) 
results = results_prediction(pred)
plot_results(accur, results)
results
confusionMatrix(pred$svmFit, pred$real)
confusionMatrix(pred$nnetFit, pred$real)
write.table(pred, file="pred.txt")
write.table(pred2, file="pred2.txt")
write.table(results, file="results.txt")  #### CONCLUSSION ####  --> BEST ARE SVM AND NNET WITH 85%

# Doing two times the same model (3 classes) --> [2](2 classes: HR & TR)
same_pred = pred
for (i in 1:length(pred[,1])) {
  for (j in 1:5) {
    if (pred[i,j] == 2 ) {
      same_pred[i,j] = pred2[i,j] + 1
    }
  }
}
same_results = results_prediction(same_pred)
same_results
confusionMatrix(same_pred$svmFit, pred$real)
write.table(same_pred, file="same_pred.txt")
write.table(same_results, file="same_results.txt")  #### CONCLUSSION #### --> SVM WITH 87 % 

# Doing two times the same model [1](3 classes) --> [2](3 classes) --> [2](2 classes: HR & TR)
same3_pred = pred
for (i in 1:length(pred[,1])) {
  for (j in 1:5) {
    if (pred[i,j] == 1 ) {
      same3_pred[i,j] = pred[i,j]
    }
    if (pred[i,j] == 2 ) {
      same3_pred[i,j] = pred2[i,j] + 1
    }
    if (pred[i,j] == 3 ) {
      same3_pred[i,j] = pred2[i,j] + 1
    }
  }
}
same3_results = results_prediction(same3_pred)
same3_results
confusionMatrix(same3_pred$svmFit, pred$real)
write.table(same3_pred, file="same3_pred.txt")
write.table(same_results, file="same3_results.txt") #### CONCLUSSION #### --> SVM WITH 88 %

# Doing combinations of models

# rf (3 classes) --> svm [2](2 classes: HR & TR)
rf_svm_pred = pred$rfFit
for (i in 1:length(rf_svm_pred)) {
  if (pred$rfFit[i] == 2 ) {
    rf_svm_pred[i] = pred2$svmFit[i] + 1
  }
}
confusionMatrix(rf_svm_pred, pred$real)
write.table(rf_svm_pred, file="rf_svm_pred.txt") ## 86 %

# nnet (3 classes) --> svm [2](2 classes: HR & TR)
nnet_svm_pred = pred$nnetFit
for (i in 1:length(nnet_svm_pred)) {
  if (pred$nnetFit[i] == 2 ) {
    nnet_svm_pred[i] = pred2$svmFit[i] + 1
  }
}
confusionMatrix(nnet_svm_pred, pred$real)
write.table(nnet_svm_pred, file="nnet_svm_pred.txt") ## 86 %

# svm (3 classes) --> nnet [2](2 classes: HR & TR)
svm_nnet_pred = pred$svmFit
for (i in 1:length(svm_nnet_pred)) {
  if (pred$svmFit[i] == 2 ) {
    svm_nnet_pred[i] = pred2$nnetFit[i] + 1
  }
}
confusionMatrix(svm_nnet_pred, pred$real)
write.table(svm_nnet_pred, file="svm_nnet_pred.txt")

# Doing combinations of models with more combinations

# rf (3 classes) --> svm [2](2 classes: HR & TR) --> svm [3](2 classes: HR & TR)
rf_2svm_pred = pred$rfFit
for (i in 1:length(rf_svm_pred)) {
  if (pred$rfFit[i] == 2 ) {
    rf_2svm_pred[i] = pred2$svmFit[i] + 1
  }
  if (pred$rfFit[i] == 3 ) {
    rf_2svm_pred[i] = pred2$svmFit[i] + 1
  }
}
confusionMatrix(rf_2svm_pred, pred$real)
write.table(rf_2svm_pred, file="rf_2svm_pred.txt") ## 93 %

# nnet (3 classes) --> svm [2](2 classes: HR & TR) --> svm [3](2 classes: HR & TR)
nnet_2svm_pred = pred$nnetFit
for (i in 1:length(nnet_2svm_pred)) {
  if (pred$nnetFit[i] == 2 ) {
    nnet_2svm_pred[i] = pred2$svmFit[i] + 1
  }
  if (pred$nnetFit[i] == 3 ) {
    nnet_2svm_pred[i] = pred2$svmFit[i] + 1
  }
}
confusionMatrix(nnet_2svm_pred, pred$real)
write.table(nnet_2svm_pred, file="nnet_2svm_pred.txt") ## 91 %

# most complex models !!!!!!!!!!!!!!!!!!!

# nnet (3 classes) --> rf [1](3 classes) --> [2]svm (2 classes: HR & TR) --> [3]svm (2 classes: HR & TR)
nnet_2svm_rf_pred = pred$nnetFit
for (i in 1:length(nnet_2svm_rf_pred)) {
  if (pred$nnetFit[i] == 1 ) {
    nnet_2svm_rf_pred[i] = pred$rfFit[i]
  }
  if (pred$nnetFit[i] == 2 ) {
    nnet_2svm_rf_pred[i] = pred2$svmFit[i] + 1
  }
  if (pred$nnetFit[i] == 3 ) {
    nnet_2svm_rf_pred[i] = pred2$svmFit[i] + 1
  }
}
confusionMatrix(nnet_2svm_rf_pred, pred$real)
write.table(nnet_2svm_rf_pred, file="nnet_2svm_rf_pred.txt") ## 91 %


# svm (3 classes) --> rf [1](3 classes) --> nnet [2](2 classes: HR & TR) --> nnet [3](2 classes: HR & TR)
svm_2nnet_rf_pred = pred$svmFit
for (i in 1:length(svm_nnet_pred)) {
  if (pred$svmFit[i] == 1 ) {
    svm_2nnet_rf_pred[i] = pred$rfFit[i]
  }
  if (pred$svmFit[i] == 2 ) {
    svm_2nnet_rf_pred[i] = pred2$nnetFit[i] + 1
  }
  if (pred$svmFit[i] == 3 ) {
    svm_2nnet_rf_pred[i] = pred2$nnetFit[i] + 1
  }
}
confusionMatrix(svm_2nnet_rf_pred, pred$real)
write.table(svm_2nnet_rf_pred, file="svm_2nnet_rf_pred.txt") ## 84 %

# svm (3 classes) --> rf [1](3 classes) --> svm [2](2 classes: HR & TR) --> [3]svm (2 classes: HR & TR)
svm_2svmt_rf_pred = pred$svmFit
for (i in 1:length(svm_nnet_pred)) {
  if (pred$svmFit[i] == 1 ) {
    svm_2svmt_rf_pred[i] = pred$rfFit[i]
  }
  if (pred$svmFit[i] == 2 ) {
    svm_2svmt_rf_pred[i] = pred2$svmFit[i] + 1
  }
  if (pred$svmFit[i] == 3 ) {
    svm_2svmt_rf_pred[i] = pred2$svmFit[i] + 1
  }
}
confusionMatrix(svm_2svmt_rf_pred, pred$real)
write.table(svm_2svmt_rf_pred, file="svm_2svmt_rf_pred.txt") ## 91 %


comb_pred = data.frame(rf_svm_pred, rf_2svm_pred, nnet_2svm_pred, nnet_2svm_rf_pred, svm_2svmt_rf_pred, pred$real)
colnames(comb_pred) <- c("gbmFit","svmFit","nnetFit","mrFit","rfFit","real")
comb_pred
comb_results = results_prediction(comb_pred)
colnames(comb_results) <- c("rfsvm","rf2svm","nnet2svm","nnet2svmrf","svm2svmrf","real")
comb_results
confusionMatrix(rf_2svm_pred, pred$real)
write.table(comb_results, file="comb_results.txt")

#### SO FINAL MODEL IS rf2Svm :D 

######################################## TRAIN THE FINAL MODEL ON ALL DATA #############################

# get data (all):
train_set = combine 
train_set2 = train_set[combine$Subgroup != "HER2Plus",]
train_set2["Subgroup"] = droplevels(train_set2$Subgroup)  # 60%

# set features
features = Huge_Unic
features2 = HRTR_features

# train first split
rfFit = rf_tuning(train_set,features)
rfFit_features = feature_var_imp(rfFit, 10)
rfFit = rf_tuning(train_set,rfFit_features)

# train second split
svmFit = svm_tuning(train_set2,features2)

# testing on the trainset 
pred1test = predict(rfFit,combine)
pred2test = predict(svmFit,combine)
predtest = pred1test
for (i in 1:length(predtest)) {
  if (pred1test[i] == 2 ) {
    predtest[i] = pred2test[i] + 1
  }
  if (pred1test[i] == 3 ) {
    predtest[i] = pred2test[i] + 1
  }
}
confusionMatrix(predtest, combine$Subgroup)


#################################### DO THE FINAL PREDICTIONS ##########################################


save(rfFit, file= "final_model_rfFit.pkl")
save(svmFit, file= "final_model_svmFit.pkl")
Unlabelled_data <- Load_Unlabeled_Data("Data/unlabelled_sample.txt")
pred1 = predict(rfFit,Unlabelled_data)
pred2 = predict(svmFit,Unlabelled_data)
predfinal = pred1
for (i in 1:length(predfinal)) {
  if (pred1[i] == 2 ) {
    predfinal[i] = pred2[i] + 1
  }
  if (pred1[i] == 3 ) {
    predfinal[i] = pred2[i] + 1
  }
}
predfinal


################################## Features Inspection ################################################
Features_RF_Final <- rfFit$coefnames

svmFit$bestTune
Features_HERTN_Final <- svmFit$coefnames
svmFit$modelInfo
svmFit$results
svmFit$modelInfo
svmFit$maximize
svmFit
