
###Shrunken centroid based single sample classifier based on signature identified in the script Network_MIBC

#Load R objects
load("category_PURE.RData") ##classification of immune categories based on median stratification
load("pure01.metadata.RData")###metadata for PURE01
load("PURE0_immune.RData")##PURE01 expression object with signature genes

load("clinicaldata_ABACUS_compact") ## ABACUS1 clinical data
load("expression_ABACUS1_immune1") ## ABACUS1 expression object with signature genes


###Look at pamr


data<-list(x=as.matrix(PURE1_immune),y=factor(category))

set.seed(123)
pam_train <- pamr.train(data)
myresults2 <- pamr.cv(pam_train, data,nfold=10)
best_index<-max(which(myresults2$error==min(myresults2$error)))
train_threshold<-myresults2$threshold[best_index]
pamr.confusion(myresults2,train_threshold)
yhat_train<-pam_train$yhat[,best_index]
names(yhat_train)<-colnames(PURE1_immune)

association_pam<-merge(PURE1_clinical,yhat_train,by="row.names")
res.cox<- coxph(Surv(censored_time, censored_status) ~y, data = association_pam)
summary(res.cox)

PURE1_merge_pam<-merge(pure01.metadata.compact,yhat_train,by="row.names")
colnames(PURE1_merge_pam)[ncol(PURE1_merge_pam)]<-"category"

PURE1_merge_high_pam<-PURE1_merge_pam[PURE1_merge_pam$category=="high",]
table(PURE1_merge_high_pam[,3]) 
                           
PURE1_merge_low_pam<-PURE1_merge_pam[PURE1_merge_pam$category=="low",]
table(PURE1_merge_low_pam[,3]) 



###Look at the single sample classifier in ABACUS

yhat_test <- pamr.predict(pam_train,expression_ABACUS1_immune1,train_threshold)
names(yhat_test)<-colnames(expression_ABACUS1_immune1)

ABACUS_category_pam<-merge(clinicaldata_ABACUS_compact,yhat_test,by="row.names")
rownames(ABACUS_category_pam)<-ABACUS_category_pam$Row.names
colnames(ABACUS_category_pam)[ncol(ABACUS_category_pam)]<-"category"

ABACUS_category_high_pam<-ABACUS_category_pam[ABACUS_category_pam$category=="high",]
table(ABACUS_category_high_pam[,4]) 
ABACUS_category_low_pam<-ABACUS_category_pam[ABACUS_category_pam$category=="low",]

table(ABACUS_category_low_pam[,4]) 	


 