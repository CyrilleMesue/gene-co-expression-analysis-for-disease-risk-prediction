#Load required packages
library(ggplot2)
library(caret)
library(e1071) 
library(ROCR) 
library(Boruta)
library(glmnet)
library(ggpubr)
library(MLeval)
library(randomForest)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load differential Co-expressed genes (DCEGs) including asthma type ( asthmatic vs controls)
# Example 
dat<-read.csv("data/GSE152004_train.csv", row.names = 1)
# Logistic regression with Lasso penalty
y<-dat$Type
x<-model.matrix(y~., dat[,-1])
cv.fit <- cv.glmnet(x,y, family="binomial", alpha= 1)
plot(cv.fit)
fit <- glmnet(x,y, family = "binomial", alpha = 1)
plot(fit)
cv.fit$lambda.min

Coefficients <- coef(fit, s = fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
coef(cv.fit, s = "lambda.min")
#************************************
# Recursive Feature Elimination (RFE) for Feature Selection
# Define the control using a random forest selection function
control <- rfeControl(functions = rfFuncs, # random forest
                      number = 10) # number of folds

# run the RFE algorithm
results <- rfe(dat[,2:86], dat[,1], sizes=c(1:86), rfeControl=control) # column 1 target var
# summarize the results

# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))


varimp_RFE <- data.frame(feature = row.names(varImp(results))[1:854],importance = varImp(results)[1:854, 1])

write.csv(varimp_RFE, "varimp_data_rfe2.csv")

#Random Forest 
Control <- trainControl(method = "cv",number = 5)
rf <- train(Type~., dat,trCtrl= Control,method = "rf")
#VIP
varImp(rf)
importance <- varImp(rf, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)
df<-importance$importance
plot(importance, top=30)
write.csv(df,"predictor._rf.csv")

# Run Boruta Algorithm
boruta <- Boruta(Type~., data = dat, doTrace = 2)
print(boruta)
plot(boruta)
xx<-attStats(boruta)
# Select final confirmed final vars
finalvars = getSelectedAttributes(boruta, withTentative = F)
finalvars
write.csv(xx, "finalvars_boruta_30.csv")

# The bar plot to show rank of selected genes using machine learning methods  

ggbarplot(dt, x = "ID", y = "Rank",
          fill = "Method",           
          color = "white",            
          palette = c("brown"), 
          sort.val = "desc",          
          sort.by.groups = FALSE,     
          x.text.angle = 90,          
          ylab = "Rank",
          legend.title = " ",
          rotate = TRUE,
          font.caption = "bold",
          font.x = c(10, "bold", "black"),
          font.y = c(10, "bold", "black"),
          font.tickslab = c(5,"bold", "black"),
          
          ggtheme = theme_minimal()
)


# machine learning analysis based RF and SVM classifier comparison 
# load selected Genes for asthma classification model using ML methods.

metric <- "Accuracy"
control <- trainControl(method="cv", 
                        number=5, 
                        repeats = 5, 
                        savePredictions = "final", 
                        search = "grid", 
                        classProbs = TRUE)

datM<-read.csv("GSE67472.modulegene854.csv", row.names = 1)
datD<-read.csv("GSE67472.DEGgene3665.csv", row.names = 1)
datB<-read.csv("GSE67472_boruta30genes.csv", row.names = 1)
datL<-read.csv("GSE67472_Lasso30genes.csv", row.names = 1)
datRF<-read.csv("GSE67472_RFOrder30genes.csv", row.names = 1)
datRFE<-read.csv("GSE67472_RFE30genes.csv", row.names = 1)
datM$Type<-as.factor(datM$Type)
datD$Type<-as.factor(datD$Type)
datB$Type<-as.factor(datB$Type)
datL$Type<-as.factor(datL$Type)
datRF$Type<-as.factor(datRF$Type)
datRFE$Type<-as.factor(datRFE$Type)


datM1<-datM[,1:854]
datD1<-datD[,1:3565]
datB1<-datB[,1:31]
datL1<-datL[,1:31]
datRF1<-datRF[,1:31]
datRFE1<-datRFE[,1:31]

set.seed(7)
fit.rfM1 <- train(Type ~ .,data=datM1, method="rf", metric=metric, trControl=control, verbose = FALSE)
fit.rfD1 <- train(Type ~ .,data=datD1, method="rf", metric=metric, trControl=control, verbose = FALSE)
fit.rfB1 <- train(Type ~ .,data=datB1, method="rf", metric=metric, trControl=control, verbose = FALSE)
fit.rfL1 <- train(Type ~ .,data=datL1, method="rf", metric=metric, trControl=control, verbose = FALSE)
fit.rfRF1 <- train(Type ~ .,data=datRF1, method="rf", metric=metric, trControl=control, verbose = FALSE)
fit.rfRFE1 <- train(Type ~ .,data=datRFE1, method="rf", metric=metric, trControl=control, verbose = FALSE)


set.seed(7)
fit.rfM1 <- train(Type ~ .,data=datM1, method="svmRadial", metric=metric, trControl=control, verbose = FALSE)
fit.rfD1 <- train(Type ~ .,data=datD1, method="svmRadial", metric=metric, trControl=control, verbose = FALSE)
fit.rfB1 <- train(Type ~ .,data=datB1, method="svmRadial", metric=metric, trControl=control, verbose = FALSE)
fit.rfL1 <- train(Type ~ .,data=datL1, method="svmRadial", metric=metric, trControl=control, verbose = FALSE)
fit.rfRF1 <- train(Type ~ .,data=datRF1, method="svmRadial", metric=metric, trControl=control, verbose = FALSE)
fit.rfRFE1 <- train(Type ~ .,data=datRFE1, method="svmRadial", metric=metric, trControl=control, verbose = FALSE)


models <- list( Allmodulegenes = fit.rfM1,
                AllDEGs = fit.rfD1,
                Borutagenes = fit.rfB1,
                Lassogenes = fit.rfL1,
                RFgenenes = fit.rfRF1,
                RFEgenes = fit.rfRFE1 )
models <- models[order(names(models))]

evalm(
  models,
  gnames = names(models),
  rlinethick = 0.8, fsize = 8, plots = "r"
)


