library(glmnet)
library(caret)

rm(list = ls())
library(MASS)
data <- Boston

ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.6, 0.3))
#ind <- sample(3, nrow(data), replace = TRUE, prob = c(0.6, 0.3,0.1))

train <- data[ind==1,]
test <- data[ind==2,]

custom <- trainControl(method = "repeatedcv",
                       
                       number = 10,
                       
                       repeats = 5,
                       
                       verboseIter = TRUE)


en <- train(medv~.,
            train,
            method='glmnet',
            tuneGrid =expand.grid(alpha=seq(0,1,length=10),
                                  lambda = seq(0.0001,0.2,length=5)),
            trControl=custom)
en


en$bestTune
coef(en$finalModel, en$bestTune$lambda)




