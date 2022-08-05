# glmnet package
library(glmnet)

L = seq(0.001,0.1,by = 0.02)
las.glm <- cv.glmnet(x=Xtrain, y=y,alpha=1,type.measure="mse",
                     nfolds = 5, lambda = L,standardize=FALSE)

