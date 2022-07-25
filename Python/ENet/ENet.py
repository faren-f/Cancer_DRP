import os
import pandas as pd
import numpy as np
from matplotlib import pyplot
import sklearn as sk


# load dataset
root= 'Raw_Data_From_R/ENet'

path = os.path.join(root, 'GE.csv')
GE = pd.read_csv(path, header = None, sep = ',' )


path = os.path.join(root, 'FP.csv')
FP = pd.read_csv(path, header = None, sep = ',' )

path = os.path.join(root, 'sen.csv')
sen = pd.read_csv(path, header = None, sep = ',' )

path = os.path.join(root, 'cellline_drug_index.csv')
cellline_drug_index = pd.read_csv(path, header = None, sep = ',' )

row1, row2 = cellline_drug_index
Input = pd.concat([GE[row1],FP[row2]],ignore_index=True)

#Split train and test data
X_train, X_test, y_train, y_test = sk.model_selection.train_test_split(
    Input, sen, test_size=0.2, random_state=0)

# Find hyperparemeters
cv = sk.model_selection.RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
ratios = np.arange(0, 1, 0.01)
alphas = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]
model = sk.linear_model.ElasticNetCV(l1_ratio=ratios, alphas=alphas, cv=cv, n_jobs=-1)
# fit model
model.fit(X_train, y_train)
# summarize chosen configuration
print('alpha: %f' % model.alpha_)
print('l1_ratio_: %f' % model.l1_ratio_)


# define model
model = sk.linear_model.ElasticNet(alpha=1.0, l1_ratio=0.5)
# fit model
model.fit(X_train, y_train)

# make a prediction
y_pred = model.predict([X_test])
print('Predicted: %.3f' % y_pred)

score = model.score(X_test, y_test)
mse = sk.metrics.mean_squared_error(y_test, y_pred)
print("R2:{0:.3f}, MSE:{1:.2f}, RMSE:{2:.2f}"
      .format(score, mse, np.sqrt(mse)))
 

x_ax = range(len(X_test))
pyplot.plt.scatter(x_ax, y_test, s=5, color="blue", label="original")
pyplot.plt.plot(x_ax, y_pred, lw=0.8, color="red", label="predicted")
pyplot.plt.legend()
pyplot.plt.show()













