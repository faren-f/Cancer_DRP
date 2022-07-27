import numpy as np
import pandas as pd
import numpy.ma as ma

def compute_loss(U, V, W, lambda_l,lambda_d,lambda_c, intMat, Sd, Sc):
    loss = np.sum(np.power(np.multiply(W, (intMat - np.matmul(U,V.T))),2))
    loss = loss + lambda_l * (np.sum(np.power(U,2)) + np.sum(np.power(V,2)))
    loss = loss + lambda_d * np.sum(np.power((Sd-np.matmul(U,U.T)),2)) 
    loss = loss + lambda_c * np.sum(np.power((Sc-np.matmul(V,V.T)),2))
    return(loss)


def algorithm_update(U, V, W, Y_na_zero, S, lambda_l, lambda_):
    A = (np.matmul(Y_na_zero,V)) + (2*lambda_*(np.matmul(S,U)))
    B = 2*lambda_*(np.matmul(U.T,U))
    U0 = np.zeros(np.shape(U))
    nu = np.shape(U)[1]
    C = np.matmul(V.T,V)
    m,n = np.shape(W)
    for i in range(m):
        ii = np.where(W[i,]>0)
        ii = ii[0]
        if len(ii) == 0:
            D = B + lambda_l*np.eye(nu)
        elif len(ii) == n:
            D = C + B + lambda_l*np.eye(nu)
        else:
            E = np.matmul(V[ii,].T,V[ii,])
            D = E + B + lambda_l*np.eye(nu)
        U0[i,] = np.matmul(A[i,],np.linalg.pinv(D))                       ##see mldivide
    return(U0)


def cmf(W, Y_na_zero, Sd, Sc, lambda_l,lambda_d,lambda_c,K,max_iter):
    m,n = np.shape(W)
 
    U0 = np.sqrt(1/K)*np.random.randn(m, K)
    V0 = np.sqrt(1/K)*np.random.randn(n, K)
    bestU = U0.copy()
    bestV = V0.copy()
    last_loss = compute_loss(U0, V0, W, lambda_l,lambda_d,lambda_c, Y_na_zero, Sd, Sc)
    bestloss = last_loss.copy()
    for t in (range(1,max_iter+1)):
        U = algorithm_update(U0, V0, W, Y_na_zero, Sd, lambda_l, lambda_d)
        V = algorithm_update(V0, U, W.T, Y_na_zero.T, Sc, lambda_l, lambda_c)
        curr_loss = compute_loss(U, V, W, lambda_l,lambda_d,lambda_c, Y_na_zero, Sd, Sc)
        if curr_loss < bestloss:
            bestU = U.copy()
            bestV = V.copy()
            bestloss = curr_loss.copy()
        delta_loss = (curr_loss-last_loss)/last_loss
        if abs(delta_loss) < 1e-5:
           break
        last_loss = curr_loss
        U0 = U.copy()
        V0 = V.copy()
    return (bestU,bestV)


Y = pd.read_csv("Data/sensitivity_matrix.csv" , header = 0, index_col=0, sep = ',')

Y_rownames = Y.index     #after transpose it would be colnames
Y_colnames = Y.columns   #after transpose it would be rownames

Y = np.array(Y)
Y = Y.T

Y_all_num = Y[~np.isnan(Y)]

#Y_norm = Y/np.abs(Y_all_num).max()  
Y_norm = Y/max(max(Y_all_num),abs(min(Y_all_num))) # Y is normalized to be prepared for finding U, V

Sd = pd.read_csv("Data/Fingerprints_sim.csv", header = 0, index_col=0, sep= ',')
Sd = np.array(Sd)


Sc = pd.read_csv("Data/Sample_Sim_Exp.csv", header = 0, index_col=0, sep = ',' )
Sc = np.array(Sc)

Corr_di = np.empty((np.shape(Y_norm)[0],1))
Corr_di[:] = np.NaN

MSE_di = np.empty((np.shape(Y_norm)[0],1))
MSE_di[:] = np.NaN

i1 = -2
i3 = -2
K = 45
lambda_l = 2 ** i1 
lambda_d = 0  
lambda_c = 2 ** i3 
max_iter=50
Y_na_zero = Y_norm.copy()
W = ~np.isnan(Y_na_zero)
Y_na_zero[np.isnan(Y_na_zero)] = 0

[U,V] = cmf(W,Y_na_zero,Sd,Sc,lambda_l,lambda_d,lambda_c,K,max_iter)

Y1 = Y_norm *max(max(Y_all_num),abs(min(Y_all_num)))              #returning Y_norm to Y
#Y1 = Y_na_zero *max(max(Y_all_num),abs(min(Y_all_num)))              #returning Y_norm to Y
Y_pred = np.matmul(U , V.T) * max(max(Y_all_num),abs(min(Y_all_num)))

for d in range(0,(np.shape(Y1)[0])):
    
    Y_di = Y1[d,]
    Y_pred_di = Y_pred[d,]
        
    Corr_di[d] = ma.corrcoef(ma.masked_invalid(Y_di.T), ma.masked_invalid(Y_pred_di.T))[0,1]
    MSE_di[d] = np.sqrt(np.nanmean((Y_di-Y_pred_di)**2))
    
    Corr_all_drugs = np.mean(Corr_di)
    MSE_all_drugs = np.mean(MSE_di)












