
import numpy as np
import pandas as pd

def compute_loss(U, V, W, lambda_l,lambda_d,lambda_c, intMat, Fingerprints_sim, Sample_Sim_Exp):
    loss = np.sum(np.power(np.multiply(W, (intMat - np.matmul(U,V.T))),2))
    loss = loss + lambda_l * (np.sum(np.power(U,2)) + np.sum(np.power(V,2)))
    loss = loss + lambda_d * np.sum(np.power((Fingerprints_sim-np.matmul(U,U.T)),2)) 
    loss = loss + lambda_c * np.sum(np.power((Sample_Sim_Exp-np.matmul(V,V.T)),2))
    return(loss)


def algorithm_update(U, V, W, WR, S, lambda_l, lambda_d_c):
    X = (np.matmul(WR,V)) + (2*lambda_d_c*(np.matmul(S,U)))
    Y = 2*lambda_d_c*(np.matmul(U.T,U))
    U0 = np.zeros(np.shape(U))
    nu = np.shape(U)[1]
    D = np.matmul(V.T,V)
    m,n = np.shape(W)
    for i in range(m):
        ii = np.where(W[i,]>0)
        ii = ii[0]
        if len(ii) == 0:
            B = Y + lambda_l*np.eye(nu)
        elif len(ii) == n:
            B = D + Y + lambda_l*np.eye(nu)
        else:
            A = np.matmul(V[ii,].T,V[ii,])
            B = A + Y + lambda_l*np.eye(nu)
    U0[i,] = np.matmul(X[i,],np.linalg.pinv(B))                       ##see mldivide
    return(U0)




def cmf(W, sen_na_zero, Fingerprints_sim, Sample_Sim_Exp, lambda_l,lambda_d,lambda_c,K,max_iter):
    m,n = np.shape(W)
 
    U0 = np.sqrt(1/K)*np.random.randn(m, K)
    V0 = np.sqrt(1/K)*np.random.randn(n, K)
    bestU = U0.copy()
    bestV = V0.copy()
    last_loss = compute_loss(U0, V0, W, lambda_l,lambda_d,lambda_c, sen_na_zero, Fingerprints_sim, Sample_Sim_Exp)
    bestloss = last_loss.copy()
    WR = np.multiply(W ,sen_na_zero)
    for t in (range(1,max_iter+1)):
        U = algorithm_update(U0, V0, W, WR, Fingerprints_sim, lambda_l, lambda_d)
        V = algorithm_update(V0, U, W.T, WR.T, Sample_Sim_Exp, lambda_l, lambda_c)
        curr_loss = compute_loss(U, V, W, lambda_l,lambda_d,lambda_c, sen_na_zero, Fingerprints_sim, Sample_Sim_Exp)
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


sen = pd.read_csv("Data/sensitivity_matrix.csv" , header = 0, index_col=0, sep = ',')

sen_rownames = sen.index     #after transpose it would be colnames
sen_colnames = sen.columns   #after transpose it would be rownames

sen = np.array(sen)
sen = sen.T

sen_all_num = sen[~np.isnan(sen)]

#sen_norm = sen/np.abs(sen_all_num).max()  
sen_norm = sen/max(max(sen_all_num),abs(min(sen_all_num))) # sen is normalized to be prepared for finding U, V

Fingerprints_sim = pd.read_csv("Data/Fingerprints_sim.csv", header = 0, index_col=0, sep= ',')
Fingerprints_sim = np.array(Fingerprints_sim)


Sample_Sim_Exp = pd.read_csv("Data/Sample_Sim_Exp.csv", header = 0, index_col=0, sep = ',' )
Sample_Sim_Exp = np.array(Sample_Sim_Exp)

drugwisecorr = np.empty((np.shape(sen_norm)[0],1))
drugwisecorr[:] = np.NaN

drugwise_qt = np.empty((np.shape(sen_norm)[0],1))
drugwise_qt[:] = np.NAN

drugwiseerr = np.empty((np.shape(sen_norm)[0],1))
drugwiseerr[:] = np.NaN

drugwiseerr_qt = np.empty((np.shape(sen_norm)[0],1))
drugwiseerr_qt[:] = np.NAN

drugwiserepn = np.empty((np.shape(sen_norm)[0],1))
drugwiserepn[:] = np.NaN

i1 = -2
i3 = -2
K = 45
lambda_l = 2 ** i1 
lambda_d = 0  
lambda_c = 2 ** i3 
max_iter=50
sen_na_zero = sen_norm.copy()
W = ~np.isnan(sen_na_zero)
sen_na_zero[np.isnan(sen_na_zero)] = 0
[U,V] = cmf(W,sen_na_zero,Fingerprints_sim,Sample_Sim_Exp,lambda_l,lambda_d,lambda_c,K,max_iter)

#sen1 = sen_norm *max(max(sen_all_num),abs(min(sen_all_num)))              #returning sen_norm to sen
sen_na_zero = sen_na_zero *max(max(sen_all_num),abs(min(sen_all_num)))              #returning sen_norm to sen
sen_pred = np.matmul(U , V.T) * max(max(sen_all_num),abs(min(sen_all_num)))

for d in range(0,(np.shape(sen_na_zero)[0])):
    
    curtemp1 = sen_na_zero[d,]
    y1 = np.percentile(curtemp1,75)
    xia1 = np.where(curtemp1 >= y1)
    xia1 = np.array(xia1[0])
    y2 = np.percentile(curtemp1,25)
    xia2 = np.where(curtemp1 <= y2)
    xia2 = np.array(xia2[0])

    xia = np.concatenate((xia1,xia2))
    
    drugwise_qt[d] = np.corrcoef(curtemp1[xia].T,sen_pred[d,xia].T)[0,1]
    drugwiseerr_qt[d] = np.sqrt(sum(np.power((curtemp1[xia]-sen_pred[d,xia]),2))/sum(~np.isnan(curtemp1[xia])))
    curtemp2 = sen_pred[d,]
    
    
    curtemp2[np.isnan(curtemp1)] = []
    curtemp1[np.isnan(curtemp1)] = []
    
    
    drugwiserepn[d] = len(curtemp1) 
    drugwisecorr[d] = np.corrcoef(curtemp1.T,curtemp2.T)[0,1]
    drugwiseerr[d] = np.sqrt(sum(np.power((curtemp1-curtemp2),2))/sum(~np.isnan(curtemp1)))

#save('drugwise_predict1.mat')












