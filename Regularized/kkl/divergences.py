import sys
sys.path.append('C:/Users/cleme/DOC/Annee_2023_2024/code/Regularized/kkl')
import numpy as np
import scipy.stats as scs
import kernels as kl



def log_ou_0(t):
    return np.where(t > 0, np.log(t), 0.)



####### KKL ########

def KKL(x,y,k,Ky,alpha):
    n = len(x)
    m = len(y) 
    Kx = 1/n * k(x,x)
    
     
    Kxy = k(x,y) 
    Lx,U = np.linalg.eig(Kx)
    U = np.real(U).T
    Lx = np.real(Lx)
    Trxy = 0
    
    K = np.concatenate([np.concatenate([alpha * Kx, np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy],axis = 1),np.concatenate([np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy.transpose(),(1-alpha)* Ky],axis = 1)],axis = 0)
    Lz,W = np.linalg.eig(K)
    W = W.T
    
    A = np.concatenate([np.concatenate([1/alpha * np.identity(n), np.zeros((n,m))],axis = 1),np.concatenate([np.zeros((m,n)),np.zeros((m,m))],axis = 1)],axis = 0)
                
    Trxx = np.sum(Lx * log_ou_0(Lx))
    Trxy = np.trace(A @ W.T @ np.diag(Lz * log_ou_0(Lz)) @ W) 
    return Trxx - Trxy 


    
def WGrad_KKL(x,y,k,dk,Ky,alpha,sigma):
    n = len(x)
    m = len(y)
    
    Kx = 1/n * k(x,x)
    Lx,U = np.linalg.eig(Kx)
    Lx = np.real(Lx)
    U = np.real(U)
    U = U.T
    logLx = log_ou_0(Lx)
    
    Kxy = k(x,y)
    
    K = np.concatenate([np.concatenate([alpha * Kx, np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy],axis = 1),np.concatenate([np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy.transpose(),(1-alpha)* Ky],axis = 1)],axis = 0)
    K = K 
    Lz,W = np.linalg.eig(K)
    Lz = np.real(Lz)
    W = np.real(W)
    W = W.transpose()
    logLz = log_ou_0(Lz)
    U_w = 1/np.sqrt(n) * k(x,x)
    DU_w = (x[None,:,:] - x[:,None,:])/sigma**2 * U_w[:,:,None] # ONLY FOR THE GAUSSIAN KERNEL
    #DU_w = 1/ np.sqrt(n) * dk(w[None,:],x)[0]
    
    V_w = np.concatenate([np.sqrt(alpha) * U_w,np.sqrt((1-alpha)/m) * k(x,y)],axis = 1) 
    DV_w = (np.concatenate([x,y])[None,:,:] - x[:,None,:])/sigma**2 * V_w[:,:,None] #only for gaussian kernel
    #DV_w = np.concatenate([np.sqrt(alpha) * DU_w,np.sqrt((1-alpha)/m) * dk(w[None,:],y)[0]]) #np.concatenate([[np.sqrt(alpha/n) * dk(w,x[i]) for i in range(n)],[np.sqrt((1-alpha)/m) * dk(w,y[j]) for j in range(m)]])
    
    VW = V_w @ W.T 
    DVW = np.einsum('ik,nkj->nij', W , DV_w )
    
    Trx = 2 * np.einsum('nk,nkj->nj',(U_w @ U.T @ np.diag(logLx/Lx) @ U),DU_w )
    Trz = 2 * np.einsum('nk,nkj->nj',(VW @ (np.diag(logLz/(Lz+1e-9)) + (logLz[:,None] - logLz[None,:]) / (Lz[:,None] - Lz[None,:] + np.identity(n+m) + 1e-9 * np.ones((n+m,n+m))) * (W[:, :n] @ W.T[:n,:]) + np.diag(np.linalg.norm(W[:,:n],axis =1)**2 / Lz))), DVW) # + np.diag(np.linalg.norm(W[:,:n],axis =1)**2 / Lz) + ((logLz[:,None] - logLz[None,:]) / (Lz[:,None] - Lz[None,:] + np.identity(n+m)) * (W[:, :n] @ W.T[:n,:]))) @ DVW
    T = Trx - Trz
    return  T 

      

def energy_distance(x,y):
    XX = np.linalg.norm((x[:,None] - x[None,:]),axis = 2)
    YY = np.linalg.norm((y[:,None] - y[None,:]),axis = 2)
    XY = np.linalg.norm((x[:,None] - y[None,:]),axis = 2)
    return -np.mean(XX) - np.mean(YY) + 2 * np.mean(XY)

def KL_gauss(mux,Sigmax,muy,Sigmay):
    d = len(mux)
    detx = np.linalg.det(Sigmax)
    dety = np.linalg.det(Sigmay)
    print((mux-muy).T @ np.linalg.inv(Sigmay) @ (mux-muy))
    print(np.trace(np.linalg.inv(Sigmay) @ Sigmax))
    return 1/2 * (np.log(dety/detx) - d + (mux-muy) @ np.linalg.inv(Sigmay) @ (mux-muy) + np.trace(np.linalg.inv(Sigmay) @ Sigmax))
    
    
    
def MMD(x,y,k):
    n = len(x)
    Kxx = k(x,x)
    Kyy = k(y,y)
    Kxy = k(x,y)
    A = 1/((n-1)*n) * (np.sum(Kxx) - np.sum(np.diag(Kxx)))
    C = 1/((n-1)*n) * (np.sum(Kyy) - np.sum(np.diag(Kyy)))
    B = 1/n**2* np.sum(Kxy)
    return A - B + C


