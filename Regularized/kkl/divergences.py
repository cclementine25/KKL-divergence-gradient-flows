import numpy as np
import scipy.stats as scs
import kkl.kernels as kl


######## MMD ##########

def MMD(x,y,k):
    n = len(x)
    Kxx = np.array([[k(x[i],x[j]) for i in range(n)] for j in range(n)])
    Kyy = np.array([[k(y[i],y[j]) for i in range(n)] for j in range(n)])
    Kxy = np.array([[k(x[i],y[j]) for i in range(n)] for j in range(n)])
    A = 1/((n-1)*n) * (np.sum(Kxx) - np.sum(np.diag(Kxx)))
    C = 1/((n-1)*n) * (np.sum(Kyy) - np.sum(np.diag(Kyy)))
    B = 1/n**2* np.sum(Kxy)
    return A - B + C


#gradient in x of MMD 
def grad_MMD(x,y,k,dk):
    d = len(x[0])
    n = len(x)
    m = len(y)
    dKx = np.array([[dk(x[i],x[j]) for j in range(n)] for i in range(n)])
    dKx[:,:,0] = dKx[:,:,0] - np.diag(np.diag(dKx[:,:,0]))
    dKx[:,:,1] = dKx[:,:,1] - np.diag(np.diag(dKx[:,:,1]))
    dKy = np.array([[dk(x[i],y[j]) for j in range(m)] for i in range(n)])
    R = np.zeros((n,d))
    R[:,0] = 2/(n * (n-1)) * dKx[:,:,0] @ np.ones(n) - 2/m**2 * dKy[:,:,0] @ np.ones(m)
    R[:,1] = 2/(n * (n-1)) * dKx[:,:,1] @ np.ones(n) - 2/m**2 * dKy[:,:,1] @ np.ones(m)
    return R



def log_ou_0(t):
    return np.where(t > 0, np.log(t), 0.)
    # t_log = np.zeros(len(t))
    # for i in range(len(t)):
    #     if t[i] > 0:
    #         t_log[i] = np.log(t[i])
    # return t_log




####### KKL ########


def KKL(x,y,k,Packy,alpha):
    n = len(x)
    m = len(y) 
    Kx = 1/n * k(x,x)
    
    
    Ky = Packy[0] 
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
    Trxy = np.trace(A @ W.T @ np.diag(Lz * log_ou_0(Lz)) @ W) #K @ Klog)
    return Trxx - Trxy 


    
def WGrad_KKL(x,y,k,dk,Packy,alpha,sigma):
    n = len(x)
    m = len(y)
    
    Kx = 1/n * k(x,x)
    Lx,U = np.linalg.eig(Kx)
    Lx = np.real(Lx)
    U = np.real(U)
    U = U.T
    logLx = log_ou_0(Lx)
    
    Ky = Packy[0] 
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
    return  T #- Tr3 - Tr4 

      

def energy_distance(x,y):
    XX = np.linalg.norm((x[:,None] - x[None,:]),axis = 2)
    YY = np.linalg.norm((y[:,None] - y[None,:]),axis = 2)
    XY = np.linalg.norm((x[:,None] - y[None,:]),axis = 2)
    return -np.mean(XX) - np.mean(YY) + 2 * np.mean(XY)
    
    
    
    


