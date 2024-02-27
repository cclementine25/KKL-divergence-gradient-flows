import numpy as np
import scipy.stats as scs
import kernels as kl


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
    return np.where(t > 0., np.log(t), 0.)
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

#(log_ou_0([Lz[t]])[0]
#Wasserstein Gradient of KKL


def WGrad_KKL(w,x,y,k,dk,Packy,alpha,sigma):
    n = len(x)
    m = len(y)
    
    
    Kx = 1/n * k(x,x) 
    Lx,U = np.linalg.eig(Kx)
    U = U.T
    logLx = log_ou_0(Lx)
    
    U_w = 1/np.sqrt(n) * k(w[None,:],x)[0] 
    DU_w = (x - w[None,:])/sigma**2 * U_w[:,None] # ONLY FOR THE GAUSSIAN KERNEL
    #DU_w = 1/ np.sqrt(n) * dk(w[None,:],x)[0]
    
    Ky = Packy[0] 
    Kxy = k(x,y)
    
    K = np.concatenate([np.concatenate([alpha * Kx, np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy],axis = 1),np.concatenate([np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy.transpose(),(1-alpha)* Ky],axis = 1)],axis = 0)
    Lz,W = np.linalg.eig(K)
    W = W.transpose()
    logLz = log_ou_0(Lz)
    
    V_w = np.concatenate([np.sqrt(alpha) * U_w,np.sqrt((1-alpha)/m) * k(w[None,:],y)[0]]) 
    DV_w = (np.concatenate([x,y]) - w[None,:])/sigma**2 * V_w[:,None] #only for gaussian kernel
    #DV_w = np.concatenate([np.sqrt(alpha) * DU_w,np.sqrt((1-alpha)/m) * dk(w[None,:],y)[0]]) #np.concatenate([[np.sqrt(alpha/n) * dk(w,x[i]) for i in range(n)],[np.sqrt((1-alpha)/m) * dk(w,y[j]) for j in range(m)]])
    
    VW = V_w @ W.T 
    DVW = W @ DV_w
    
    Trx = 2 * U_w @ U.T @ np.diag(logLx/Lx)  @ U @ DU_w 
    Trz = 2 * VW @ np.diag(logLz/Lz) @ DVW 
    Tr3 = 2 * VW @ np.diag(np.linalg.norm(W[:,:n],axis =1)**2 / Lz) @ DVW
    Tr4 = 2 * VW @ ((logLz[:,None] - logLz[None,:]) / (Lz[:,None] - Lz[None,:] + np.identity(n+m)) * (W[:, :n] @ W.T[:n,:])) @ DVW
    
    return  Trx - Trz - Tr3 - Tr4 # 
    

      


######## Kernel density estimation ###############

#base distribution sample
x_tau = scs.multivariate_normal.rvs(np.zeros(2),np.identity(2),100)    


def h(x,y,k):
    return np.mean(np.array([k(x,x_tau[i]) * k(y,x_tau[i]) * np.exp(np.linalg.norm(x_tau[i])) /(np.sqrt(2 * np.pi)) for i in range(len(x_tau))]))
    
    

def DE(x,k,y):
    n = len(x)
    return 1/n * np.sum(np.array([h(x[i],y,k) for i in range(n)]))

def KDE(x, y, k):
    n = len(x)
    Q = np.array([DE(x,k,x[i]) for i in range(n)])
    P = np.array([DE(y,k,x[i]) for i in range(n)])
    return 1/n * np.sum(np.log(Q) * Q - np.log(P) * Q)
    
    
    
######### TRACE #######################

def K_trace(x,k):
    n = len(x)
    Kx = 1/n * np.array([[k(x[i],x[j]) for i in range(n)] for j in range(n)])
    Lambdx,_ = np.linalg.eig(Kx)
    return np.sum(Lambdx)


#######################OTHER SHAPES############


