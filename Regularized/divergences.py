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
    t_log = np.zeros(len(t))
    for i in range(len(t)):
        if t[i] > 0:
            t_log[i] = np.log(t[i])
    return t_log




####### KKL ########


def KKL(x,y,k,Packy,alpha):
    n = len(x)
    m = len(y) 
    Kx = 1/n * np.array([[k(x[i],x[j]) for i in range(n)] for j in range(n)])
    Ky = Packy[0] 
    Kxy = np.array([[k(x[i],y[j]) for j in range(m)] for i in range(n)])
    Lx,U = np.linalg.eig(Kx)
    U = np.real(U).transpose()
    Lx = np.real(Lx)
    Trxy = 0
    
    K = np.concatenate([np.concatenate([alpha * Kx, np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy],axis = 1),np.concatenate([np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy.transpose(),(1-alpha)* Ky],axis = 1)],axis = 0)
    Lw,W = np.linalg.eig(K)
    W = W.transpose()
    
    A = np.concatenate([np.concatenate([1/alpha * np.identity(n), np.zeros((n,m))],axis = 1),np.concatenate([np.zeros((m,n)),np.zeros((m,m))],axis = 1)],axis = 0)
    Klog = W.transpose() @ np.diag(log_ou_0(Lw)) @ W
                
    Trxx = np.sum(Lx * log_ou_0(Lx))
    Trxy = np.trace(A @ K @ Klog)
    return Trxx - Trxy # 

#(log_ou_0([Lz[t]])[0]
#Wasserstein Gradient of KKL


def WGrad_KKL(w,x,y,k,dk,Packy,alpha):
    n = len(x)
    m = len(y)
    
    Kx = 1/n * np.array([[k(x[i],x[j]) for i in range(n)] for j in range(n)])
    Lx,U = np.linalg.eig(Kx)
    U = U.transpose()
    Kx_log = U.transpose() @ np.diag(log_ou_0(Lx)) @ U
    U_w = np.array([1/np.sqrt(n) * k(w,x[i]) for i in range(n)])
    DU_w = np.array([1/ np.sqrt(n) * dk(w,x[i]) for i in range(n)])
    
    Ky = Packy[0] 
    Kxy = np.array([[k(x[i],y[j]) for j in range(m)] for i in range(n)])
    K = np.concatenate([np.concatenate([alpha * Kx, np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy],axis = 1),np.concatenate([np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy.transpose(),(1-alpha)* Ky],axis = 1)],axis = 0)
    Lz,W = np.linalg.eig(K)
    W = W.transpose()
    K_log = W.transpose() @ np.diag(log_ou_0(Lz)) @ W
    V_w = np.concatenate([[np.sqrt(alpha/n) * k(w,x[i]) for i in range(n)], [np.sqrt((1-alpha)/m) *k(w,y[j]) for j in range(m)]])
    DV_w = np.concatenate([[np.sqrt(alpha/n) * dk(w,x[i]) for i in range(n)],[np.sqrt((1-alpha)/m) * dk(w,y[j]) for j in range(m)]])

    Trx = (U_w @ U.transpose() @ np.diag(log_ou_0(Lx)/Lx)  @ U @ DU_w) + (DU_w.transpose() @ U.transpose() @ np.diag(log_ou_0(Lx)/Lx)  @ U @ U_w)
    Trz = (V_w @ W.transpose() @ np.diag(log_ou_0(Lz)/Lz) @ W @ DV_w) + (DV_w.transpose() @ W.transpose() @ np.diag(log_ou_0(Lz)/Lz) @ W @ V_w)
    Tr3 = 0
    for j in range(n+m):
        Tr3 = Tr3 + np.linalg.norm(W[j])**2/Lz[j] *( V_w @ W[j].reshape(n+m,1) @ W[j].reshape(1,n+m) @ DV_w + DV_w.transpose() @ W[j].reshape(n+m,1) @ W[j].reshape(1,n+m) @ V_w)
        #Tr3 = Tr3 + np.linalg.norm(W[j])**2/Lz[j] * np.array([np.trace(W[j].reshape(n+m,1) @ W[j].reshape(1,n+m) @ (V_w.reshape(n+m,1) @ DV_w[:,0].reshape(1,n+m) + DV_w[:,0].reshape(n+m,1) @ V_w.reshape(1,n+m))),np.trace(W[j].reshape(n+m,1) @ W[j].reshape(1,n+m) @ (V_w.reshape(n+m,1) @ DV_w[:,1].reshape(1,n+m) + DV_w[:,1].reshape(n+m,1) @ V_w.reshape(1,n+m)))])
        for k in range(n+m):
            if k != j:
                Tr3 = Tr3 + (log_ou_0([Lz[j]])[0] - log_ou_0([Lz[k]])[0])/(Lz[j] - Lz[k]) * W[j,:n] @ W[k,:n] * ( V_w @ W[j].reshape(n+m,1) @ W[k].reshape(1,n+m)  @ DV_w + DV_w.transpose() @ W[j].reshape(n+m,1) @ W[k].reshape(1,n+m)  @ V_w)
    print(Trx)
    print(Trz)
    print(Tr3)
    print()
    return  Trx - Trz - Tr3 # 
    
#U[s] @ Kxz @ V[t]        @ np.diag(np.log(Lz)/Lz)                           
    

# def WGrad_KKL(w,x,y,k,dk,Packy,alpha):
#     n = len(x)
#     m = len(y)
#     Kx = 1/n * np.array([[k(x[i],x[j]) for i in range(n)] for j in range(n)])
#     Ky = Packy[0] #1/m * np.array([[k(y[i],y[j]) for i in range(m)] for j in range(m)])
#     Kxy = np.array([[k(x[i],y[j]) for j in range(m)] for i in range(n)])
#     Kz = np.concatenate([np.concatenate([alpha * Kx, alpha/n *Kxy],axis = 1),np.concatenate([(1-alpha)/m *Kxy.transpose(),(1-alpha)* Ky],axis = 1)],axis = 0)
#     Kxz = np.concatenate([n*Kx, Kxy],axis = 1)
#     Lx,U = np.linalg.eig(Kx)
#     U = U.transpose()
#     Lx = np.real(Lx)
#     Lz,V = np.linalg.eig(Kz) 
#     V = V.transpose()
#     Lz = np.real(Lz)
#     Kwx = np.array([k(w,x[i]) for i in range(n)]).transpose()
#     Kwy = np.array([k(w,y[j]) for j in range(m)]).transpose()
#     Kwz = np.concatenate([Kwx,Kwy])
#     DKx = np.array([dk(w,x[i]) for i in range(n)]).transpose()
#     DKy = np.array([dk(w,y[j]) for j in range(m)]).transpose()
#     DKz = np.concatenate([DKx.transpose(),DKy.transpose()]).transpose()
#     Trwx = 0
#     Trwz = 0 
#     Tr3 = 0
#     for s in range(n):
#         Trwx = Trwx + log_ou_0([Lx[s]])[0] / Lx[s] * 2 * (U[s] @ Kwx)* (DKx @ U[s]) 
#         #print(U[s] @ (n * Kx) @ U[s])
#     for t in range(m):
#         Trwz = Trwz + log_ou_0([Lz[t]])[0] / (Lz[t] *(n/alpha * np.linalg.norm(V[t][:n])**2 + m/(1-alpha) * np.linalg.norm(V[t][n:])**2)) * 2 * (V[t] @ Kwz)* (DKz @ V[t]) 
        
#     for s in range(n):
#         for t in range(n+m):
#             if Lz[t] > 1e-15 and Lx[s] > 0:
#                 Tr3 = Tr3 + 1/(Lz[t]**2 * (n/alpha * np.linalg.norm(V[t][:n])**2 + m/(1-alpha) * np.linalg.norm(V[t][n:])**2)) *  (n/alpha *Lz[t]* U[s] @ V[t][:n]) * ((U[s] @ Kwx) * (DKz @ V[t]) + (DKx @ U[s]) * (V[t] @ Kwz))
#                 # if np.linalg.norm(1/(Lz[t]**2 * (n/alpha * np.linalg.norm(V[t][:n])**2 + m/(1-alpha) * np.linalg.norm(V[t][n:])**2)) *  (U[s] @ Kxz @ V[t]) * ((U[s] @ Kwx) * (DKz @ V[t]) + (DKx @ U[s]) * (V[t] @ Kwz))) > 50:
#                 #     print(U[s] @ Kxz @ V[t])
#                 #     print((n/alpha * Lz[t] * U[s] @ V[t][:n]))
#                 #     print(np.linalg.norm(((U[s] @ Kwx) * (DKz @ V[t]) + (DKx @ U[s]) * (V[t] @ Kwz))))
#                 #     print(Lz[t])
#                 #     print()
#     return 1/n * Trwx  -  Trwz #- alpha/n * Tr3
      


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


