


000000000000000000000000000
# def KKL(x,y,k):
#     n = len(x)
#     m = len(y)
#     Kx = 1/n * np.array([[k(x[i],x[j]) for i in range(n)] for j in range(n)])
#     Ky = 1/m * np.array([[k(y[i],y[j]) for i in range(m)] for j in range(m)])
#     regx = 1e-9*np.eye(n)
#     regy = 1e-9*np.eye(m)
#     Kx = Kx +regx
#     Ky = Ky+regy
#     Lx,U = np.linalg.eig(Kx)
#     U = np.real(U).transpose()
#     Lx = np.real(Lx)
#     Ly,V = np.linalg.eig(Ky)
#     V = np.real(V).transpose()
#     Ly = np.real(Ly)
#     Trxy = 0
#     Kxy = np.array([[k(x[i],y[j]) for j in range(m)] for i in range(n)])
#     Trxx = np.sum(Lx * log_ou_0(Lx))
#     for s in range(n):
#         for t in range(m):
#             Trxy = Trxy + (Lx[s] * log_ou_0([Ly[t]])[0] * (U[s] @ Kxy @ V[t])**2) / ((U[s] @ (n * Kx) @ U[s]) * (V[t] @ (m*Ky) @ V[t]))
#             # if n*Lx[s] != (U[s] @ (n * Kx) @ U[s]):
#             #     print((U[s] @ (n * Kx) @ U[s]) / (n*Lx[s]))
#                 #print(U[s])#np.abs((U[s] @ Kxy @ V[t])**2 / ((U[s] @ (n * Kx) @ U[s]) * (V[t] @ (m * Ky) @ V[t]))) > 0: #(n*Lx[s]*m*Ly[t])) > 1:#
#                 # print(" n x lambda_x = " + str(n*Lx[s]))
#                 # print("m x lambda_y = " + str(m*Ly[t]))
#                 # print(" produit scalire = " + str((U[s] @ Kxy @ V[t])**2))
#                 # print("norme^2 de f = " + str((U[s] @ (n * Kx) @ U[s])))
#                 # print("norme^2 de g = " + str((V[t] @ (m * Ky) @ V[t])))
#             #print((U[s] @ Kxy @ V[t])**2 / ((U[s] @ (n * Kx) @ U[s]) * (V[t] @ (m*Ky) @ V[t])))
            
#     Trxx = np.sum(Lx * log_ou_0(Lx))
#     #print(UU)
#     return Trxx - Trxy


# def WGrad_KKL(w,x,y,k,dk):
#     n = len(x)
#     m = len(y)
#     Kx = 1/n * np.array([[k(x[i],x[j]) for i in range(n)] for j in range(n)])
#     Ky = 1/m * np.array([[k(y[i],y[j]) for i in range(m)] for j in range(m)])
#     Lx,U = np.linalg.eig(Kx)
#     U = U.transpose()
#     Lx = np.real(Lx)
#     Ly,V = np.linalg.eig(Ky)
#     V = V.transpose()
#     Ly = np.real(Ly)
#     Kwx = np.array([k(w,x[i]) for i in range(n)]).transpose()
#     Kwy = np.array([k(w,y[j]) for j in range(m)]).transpose()
#     DKx = np.array([dk(w,x[i]) for i in range(n)]).transpose()
#     DKy = np.array([dk(w,y[j]) for j in range(m)]).transpose()
#     Trwx = 0
#     Trwy = 0 
#     for s in range(n):
#         Trwx = Trwx + log_ou_0([Lx[s]])[0] * 2 * (U[s] @ Kwx)* (DKx @ U[s]) / (U[s] @ (n * Kx) @ U[s])
#         #print(U[s] @ (n * Kx) @ U[s])
#     for t in range(m):
#         Trwy = Trwy + log_ou_0([Ly[t]])[0] * 2 * (V[t] @ Kwy)* (DKy @ V[t]) / (V[t] @ (n * Ky) @ V[t])
#     return Trwx - Trwy




000000000000000000
# ###### GAUSSIAN ######
# #initial distribution
# mux = np.array([2,5])
# Lx = np.array([[1/2,1/3],[1/4,-2]])
# Sigmax = Lx @ Lx.transpose()
# #x0 = scs.multivariate_normal.rvs(mux,Sigmax,n)

# #Simulation of (Y_i)_i<n ~ p -> objective distribution
# muy = np.array([0,0])
# Ly = np.array([[1/5, -1],[1/2,1/2]])
# Sigmay = Ly @ Ly.transpose()
# #y = scs.multivariate_normal.rvs(muy,Sigmay,m)



# # ########### Mixture de gaussienne ###########
# MU = np.array([[-2,-1],[5,0]])
# Z = np.random.choice([0,1],m,p=[1/2,1/2])
# #y = np.array([scs.multivariate_normal.rvs(MU[Z[i]],0.5 * np.identity(2)) for i in range(m)])


# ########## RINGS ##########
# y,x0 = rg.generate_rings(n,m, 0.5,1, 0.5,1)

# ########################

# x0 = scs.multivariate_normal.rvs(muy,0.2 * np.diag([1,2]),n)


0000000000000000000
# def GD_k_gauss(J,dJ,x0,h,eps,n_it_max):
#     x = x0  
#     grad = dJ(x)
#     X = [x0]
#     i = 0 
#     liste_J = []
#     Grad = []
#     while np.linalg.norm(grad) > eps and i < n_it_max:
#         liste_J.append(J(x))
#         grad = dJ(x)
#         Grad.append(np.linalg.norm(grad))
#         x = x - h * grad
#         X.append(x)
#         i = i + 1
#     return np.array(X), liste_J,Grad


0000000000000000000000000


# def KKL(x,y,k,Packy,alpha):
#     n = len(x)
#     m = len(y) 
#     Kx = 1/n * np.array([[k(x[i],x[j]) for i in range(n)] for j in range(n)])
#     Ky = Packy[0] 
#     Kxy = np.array([[k(x[i],y[j]) for j in range(m)] for i in range(n)])
#     #Kz = np.concatenate([np.concatenate([alpha * Kx, alpha/n *Kxy],axis = 1),np.concatenate([(1-alpha)/m *Kxy.transpose(),(1-alpha)* Ky],axis = 1)],axis = 0)
#     # regx = 1e-9*np.eye(n)
#     # regy = 1e-9*np.eye(m)
#     # Kx = Kx +regx
#     # Ky = Ky+regy
#     Lx,U = np.linalg.eig(Kx)
#     U = np.real(U).transpose()
#     Lx = np.real(Lx)
#     #Lz,V = np.linalg.eig(Kz)
#     #V = np.real(V).transpose()
#     #Lz = np.real(Lz)
#     #Kxz = np.concatenate([n*Kx, Kxy],axis = 1)
#     Trxy = 0
    
#     K = np.concatenate([np.concatenate([alpha * Kx, np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy],axis = 1),np.concatenate([np.sqrt(alpha*(1-alpha)/(n*m)) *Kxy.transpose(),(1-alpha)* Ky],axis = 1)],axis = 0)
#     Lw,W = np.linalg.eig(K)
#     W = W.transpose()
    
#     A = np.concatenate([np.concatenate([1/alpha * np.identity(n), np.zeros((n,m))],axis = 1),np.concatenate([np.zeros((m,n)),np.zeros((m,m))],axis = 1)],axis = 0)
#     Klog = W.transpose() @ np.diag(log_ou_0(Lw)) @ W
#     # for s in range(n):
#     #     for t in range(n+m):
#     #         if Lz[t] > 1e-15:
#     #             Tr = Tr + (1 /(n/alpha * np.linalg.norm(V[t][:n])**2 + m/(1-alpha) * np.linalg.norm(V[t][n:])**2)) * (np.log(Lz[t]) / Lz[t]) * (U[s] @ Kxz @ V[t])**2 
#     #             D = np.array([np.sqrt(alpha/n) for i in range(n)] + [np.sqrt((1-alpha)/m) for j in range(m)]) * W[t]
#     #             Trxy = Trxy + (np.log(Lw[t]) / Lw[t]) * (U[s] @ Kxz @ D)**2
                
#     Trxx = np.sum(Lx * log_ou_0(Lx))
#     Trxy = np.trace(A @ K @ Klog)
#     # print(Trxx)
#     # print(1/n*Trxy)
#     # print(np.trace(A @ K @ Klog))
#     return Trxx - Trxy # 
# #(log_ou_0([Lz[t]])[0]
# #Wasserstein Gradient of KKL


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



000000000000000000000000000000000000
# Tr5 = 0
# Tr6 = 0
# Tr55 = 0
# for j in range(n+m):
#     Tr5 = Tr5 + np.linalg.norm(W[j][:n])**2/Lz[j] *( (V_w @ W[j]) * (W[j].reshape(1,n+m) @ DV_w)[0] + (DV_w.T @ W[j]) * (W[j] @ V_w))
#     Tr55 = Tr55 + 2*np.linalg.norm(W[j][:n])**2/Lz[j] *( (V_w @ W[j]) * (W[j].reshape(1,n+m) @ DV_w)[0])
#     #Tr5 = Tr5 + np.linalg.norm(W[j])**2/Lz[j] * np.array([np.trace(W[j].reshape(n+m,1) @ W[j].reshape(1,n+m) @ (V_w.reshape(n+m,1) @ DV_w[:,0].reshape(1,n+m) + DV_w[:,0].reshape(n+m,1) @ V_w.reshape(1,n+m))),np.trace(W[j].reshape(n+m,1) @ W[j].reshape(1,n+m) @ (V_w.reshape(n+m,1) @ DV_w[:,1].reshape(1,n+m) + DV_w[:,1].reshape(n+m,1) @ V_w.reshape(1,n+m)))])
#     for k in range(n+m):
#         if k != j:
#             Tr6 = Tr6 + (log_ou_0(np.array([Lz[j]]))[0] - log_ou_0(np.array([Lz[k]]))[0])/(Lz[j] - Lz[k]) * W[j,:n] @ W[k,:n] * ( V_w @ W[j].reshape(n+m,1) @ W[k].reshape(1,n+m)  @ DV_w + DV_w.transpose() @ W[j].reshape(n+m,1) @ W[k].reshape(1,n+m)  @ V_w)
# print(Tr3)
# print(Tr5)
# print(Tr55)
# print()
    
000000000000000000000000000000000000
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