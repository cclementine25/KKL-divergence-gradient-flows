import numpy as np
import scipy.stats as scs
import matplotlib.pyplot as plt
import scipy.optimize as sco 

import divergences as dv
import kernels as kl
import generate_y as gy 

import cProfile

d = 2 #dimension of the particles 
n = 30 #nombre de particules pour p
m = 30 #nombre de particules pour q
alpha = 0.1

########## EXPE ##########
expe_var_mu = False
expe_var_n = True
linesearch = False
test_lipschitz = False
expe_var_alpha = False


## KERNEL ###
sigma = 1
k = lambda x,y :  kl.k_gauss(x,y,sigma)
dk = lambda x,y : kl.dk_gauss(x, y, sigma)



#######SIMULATION OF THE DATA#########

#Simulation of (X_i)_i<n ~ p
mux = np.array([10,10])
Mux = np.array([1/k *mux for k in range(1,17)])
Lx = np.array([[1/2,1/3],[1/4,-2]])
Sigmax = Lx @ Lx.transpose()
SS = np.array([Sigmax * k for k in [6,5,4,3,2,1]])

#Simulation of (Y_i)_i<n ~ q 
muy = np.array([0,0])
Ly = np.array([[1/5, -1],[1/2,1/2]])
Sigmay = Ly @ Ly.transpose() #Lx @ Lx.transpose()
#Y = scs.multivariate_normal.rvs(muy,Sigmay,m)

#######################################################

###Generation of y ###
y = gy.rings(1, 1, 0.5, 1, m)
#x = scs.multivariate_normal.rvs(mux,Sigmax,n) 
x = scs.multivariate_normal.rvs(0.5 * np.array([np.cos(-3 * np.pi/4),np.sin(-3* np.pi/4)]),0.01 * np.identity(2),n)


###### Ky ####
Ky = 1/m * k(y,y) #1/m * np.array([[k(y[i],y[j],sigm(x,y)) for i in range(m)] for j in range(m)])
Ly,V = np.linalg.eig(Ky)
V = V.transpose()
Ly = np.real(Ly)
Packy = [Ky,Ly,V]


########PLOTS#########


J = lambda x : dv.KKL(x, y, k,Packy,alpha) 
dJ = lambda x : np.array([dv.WGrad_KKL(x[i],x, y,k, dk,Packy,alpha,sigma) for i in range(n)])


""" Here we plot two gaussian distributions with same variance and
different mean and we draw the evolution of MMD, KKL and KDE when the
 means of one of the  distribution become closer to the other one"""


if expe_var_mu:
    mmd = []
    kkl = []
    kde = []
    k_trace = []
    trxx = []
    dkkl = []
    
    #sig = [0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1,1.5,2,3,4]
    
    fig, axs = plt.subplots(4, 4, figsize=(20,20))
    for j in range(16):
        axs[j//4,j%4].axis([-3,15,-10,20])
        X = scs.multivariate_normal.rvs(Mux[j],Sigmax,n)
        #mmd.append(dv.MMD(X,y,k))
        kkl.append(dv.KKL(X,y,k,Packy,alpha))
        dkkl.append(np.linalg.norm(dJ(X)))
        cProfile.run('dJ(X)')
        axs[j//4,j%4].scatter(y[:,0],y[:,1],color = "blue")
        axs[j//4,j%4].scatter(X[:,0],X[:,1],color = "red")
        print(j)
    
        
    # 
        
    plt.figure()
    #plt.plot(mmd,label = "mmd")
    plt.plot(kkl,label = "kkl")
    #plt.title("evolution of mmd for 2 distribution of same variance when their means get closer ")
    
    plt.figure()
    plt.plot(kkl,label = "kkl")
    #plt.plot(trxx)
    #plt.plot(k_trace,label = "k_trace")
    plt.legend()
    plt.title("kkl for sigma = " + str(sigma))
    #plt.title("evolution of kkl for 2 distribution of same variance when their means get closer ")
    
    # plt.figure()
    # plt.plot(kde,label = "kde")
    # plt.title("evolution of kde for 2 distribution of same variance when their means get closer ")
    
    plt.figure()
    plt.plot(dkkl,label = "dkkl")
    """ same experience switching the roles of the mean and variance"""


if expe_var_n:
    
    N = [5,10,15,20,30,40,50,60,80,100]#,120,150]#,200,300,400,500]#,600,700,800,900]#,1000,1100,1300,1400,1500,1600,1800,2000]
    nb_it = 100
    #fig, axs = plt.subplots(3, 4, figsize=(20,20))
    for alpha in [0.0001,0.001,0.01,0.1,0.5]:
        KKL_n = np.zeros(len(N))
        I_KKL_n = np.zeros(len(N))
        #KKL_inv_n = []
        for j in range(len(N)):
            F = np.zeros(nb_it)
            print(j)
            for i in range(nb_it):
                y = gy.mixt_gauss(gy.MU,[np.identity(2),np.identity(2)], [1/2,1/2], N[j])
                x = gy.gaussian(1/5 *mux,Sigmax,N[j])
                
                Ky = 1/N[j] * k(y,y) 
                Ly,V = np.linalg.eig(Ky)
                V = V.transpose()
                Ly = np.real(Ly)
                Packy = [Ky,Ly,V]
                
                J = lambda x,y : dv.KKL(x, y,k,Packy,alpha) 
                dJ = lambda x : np.array([dv.WGrad_KKL(x[i],x, y,k, dk,Packy,alpha,sigma) for i in range(N[j])])
                
                F[i] = J(x,y)
            KKL_n[j] = np.mean(F)
            I_KKL_n[j] = np.std(F)
            
            #KKL_inv_n.append(J(y,x))
            
            #axs[j//4,j%4].scatter(x[:,0],x[:,1])
            #axs[j//4,j%4].scatter(y[:,0],y[:,1])
            
            
        #plt.figure()
        plt.plot(N,KKL_n,label = r"$\alpha =$ " + str(alpha))
        plt.fill_between(N,KKL_n - I_KKL_n,KKL_n + I_KKL_n,alpha = 0.5)
        #plt.plot(np.linspace(0,150,100), 0.2/ np.sqrt(np.linspace(0,150,100))+ 1.54)
        plt.title("Evolution of " + r"$KKL_{\alpha}$" + " for sets of points from gaussians distribution and mixture of gaussian \nwith increasing number of particules " + r"$n$." + " Parameters : " + r"$\sigma =$ " + str(sigma) + r" $\alpha = $" + str(alpha) ,fontsize = 16)
        plt.xlabel('number of particules')
        plt.ylabel(r'$KKL_{\alpha}$')
        plt.legend()

    # plt.figure()
    # plt.plot(N,KKL_inv_n)
if linesearch:
    X= []
    l_J = []
    ED = []
    
    def callback(x):
        X.append(np.array([x[:n],x[n:]]).T)
        l_J.append(J1(x))
        ED.append(dv.energy_distance(np.array([x[:n],x[n:]]).T,y))
        
    options = {'maxiter': 100}
        
        
        
    x0 = np.hstack([x[:,0],x[:,1]])
        
        
    J1 = lambda x : dv.KKL(np.array([x[:n],x[n:]]).T, y, k,Packy,alpha) 
    dJ1 = lambda x : np.hstack([np.array([dv.WGrad_KKL(np.array([x[i],x[i+n]]),np.array([x[:n],x[n:]]).T, y,k, dk,Packy,alpha,sigma) for i in range(n)])[:,0], np.array([dv.WGrad_KKL(np.array([x[i],x[i+n]]),np.array([x[:n],x[n:]]).T, y,k, dk,Packy,alpha,sigma) for i in range(n)])[:,1]]) 
    
    JJ = lambda x : (J1(x), dJ1(x))
    result = sco.minimize(J1,x0,jac = dJ1,options= options, callback = callback)
    cProfile.run('sco.minimize(J1,x0,jac = dJ1,options= options, callback = callback)')
    
    x_fin = np.array([result.x[:n],result.x[n:]]).T
    
    X = np.array(X)
    l_J = np.array(l_J)
    
    plt.scatter(x_fin[:,0],x_fin[:,1])
    plt.scatter(y[:,0],y[:,1])
    
    
    fig, axs = plt.subplots(5, 4, figsize=(20,20))
    for i in range(0,len(X)-1,len(X)//20):
        j = i//(len(X)//20)
        #axs[j//4,j%4].axis([-3,3.5,-4,1])
        axs[j//4,j%4].scatter(y[:,0],y[:,1],color = "orange")
        axs[j//4,j%4].scatter(X[i,:,0], X[i,:,1], color = "blue")
    
    
    plt.figure()    
    plt.plot(l_J)
    plt.title(r"values of $KKL_{\alpha}")


if test_lipschitz:
    kkl = []
    dkkl = []
    T = []

    
    fig, axs = plt.subplots(4, 4, figsize=(20,20))
    for j in range(16):
        axs[j//4,j%4].axis([-3,15,-10,20])
        X = scs.multivariate_normal.rvs(Mux[j],Sigmax,n)
        T.append(np.linalg.norm(dJ(X)) / scs.wasserstein_distance(X, y))
        #mmd.append(dv.MMD(X,y,k))
        kkl.append(dv.KKL(X,y,k,Packy,alpha))
        dkkl.append(np.linalg.norm(dJ(X)))
        axs[j//4,j%4].scatter(y[:,0],y[:,1],color = "blue")
        axs[j//4,j%4].scatter(X[:,0],X[:,1],color = "red")
        print(j)
    
        
    plt.plot(T)
    
    
if expe_var_alpha:
    #Alpha = [0.01,0.005,0.001,0.0001,0.00001,0.000001]
    Alpha = [0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0.01,0.005,0.001,0.0001,0.00001,0.000001]
    nb_it = 100
    KKL_a = np.zeros(len(Alpha))
    I_KKL_a = np.zeros(len(Alpha))
    for j in range(len(Alpha)):
        F = np.zeros(nb_it)
        print(j)
        for i in range(nb_it):
            y = gy.mixt_gauss(gy.MU,[np.identity(2),np.identity(2)], [1/2,1/2], m)
            x = gy.gaussian(1/5 *mux,Sigmax,n)
            
            Ky = 1/m * k(y,y) 
            Ly,V = np.linalg.eig(Ky)
            V = V.transpose()
            Ly = np.real(Ly)
            Packy = [Ky,Ly,V]
            
            J = lambda x,y : dv.KKL(x, y,k,Packy,Alpha[j]) 
            dJ = lambda x : np.array([dv.WGrad_KKL(x[i],x, y,k, dk,Packy,Alpha[j],sigma) for i in range(n)])
            
            F[i] = J(x,y)
        KKL_a[j] = np.mean(F)
        I_KKL_a[j] = np.std(F)
        
        #KKL_inv_n.append(J(y,x))
        
        #axs[j//4,j%4].scatter(x[:,0],x[:,1])
        #axs[j//4,j%4].scatter(y[:,0],y[:,1])
        
        
    #plt.figure()
    plt.plot(np.array(Alpha),KKL_a)
    plt.fill_between(np.array(Alpha),KKL_a - I_KKL_a,KKL_a + I_KKL_a,alpha = 0.5)
    #plt.plot(np.linspace(0,150,100), 0.2/ np.sqrt(np.linspace(0,150,100))+ 1.54)
    plt.title("Evolution of " + r"$KKL_{\alpha}$" + " with " + r"$1- \alpha$" ,fontsize = 16)
    plt.xlabel(r"$1-\alpha$")
    plt.ylabel(r'$KKL_{\alpha}$')
    plt.legend()

    