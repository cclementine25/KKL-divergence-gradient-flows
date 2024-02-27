import numpy as np
import scipy.stats as scs
import matplotlib.pyplot as plt
import scipy

import divergences as dv
import kernels as kl
import gradient_descent as gd
import generate_y as gy

##############################
######## PARAMETERS ##########
##############################

d = 2 #dimension of the particles 
n = 10 # nombre de particules pour q
m = 10 # nombre de particules pour p
T = 60 # nombre d'itérations
h = 0.1 # stepsize gradient descent
eps = 0.0001
alpha = 0.01
sigm = 1

Alpha = [0.00001,0.001,0.01,0.1,0.5]
N = [5,10,30,40]
Sigma = [0.5,1,3,5,10,20]

config_y = lambda : gy.gaussian(gy.muy, gy.Sigmay, m)

####### INITIAL DISTRIBUTIONS P AND Q  ########

x0 = scs.multivariate_normal.rvs(gy.mux,gy.Sigmax,n)
#x0 = scs.multivariate_normal.rvs(0.1 *gy.mux,0.02*np.identity(2),n) 
#x0 = scs.multivariate_normal.rvs(np.array([np.cos(-3 * np.pi/4),np.sin(-3* np.pi/4)]),0.01 * np.identity(2),n)
y = config_y()


### KERNEL ###
 #np.max(np.linalg.norm(X-Y,axis = 1)) / (np.sqrt(1000 * np.log(10))) # np.abs(np.mean(np.linalg.norm(X,axis = 1)) - np.mean(np.linalg.norm(Y,axis = 1)))#max(2,np.linalg.norm(np.mean(x) - np.mean(y)))
k = lambda x,y :  kl.k_gauss(x,y,sigm)
dk = lambda x,y : kl.dk_gauss(x,y,sigm)


#### Matrice Ky, eigenvalues and eigenvectors ####
Ky = 1/m * np.array([[k(y[i],y[j]) for i in range(m)] for j in range(m)])
Ly,V = np.linalg.eig(Ky)
V = V.transpose()
Ly = np.real(Ly)
Packy = [Ky,Ly,V]

#################################################
############# PLOTS #############################
#################################################


######### Variation alpha #############

fig, axs = plt.subplots(2,figsize = (20,20))
fig.suptitle(r"$p$ and $q$ gaussians with parameters " + r"$\sigma = $ " + str(sigm) + r"$\gamma = $" + str(h))
for alpha in Alpha:
    J = lambda x : dv.KKL(x, y, k,Packy,alpha) 
    dJ = lambda x : np.array([dv.WGrad_KKL(x[i],x, y, k, dk,Packy,alpha) for i in range(n)])
    X,l_J,Grad = gd.gradient_descent(J, dJ, x0, h, eps, T)
    axs[0].plot(l_J,label = r"$\alpha = $ " + str(alpha))
    axs[0].legend()
    axs[1].plot(Grad,label = r"$\alpha = $ " + str(alpha))
    axs[1].legend()
    plt.figure()
    plt.scatter(x0[:,0],x0[:,1],color = "blue",label = r"$p_0$")
    plt.scatter(X[-1,:,0],X[-1,:,1],color = "cyan",label = r"$p_{fin}$")
    plt.scatter(y[:,0],y[:,1],color = "red",label =r"$q$")
    for i in range(n):
        plt.plot(X[:,i,0],X[:,i,1],color = "black",linewidth = 0.2)
    plt.title(r"$\alpha = $" + str(alpha))
    plt.legend()
axs[0].set_title("Values of "+r" $KKL_{\alpha}$" +" for different values of " +r"$\alpha$")
axs[1].set_title("Norm of the gradient of "r"$KKL_{\alpha}$" + "for different values of"+ r"$\alpha$")

    
########## Variation of sigma ##################



# fig, axs = plt.subplots(2,figsize = (20,20))
# fig.suptitle(r"$p$ and $q$ gaussians with parameters " + r"$n = $ " + str(n) + r"$\gamma = $" + str(h))
# for sigm in Sigma:
#     k = lambda x,y :  kl.k_gauss(x,y,sigm)
#     dk = lambda x,y : kl.dk_gauss(x,y,sigm)
#     Ky = 1/m * np.array([[k(y[i],y[j]) for i in range(m)] for j in range(m)])
#     Ly,V = np.linalg.eig(Ky)
#     V = V.transpose()
#     Ly = np.real(Ly)
#     Packy = [Ky,Ly,V]


#     J = lambda x : dv.KKL(x, y, k,Packy,alpha) 
#     dJ = lambda x : np.array([dv.WGrad_KKL(x[i],x, y, k, dk,Packy,alpha) for i in range(n)])
#     X,l_J,Grad = gd.gradient_descent(J, dJ, x0, h, eps, T)
#     axs[0].plot(l_J,label = r"$\sigma = $ " + str(sigm))
#     axs[0].legend()
#     axs[1].plot(Grad,label = r"$\sigma = $ " + str(sigm))
#     axs[1].legend()
#     plt.figure()
#     plt.scatter(x0[:,0],x0[:,1],color = "blue",label = r"$p_0$")
#     plt.scatter(X[-1,:,0],X[-1,:,1],color = "cyan",label = r"$p_{fin}$")
#     plt.scatter(y[:,0],y[:,1],color = "red",label =r"$q$")
#     for i in range(n):
#         plt.plot(X[:,i,0],X[:,i,1],color = "black",linewidth = 0.2)
#     plt.title(r"$\sigma = $ " + str(sigm))
#     plt.legend()
# axs[0].set_title("Values of "+r" $KKL_{\alpha}$" +" for different values of " +r"$\sigma$")
# axs[1].set_title("Norm of the gradient of "r"$KKL_{\alpha}$" + "for different values of"+ r"$\sigma$")