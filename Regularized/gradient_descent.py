import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as scs
import scipy.optimize as sco


def gradient_descent(J,dJ,x0,h,eps,n_it_max):
    x = x0 
    grad = dJ(x)
    J_val = J(x)
    X = [x0]
    i = 0 
    liste_J = []
    Grad = []
    while np.linalg.norm(grad) > eps and i < n_it_max:
        print("hop")
        # print(np.linalg.norm(grad))
        # print(J_val)
        # print()
        #y = x #+ 0 * scs.multivariate_normal.rvs(np.zeros(2),np.identity(2),len(x))
        Grad.append(np.linalg.norm(grad))
        liste_J.append(J_val)
        grad = dJ(x)
        #h = sco.line_search(J, dJ, x, -grad)
        x = x - h * grad
        X.append(x)
        J_val = J(x)
        i = i + 1
        #print(x)
        #plt.figure()
        #plt.scatter(x[:,0],x[:,1])
        #plt.scatter(y_0[:,0],y_0[:,1])
    return np.array(X),liste_J,Grad
        
def GD(J,dJ,x0,h,eps,n_it_max):
    x = x0 
    grad = dJ(x)
    X = [x0]
    i = 0 
    liste_J = []
    Grad = []
    while np.linalg.norm(grad) > eps and i < n_it_max:
        liste_J.append(J(x))
        y = x #+ 0 * scs.multivariate_normal.rvs(np.zeros(2),np.identity(2),len(x))
        grad = dJ(y)
        # grad = np.real(grad)
        Grad.append(np.linalg.norm(grad))
        # print(np.real(sco.line_search(J, dJ, x, -h*grad,gfk = grad)))
        # print(h)
        x = x - h * grad
        X.append(x)
        i = i + 1
    return np.array(X),liste_J,Grad
        