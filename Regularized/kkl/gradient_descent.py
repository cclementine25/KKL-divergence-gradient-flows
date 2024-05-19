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
        #if i%10 == 0:
           # print("T = " + str(i))
        Grad.append(np.linalg.norm(grad))
        liste_J.append(J_val)
        grad = dJ(x)
        x = x - h * grad
        X.append(x)
        J_val = J(x)
        i = i + 1
    return np.array(X),liste_J,Grad
        
