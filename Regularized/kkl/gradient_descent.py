import numpy as np

def gradient_descent(J,dJ,x0,h,eps,n_it_max):
    x = x0 
    grad = dJ(x)
    J_val = J(x)
    X = [x0]
    i = 0 
    while np.linalg.norm(grad) > eps and i < n_it_max:
        grad = dJ(x)
        x = x - h * grad
        X.append(x)
        i = i + 1
    return np.array(X)
        
