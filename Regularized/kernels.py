import numpy as np



###################################

def k_gauss(x,y,sigma):
    X_norm = np.sum(x**2, axis=1)
    Y_norm = np.sum(y**2, axis=1)
    xy = np.dot(x, y.T)
    dist = X_norm[:,None] + Y_norm[None,:] - 2 * xy
    return np.exp(-dist/(2*sigma**2))

def dk_gauss(x,y,sigma):
    x_mat = np.tile(x, (len(y), 1, 1))
    y_mat = np.tile(y, (len(x), 1, 1))
    x_mat = np.transpose(x_mat, axes=(1, 0, 2))
    diff = y_mat - x_mat
    return diff/(sigma**2) * k_gauss(x,y,sigma)[:,:,None]
    
    
def k_gauss1(x,y,sigma):
    return np.exp(-np.dot(x-y,x-y) / (2 * sigma**2))

def dk_gauss1(x,y,sigma):
    return -1/sigma**2 * np.exp(-np.dot(x-y,x-y)/(2 * sigma**2)) * (x-y)

# def dkk_gauss(x,y,sigma):
#     return np.zeros(2)

####################################

def k_RBF(x,y,beta):
    return 1/(1 + np.linalg.norm(x-y))**beta

def dk_RBF(x,y,beta):
    if x[0] == y[0]:
        return -2 * beta * (x-y)/ np.linalg.norm(x-y) * k_RBF(x,y,beta + 1)
    return 0


##################################

def k_lin(x,y):
    return np.dot(x,y)

def dk_lin(x,y):
    return y

##################################

def k_poly2(x,y):
    return np.dot(x,x) * np.dot(y,y) + 0.1*np.dot(x,y)

def dk_poly2(x,y):
    return 2 * np.dot(y,y) * x

####################################

def k_(x,y,b):
    return 1 / (1 + np.linalg.norm(x-y)**2)**b
    