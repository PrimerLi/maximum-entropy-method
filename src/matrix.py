import numpy as np

def isIdentity(matrix):
    shape = matrix.shape
    if (shape[0] != shape[1]):
        return False
    dimension = shape[0]
    identity = np.zeros(shape)
    for i in range(dimension):
        identity[i,i] = 1
    diff = matrix - identity
    eps = 1.0e-12
    return np.linalg.norm(diff) < eps

def isOrthogonal(matrix): 
    shape = matrix.shape
    if (shape[0] != shape[1]):
        return False
    return isIdentity(matrix.dot(np.transpose(matrix)))

def isSymmetric(matrix):
    shape = matrix.shape
    if (shape[0] != shape[1]):
        return False
    dimension = shape[0]
    diff = matrix - np.trasnspose(matrix)
    eps = 1.0e-12
    return np.linalg.norm(diff) < eps
