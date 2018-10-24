import numpy as np
import sys

def createPD(dimension):
    import random
    shape = (dimension, dimension)
    A = np.zeros(shape)
    for i in range(dimension):
        for j in range(dimension):
            A[i, j] = random.randint(-2*dimension, 2*dimension)
    return np.transpose(A).dot(A)

def norm(vector):
    return np.sqrt(vector.dot(vector))

def conjugateGradient(A, b):#A must be positive definite. 
    length = len(b)
    shape = A.shape               
    if (shape[0] != shape[1] or shape[1] != length):
        sys.exit("Incompatible dimension in conjugate gradient method. ")
    x0 = np.zeros(length)
    r0 = b - A.dot(x0) 
    p0 = r0
    iterationMax = 100
    counter = 0
    eps = 1.0e-12
    while(True):
        counter = counter + 1
        if (counter > iterationMax):
            break
        alpha = r0.dot(r0)/(p0.dot(A.dot(p0)))
        x1 = x0 + alpha*p0
        r1 = r0 - alpha*A.dot(p0)
        residual = norm(r1)
        if (residual < eps):
            break
        beta = r1.dot(r1)/(r0.dot(r0))
        p1 = r1 + beta*p0
        x0 = x1
        r0 = r1
        p0 = p1
    return x0

def main():
    import sys
    import random

    dimension = 101
    A = createPD(dimension)
    b = np.zeros(dimension)
    for i in range(dimension):
        b[i] = random.randint(0, dimension)
    solution = conjugateGradient(A, b)
    test = A.dot(solution)
    for i in range(len(b)):
        print b[i], "    ", test[i]
    print norm(b - test)
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
