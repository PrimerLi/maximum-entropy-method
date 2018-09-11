def entropy(omega, A, D):
    import numpy as np
    s = 0.0
    for i in range(len(omega)-1):
        s = s + (omega[i+1] - omega[i])*A[i]*np.log(A[i]/D[i])
    return -s
