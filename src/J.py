def J(alpha, A, omega, K_real_rotated, K_imag_rotated, LambdaInverse, lambd):
    import numpy as np
    
    Nomega = len(omega)
    domega = omega[1] - omega[0]
    result = np.zeros((Nomega, Nomega))
    for i in range(Nomega-2):
        result[i, i+2] = 0.25*lambd/domega
    for i in range(2, Nomega):
        result[i-2, i] = 0.25*lambd/domega
    
    result = -np.transpose(K_real_rotated).dot(LambdaInverse).dot(K_real_rotated)
    result = result -np.transpose(K_imag_rotated).dot(LambdaInverse).dot(K_imag_rotated)
    for i in range(Nomega):
        result[i,i] = result[i,i] - alpha*domega/A[i] - 0.5*lambd/domega
    return result
