import J
import f
import nan

def diff(a, b):
    import numpy as np
    s = 0.0
    for i in range(len(a)):
        s = s + (a[i] - b[i])**2
    return np.sqrt(s)

def newton(alpha, G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, omega, A_initial, D, LambdaInverse):
    import numpy as np
    
    Nomega = len(omega)
    
    iterationMax = 30
    counter = 0
    
    eps = 0.00001
    while(True):
        counter = counter + 1
        if (counter > iterationMax):
            break
        function = f.f(alpha, G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, A_initial, D, omega, LambdaInverse)
        Jacobian = J.J(alpha, A_initial, omega, K_real_rotated, K_imag_rotated, LambdaInverse)
        A_updated = A_initial - np.linalg.inv(Jacobian).dot(function)
        for i in range(len(A_updated)):
            if (A_updated[i] < 0):
                A_updated[i] = 1.0e-14
        if (nan.array_isnan(A_updated)):
            A_initial = A_updated
            break
        error = diff(A_initial, A_updated)
        if (error < eps):
            break
        print "counter = ", counter, ", diff = ", error
        A_initial = A_updated
    
    return A_initial
