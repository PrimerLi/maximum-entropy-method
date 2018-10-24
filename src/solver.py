import conjugateGradient
import f
import J
import numpy as np

def solver(alpha, G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, omega, A_initial, D, LambdaInverse, lambd):
    eps = 1.0e-15
    counter = 0
    iterationMax = 30

    while(True):
        counter = counter + 1
        if (counter > iterationMax):
            break
        Jacobian = J.J(alpha, A_initial, omega, K_real_rotated, K_imag_rotated, LambdaInverse, lambd)
        function = f.f(alpha, G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, A_initial, D, omega, LambdaInverse, lambd)
        det = np.linalg.det(Jacobian)
        if (abs(det) > 1.0e-3):
            solution = np.linalg.inv(-Jacobian).dot(function)
        else:
            solution = conjugateGradient.conjugateGradient(-Jacobian, function)
        error = conjugateGradient.norm(solution)
        print "counter = ", counter, ", norm of gradient = ", np.linalg.norm(function)
        #print "counter = ", counter, ", error = ", error, ", det of Hessian matrix = ", det, ", norm of gradient = ", np.linalg.norm(function)
        if (error < eps):
            break
        A_updated = A_initial + solution
        for i in range(len(A_updated)):
            if (A_updated[i] < 0):
                A_updated[i] = 1.0e-10
        A_initial = A_updated
    return A_initial
