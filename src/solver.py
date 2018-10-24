import conjugateGradient
import f
import J
import numpy as np

def adam(alpha, G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, omega, A_initial, D, LambdaInverse, lambd):
    eps = 1.0e-8
    counter = 0
    iterationMax = 1000

    rate = -0.01
    beta1 = 0.9
    beta2 = 0.99
    m = np.zeros(len(omega))
    v = np.zeros(len(omega))

    while(True):
        counter += 1
        if (counter > iterationMax):
            error = np.linalg.norm(gradient)
            print "counter = " + str(counter) + ", norm of gradient = " + str(error)
            break
        gradient = f.f(alpha, G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, A_initial, D, omega, LambdaInverse, lambd)
        m = beta1*m + (1 - beta1)*gradient
        v = beta2*v + (1 - beta2)*gradient**2
        m = m/(1 - beta1**counter)
        v = v/(1 - beta2**counter)
        A_initial = A_initial - rate*m/(np.sqrt(v) + eps)
        for i in range(len(A_initial)):
            if (A_initial[i] < 0):
                A_initial[i] = 1.0e-10
        #error = np.linalg.norm(gradient)
        #print "counter = " + str(counter) + ", norm of gradient = " + str(error)
        #if (error < 1.0e-3):
        #    break
    return A_initial

def solver(alpha, G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, omega, A_initial, D, LambdaInverse, lambd):
    eps = 1.0e-7
    counter = 0
    iterationMax = 30

    while(True):
        counter = counter + 1
        if (counter > iterationMax):
            #result = adam(alpha, G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, omega, A_initial, D, LambdaInverse, lambd)
            #return result
            break
        Jacobian = J.J(alpha, A_initial, omega, K_real_rotated, K_imag_rotated, LambdaInverse, lambd)
        function = f.f(alpha, G_real_rotated, G_imag_rotated, K_real_rotated, K_imag_rotated, A_initial, D, omega, LambdaInverse, lambd)
        #solution = conjugateGradient.conjugateGradient(-Jacobian, function)
        eigenvalues = np.linalg.eigvals(-Jacobian)
        eigen_max = max(eigenvalues)
        eigen_min = min(eigenvalues)
        print "max and min of eigenvalues of Hessian matrix are " + str(eigen_max) + ", " + str(eigen_min)
        if (abs(eigen_max) > 1.0e-3 and abs(eigen_min) > 1.0e-3):
            solution = np.linalg.inv(-Jacobian).dot(function)
        else:
            solution = conjugateGradient.conjugateGradient(-Jacobian, function)
        #error = np.linalg.norm(solution)
        error = np.linalg.norm(function)
        print "counter = ", counter, ", norm of gradient = ", np.linalg.norm(function) #, ", error = ", error
        #print "counter = ", counter, ", error = ", error, ", det of Hessian matrix = ", det, ", norm of gradient = ", np.linalg.norm(function)
        if (error < eps):
            break
        A_updated = A_initial + solution
        for i in range(len(A_updated)):
            if (A_updated[i] < 0):
                A_updated[i] = 1.0e-10
        A_initial = A_updated
    return A_initial
