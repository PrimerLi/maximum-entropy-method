import kernel 
import numpy as np
import solver
import generateSpectral
import generateGreenFunction
import printFile
import chi
import entropy

def normalize(omega, A_original):
    s = 0.0
    for i in range(len(omega) - 1):
        s += 0.5*(omega[i+1] - omega[i])*(A_original[i+1] + A_original[i])
    return A_original*(1.0/s)

def generateModel(Nomega, modelFileName):
    omega = np.zeros(Nomega)
    A = np.zeros(Nomega)
    omegaLower = -5
    omegaUpper = 5
    domega = (omegaUpper - omegaLower)/float(Nomega - 1)
    for i in range(Nomega):
        omega[i] = omegaLower + i*domega
        A[i] = generateSpectral.gauss(omega[i], 0, 1)
    printFile.printFile(omega, A, modelFileName)
    return omega, A

def main():
    import sys
    if (len(sys.argv) != 3):
        print "alphaInitial = sys.argv[1], lambd = sys.argv[2]. "
        return -1

    alphaInitial = float(sys.argv[1])
    omega, temp = generateGreenFunction.readA("A.txt")
    model = map(lambda ele: generateGreenFunction.gauss(ele, 0, 1), omega)
    printFile.printFile(omega, model, "model.txt")

    omega_n, GReal = generateGreenFunction.readA("G_real.txt")
    omega_n, GImag = generateGreenFunction.readA("G_imag.txt")
    Niom = len(omega_n)

    KReal = kernel.K_matrix_real(omega_n, omega)
    KImag = kernel.K_matrix_imag(omega_n, omega)

    LambdaInverse = np.zeros((Niom, Niom))
    for i in range(Niom):
        LambdaInverse[i, i] = 1.0
    
    alpha = alphaInitial
    A = model
    eps = 1.0e-4

    lambd = float(sys.argv[2])

    A_averaged = np.zeros(len(A))

    while(alpha > 1.0e-7):
        print "alpha = " + str(alpha)
        AUpdated = solver.solver(alpha, GReal, GImag, KReal, KImag, omega, A, model, LambdaInverse, lambd)
        error = np.linalg.norm(AUpdated - A)
        print "error = " + str(error)
        printFile.printFile(omega, A, "A_final_" + str(lambd) + ".txt")
        alpha = alpha/1.1
        #lambd = max(lambd/1.03, 1.0e-2)
        A = normalize(omega, AUpdated)
        deviation = chi.chi(GReal, GImag, KReal, KImag, A, omega, LambdaInverse)
        S = entropy.entropy(omega, A, model)
        A_averaged = A_averaged + A*np.exp(alpha*S - 0.5*deviation)
        A_averaged = normalize(omega, A_averaged)
        #if (error < eps):
        #    break
    printFile.printFile(omega, A_averaged, "A_averaged_" + str(lambd) + ".txt")
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
