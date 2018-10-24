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
    import os

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
    s = []
    if (os.path.exists("s.txt")):
        ifile = open("s.txt", "r")
        for (index, string) in enumerate(ifile):
            s.append(float(string))
        ifile.close()

    if (len(s) == 0):
        for i in range(Niom):
            LambdaInverse[i, i] = 1.0
    else:
        assert(len(s) == Niom)
        for i in range(len(s)):
            LambdaInverse[i, i] = 1.0/s[i]
    
    alpha = alphaInitial
    A = model
    eps = 1.0e-4

    lambd = float(sys.argv[2])

    #A_averaged = np.zeros(len(A))

    alpha_values = []
    chi_values = []
    while(alpha >= 1.0e-8):
        print "alpha = " + str(alpha)
        alpha_values.append(alpha)
        AUpdated = solver.solver(alpha, GReal, GImag, KReal, KImag, omega, A, model, LambdaInverse, lambd)
        AUpdated = normalize(omega, AUpdated)
        error = np.linalg.norm(AUpdated - A)
        print "error = " + str(error)
        A = AUpdated
        printFile.printFile(omega, A, "A_alpha_" + str(alpha) + ".txt")
        deviation = chi.chi(GReal, GImag, KReal, KImag, A, omega, LambdaInverse)
        chi_values.append(deviation)
        S = entropy.entropy(omega, A, model)
        #A_averaged = A_averaged + A*np.exp(alpha*S - 0.5*deviation)
        #A_averaged = normalize(omega, A_averaged)
        alpha = alpha/3.0
        #if (error < eps):
        #    break
    #printFile.printFile(omega, A_averaged, "A_averaged_" + str(lambd) + ".txt")
    printFile.printFile(alpha_values, chi_values, "alpha_chi_" + str(lambd) + ".txt")
    printFile.printFile(omega, A, "A_lambda_" + str(lambd) + ".txt")
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
