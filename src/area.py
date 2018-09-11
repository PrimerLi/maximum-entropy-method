#!/usr/bin/env python
def printFile(x, y, output):
    ofile = open(output, "w")
    for i in range(min(len(x), len(y))):
        ofile.write(str(x[i]) + "    " + str(y[i]) + "\n")
    ofile.close()

def integrate(omega, A, N):
    import sys
    if (N > len(omega)-1):
        sys.exit("N is too large. ")
    s = 0.0
    for i in range(N):
        s = s + 0.5*(A[i+1] + A[i])*(omega[i+1] - omega[i])
    return s

def main():
    import os
    import sys

    if (len(sys.argv) != 2):
        print "fileName = sys.argv[1]. "
        return -1

    fileName = sys.argv[1]
    omega = []
    A = []
    ifile = open(fileName, "r")
    for i,string in enumerate(ifile):
        a = string.split()
        omega.append(float(a[0]))
        A.append(float(a[1]))
    ifile.close()

    area = []
    for i in range(len(omega)-1): 
        area.append(integrate(omega, A, i))
    printFile(omega, area, "area_omega.txt") 

    return 0

main()
