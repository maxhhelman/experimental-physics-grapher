import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math

bSigmas = []
BigI = []
bList = []

def func(x, a, b, c):
    return a + b*x


def SigmaB(i, m):
    L = 0.101
    g = 9.8
    sigmaM = 0.00001
    sigmaSmallI = 0.04
    sigmaL = 0.01 
    partialM = (g*sigmaL/(i*L))**2
    partialI = (-m*g*sigmaSmallI/(i*i*L))**2
    partialL = (-m*g*sigmaSmallI/(i*L*L))**2
    sigmaB = partialM + partialI + partialL
    return math.sqrt(sigmaB)


def graph(filename, color):
    x0 = np.array([0.0, 0.0, 0.0])
    L = 0.101
    g = 9.8
    data = pd.read_csv(filename)
    mass = np.transpose(data.values)[0]
    littleI = np.transpose(data.values)[1]
    bigI = np.transpose(data.values)[2]
    x_error = [0.04*0.01] * len(littleI)
    y_error = [0.00001*g] * len(mass)
    iL = littleI * L
    mg = mass * g
    w, _ = opt.curve_fit(func, iL, mg, x0, y_error) 
    bList.append(w[1])
    BigI.append(bigI)

    print("Estimated parameters", w)
    yfit = func(iL,*w)
    plt.errorbar(iL, mg, yerr=y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
    plt.plot(iL,yfit,label="Fit")

graph('lab3_5.csv', '.r')
graph('lab3_4.csv', '.g')
graph('lab3_3.csv', '.b')
graph('lab3_2.csv', '.y')


plt.xlabel('iL')
plt.ylabel('mg')
plt.title('Measured Values of iL and mg for All Trials')
plt.show()

plt.title('B vs. I for all Trials')
plt.xlabel('I (A)')
plt.ylabel('B (T)')
plt.plot(BigI,bList,'.r-', label="Graph 2")
plt.show()

graph('lab3_5.csv', '.r')
plt.xlabel('iL')
plt.ylabel('mg')
plt.title('Measured Values of iL and mg for I = 5A')
plt.show()

graph('lab3_4.csv', '.g')
plt.xlabel('iL')
plt.ylabel('mg')
plt.title('Measured Values of iL and mg for I = 4A')
plt.show()

graph('lab3_3.csv', '.b')
plt.xlabel('iL')
plt.ylabel('mg')
plt.title('Measured Values of iL and mg for I = 3A')
plt.show()

graph('lab3_2.csv', '.y')
plt.xlabel('iL')
plt.ylabel('mg')
plt.title('Measured Values of iL and mg for I = 2A')
plt.show()
