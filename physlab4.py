import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math
radii = [
	5.77 * 10**(-2),
	5.15 * 10**(-2),
	4.51 * 10**(-2),
	3.88 * 10**(-2),
	3.24 * 10**(-2),
]
uncertCurrent = 0.2
uncertVoltage = 2
C = 0.0001961838899

def func(x, a, b,c):
    return (a*x) + b

def bestFit(x, y):
    x0 = np.array([0.0, 0.0, 0.0])
    yerr = [uncertCurrent] * len(y)
    w, er = opt.curve_fit(func, x, y, x0, sigma=yerr,  bounds=(-200, [3000., 2000., 1000.])) 
    return w, er

# popt, pcov = curve_fit(func, xdata, ydata)

# print(pcov)


def getData(filename):
	raw = pd.read_csv(filename)
	v40 = raw.values[0]
	v60 = raw.values[1]
	v80 = raw.values[2]
	v100 = raw.values[3]
	data = np.transpose(raw.values)
	voltage = data[0]
	bar1 = data[1]
	bar2 = data[2]
	bar3 = data[3]
	bar4 = data[4]
	bar5 = data[5]
	return voltage, bar1, bar2, bar3, bar4, bar5, v40[1:], v60[1:], v80[1:], v100[1:]

def plotVoltageCurrent(x,y):
	#add erorrs from sheet
	y_error =  0.2#Current
	x_error = 2 #Voltage
	w, er = bestFit(x,y)
	a = w[0]
	b = w[1]
	yfit = (a*x) + b
	plt.plot(x,yfit)
	plt.errorbar(x, y, yerr=y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x,y,'.')

def sigmaAandD(x, sigma):
	xSquared = [p**2 for p in x]
	D = len(x)*(sum(xSquared)) - (sum(x))**2
	sigmaA = sigma*math.sqrt((len(x)/D))
	sigmaD = sigma*math.sqrt(sum(xSquared)/D)
	return sigmaA, sigmaD


def findSigmaBandSigmaEM(sigmaA, sigmaD, A, D, C, uncertV, V):
	sigmaBe = sigmaD*(C)
	sigmaEM = math.sqrt(
		(uncertV**2)*((2/((A**2)*(C**2)))**2) + 
		(sigmaA**2)*(((-4*V)/((A**3)*(C**2)))**2))
	return sigmaBe, sigmaEM

def findEMandBE(C, A, D, V):
	Be = D * C
	EM = (2*V)/((A**2)*(C**2))
	return EM, Be


def plotCurrentAndRadius(x, y, voltage):
	y_error = 0.2
	w, er = bestFit(x,y)
	a = w[0]
	b = w[1]
	A = a
	D = b 
	errorArr = np.sqrt(np.diag(er))
	sigmaA, sigmaD = errorArr[0], errorArr[1]
	sigmaBe, sigmaEM = findSigmaBandSigmaEM(sigmaA, sigmaD, A, D, C, uncertVoltage, voltage)
	EM, Be = findEMandBE(C, A, D, voltage)
	print("====== Voltage ",voltage, "======")
	print("EM 10^-11:", EM *(10**(-11)))
	print("BE: ", Be )
	print("Sigma EM 10^-11: ",sigmaEM*(10**(-11)))
	print("Sigma BE: ", sigmaBe)
	print("A: ", A)
	print("sigma A: ", sigmaA)
	print("D: ", D)
	print("sigma D: ", sigmaD)
	yfit = ([a*p for p in x]) + b
	plt.errorbar(x, y, yerr=y_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x,yfit)
	plt.plot(x,y,'.')
	plt.xlabel('1/r (m^-1)') #latex it
	plt.ylabel('Current (A)')
	plt.title("Current and Radii Plot For Different Voltages")
	return EM, Be, sigmaEM, sigmaBe

def findWeighedMeandAndWeighedError(x, sigma):
	xMean = sum((xi/(sigmai**2)) for xi,sigmai in zip(x,sigma))/sum([1/(p**2) for p in sigma])	
	sigmaMean = math.sqrt(1/(sum([1/(p**2) for p in sigma])))
	return xMean, sigmaMean

voltage, bar1, bar2, bar3, bar4, bar5, v40, v60, v80, v100 = getData('data.csv')

#V and I plot
#plotVoltageCurrent(voltage, bar1)
#plotVoltageCurrent(voltage, bar2)
#plotVoltageCurrent(voltage, bar3)
#plotVoltageCurrent(voltage, bar4)
#plotVoltageCurrent(voltage, bar5)
#plt.ylabel('Current (A)')
#plt.xlabel('Voltage (V)')
#plt.title("Voltage and Current Measurements")
# plt.show()

invRadii =  [1/x for x in radii]

EM1, Be1, sigmaEM1, sigmaBe1 = plotCurrentAndRadius(invRadii, v40, 40)
EM2, Be2, sigmaEM2, sigmaBe2 = plotCurrentAndRadius(invRadii, v60, 60)
EM3, Be3, sigmaEM3, sigmaBe3 = plotCurrentAndRadius(invRadii, v80, 80)
EM4, Be4, sigmaEM4, sigmaBe4 = plotCurrentAndRadius(invRadii, v100, 100)
#plt.show()

EMList = [EM1, EM2, EM3, EM4]
BeList = [Be1, Be2, Be3, Be4]
sigmaEMList = [sigmaEM1, sigmaEM2, sigmaEM3, sigmaEM4]
sigmaBeList = [sigmaBe1, sigmaBe2, sigmaBe3, sigmaBe4]

meanBe, sigmaBe = findWeighedMeandAndWeighedError(BeList, sigmaBeList)
meanEM, sigmaEM = findWeighedMeandAndWeighedError(EMList, sigmaEMList)
print("=================")
print("MeanBE: ", meanBe)
print("SigmaBe: ", sigmaBe)
print("MeanEM * 10^-11: ", meanEM*(10**(-11)))
print("SigmaEM * 10^-11: ", sigmaEM*(10**(-11)))


