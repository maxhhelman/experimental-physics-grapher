import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math


def compute_n_glass(n_air, d, lbd_0, m, theta):
	theta = math.radians(theta)
	up = (n_air*(2*d*n_air -m*lbd_0)*(1-math.cos(theta)))
	low = (2*d*n_air*(1-math.cos(theta)) - m*lbd_0)
	return up/low	

def compute_n_glass_error(err_m,err_th,err_d,n,d,l,m,th):
	th = math.radians(th)
	pi = 3.14159265359
	err_th = (pi*err_th)/(180)
	partial_m = (2*(n**2)*l*d*math.cos(th)*(1-math.cos(th)))/((2*n*d*(1-math.cos(th)) - m*l)**2)
	partial_th = - (n*m*l*math.sin(th)*(2*n*d-m*l))/((2*n*d*(1-math.cos(th)) -m*l)**2)
	partial_d = - (2*l*(n**2)*m*math.cos(th)*(1-math.cos(th)))/((2*n*d*(1-math.cos(th)) -m*l)**2)
	err = math.sqrt( (partial_m**2)*(err_m**2) + (partial_d**2)*(err_d**2) + (partial_th**2)*(err_th**2))
	return err


def valsToWavelength(m, dm):
	return (2/m)*dm

def lambdasErrors(errM, errDM, dmVals):
	result = [math.sqrt( (((2*dm*errM)/(m**2))**2) + (((2*errDM)/(m))**2) ) for dm in dmVals]
	return result

def meanAndErr(x, sigma):
	#weithed mean and error, manual pg 125
	xs = zip(x, sigma)
	x_mean = sum( [i[0]/(i[1]**2) for i in xs] )/sum([1/(sig**2) for sig in sigma])
	x_mean_err = math.sqrt(1/sum([1/(sig**2) for sig in sigma]))
	return x_mean, x_mean_err

def func(x, a, b,c):
    return (a*x) + b

def bestFit(x, y, yerr):
    x0 = np.array([0.0, 0.0, 0.0])
    w, er = opt.curve_fit(func, x, y, x0, sigma=yerr,  bounds=(-200, [3000., 2000., 1000.])) 
    return w, er

def airPlot(x,y, x_error, y_error):
	w, err = bestFit(x,y,y_error)
	A = w[0]
	B = w[1]
	errorArr = np.sqrt(np.diag(err))
	Aerr = errorArr[0]
	Berr = errorArr[1]
	y_fit =  ([A*p for p in x]) + B
	plt.errorbar(x, y, yerr=y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x, y_fit, '-')
	plt.xlabel(r"$P - P_{atm} (cmHg)$", fontsize=12) 
	plt.ylabel(r"$n-n_{atm}$", rotation=90, fontsize=12)
	plt.title(r"$n-n_{atm}$ versus $P - P_{atm}$ Plot")
	plt.show()
	return A, B, Aerr, Berr
	
#Part 1 Errors
errM = 1
errX = 0.5
errDm = math.sqrt(errX**2 + errX**2)
m = 20

#Part 1 Read Data
fpRaw = pd.read_csv('dataFarbyPerot.csv')
fpVals = np.transpose(fpRaw.values)
fpXi = fpVals[1]
fpXf = fpVals[2]
fpDeltaX = fpVals[3]
fp_err_lambdas = lambdasErrors(errM, errDm, fpDeltaX)
fp_lambdas = valsToWavelength(m, fpDeltaX)
fp_lambdaMean, fp_err_lambdaMean = meanAndErr(fp_lambdas, fp_err_lambdas)

miRaw = pd.read_csv('dataMichelson.csv')
miVals = np.transpose(miRaw.values)
miXi = miVals[1]
miXf = miVals[2]
miDeltaX = miVals[3]
mi_err_lambdas = lambdasErrors(errM, errDm, miDeltaX)
mi_lambdas = valsToWavelength(m, miDeltaX)
mi_lambdaMean, mi_err_lambdaMean = meanAndErr(mi_lambdas, mi_err_lambdas)

print("=== Part 1 ===")
print("Farby Perot")
print("Wavelength: ", fp_lambdaMean, " Error: ", fp_err_lambdaMean, " micro meters")
print("Michelson")
print("Wavelength: ", mi_lambdaMean, " Error: ", mi_err_lambdaMean, " micro meters")

#Part 2 Values
d = 3 #cm
P_atm = 76 #cmHg
#https://en.wikipedia.org/wiki/Helium%E2%80%93neon_laser
lambda_0 = 632.991 #nm #POTENTIAL BS
err_p = 1
err_m_air = 1 
airRaw = pd.read_csv('air.csv')
airVals = np.transpose(airRaw.values)
deltaP = airVals[1]
m_air = airVals[2]
nf_ni = (m_air*lambda_0*(10**(-7)))/(2*(d))
err_nf_ni = ((lambda_0*(10**(-7))*err_m_air)/(2*d))

air_A, air_B, air_Aerr, air_Berr = airPlot(deltaP,nf_ni,[err_p] * len(deltaP),[err_nf_ni] * len(nf_ni))
nAir = 1 - air_A*deltaP[0]

print("=== Part 2 ===")
print("C: ", air_A, " Error: ",air_Aerr)
print("n air: ", nAir, " error: ", math.sqrt( ((deltaP[0]*air_Aerr)**2) + ((air_A*err_p)**2) ) )


#Part 3
err_glass_th = 0.5 * 10**(-3)
err_theta = 0.05
err_m_glass = 2

glassRAW = pd.read_csv('glass.csv')
glassVals = np.transpose(glassRAW.values)
thetaI = glassVals[0] #this is always 0
thetaF = glassVals[1] 
m_glass = glassVals[2]
n_air_val = 1.000293
d_glass = 5 * 10**(-3)
lambda_laser = 632.8 * 10**(-9) #from wikipedia

n_glass = [compute_n_glass(n_air_val, d_glass, lambda_laser, m_gl, theta) for theta,m_gl in zip(thetaF,m_glass)]
err_n_glass = [compute_n_glass_error(err_m_glass, err_theta, err_glass_th, n_air_val, d_glass, lambda_laser, mx, thx) for thx,mx in zip(thetaF,m_glass)]
n_glass_mean, err_n_glass_mean = meanAndErr(n_glass, err_n_glass)
print("=== Part 3 ===")
print("n glass mean: ", n_glass_mean, " Error: ",err_n_glass_mean)

