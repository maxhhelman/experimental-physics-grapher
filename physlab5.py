import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math

#Part 1 Error Values
intensityValues_err = 0.1
degrees_err = 1
xm_err = 0.05
D_err = 0.5 
D = 110.3
d = 0.250
#Helper Methods
def func(x, a, b,c):
    return (a*x) + b

def bestFit(x, y, yerr):
    x0 = np.array([0.0, 0.0, 0.0])
    w, er = opt.curve_fit(func, x, y, x0, sigma=yerr,  bounds=(-200, [3000., 2000., 1000.])) 
    return w, er

def sigmaAandD(x, sigma):
	xSquared = [p**2 for p in x]
	D = len(x)*(sum(xSquared)) - (sum(x))**2
	sigmaA = sigma*math.sqrt((len(x)/D))
	sigmaD = sigma*math.sqrt(sum(xSquared)/D)
	return sigmaA, sigmaD

def part1Plot(x,y,sigmaX,sigmaY):
	#x - degrees, y - intensity
	x_error = [abs(2*(math.sin(i)*math.cos(i))*sigmaX) for i in x] #2sinXcosX
	y_error = [ math.sqrt(((1/y[0])**2)*(sigmaY**2) + (i/(y[0]**2))*(sigmaY**2)) for i in y] #partials formula
	x = [(math.cos(i))**2 for i in x] #cos^2
	y = [i/y[0] for i in y] #I/I0
	w, err = bestFit(x,y,y_error)
	A = w[0]
	B = w[1]
	errorArr = np.sqrt(np.diag(err))
	Aerr = errorArr[0]
	Berr = errorArr[1]
	y_fit =  ([A*p for p in x]) + B
	plt.errorbar(x, y, yerr=y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x, y_fit, '-')
	plt.xlabel(r"$cos(\theta)^{2}$", fontsize=17) 
	plt.ylabel(r"$\frac{I}{I_0}$", rotation=0, fontsize=17)
	plt.title("Relative Intensity versus Squared Cosine of Angle Difference")
	plt.show()
	return A, B, Aerr, Berr

def part2Plot(x,y, x_error, y_error):
	#x - m (order number),y - xm values
	w, err = bestFit(x,y,y_error)
	A = w[0]
	B = w[1]
	errorArr = np.sqrt(np.diag(err))
	Aerr = errorArr[0]
	Berr = errorArr[1]
	y_fit =  ([A*p for p in x]) + B
	plt.errorbar(x, y, yerr=y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x, y_fit, '-')
	plt.xlabel(r"$m$", fontsize=12) 
	plt.ylabel(r"$x_{m} (cm)$", rotation=0, fontsize=12)
	plt.title("Maxima Distances versus Peak Order")
	plt.show()
	return A, B, Aerr, Berr

def part3Plot(x, y, x_error, y_error):
	plt.errorbar(x, y, yerr=y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x, y, '-')
	plt.xlabel(r"$x_{m} (cm)$", fontsize=17) 
	plt.ylabel(r"$\frac{I}{I_0}$", rotation=0, fontsize=17)
	plt.title("Relative Intensity versus Peak Distance")
	plt.show()

def part3Plot2(x,y,x_error, y_error):
	w, err = bestFit(x,y,y_error)
	A = w[0]
	B = w[1]
	errorArr = np.sqrt(np.diag(err))
	Aerr = errorArr[0]
	Berr = errorArr[1]
	y_fit =  ([A*p for p in x]) + B
	plt.errorbar(x, y, yerr  = y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x, y_fit, '-')
	plt.xlabel(r"$n$", fontsize=17) 
	plt.ylabel(r"$x_{n} (cm)$", rotation=0, fontsize=10)
	plt.title("Single Slit Minima Distance and Order")
	plt.show()
	return A, B, Aerr, Berr

def findLambda(slope2, slope2_err, d,D, D_err):
	#slope and slope error are in cm, d in mm, D, D_err are in cm
	#we assume d has no error

	#convert to meters
	slope2 = slope2*(10**(-2))
	slope2_err = slope2_err*(10**(-2))
	d = d*(10**(-3))
	D = D*(10**(-2))
	D_err = D_err*(10**(-2))

	lambdaVal = (slope2*d)/D #10^-7
	lambda_err = math.sqrt( ((d/D)**2)*(slope2_err**2) + ((slope2*d/D)**2)*(D_err**2) ) #10^-8

	return lambdaVal, lambda_err

def findSmallA(slope3, slope3_err, lambdaVal, lambda_err, D, D_err):
	#slope and D in cm, lambda is in m
	slope3 = slope3*(10**(-2))
	slope3_err = slope3_err*(10**(-2))
	D = D*(10**(-2))
	D_err = D_err*(10**(-2))

	#just to make it faster to type
	s = slope3
	s_err = slope3_err
	l = lambdaVal
	l_err = lambda_err

	a = (lambdaVal*D/slope3)
	a_err = math.sqrt(  ((l/s)**2)*(D_err**2) + ((D/s)**2)*(l_err**2) + (((l*d)/(s**2))**2)*(s_err**2) )
	return a, a_err

#Get Intensities and Degrees
print("==== PART 1 ====")
intensities = pd.read_csv("intensities.csv")
intensities = np.transpose(intensities.values)

degreeDifference = intensities[0]
intensityValues = intensities[1]
radiansDifference = [math.radians(i) for i in degreeDifference]
radian_err = math.radians(degrees_err)

slope, interc, slope_err, interc_err = part1Plot(radiansDifference, intensityValues, radian_err, intensityValues_err)
print("Slope: ",slope," Error: ",slope_err)
print("Intercept: ",interc, " Error: ",interc_err)

#Get Maxima Values
print("==== PART 2 ====")
maxima = pd.read_csv("maxima.csv")
maxima = np.transpose(maxima.values)
xm = maxima[1]
#this is an example of the deadly sin of hardcoding
m = [-4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
intensity = maxima[2]
intensity_err= maxima[4]
realIntensities = maxima[5]

slope2, interc2, slope2_err, interc2_err = part2Plot(m, xm, 0, [xm_err for i in xm])
lambdaValue, lambda_err = findLambda(slope2, slope2_err, d,D, D_err)

print("Slope: ",slope2," Error: ",slope2_err)
print("Intercept: ",interc2, " Error: ",interc2_err)
print("Lambda Estimate (m): ", lambdaValue, " Error (m): ", lambda_err)

print("==== PART 3 ====")
relativeIntensityValuesAdjusted = [i/max(intensity) for i in intensity]
part3Plot(xm, relativeIntensityValuesAdjusted, [xm_err for i in xm], intensity_err)
print("Just the graph here")


print("==== PART 3.2 ====")
#estimates for xn
xn = [intensity[2], intensity[8]]
n = [-1, 1]
slope3, interc3, slope3_err, interc3_err = part3Plot2(n, xn, [0,0], [xm_err, xm_err])
slope3_err, interc3_err = sigmaAandD(xn, 0.05)
a, a_err = findSmallA(slope3, slope3_err, lambdaValue, lambda_err, D, D_err)

#ERRORS ARE BS, BC IT CANT CALCULATE THEM
print("Slope: ",slope3," Error: ", slope3_err) 
print("Intercept: ",interc3, " Error: ",interc3_err)
print("a estimate (m): ", a, " Error (m): ", a_err)



