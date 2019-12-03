import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math

def func(x, a, b,c):
    return (a*x) + b

def bestFit(x, y, yerr):
    x0 = np.array([0.0, 0.0, 0.0])
    w, er = opt.curve_fit(func, x, y, x0, sigma=yerr,  bounds=(-200, [3000., 2000., 1000.])) 
    return w, er

def func_const(x, a, b,c):
    return a

def bestFit_const(x, y, yerr):
	x0 = np.array([0.0, 0.0, 0.0])
	w, er = opt.curve_fit(func_const, x, y, x0, sigma=yerr,  bounds=(-200, [3000., 2000., 1000.])) 
	return w, er

def backgroundPlot(x,y, x_error, y_error):
	plt.errorbar(x, y, yerr=y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x, y, '-')
	plt.xlabel(r"Voltage V", fontsize=12) 
	plt.ylabel(r"Particle Count", rotation=90, fontsize=12)
	plt.title(r"Background Measurements")
	plt.show()	

def betaPlot(x,y, x_error, y_error):
	w, err = bestFit_const(x[-4:],y[-4:],y_error[-4:])
	A = w[0]
	errorArr = np.sqrt(np.diag(err))
	Aerr = errorArr[0]
	y_fit = [A] * len(x)

	plt.errorbar(x, y, yerr=y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x, y, '-')
	plt.plot(x, y_fit, ':')
	plt.xlabel(r"Thickness cm", fontsize=12) 
	plt.ylabel(r"Particle Count (log scale)", rotation=90, fontsize=12)
	plt.title(r"$\beta$ Particle Count")
	plt.show()

def gammaPlot(x,y, x_error, y_error):
	w, err = bestFit(x,y,y_error)
	A = w[0]
	B = w[1]
	errorArr = np.sqrt(np.diag(err))
	Aerr = errorArr[0]
	y_fit = [A*p for p in x] + B
	plt.errorbar(x, y, yerr=y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x, y, '-')
	plt.plot(x, y_fit, ':')
	plt.xlabel(r"Thickness cm", fontsize=12) 
	plt.ylabel(r"Particle Count (log scale)", rotation=90, fontsize=12)
	plt.title(r"$\gamma$ Particle Count")
	plt.show()
	return A, Aerr

#Background Measurements
print("=== Background ===\n Just graph") 
backgroundRAW = pd.read_csv('background.csv')
background_data = np.transpose(backgroundRAW.values)
bg_voltage = background_data[0] 
bg_counts = background_data[1]
err_bg_voltage = [5] * len(bg_voltage)
err_bg_count = [math.sqrt(x) for x in bg_counts]
backgroundPlot(bg_voltage, bg_counts, err_bg_voltage, err_bg_count)

#Beta Rays
energy = 0.765 #MeV
desnsityAl = 2.702 #g cm^-3
rng_exp = (0.412/desnsityAl)*(energy**(1.29)) #cm
betaRAW = pd.read_csv('beta.csv')
beta_data = np.transpose(betaRAW.values)
beta_thickness = 2*beta_data[0] #milis
beta_count = beta_data[1] #count
beta_log_count = [math.log(x) for x in beta_count]
beta_cm_thickness = [x*2.54*(10**(-3)) for x in beta_thickness]
err_beta_count = [math.sqrt(x)/x for x in beta_count] #1/x * err_x
err_beta_thickness = [0] * len(beta_thickness)
print("=== Beta ===")
print("Max Range estimated by taking the last 4 (~28.5%) datapoints and doing a best fit.")

betaPlot(beta_cm_thickness, beta_log_count, err_beta_thickness, err_beta_count)
print("Range Obtained: ", 0.05, " cm, Obtained by observation" )
print("Range Expected: ", rng_exp)

print("=== Gamma ===")
background_counts = 19;
plate_thickness_cm = 0.15748
gammaRAW = pd.read_csv('gamma.csv')
gamma_data = np.transpose(gammaRAW.values)
gamma_plates = gamma_data[0] #milis
gamma_count = gamma_data[1] #count
gamma_log_count = [math.log(x-background_counts) for x in gamma_count]
gamma_cm_thickness = [x*plate_thickness_cm for x in gamma_plates]
err_gamma_count = [math.sqrt(x)/x for x in gamma_count] #1/x * err_x
err_gamma_thickness = [0] * len(gamma_plates)
mu, mu_err = gammaPlot(gamma_cm_thickness, gamma_log_count, err_gamma_thickness, err_gamma_count)
print("Absorption coeff: ", -1*mu, " Error: ", mu_err)
print("Ray Energy: ", 0.94, " MeV, Falls in accpeted range (10^4 to 10^7 eV)")


