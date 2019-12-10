import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math

def e_bar_and_err(e):
	N = len(e)
	e_bar = sum(e)/N
	s = math.sqrt(sum( (e-e_bar)**2 )/(N - 1))
	e_bar_err = (s)/math.sqrt(N)
	return e_bar, e_bar_err

def e_weighted_and_err(e, vf, vi, err_vf, err_vi):
	#e = vf/vi 
	#find errors
	partial_vi = np.asarray(-vf/((vi)**2)) 
	partial_vf = np.asarray([1/x for x in vi])
	part1 = (partial_vf**2)*(err_vf**2)
	part2 = (partial_vi**2)*(err_vi**2)
	err_e_sq = part1+part2
	err_e = np.asarray([math.sqrt(x) for x in err_e_sq])

	#find weighted e
	e_w_bar = sum( [e/(sg**2) for e,sg in zip(e,err_e)] )/sum( [1/(sg**2) for sg in err_e] )
	err_e_w_bar = 1/(math.sqrt( sum( [1/(sg**2) for sg in err_e] )))
	return e_w_bar, err_e_w_bar, err_e

def func(x, a, b,c):
    return (a*x) + b

def bestFit(x, y, yerr):
    x0 = np.array([0.0, 0.0, 0.0])
    w, er = opt.curve_fit(func, x, y, x0, sigma=yerr,  bounds=(-200, [3000., 2000., 1000.])) 
    return w, er

def part1Plot1(x,y,sigmaX,sigmaY):
	plt.errorbar(x, y, yerr=sigmaY, xerr = sigmaX, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.xlabel(r"$v_i$ m/s", fontsize=17) 
	plt.ylabel(r"e", rotation=0, fontsize=17)
	plt.title(r"Plot of e and $v_i$")
	plt.show()

def part1Plot2(x,y,sigmaX, sigmaY, x_2, y_2, sigmaX_2, sigmaY_2):
	#6 is e_bar, 7 is e_w_bar
	plt.errorbar(x, y, xerr = sigmaX, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.errorbar(x_2, y_2, xerr = sigmaX_2, marker='^', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")	
	plt.xlabel(r"e", fontsize=17) 
	plt.ylabel(r"Measurement", rotation=90, fontsize=17)
	plt.title(r"Plot of the values of e")
	plt.show()

def part2Plot(x,y, x_error, y_error):
	w, err = bestFit(x,y,y_error)
	A = w[0]
	B = w[1]
	errorArr = np.sqrt(np.diag(err))
	Aerr = errorArr[0]
	Berr = errorArr[1]
	y_fit =  ([A*p for p in x]) + B
	plt.errorbar(x, y, yerr=y_error, xerr=x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x, y_fit, '-')
	plt.xlabel(r"Height m", fontsize=12) 
	plt.ylabel(r"$\bar{a_{x}} \ m/s^{2}$", rotation=90, fontsize=12)
	plt.title(r"Plot of $\bar{a_{x}}$ against the added height")
	plt.show()
	return A, B, Aerr, Berr


#Part 1
print("=== PART 1 ===")
part1_data_raw = pd.read_csv("part1.csv")
part1_data = np.transpose(part1_data_raw.values)

v_i = part1_data[0]
sigma_v_i = part1_data[1]
v_f = part1_data[2]
sigma_v_f = part1_data[3]

# unweighted avg
e = abs(v_f)/abs(v_i)
e_bar, err_e_bar = e_bar_and_err(e)
#weighted avg
#get sigma e from partial derivative formula
e_w_bar, err_e_w_bar, err_e = e_weighted_and_err(e, v_f, v_i, sigma_v_f, sigma_v_i)
print("Unweighted: e_bar = ", e_bar, " Error: ", err_e_bar)
print("Weighted: e_w_bar = ", e_w_bar, " Error: ", err_e_w_bar)

#e and v_i plot
part1Plot1(v_i,e,sigma_v_i,err_e)

#construct vectors for second plot
indices = [i for i in range(1, len(e)+1)]
indices_2 = [i for i in range(len(e)+1, len(e)+3)]
points_2 = [e_bar, e_w_bar]
err_pts_2 = [err_e_bar, err_e_w_bar]

#plot spread of points around e_bar and e_w_bar
part1Plot2(e, indices, err_e, None, points_2, indices_2, err_pts_2, None)


#Part 2
print("=== PART 2 ===")
part2_data_raw = pd.read_csv("part2.csv")
part2_data = np.transpose(part2_data_raw.values)
shim_th = 1.24 * 10**(-3)
err_shim_th = 0.05 * 10**(-3)
shims = part2_data[0] 
a_x = part2_data[1]
sigma_a_x = part2_data[2]
l_2 = part2_data[3]
v_1 = part2_data[4]
sigma_v_1 = part2_data[5]
coeff_a = part2_data[6]
sigma_coeff_a = part2_data[7]

#Get avg for each h
#so we have the errors in a but they want unweighted
# aka the formula from the manual
# i.e we disregard the ones we measured
meas = 3
L = 1 #m
heights = sorted(set(shims))
a_x_bar = []
err_a_x_bar = []
for i in heights:
	start = (int)((i-1)*meas)
	end = (int)(i*meas)
	mean, err_mean = e_bar_and_err( a_x[ start:end ]  )
	a_x_bar.append(mean)
	err_a_x_bar.append(err_mean)

m, m_err, b, b_err = part2Plot([i*shim_th for i in heights], a_x_bar, [i*err_shim_th for i in heights], err_a_x_bar)
 
print("Slope: ",m, " Error: ", m_err, "\nB intercept: ", b, " Error: ", b_err)
print("g estimate: ", m, " Error: ", m_err) #L is 1m so doesnt really have an effect

#Part 3
for i in heights:
	print("HEIGHT: ", (int)(i))
	start = (int)((i-1)*meas)
	end = (int)(i*meas)
	ax = a_x[ start:end ]
	l = l_2[ start:end ]
	v1 = v_1[ start:end ]
	d = (v1**2)/(2*ax*l) - 1
	print(d)
