import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math


err_angle = 1/60
pi = 3.14159265359

def clearNan(items):
	return [x for x in items if x == x]

def convertToDegrees(degrees, minutes):
	minutes = [(1/60)*x for x in minutes]
	degrees = [x+y for x,y in zip(degrees,minutes)]
	return degrees

def avg(a,b):
	assert len(a) == len(b)
	sums = [x-y for x,y in zip(a,b)]
	sums = [x/2 for x in sums]
	return sums

def computeNi(R, nf, m, latConst, theta_rad):
	ni = math.sqrt( abs((R)/((R/4) - (m/(latConst*math.sin(theta_rad))))) )
	return ni


def errorNi(R,d,m,theta_rad, err_d, err_th):
	partial_theta_squared = ( 16*R*(m**2)*d*(1/math.tan(theta_rad)**2)*math.sin(theta_rad) )/( (R*d*math.sin(theta_rad) -4*m)**3)
	partial_d_squared = ( 16*R*(m**2)*math.sin(theta_rad) )/(d*(d*R*math.sin(theta_rad) - 4*m)**3)
	err = math.sqrt( partial_d_squared*(err_d**2) + partial_theta_squared*(err_th**2)) 
	return err

#Lattice constant
m = 1
left_yellow = 159 +  33/60
right_yellow = 201 + 22/60
avgDegLC = (left_yellow - right_yellow)/2
rads = np.radians(avgDegLC)
wavelength = 5.8756 * 10**(-7)
latticeConst = abs(m*wavelength/(math.sin(rads))) 

error_in_degrees_when_we_average = math.sqrt( (err_angle**2)*( ((1/2)**2) + ((1/2)**2) ) )
err_rad_average_LC = math.sqrt( ((pi/180)**2)*((error_in_degrees_when_we_average)**2) )
err_d = math.sqrt( (((m*wavelength*math.cos(rads))/(math.sin(rads)**2))**2)*(err_rad_average_LC**2) )


print("=== Lattice Constant d  (meters) ===")
print(latticeConst, " Error: ", err_d)


dataRAW = pd.read_csv('data_p2.csv')
data = np.transpose(dataRAW.values)

#Dont need to read this
Left1Order_col = data[0]
Left1Order_deg = convertToDegrees(data[1], data[2])
Right1Order_col = data[3]
Right1Order_deg = convertToDegrees(data[4], data[5])
Left2Order_col = clearNan(data[6])
Left2Order_deg = convertToDegrees(clearNan(data[7]),clearNan(data[8]))
Right2Order_col = clearNan(data[9])
Right2Order_deg = convertToDegrees(clearNan(data[10]), clearNan(data[11]))


Order1_col = Left1Order_col
Order1_avg_deg = avg(Right1Order_deg, Left1Order_deg)
Order2_col = Left2Order_col
Order2_avg_deg = avg(Right2Order_deg, Left2Order_deg)

#Average Degrees
print("=== Averages ===")
print("Order 1")
print(Order1_col)
print(Order1_avg_deg, " Errors: ", error_in_degrees_when_we_average)

print("Order 2")
print(Order2_col)
print(Order2_avg_deg, " Errors: ", error_in_degrees_when_we_average)


#ni
R = 1.0974 * 10**(7)
ni_vals_ord1 = [computeNi(R, 2, 1, latticeConst, np.radians(x)) for x in Order1_avg_deg]
ni_vals_ord2 = [computeNi(R, 2, 2, latticeConst, np.radians(x)) for x in Order2_avg_deg]
err_ni_o1 = [errorNi(R, latticeConst, 1, np.radians(x), err_d, err_rad_average_LC) for x in Order1_avg_deg]
err_ni_o2 = [errorNi(R, latticeConst, 2, np.radians(x), err_d, err_rad_average_LC) for x in Order2_avg_deg]



print("=== n_i ===")
print("Order 1")
print(Order1_col)
print(ni_vals_ord1, " Errors: ", err_ni_o1)

print("Order 2")
print(Order2_col)
print(ni_vals_ord2, " Errors: ", err_ni_o2)
