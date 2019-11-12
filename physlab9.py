import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math
pi = 3.14159265359

def normalize(data):
	biggest = max(data)
	smallest = min(data)
	result = [(x - smallest)/(biggest - smallest) for x in data]
	return result

def getNormalizeError(data, normalized):
	data = data + 0.04
	biggest = max(data)
	smallest = min(data)
	result = [(x - smallest)/(biggest - smallest) for x in data]
	result = [(x-y) for x,y in zip(result, normalized)]
	return result

def getData(filename):
	raw = pd.read_csv(filename).sort_values(by=["Frequency"])
	data = np.transpose(raw.values)
	frequency = data[0] * 2 * pi
	voltage = data[1]
	voltage = normalize(voltage)
	voltageErr = getNormalizeError(data[1], voltage)
	return frequency, voltage, voltageErr


def plotData(x,y):
    # x_error = [*2*pi for p in x]
    # plt.errorbar(x, y, yerr=y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
	plt.plot(x,y)
	plt.xlabel('w (Hz)')
	plt.ylabel('Normalized voltage (V)')
	plt.title("Measurements")


w10, v10, err10 = getData('10ohm.csv')
w50, v50, err50 = getData('50ohm.csv')
w500, v500, err500 = getData('500ohm.csv')

print(err10)

plotData(w10/6.28, v10)
plotData(w50/6.28, v50)
plotData(w500/6.28, v500)

plt.show()

w10, v10, err10 = getData('10ohm.csv')
w50, v50, err50 = getData('50ohm.csv')
w500, v500, err500 = getData('500ohm.csv')
bcw, bcv, bcerr = getData('bigCoil.csv')

# plotData(w10, v10)
# plotData(w50, v50)
# plotData(w500, v500)
plotData(bcw, bcv)