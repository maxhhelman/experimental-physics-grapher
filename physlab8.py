import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math
def cleanNaN(data):
    x = [x for x in data if str(x) != 'nan']
    return x

def func(x, a, b, c):
    return (a*x) + b

def bestFit(x, y, y_error):
    x0 = np.array([0.0, 0.0, 0.0])
    w, _ = opt.curve_fit(func, x, y, x0, y_error) 
    return w


def grapLog(x, y, mode):
    k = y
    y = np.log(y)   
    #sigma x over x
    y_error = [0.5/p for p in k]
    x_error = [0.5/p for p in x]
    plt.errorbar(x, y, yerr=y_error, xerr = x_error, marker='x', ecolor='black', mec='red', linestyle='None',ms=4, mew=4, label="Data")
    w = bestFit(x, y, y_error)
    x = np.asarray(x)
    w = np.asarray(w)
    yfit = func(x ,w[0], w[1], w[2])
    plt.plot(x,yfit)
    plt.xlabel('Time (s)')
    plt.ylabel('Log(Current (Micro Amperes))')
    plt.title("Logarithmic Graph of " + mode + " muF Capacitor")
    plt.show() 

def plotNormal(x,y,mode):
    plt.plot(x,y)
    plt.xlabel('Time (s)')
    plt.ylabel('Current (Micro Amperes)')
    plt.title("Graph of "+mode+ " muF Capacitor")
    plt.show()

def readData(filename, mode):
    data = pd.read_csv(filename)
    dataIn = np.transpose(data.values)

    time10 = cleanNaN(dataIn[0])
    I10 = cleanNaN(dataIn[1])
    plotNormal(time10, I10, mode + ' 10')
    grapLog(time10, I10, mode+ ' 10')

    time20 = cleanNaN(dataIn[2])
    I20 = cleanNaN(dataIn[3])
    plotNormal(time20, I20, mode+ ' 20')
    grapLog(time20, I20, mode+ ' 20')

    time30 = cleanNaN(dataIn[4])
    I30 = cleanNaN(dataIn[5])
    plotNormal(time20, I20, mode+ ' 30')
    grapLog(time30, I30, mode + ' 30')

readData('data_charge.csv', 'charge')
readData('data_discharge.csv', 'discharge')
