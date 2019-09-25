import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import math
import sys

def graph(filename, color):
    data = pd.read_csv(filename)
    mass = np.transpose(data.values)[0]
    littleI = np.transpose(data.values)[1]
    print("Mass: {}".format(mass))
    print("Little I: {}".format(littleI))
    plt.plot(littleI, mass, color)
graph('lab3_5.csv', '.r-')
graph('lab3_4.csv', '.b-')
graph('lab3_3.csv', '.g-')
graph('lab3_2.csv', '.y-')


plt.xlabel('i')
plt.ylabel('m')
plt.title('Measurements')
plt.show()

