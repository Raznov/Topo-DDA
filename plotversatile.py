import numpy as np
import scipy as sp 
import math
import scipy.spatial.distance as dt
import scipy.sparse.linalg as la
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from mpl_toolkits.mplot3d import Axes3D
import time



time1 = np.genfromtxt("DDATime1.txt")
time6 = np.genfromtxt("DDATime6.txt")

fig = plt.figure()
plt.plot(time1[:,0], time1[:,1])
plt.plot(time6[:,0], time6[:,1])
plt.ylabel("Time/s")
plt.xlabel("Number of dipoles/s")
plt.legend(["Threads = 1", "Threads = 6"])
plt.show()



