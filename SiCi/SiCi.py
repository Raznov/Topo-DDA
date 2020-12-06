import numpy as np
import scipy as sp
import scipy.special
import math
import matplotlib
import matplotlib.pyplot as plt



number = 1000000
dis = 0.1
gamma = 0.577216

pos = np.linspace(dis, dis*number, number)
pos_plot = np.linspace(dis, dis*(number-1), number-1)
Si_plot, Ci_plot = sp.special.sici(pos_plot)

fSi = open("Si.txt", "w")
fCi = open("Ci.txt", "w")
fSi.write(str(number-1)+"\n"+str(dis)+"\n")
fCi.write(str(number-1)+"\n"+str(dis)+"\n")
for i in range(number-1):
    fSi.write(str(Si_plot[i])+"\n")
    fCi.write(str(Ci_plot[i])+"\n")
fSi.close()
fCi.close()


plt.plot(pos_plot,Si_plot)
plt.plot(pos_plot,Ci_plot)
plt.show()
