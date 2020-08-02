import numpy as np
import scipy as sp 
import math
import cmath
import scipy.spatial.distance as dt
import scipy.sparse.linalg as la
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from mpl_toolkits.mplot3d import Axes3D
import time
import os

def Get_material(stwl,edwl,step,TYPE,DIEL=None,BUILTIN=None,PATH=None,UNIT="nm"):
    """Calculate material dielectric function (epsilon) basing on wavelength.
    #Input
    # --stwl,edwl,step                                      float
    #   Start wavelength, end wavelength, wavelength step.
    # --TYPE                                                char
    #   Type of material data source. 
    #         "DUMMY": constant complex dielectric function from input.
    #         "BUILTIN": material data from code library. 
    #         "FROMFILE": material data from data file.(Format: wavelength epsilon_real, epsilon_imag)
    # --DIEL                                                list of complex number
    #   Complex epsilon for DUMMY method.
    # --BUILTIN                                             list of char
    #   Name of built-in material data.
    #         Supported material: none.
    # --PATH                                                list of char
    #   Path to the target material file.
    # --UNIT                                                char
    #   Unit of wavelength, used for built-in material type.
    #Output:
    #  --wl                                                 list of float
    #    Target wavelengths.
    #  --material                                           list of ndarray         complex
    #    List of material used. Wavelength and every array in material should have one-to-one correspondance.
    """
    wl = np.arange(stwl,edwl+0.1*step,step,dtype=float)    
    if TYPE == "DUMMY":
        material = []
        for i0 in range(len(DIEL)):
            material.append(DIEL[i0]*np.ones(len(wl),dtype=complex))
        return wl.tolist(), material
    elif TYPE == "BUILTIN":
        unit_dic = {
            "nm"    : 1E9,
            "um"    : 1E6,
            "m"     : 1
        }
        diel_dic = {
            "Ag"    : "./diel/Ag (Silver) - CRC (raw)",
            "Al"    : "./diel/Al (Aluminium) - Palik (raw)",
            "Au"    : "./diel/Au (Gold) - Johnson and Christy (raw)",
            "Si"    : "./diel/Si (Silicon) - Palik (raw)",
            "SiO2"  : "./diel/SiO2 (Glass) - Palik (raw)",
            "Air"   : "./diel/Air"
        }
        unit_factor = unit_dic.get(UNIT)
        if unit_factor == None:
            raise Exception("In PreProcessing.py/Get_Material: Unit {} not supported.".format(UNIT))
        material = []
        for i0 in range(len(BUILTIN)):
            Dirname = os.path.dirname(__file__)
            Filename = diel_dic.get(BUILTIN[i0])
            if Filename == None:
                raise Exception("In PreProcessing.py/Get_Material: {} not supported in builtin material .".format(BUILTIN[i0]))            
            Im_eps = np.loadtxt(Filename+"Im_eps.txt")
            Re_eps = np.loadtxt(Filename+"Re_eps.txt")
            Im_eps = Im_eps[Im_eps[:,0].argsort()]
            Re_eps = Re_eps[Re_eps[:,0].argsort()]
            if stwl<Im_eps[0,0]*unit_factor or edwl<Im_eps[0,0]*unit_factor or \
               stwl>Im_eps[-1,0]*unit_factor or edwl>Im_eps[-1,0]*unit_factor:
                raise Exception("In PreProcessing.py/Get_Material: Target wavelength out of bound.")
            Im_interp = np.interp(wl,Im_eps[:,0]*unit_factor,Im_eps[:,1])
            Re_interp = np.interp(wl,Re_eps[:,0]*unit_factor,Re_eps[:,1])
            material.append(Re_interp+Im_interp*1j)
        return wl.tolist(), material
    else:
        raise Exception("Not supported material type. Please choose one of following:\n\
                        \tDUMMY,\tBUILTIN,\tFROMFILE.")

def solve_h(n, f, H, r):
    a = n
    b = -2*(f+(n-1)*H)
    c = (n-2)*H**2+2*H*f-(r**2)/(n-1)
    
    delta = b**2-4*a*c
    root1 = np.real((-b - cmath.sqrt(delta))/(2*a))           #root1 is the correct one
    #root2 = (-b + cmath.sqrt(delta))/(2*a)
    return root1
    


wavelength = 500
epsilon = Get_material(wavelength, wavelength, 10, "BUILTIN", BUILTIN=["SiO2"])[1][0][0]
epsilon1 = np.real(epsilon)
epsilon2 = np.imag(epsilon)

n2 = math.sqrt((abs(epsilon1)+epsilon1)/2)        #actual refractive index of the high index material,  the output would be normalized between the low and high value and would be between 0~1
n2 = np.real(n2)                                  #now use only the real part for dielectric material


focus = 1600                                     #distance from the bottom (including thickness)
#thickness = 275
d = 25         #discrete interval
n1 = 1         #actual refractive index of the low index material
         
ncenter = n2  #actual refractive index in the middle of the lens


# For fixed n len and n1 must be 1
Dmax = 1975
#Dmax = 2*math.sqrt((2*n2-2)*focus*thickness + (n2**2-3*n2+2)*(thickness**2))
print("Diameter of the len is: {}".format(Dmax))
thickness = ((2-2*n2)*focus+math.sqrt(((2*n2-2)**2)*(focus**2)+4*(n2**2-3*n2+2)*(Dmax/2)**2))/(2*(n2**2-3*n2+2))
print("thickness of the len is {}.".format(thickness))
print("focus of the len(from the buttom) is {}".format(focus))

D = Dmax
xscale = int(np.ceil(D/d) + 1)
yscale = xscale
zscale = int(np.ceil(thickness/d) + 1)
zscale_all = int(np.ceil(focus/d) + 1)

X, Y, Z = np.meshgrid(range(xscale), range(yscale), range(zscale_all), indexing = 'ij')
X = np.reshape(X, (X.size, 1), order='C')
Y = np.reshape(Y, (Y.size, 1), order='C')
Z = np.reshape(Z, (Z.size, 1), order='C')
Geometry = np.concatenate((X, Y, Z), axis = 1)
Geometry = np.reshape(Geometry, (Geometry.size, 1), order='C')

numincircle = 0
Geometry_circle = []
for i in range(int(round(Geometry.size/3))):                                         ##Get a circular region
    r = math.sqrt((Geometry[3*i]*d - D/2)**2 + (Geometry[3*i+1]*d - D/2)**2)
    if(r<=Dmax/2):
        numincircle = numincircle + 1
        Geometry_circle.append(Geometry[3*i])
        Geometry_circle.append(Geometry[3*i+1])
        Geometry_circle.append(Geometry[3*i+2])

print("numincircle{}".format(numincircle))
Geometry_circle = np.array(Geometry_circle)
print(Geometry_circle.shape)
print("xmax{}".format(np.max([Geometry_circle[3*i] for i in range(int(round(Geometry_circle.size/3)))])))     
print("ymax{}".format(np.max([Geometry_circle[3*i+1] for i in range(int(round(Geometry_circle.size/3)))])))     
print("zmaxforwhole{}".format(np.max([Geometry_circle[3*i+2] for i in range(int(round(Geometry_circle.size/3)))])))    
print("xmin{}".format(np.min([Geometry_circle[3*i] for i in range(int(round(Geometry_circle.size/3)))])))     
print("ymin{}".format(np.min([Geometry_circle[3*i+1] for i in range(int(round(Geometry_circle.size/3)))])))     
print("zminforwhole{}".format(np.min([Geometry_circle[3*i+2] for i in range(int(round(Geometry_circle.size/3)))])))    
print("zmaxforwhole-ref{}".format(zscale_all))   
print("zmaxforlens{}".format(zscale-1))     

Diel = np.zeros(Geometry_circle.shape)
for i in range(int(round(Geometry_circle.size/3))):
    r = math.sqrt((Geometry_circle[3*i]*d - D/2)**2 + (Geometry_circle[3*i+1]*d - D/2)**2)
    h = solve_h(n2, focus, thickness, r)
    if Geometry_circle[3*i+2]*d >= (thickness-h)/2 and Geometry_circle[3*i+2]*d <= (thickness+h)/2:
        Diel[3*i] = n2
    else:
        Diel[3*i] = n1
    Diel[3*i+1] = Diel[3*i]
    Diel[3*i+2] = Diel[3*i]

Diel = (Diel - n1)/(n2-n1)

Geometry_circle_str = []
Diel_str = []
for i in range(int(round(Geometry_circle.size/3))):
    if(Geometry_circle[3*i+2]*d <= thickness):
        Geometry_circle_str.append(Geometry_circle[3*i])
        Geometry_circle_str.append(Geometry_circle[3*i+1])
        Geometry_circle_str.append(Geometry_circle[3*i+2])
        Diel_str.append(Diel[3*i])
        Diel_str.append(Diel[3*i+1])
        Diel_str.append(Diel[3*i+2])
Geometry_circle_str = np.array(Geometry_circle_str)
Diel_str = np.array(Diel_str)
print("------------for str only---------------:")
print("xmax{}".format(np.max([Geometry_circle_str[3*i] for i in range(int(round(Geometry_circle_str.size/3)))])))     
print("ymax{}".format(np.max([Geometry_circle_str[3*i+1] for i in range(int(round(Geometry_circle_str.size/3)))])))     
print("zmaxforwhole{}".format(np.max([Geometry_circle_str[3*i+2] for i in range(int(round(Geometry_circle_str.size/3)))])))    
print("xmin{}".format(np.min([Geometry_circle_str[3*i] for i in range(int(round(Geometry_circle_str.size/3)))])))     
print("ymin{}".format(np.min([Geometry_circle_str[3*i+1] for i in range(int(round(Geometry_circle_str.size/3)))])))     
print("zminforwhole{}".format(np.min([Geometry_circle_str[3*i+2] for i in range(int(round(Geometry_circle_str.size/3)))])))    


np.savetxt('LensDiel.txt', Diel_str, '%f')
np.savetxt('LensGeometry.txt', Geometry_circle_str, '%i')

Diel1num = 0
Diel0num = 0
Geo1 = []
Geo0 = []
for i in range(int(round(Geometry_circle.size/3))):
    if(Geometry_circle[3*i+2]*d <= thickness):
        Diel1num += 0
        Geo1.append(Geometry_circle[3*i])
        Geo1.append(Geometry_circle[3*i+1])
        Geo1.append(Geometry_circle[3*i+2])
    else:
        Diel0num += 1
        Geo0.append(Geometry_circle[3*i])
        Geo0.append(Geometry_circle[3*i+1])
        Geo0.append(Geometry_circle[3*i+2])
    
Geo1 = np.array(Geo1)
Geo0 = np.array(Geo0)

Diel1 = np.ones(Geo1.shape)                                       #material with diel=1 (relatively)
Diel0 = np.zeros(Geo0.shape)                                      #with diel=0 (relatively)

np.savetxt('1Diel.txt', Diel1, '%f')
np.savetxt('1Geometry.txt', Geo1, '%i')
np.savetxt('0Diel.txt', Diel0, '%f')
np.savetxt('0Geometry.txt', Geo0, '%i')






















'''
# For varying n len

Dmax = 2*math.sqrt(((ncenter - n1)*wavelength*thickness/(2*math.pi) + focus)**2-focus**2)
print("Diameter of the len is: {}".format(Dmax))

D = Dmax
xscale = int(np.ceil(D/d) + 1)
yscale = xscale
zscale = int(np.ceil(thickness/d) + 1)

X, Y, Z = np.meshgrid(range(xscale), range(yscale), range(zscale), indexing = 'ij')
X = np.reshape(X, (X.size, 1), order='C')
Y = np.reshape(Y, (Y.size, 1), order='C')
Z = np.reshape(Z, (Z.size, 1), order='C')
Geometry = np.concatenate((X, Y, Z), axis = 1)
Geometry =np.reshape(Geometry, (Geometry.size, 1), order='C')

Diel = np.zeros(Geometry.shape)
for i in range(int(round(Geometry.size/3))):
    r = math.sqrt((Geometry[3*i]*d - D/2)**2 + (Geometry[3*i+1]*d - D/2)**2)
    Diel[3*i] =  ncenter - (2*math.pi/(wavelength * thickness)) * (math.sqrt(focus**2 + r**2) - focus)
    if Diel[3*i] > n2:
        Diel[3*i] = n2
    if Diel[3*i] < n1:
        Diel[3*i] = n1
    Diel[3*i+1] = Diel[3*i]
    Diel[3*i+2] = Diel[3*i]

Diel = (Diel - n1)/(n2-n1)


np.savetxt('LensDiel.txt', Diel, '%f')
np.savetxt('LensGeometry.txt', Geometry, '%i')
'''