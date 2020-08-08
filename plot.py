import numpy as np
import scipy as sp 
import math
import sys
import os
import scipy.spatial.distance as dt
import scipy.sparse.linalg as la
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from mpl_toolkits.mplot3d import Axes3D
import time

#matplotlib.use('Agg')

def Shape(geometry,diel,d,iteration=-1,position="./",decimal=0,FullLattice=False):
    """Plot the shape of object as dot matrix.
    #Input:
    # --SC                                                         SolutionClass
    #   Solution Class.
    # --FullLattice   Boolean
    #   If true. Geometry is a full n*n*n matrix. diel=0 for point seen as air.
    """
    #d = SC.Get_d()
    N=round(np.shape(geometry)[0]/3)
    if(N!=np.shape(diel)[0]/3):
        print("size not equal!")
    geometry=np.reshape(geometry,(N,3))
    #geometry = SC.Get_geometry()
    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    [X,Y,Z] = map(int,list(np.amax(geometry,axis=0)+1))
    Axis_max = max(X,Y,Z)*1.2*d

    diel=np.reshape(diel,(N,3))
    diel=diel[:,0]
    #diel = SC.Get_diel()[:,0]


    cmaparg = 'Spectral_r'
    minn, maxx = 0, 1
    norm = matplotlib.colors.Normalize(minn, maxx)
    colorset = cm.ScalarMappable(norm=norm, cmap=cmaparg)
    colorset.set_array([])
    if FullLattice:
        index = np.where(diel>0.5)
        diel = diel[index]
        geometry = geometry[index]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    geo_dic = set()
    surf_list = []
    surf_color = []
    x_grid, y_grid, z_grid = (np.indices((X+1,Y+1,Z+1))-0.5)*d
    filled, colors = np.zeros((X,Y,Z)), np.zeros((X,Y,Z))
    for i,pos in enumerate(geometry):
        geo_dic.add((pos[0],pos[1],pos[2]))
    for i,pos in enumerate(geometry):
        if (pos[0]+1,pos[1],pos[2]) not in geo_dic or (pos[0]-1,pos[1],pos[2]) not in geo_dic or\
           (pos[0],pos[1]+1,pos[2]) not in geo_dic or (pos[0],pos[1]-1,pos[2]) not in geo_dic or\
           (pos[0],pos[1],pos[2]+1) not in geo_dic or (pos[0],pos[1],pos[2]-1) not in geo_dic:
            filled[pos[0],pos[1],pos[2]] = 1
            colors[pos[0],pos[1],pos[2]] = diel[i]
    surf_list = np.array(surf_list)
    surf_color = np.array(surf_color)
    # print(x_grid.shape,y_grid.shape,z_grid.shape,filled.shape,colors.shape)
    colors = cm.Spectral_r(norm(colors))
    ln=ax.voxels(x_grid,y_grid,z_grid,filled.astype(bool),facecolors=colors,edgecolor ='white',linewidth=0.2)
    #ax2.scatter(geometry[:,0]*d,geometry[:,1]*d,geometry[:,2]*d,c=E_tot_abs, s=15, cmap=cmaparg)
    ax.set_xlim(-(Axis_max-X*d)/2,(Axis_max+X*d)/2)
    ax.set_ylim(-(Axis_max-Y*d)/2,(Axis_max+Y*d)/2)
    ax.set_zlim(-(Axis_max-Z*d)/2,(Axis_max+Z*d)/2)
    ax.set_xlabel("x[nm]")
    ax.set_ylabel("y[nm]")
    ax.set_zlabel("z[nm]")
    ax.grid(False)
    fig.colorbar(colorset, shrink=0.9, aspect=10)
     
    #plt.show()
    #plt.savefig("./optimization_geometries/Iteration{}.png".format(it))
    if iteration==-1:
        plt.savefig("Space.png")
    else:
        fig.suptitle("iteration{}".format(iteration))
        plt.savefig(position+"{}Space.png".format(str(iteration).zfill(decimal))) 
    #plt.show()

def EField(geometry,diel,d,wl,k_dir,E_dir,E_tot,iteration=-1,position="./",decimal=0,FullLattice=False):
    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """
    #print(geometry)

    N=round(np.shape(geometry)[0]/3)
    #print(N)
    geometry=np.reshape(geometry,(N,3))
    
    print(geometry.shape)
    #diel=np.reshape(diel,(N,3))
    #diel=diel[:,0]
    

    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    
    #print(geometry)
    [X,Y,Z] = list(np.amax(geometry,axis=0)+1)
    Axis_max = max(X,Y,Z)*1.2*d
    
    E_tot = E_tot.reshape(int(E_tot.size/3),3)
    E_tot_abs = np.abs(np.sqrt(np.sum((E_tot)**2,axis=1)))
    
    """
    if FullLattice:
        index = np.where(diel!=0)                                      #disabled all codes with "diel"
        geometry = geometry[index]
        E_tot = E_tot[index]
        E_tot_abs = E_tot_abs[index]
    """
    cmaparg = 'Spectral_r'
    minn, maxx = E_tot_abs.min(), E_tot_abs.max()
    print(minn,maxx,np.argmax(E_tot_abs))
    norm = matplotlib.colors.Normalize(minn, maxx)
    colorset = cm.ScalarMappable(norm=norm, cmap=cmaparg)
    colorset.set_array([])
    
    
    
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.quiver(geometry[:,0]*d,geometry[:,1]*d,geometry[:,2]*d,np.real(E_tot[:,0]),np.real(E_tot[:,1]),np.real(E_tot[:,2]),
                length=10, lw=1)
    ax1.set_xlim(-(Axis_max-X*d)/2,(Axis_max+X*d)/2)
    ax1.set_ylim(-(Axis_max-Y*d)/2,(Axis_max+Y*d)/2)
    ax1.set_zlim(-(Axis_max-Z*d)/2,(Axis_max+Z*d)/2)
    ax1.set_xlabel("x[nm]")
    ax1.set_ylabel("y[nm]")
    ax1.set_zlabel("z[nm]")
    ax1.grid(False)
    fig1.suptitle("E field - Arrow plot\n {}nm, {}".format(wl,E_dir))
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    geo_dic = set()
    surf_list = []
    surf_color = []
    x_grid, y_grid, z_grid = (np.indices((X+1,Y+1,Z+1))-0.5)*d
    filled, colors = np.zeros((X,Y,Z)), np.zeros((X,Y,Z))
    for i,pos in enumerate(geometry):
        geo_dic.add((pos[0],pos[1],pos[2]))
    for i,pos in enumerate(geometry):
        if (pos[0]+1,pos[1],pos[2]) not in geo_dic or (pos[0]-1,pos[1],pos[2]) not in geo_dic or\
           (pos[0],pos[1]+1,pos[2]) not in geo_dic or (pos[0],pos[1]-1,pos[2]) not in geo_dic or\
           (pos[0],pos[1],pos[2]+1) not in geo_dic or (pos[0],pos[1],pos[2]-1) not in geo_dic:
            filled[pos[0],pos[1],pos[2]] = 1
            colors[pos[0],pos[1],pos[2]] = E_tot_abs[i]
    surf_list = np.array(surf_list)
    surf_color = np.array(surf_color)
   

    colors = cm.Spectral_r(norm(colors))
    ax2.voxels(x_grid,y_grid,z_grid,filled.astype(bool),facecolors=colors,linewidth=0.5)
    

    ax2.set_xlim(-(Axis_max-X*d)/2,(Axis_max+X*d)/2)
    ax2.set_ylim(-(Axis_max-Y*d)/2,(Axis_max+Y*d)/2)
    ax2.set_zlim(-(Axis_max-Z*d)/2,(Axis_max+Z*d)/2)
    ax2.set_xlabel("x[nm]")
    ax2.set_ylabel("y[nm]")
    ax2.set_zlabel("z[nm]")
    ax2.grid(False)
    fig2.colorbar(colorset, shrink=0.9, aspect=10)
    if iteration==-1:
        fig2.suptitle("E field - Scatter plot\n {}nm, {}".format(wl,k_dir)) 
        plt.savefig("E_field.png")
    else:
        fig2.suptitle("E field - Scatter plot\n {}nm, {} - iteration{}".format(wl,k_dir,iteration)) 
        plt.savefig(position+"{}E_field.png".format(str(iteration).zfill(decimal)))
    #plt.show()

"""
#Code used for single DDA calculation or if you only wants to see the final results for Evooptimization
data=np.genfromtxt("Model_results.txt",dtype=complex)
N=int(np.real(data[3]))
#print(np.real(data[(3*N+4):(6*N+4)]))
geometry=np.real(data[4:(3*N+4)]).astype(int)
polarization=data[(3*N+4):(6*N+4)]
diel=np.real(data[(6*N+4):(9*N+4)])
d=np.real(data[9*N+4])
wl=np.real(data[9*N+5])
k_dir=np.real(data[(9*N+6):(9*N+9)])
E_dir=np.real(data[(9*N+9):(9*N+12)])
#E_tot=(data[(6*N+12):(9*N+12)])
N_plot=int(np.real(data[(9*N+12)]))
geometry_plot=np.real(data[(9*N+13):(9*N+13+N_plot)]).astype(int)
E_tot=(data[(9*N+13+N_plot):(9*N+13+2*N_plot)])
Shape(geometry, diel, d)
EField(geometry_plot, diel, d, wl, k_dir, E_dir, E_tot)
"""



#Code used for taking a look at the geometry you built
                                   
data=np.genfromtxt("Space.txt",dtype=complex)
N=int(np.real(data[3]))
geometry=np.real(data[4:(3*N+4)]).astype(int)
diel=np.real(data[(3*N+4):(6*N+4)])
para=0.5*np.real(data[(6*N+4):(9*N+4)])
d=1
Shape(geometry,para, d)


"""
#Code used for optimization
#final, pos=input().split()
if __name__ == "__main__":
    objective_number = 1
    pos="./" + sys.argv[1] + "/"
    it = 0
    dec = 5
    
    if objective_number == 1:
        convergence=np.genfromtxt(pos+"convergence.txt")
        fig=plt.figure()
        plt.plot(np.arange(len(convergence)),convergence)
        plt.xlim = (0,len(convergence)-1)
        plt.savefig(pos+"convergence.png")
    else:
        convergence=np.genfromtxt(pos+"convergence.txt")
        for i in range(objective_number):
            fig=plt.figure()
            plt.plot(np.arange(len(convergence)),convergence[:,i])
            plt.xlim = (0,len(convergence)-1)
            plt.savefig(pos+"convergence_{}.png".format(i))
    for filename in sorted(os.listdir(pos+"Model_output"), key = lambda x: int(x[13:x.index(".txt")])):
        if filename.endswith(".txt"):
            print(filename)
            data=np.genfromtxt(os.path.join(pos+"Model_output",filename),dtype=complex)
            N=int(np.real(data[3]))
            #print(np.real(data[(3*N+4):(6*N+4)]))
            geometry=np.real(data[4:(3*N+4)]).astype(int)
            polarization=data[(3*N+4):(6*N+4)]
            diel=np.real(data[(6*N+4):(9*N+4)])
            d=np.real(data[9*N+4])
            wl=np.real(data[9*N+5])
            k_dir=np.real(data[(9*N+6):(9*N+9)])
            E_dir=np.real(data[(9*N+9):(9*N+12)])
            N_plot=int(np.real(data[(9*N+12)]))
            geometry_plot=np.real(data[(9*N+13):(9*N+13+N_plot)]).astype(int)
            E_tot=data[(9*N+13+N_plot):(9*N+13+2*N_plot)]
            ##print(6*N,N_plot,data.shape,E_tot.shape)
            Shape(geometry, diel, d, iteration=it, position=pos+"Shape/", decimal=dec, FullLattice=False)
            Shape(geometry, diel, d, iteration=it, position=pos+"ShapeSolid/", decimal=dec, FullLattice=True)
            if(it==99):
                EField(geometry_plot, diel, d, wl, k_dir, E_dir, E_tot, iteration=it, position=pos+"E-field/", decimal=dec)
            it += 1
"""



