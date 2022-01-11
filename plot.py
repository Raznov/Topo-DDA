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
from scipy import ndimage
import time

#matplotlib.use('Agg')
plotdpi=100
shapebarrier=0.5
#plotfor="_verify"
plotfor=""
plotlimit=None
Elimitlow=0.7
Elimithigh=2.1
colormax=2

def deleteindice(object, target, col):
    result=[]
    print(object.shape)
    for i in range(object.shape[0]):
        if object[i][col]!=target:
            result+=[i]
    
    return result

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
    minn, maxx = 0, colormax
    norm = matplotlib.colors.Normalize(minn, maxx)
    colorset = cm.ScalarMappable(norm=norm, cmap=cmaparg)
    colorset.set_array([])
    if FullLattice:
        index = np.where(diel>shapebarrier)
        diel = diel[index]
        geometry = geometry[index]
    fig = plt.figure(figsize=(10,10))
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
        plt.savefig(position+"{}Space.png".format(str(iteration).zfill(decimal)),dpi=plotdpi) 
    #plt.show()

def EField(geometry,diel,d,k_dir,E_dir,E_tot,iteration=-1,position="./",decimal=0,FullLattice=False):
    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """
    #print(geometry)

    N=round(np.shape(geometry)[0]/3)
    print(N)
    geometry=np.reshape(geometry,(N,3))
    
    print(geometry.shape)
    #diel=np.reshape(diel,(N,3))
    #diel=diel[:,0]
    

    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    
    #print(geometry)
    [X,Y,Z] = list(np.amax(geometry,axis=0)+1)
    Axis_max = max(X,Y,Z)*1.2*d
    
    print("{}, {}, {}".format(X,Y,Z))

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
    fig1.suptitle("E field - Arrow plot\n {}".format(E_dir))
    plt.savefig(position+"{}E_field_arrow.png".format(str(iteration).zfill(decimal)))

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
        fig2.suptitle("E field - Scatter plot\n, {}".format(k_dir)) 
        plt.savefig("E_field.png")
    else:
        fig2.suptitle("E field - Scatter plot\n, {} - iteration{}".format(k_dir,iteration)) 
        plt.savefig(position+"{}E_field.png".format(str(iteration).zfill(decimal)))
    plt.show()


    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """

def EField_slice(geometry,diel,d,k_dir,E_dir,E_tot,Xslice=-1,Yslice=-1,Zslice=-1,iteration=-1,position="./",decimal=0,FullLattice=False):
    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """
    #print(geometry)

    N=round(np.shape(geometry)[0]/3)
    print(N)
    geometry=np.reshape(geometry,(N,3))
    
    print(geometry.shape)
    #diel=np.reshape(diel,(N,3))
    #diel=diel[:,0]
    

    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    
    #print(geometry)
    [X,Y,Z] = list(np.amax(geometry,axis=0)+1)
    Axis_max = max(X,Y,Z)*1.2*d
    
    print("{}, {}, {}".format(X,Y,Z))

    E_tot = E_tot.reshape(int(E_tot.size/3),3)
    E_tot_real=E_tot.real
    E_tot_imag=E_tot.imag
    E_tot_abs = np.sqrt(np.sum(E_tot_real**2+E_tot_imag**2,axis=1))
    
    slicedim=-1
    

    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Y, Z))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((Z, X))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((X, Y))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-2]][geometry[i][slicedim-1]]=ele

    """
    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Z, Y))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((X, Z))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((Y, X))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-1]][geometry[i][slicedim-2]]=ele
    """
    rotated_img = ndimage.rotate(Eslice, 90)
    
    fig1 = plt.figure(figsize=(10, 10))
    plt.imshow(rotated_img, cmap='jet', interpolation='bilinear')
    if plotlimit:
        plt.clim(Elimitlow, Elimithigh)
    plt.colorbar()
    plt.savefig(position+"Model{} E_slice_{}at{}.png".format(iteration, (["X", "Y", "Z"])[slicedim], slicepos), dpi=plotdpi) 

def EField_slice_dirx(geometry,diel,d,k_dir,E_dir,E_tot,Xslice=-1,Yslice=-1,Zslice=-1,iteration=-1,position="./",decimal=0,FullLattice=False):
    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """
    #print(geometry)

    N=round(np.shape(geometry)[0]/3)
    print(N)
    geometry=np.reshape(geometry,(N,3))
    
    print(geometry.shape)
    #diel=np.reshape(diel,(N,3))
    #diel=diel[:,0]
    

    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    
    #print(geometry)
    [X,Y,Z] = list(np.amax(geometry,axis=0)+1)
    Axis_max = max(X,Y,Z)*1.2*d
    
    print("{}, {}, {}".format(X,Y,Z))

    E_tot = E_tot.reshape(int(E_tot.size/3),3)
    E_tot_abs = np.real(E_tot[:, 0])
    print(E_tot_abs.shape)
    slicedim=-1
    

    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Y, Z))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((Z, X))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((X, Y))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-2]][geometry[i][slicedim-1]]=ele

    """
    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Z, Y))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((X, Z))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((Y, X))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-1]][geometry[i][slicedim-2]]=ele
    """
    rotated_img = ndimage.rotate(Eslice, 90)
    
    fig1 = plt.figure(figsize=(10, 10))
    plt.imshow(rotated_img, cmap='jet', interpolation='bilinear')
    plt.colorbar()
    plt.savefig(position+"Model{} E_slice_{}_xdir.png".format(iteration, (["X", "Y", "Z"])[slicedim]), dpi=plotdpi) 

def EField_slice_diry(geometry,diel,d,k_dir,E_dir,E_tot,Xslice=-1,Yslice=-1,Zslice=-1,iteration=-1,position="./",decimal=0,FullLattice=False):
    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """
    #print(geometry)

    N=round(np.shape(geometry)[0]/3)
    print(N)
    geometry=np.reshape(geometry,(N,3))
    
    print(geometry.shape)
    #diel=np.reshape(diel,(N,3))
    #diel=diel[:,0]
    

    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    
    #print(geometry)
    [X,Y,Z] = list(np.amax(geometry,axis=0)+1)
    Axis_max = max(X,Y,Z)*1.2*d
    
    print("{}, {}, {}".format(X,Y,Z))

    E_tot = E_tot.reshape(int(E_tot.size/3),3)
    E_tot_abs = np.real(E_tot[:, 1])
    print(E_tot_abs.shape)
    slicedim=-1
    

    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Y, Z))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((Z, X))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((X, Y))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-2]][geometry[i][slicedim-1]]=ele

    """
    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Z, Y))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((X, Z))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((Y, X))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-1]][geometry[i][slicedim-2]]=ele
    """
    rotated_img = ndimage.rotate(Eslice, 90)
    
    fig1 = plt.figure(figsize=(10, 10))
    plt.imshow(rotated_img, cmap='jet', interpolation='bilinear')
    plt.colorbar()
    plt.savefig(position+"Model{} E_slice_{}_ydir.png".format(iteration, (["X", "Y", "Z"])[slicedim]), dpi=plotdpi) 

def EField_slice_dirz(geometry,diel,d,k_dir,E_dir,E_tot,Xslice=-1,Yslice=-1,Zslice=-1,iteration=-1,position="./",decimal=0,FullLattice=False):
    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """
    #print(geometry)

    N=round(np.shape(geometry)[0]/3)
    print(N)
    geometry=np.reshape(geometry,(N,3))
    
    print(geometry.shape)
    #diel=np.reshape(diel,(N,3))
    #diel=diel[:,0]
    

    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    
    #print(geometry)
    [X,Y,Z] = list(np.amax(geometry,axis=0)+1)
    Axis_max = max(X,Y,Z)*1.2*d
    
    print("{}, {}, {}".format(X,Y,Z))

    E_tot = E_tot.reshape(int(E_tot.size/3),3)
    E_tot_abs = np.real(E_tot[:, 2])
    print(E_tot_abs.shape)
    slicedim=-1
    

    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Y, Z))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((Z, X))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((X, Y))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-2]][geometry[i][slicedim-1]]=ele

    """
    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Z, Y))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((X, Z))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((Y, X))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-1]][geometry[i][slicedim-2]]=ele
    """
    rotated_img = ndimage.rotate(Eslice, 90)
    
    fig1 = plt.figure(figsize=(10, 10))
    plt.imshow(rotated_img, cmap='jet', interpolation='bilinear')
    plt.colorbar()
    plt.savefig(position+"Model{} E_slice_{}_zdir.png".format(iteration, (["X", "Y", "Z"])[slicedim]), dpi=plotdpi) 


def EField_slice_arrow(geometry,diel,d,k_dir,E_dir,E_tot,Xslice=-1,Yslice=-1,Zslice=-1,iteration=-1,position="./",decimal=0,FullLattice=False):
    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """
    #print(geometry)

    N=round(np.shape(geometry)[0]/3)
    print(N)
    geometry=np.reshape(geometry,(N,3))
    
    print(geometry.shape)
    #diel=np.reshape(diel,(N,3))
    #diel=diel[:,0]
    

    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    
    #print(geometry)
    [X,Y,Z] = list(np.amax(geometry,axis=0)+1)
    Axis_max = max(X,Y,Z)*1.2*d
    
    print("{}, {}, {}".format(X,Y,Z))

    E_tot = E_tot.reshape(int(E_tot.size/3),3)
    E_tot_abs = np.abs(np.sqrt(np.sum((E_tot)**2,axis=1)))
    
    slicedim=-1
    
    
    slicepos= max([Xslice, Yslice, Zslice])
    ##-------------------------------arrow field---------------------
    
    E_tot_dir_real=np.real(E_tot)    




    if Xslice!=-1:
        
        deletelist=deleteindice(geometry, slicepos, 0)
        geometryslice=np.delete(geometry, deletelist, 0)
        E_tot_dir_real_slice=np.delete(E_tot_dir_real, deletelist, 0)

    if Yslice!=-1:
        
        deletelist=deleteindice(geometry, slicepos, 1)
        geometryslice=np.delete(geometry, deletelist, 0)
        E_tot_dir_real_slice=np.delete(E_tot_dir_real, deletelist, 0)
    if Zslice!=-1:
       
        deletelist=deleteindice(geometry, slicepos, 2)
        geometryslice=np.delete(geometry, deletelist, 0)
        E_tot_dir_real_slice=np.delete(E_tot_dir_real, deletelist, 0)
    
    
    fig2 = plt.figure(figsize=(10,10))
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.quiver(geometryslice[:,0]*d,geometryslice[:,1]*d,geometryslice[:,2]*d,E_tot_dir_real_slice[:,0],E_tot_dir_real_slice[:,1],E_tot_dir_real_slice[:,2], length=60, arrow_length_ratio=0.3, linewidth=0.5)
    #print(geometryslice.shape)
    ax2.set_xlim(-(Axis_max-X*d)/2,(Axis_max+X*d)/2)
    ax2.set_ylim(-(Axis_max-Y*d)/2,(Axis_max+Y*d)/2)
    ax2.set_zlim(-(Axis_max-Z*d)/2,(Axis_max+Z*d)/2)
    ax2.set_xlabel("x[nm]")
    ax2.set_ylabel("y[nm]")
    ax2.set_zlabel("z[nm]")
    ax2.grid(False)
    
    fig2.suptitle("E field - Arrow plot\n {}".format(E_dir))
    plt.show()
    ax2.view_init(elev=0, azim=0)
    plt.savefig(position+"{}E_field_arrow at {}degree.png".format(str(iteration).zfill(decimal), "0"), dpi=plotdpi)
    ax2.view_init(elev=45, azim=30)
    plt.savefig(position+"{}E_field_arrow at {}degree.png".format(str(iteration).zfill(decimal), "45"), dpi=plotdpi)
    #ax2.view_init(elev=90, azim=0)
    #plt.savefig(position+"{}E_field_arrow at {}degree.png".format(str(iteration).zfill(decimal), "90"), dpi=1200)

def P_slice(geometry,diel,d,k_dir,E_dir,E_tot,Xslice=-1,Yslice=-1,Zslice=-1,iteration=-1,position="./",decimal=0,FullLattice=False):
    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """
    #print(geometry)

    N=round(np.shape(geometry)[0]/3)
    print(N)
    geometry=np.reshape(geometry,(N,3))
    
    print(geometry.shape)
    #diel=np.reshape(diel,(N,3))
    #diel=diel[:,0]
    

    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    
    #print(geometry)
    [X,Y,Z] = list(np.amax(geometry,axis=0)+1)
    Axis_max = max(X,Y,Z)*1.2*d
    
    print("{}, {}, {}".format(X,Y,Z))

    E_tot = E_tot.reshape(int(E_tot.size/3),3)
    E_tot_real=E_tot.real
    E_tot_imag=E_tot.imag
    E_tot_abs = np.sqrt(np.sum(E_tot_real**2+E_tot_imag**2,axis=1))
    
    #P_phase=np.zeros(shape(E_tot_real))

    slicedim=-1
    

    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Y, Z))
        Exslice=np.zeros((Y, Z))
        Eyslice=np.zeros((Y, Z))
        Ezslice=np.zeros((Y, Z))
        P1=np.zeros((Y, Z))
        P2=np.zeros((Y, Z))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((Z, X))
        Exslice=np.zeros((Z, X))
        Eyslice=np.zeros((Z, X))
        Ezslice=np.zeros((Z, X))
        P1=np.zeros((Z, X))
        P2=np.zeros((Z, X))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((X, Y))
        Exslice=np.zeros((X, Y))
        Eyslice=np.zeros((X, Y))
        Ezslice=np.zeros((X, Y))
        P1=np.zeros((X, Y))
        P2=np.zeros((X, Y))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-2]][geometry[i][slicedim-1]]=ele
            Exslice[geometry[i][slicedim-2]][geometry[i][slicedim-1]]=E_tot_real[i][0]
            Eyslice[geometry[i][slicedim-2]][geometry[i][slicedim-1]]=E_tot_real[i][1]
            Ezslice[geometry[i][slicedim-2]][geometry[i][slicedim-1]]=E_tot_real[i][2]
            P1[geometry[i][slicedim-2]][geometry[i][slicedim-1]]=E_tot_real[i][slicedim-2]
            P2[geometry[i][slicedim-2]][geometry[i][slicedim-1]]=E_tot_real[i][slicedim-1]
            Exslice
    
    

    """
    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Z, Y))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((X, Z))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((Y, X))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-1]][geometry[i][slicedim-2]]=ele
    """
    rotated_img = ndimage.rotate(Eslice, 90)
    
    fig1 = plt.figure(figsize=(10, 10))
    plt.imshow(rotated_img, cmap='jet', interpolation='bilinear')
    if plotlimit:
        plt.clim(Elimitlow, Elimithigh)
    plt.colorbar()
    plt.savefig(position+"Model{} P_slice_{}at{}.png".format(iteration, (["X", "Y", "Z"])[slicedim], slicepos), dpi=plotdpi) 

    """
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    #plt.quiver(geometry[:,0]*d, geometry[:,1]*d, np.real(E_tot[:,0]),np.real(E_tot[:,1]))
    ax2.quiver(geometry[:,0]*d,geometry[:,1]*d,geometry[:,2]*d,np.real(E_tot[:,0]),np.real(E_tot[:,1]),np.real(E_tot[:,2]),
                length=10, lw=1)
    ax2.set_xlim(-(Axis_max-X*d)/2,(Axis_max+X*d)/2)
    ax2.set_ylim(-(Axis_max-Y*d)/2,(Axis_max+Y*d)/2)
    #ax2.set_zlim(-(Axis_max-Z*d)/2,(Axis_max+Z*d)/2)
    ax2.set_xlabel("x[nm]")
    ax2.set_ylabel("y[nm]")
    #ax2.set_zlabel("z[nm]")
    ax2.grid(False)
    fig2.suptitle("P field - Arrow plot\n {}".format(E_dir))
    plt.savefig(position+"{}P_field_arrow.png".format(str(iteration).zfill(decimal)))
    """
    coordinates=[X,Y,Z]
    #coord1=coordinates[([0, 1, 2])[slicedim-2]]
    coord1=coordinates[slicedim-2]
    coord2=coordinates[slicedim-1]

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.quiver([i for i in range(coord1)], [i for i in range(coord2)], P1, P2)
    #ax3.set_xlim(-(Axis_max-X*d)/2,(Axis_max+X*d)/2)
    #ax3.set_ylim(-(Axis_max-Y*d)/2,(Axis_max+Y*d)/2)
    #ax2.set_zlim(-(Axis_max-Z*d)/2,(Axis_max+Z*d)/2)
    ax3.set_xlabel("{}".format((["X", "Y", "Z"])[slicedim-2]))
    ax3.set_ylabel("{}".format((["X", "Y", "Z"])[slicedim-1]))
    #ax2.set_zlabel("z[nm]")
    ax3.grid(False)
    fig3.suptitle("P field - Arrow plot\n {}".format(E_dir))
    plt.savefig(position+"{}P_field_arrow_slice_{}at{}.png".format(str(iteration).zfill(decimal), (["X", "Y", "Z"])[slicedim], slicepos))
    
    rotated_img = ndimage.rotate(Exslice, 90)
    fig4 = plt.figure(figsize=(10, 10))
    plt.imshow(rotated_img, cmap='jet', interpolation='bilinear')
    if plotlimit:
        plt.clim(Elimitlow, Elimithigh)
    plt.colorbar()
    plt.savefig(position+"Model{} P_slicex_{}at{}.png".format(iteration, (["X", "Y", "Z"])[slicedim], slicepos), dpi=plotdpi) 

    rotated_img = ndimage.rotate(Eyslice, 90)
    fig5 = plt.figure(figsize=(10, 10))
    plt.imshow(rotated_img, cmap='jet', interpolation='bilinear')
    if plotlimit:
        plt.clim(Elimitlow, Elimithigh)
    plt.colorbar()
    plt.savefig(position+"Model{} P_slicey_{}at{}.png".format(iteration, (["X", "Y", "Z"])[slicedim], slicepos), dpi=plotdpi) 

    rotated_img = ndimage.rotate(Ezslice, 90)
    fig6 = plt.figure(figsize=(10, 10))
    plt.imshow(rotated_img, cmap='jet', interpolation='bilinear')
    if plotlimit:
        plt.clim(Elimitlow, Elimithigh)
    plt.colorbar()
    plt.savefig(position+"Model{} P_slicez_{}at{}.png".format(iteration, (["X", "Y", "Z"])[slicedim], slicepos), dpi=plotdpi) 
    

"""
#Code used for single DDA calculation or if you only wants to see the final results for Evooptimization
pos="./" + sys.argv[1] + "/"
CoreStructure=np.genfromtxt(os.path.join(pos+"Model_output","CoreStructure0.txt"),dtype=complex)
Modelresults=np.genfromtxt(os.path.join(pos+"Model_output","Model_results0it0.txt"),dtype=complex)
N=int(np.real(CoreStructure[3]))
geometry=np.real(CoreStructure[4:(3*N+4)]).astype(int)
diel=np.real(CoreStructure[(3*N+4):(6*N+4)])
d=np.real(CoreStructure[6*N+4])

k_dir=np.real(Modelresults[(3*N):(3*N+3)])
E_dir=np.real(Modelresults[(3*N+3):(3*N+6)])

E_tot=(Modelresults[(6*N+7):(9*N+7)])
#Shape(geometry, diel, d)

yslice=35

EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Yslice=yslice,position=pos)
"""


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

"""
if __name__ == "__main__":

    objective_number = 1
    pos="./" + sys.argv[1] + "/"
    it_start = sys.argv[2]
    it_end = sys.argv[3]

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
            #polarization=data[(3*N+4):(6*N+4)]
            diel=np.real(data[(3*N+4):(6*N+4)])
            d=np.real(data[6*N+4])
            #wl=np.real(data[9*N+5])
            #k_dir=np.real(data[(9*N+6):(9*N+9)])
            #E_dir=np.real(data[(9*N+9):(9*N+12)])
            #N_plot=int(np.real(data[(9*N+12)]))
            #geometry_plot=np.real(data[(9*N+13):(9*N+13+N_plot)]).astype(int)
            #E_tot=data[(9*N+13+N_plot):(9*N+13+2*N_plot)]
            ##print(6*N,N_plot,data.shape,E_tot.shape)
            if(it >= int(it_start) and it <= int(it_end)):
                Shape(geometry, diel, d, iteration=it, position=pos+"Shape/", decimal=dec, FullLattice=False)
                Shape(geometry, diel, d, iteration=it, position=pos+"ShapeSolid/", decimal=dec, FullLattice=True)

            #if(it==99):
            #    EField(geometry_plot, diel, d, wl, k_dir, E_dir, E_tot, iteration=it, position=pos+"E-field/", decimal=dec)

            it += 1
"""

"""
pos="./theta0phi0-lam500-size1000-d12d5-1-diel4-FCD/"
convergence=np.genfromtxt(pos+"convergence.txt")
fig=plt.figure()
plt.plot(np.arange(len(convergence)),convergence)
plt.xlim = (0,len(convergence)-1)
plt.savefig(pos+"convergence.png")
"""


#For Single DDA calculation using CoreSturcutre to see the geometry
"""
if __name__ == "__main__":

    objective_number = 1
    pos=".\\" + sys.argv[1] + "\\"
    filename=sys.argv[2]+".txt"

    it = 0
    dec = 5

    
    print(filename)
    data=np.genfromtxt(os.path.join(pos, filename),dtype=complex)
    N=int(np.real(data[3]))
    #print(np.real(data[(3*N+4):(6*N+4)]))
    geometry=np.real(data[4:(3*N+4)]).astype(int)
    #polarization=data[(3*N+4):(6*N+4)]
    diel=np.real(data[(3*N+4):(6*N+4)])
    d=np.real(data[6*N+4])
    #wl=np.real(data[9*N+5])
    #k_dir=np.real(data[(9*N+6):(9*N+9)])
    #E_dir=np.real(data[(9*N+9):(9*N+12)])
    #N_plot=int(np.real(data[(9*N+12)]))
    #geometry_plot=np.real(data[(9*N+13):(9*N+13+N_plot)]).astype(int)
    #E_tot=data[(9*N+13+N_plot):(9*N+13+2*N_plot)]
    ##print(6*N,N_plot,data.shape,E_tot.shape)
    Shape(geometry, diel, d, iteration=it, position=pos+"Shape/", decimal=dec, FullLattice=False)
    Shape(geometry, diel, d, iteration=it, position=pos+"ShapeSolid/", decimal=dec, FullLattice=True)

    #if(it==99):
    #    EField(geometry_plot, diel, d, wl, k_dir, E_dir, E_tot, iteration=it, position=pos+"E-field/", decimal=dec)
"""
"""
#For several DDA calculation Structure

if __name__ == "__main__":

    objective_number = 1
    pos="./" + sys.argv[1] + "/"
    it_start = sys.argv[2]
    it_end = sys.argv[3]

    it = 0
    dec = 5
    for filename in sorted(os.listdir(pos+"CoreStructure"), key = lambda x: int(x[13:x.index(".txt")])):
        if filename.endswith(".txt"):
            print(filename)
            data=np.genfromtxt(os.path.join(pos+"CoreStructure",filename),dtype=complex)
            N=int(np.real(data[3]))
            #print(np.real(data[(3*N+4):(6*N+4)]))
            geometry=np.real(data[4:(3*N+4)]).astype(int)
            #polarization=data[(3*N+4):(6*N+4)]
            diel=np.real(data[(3*N+4):(6*N+4)])
            d=np.real(data[6*N+4])
            if(it >= int(it_start) and it <= int(it_end)):
                Shape(geometry, diel, d, iteration=it, position=pos+"Shape/", decimal=dec, FullLattice=False)
                Shape(geometry, diel, d, iteration=it, position=pos+"ShapeSolid/", decimal=dec, FullLattice=True)

            it += 1
"""

"""
#For several DDA calculation E field

if __name__ == "__main__":

    objective_number = 1
    pos="./" + sys.argv[1] + "/"
    it_start = sys.argv[2]
    it_end = sys.argv[3]

    dec = 5
    for it in range(int(it_start), int(it_end)+1):
        pos="./" + sys.argv[1] + "/"
        CoreStructure=np.genfromtxt(os.path.join(pos+"CoreStructure","CoreStructure"+str(it)+".txt"),dtype=complex)
        Modelresults=np.genfromtxt(os.path.join(pos+"Model_output","Model_results"+str(it)+"it0.txt"),dtype=complex)
        N=int(np.real(CoreStructure[3]))
        geometry=np.real(CoreStructure[4:(3*N+4)]).astype(int)
        diel=np.real(CoreStructure[(3*N+4):(6*N+4)])
        d=np.real(CoreStructure[6*N+4])
        k_dir=np.real(Modelresults[(3*N):(3*N+3)])
        E_dir=np.real(Modelresults[(3*N+3):(3*N+6)])
        E_tot=(Modelresults[(6*N+7):(9*N+7)])
        zslice=10

        EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, iteration=it, Zslice=zslice,position=pos+"E-field/")
"""

"""
#For several DDA calculation Structure with simplified output

if __name__ == "__main__":
    print('fuck')
    objective_number = 1
    pos="./" + sys.argv[1] + "/"
    it_start = sys.argv[2]
    it_end = sys.argv[3]
    datacommon=np.genfromtxt(pos+"Commondata.txt")
    N=int(np.real(datacommon[3]))
    geometry=np.real(datacommon[4:(3*N+4)]).astype(int)
    d=np.real(datacommon[3*N+4])
    it = 1
    dec = 5
    for filename in sorted(os.listdir(pos+"CoreStructure"), key = lambda x: int(x[13:x.index(".txt")])):
        if filename.endswith(".txt"):
            print(filename)
            name=int((filename[13:])[:-4])
            data=np.genfromtxt(os.path.join(pos+"CoreStructure",filename),dtype=complex)
            diel=np.real(data[0:3*N])
            if(it >= int(it_start) and it <= int(it_end)):
                #Shape(geometry, diel, d, iteration=name, position=pos+"Shape/", decimal=dec, FullLattice=False)
                Shape(geometry, diel, d, iteration=name, position=pos+"ShapeSolid/", decimal=dec, FullLattice=True)
                
                it += 1
"""



#For several DDA calculation E field with simplifed

if __name__ == "__main__":

    objective_number = 1
    pos="./" + sys.argv[1] + "/"
    it_start = sys.argv[2]
    it_end = sys.argv[3]
    
    datacommon=np.genfromtxt(pos+"Commondata.txt")
    N=int(np.real(datacommon[3]))
    geometry=np.real(datacommon[4:(3*N+4)]).astype(int)
    d=np.real(datacommon[3*N+4])
    E_dir=np.real(datacommon[(3*N+5):(3*N+8)])
    k_dir=np.real(datacommon[(3*N+8):(3*N+11)])
    
    dec = 5
    cutnumber=13
    for filename in sorted(os.listdir(pos+"CoreStructure"+plotfor), key = lambda x: int(x[cutnumber:x.index(".txt")])):
        if filename.endswith(".txt"):
            #print(filename)
            nameit=int((filename[cutnumber:])[:-4])
            print(nameit)
            CoreStructure=np.genfromtxt(os.path.join(pos+"CoreStructure"+plotfor,"CoreStructure"+str(nameit)+".txt"),dtype=complex)
            Modelresults=np.genfromtxt(os.path.join(pos+"Model_output"+plotfor,"Model_results"+"it"+str(nameit)+".txt"),dtype=complex)
        
            diel=np.real(CoreStructure[(0):(3*N)])
            E_tot=(Modelresults[(0):(3*N)])
            P_tot=(Modelresults[(3*N):(6*N)])
            zslice=17
            if(nameit >= int(it_start) and nameit <= int(it_end)):
                Shape(geometry, diel, d, iteration=nameit, position=pos+"ShapeSolid"+plotfor+"/", decimal=dec, FullLattice=True)
                Shape(geometry, diel, d, iteration=nameit, position=pos+"Shape"+plotfor+"/", decimal=dec, FullLattice=False)
                EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, iteration=nameit, Zslice=zslice,position=pos+"E-field"+plotfor+"/")         #----------electric field intensity-----------
                #P_slice(geometry, diel, d, k_dir, E_dir, P_tot, iteration=nameit, Zslice=zslice,position=pos+"E-field"+plotfor+"/")         #----------For p----------------
                #EField_slice_dirx(geometry, diel, d, k_dir, E_dir, E_tot, iteration=it, Zslice=zslice,position=pos+"E-field/")     #----------real Ex-----------------
                #EField_slice_diry(geometry, diel, d, k_dir, E_dir, E_tot, iteration=it, Zslice=zslice,position=pos+"E-field/")     #----------real Ey-----------------
                #EField_slice_dirz(geometry, diel, d, k_dir, E_dir, E_tot, iteration=it, Zslice=zslice,position=pos+"E-field/")     #----------real Ez-----------------
                #EField_slice_arrow(geometry, diel, d, k_dir, E_dir, E_tot, iteration=it, Zslice=zslice,position=pos+"E-field/")   #-----------Vector field (Ex.real, Ey.real, Ez.real) in a cross section----------




pos="./" + sys.argv[1] + "/"
filename=os.path.join(pos, 'Loss.txt')
Loss=np.genfromtxt(filename)
    
fig=plt.figure()
plt.plot(Loss[:,0],Loss[:,1])
plt.xlim = (0,max(Loss[:,0]))
plt.savefig(pos+"convergence.png")





"""
if __name__ == "__main__":

    objective_number = 1
    pos="./" + sys.argv[1] + "/"
    it_start = sys.argv[2]
    it_end = sys.argv[3]
    
    datacommon=np.genfromtxt(pos+"Commondata.txt")
    N=int(np.real(datacommon[3]))
    geometry=np.real(datacommon[4:(3*N+4)]).astype(int)
    d=np.real(datacommon[3*N+4])
    E_dir=np.real(datacommon[(3*N+5):(3*N+8)])
    k_dir=np.real(datacommon[(3*N+8):(3*N+11)])
    
    dec = 5
    for it in range(int(it_start), int(it_end)+1):
        pos="./" + sys.argv[1] + "/"
        CoreStructure=np.genfromtxt(os.path.join(pos+"CoreStructure","CoreStructure"+str(it)+".txt"),dtype=complex)
        Modelresults=np.genfromtxt(os.path.join(pos+"Model_output","Model_results"+"it"+str(it)+".txt"),dtype=complex)
        
        diel=np.real(CoreStructure[(0):(3*N)])
        E_tot=(Modelresults[(0):(3*N)])
        

        EField(geometry, diel, d, k_dir, E_dir, E_tot, iteration=it, position=pos+"E-field/")
"""




"""
if __name__ == "__main__":
    objective_number = 1
    pos="./" + sys.argv[1] + "/"
    it_start = sys.argv[2]
    it_end = sys.argv[3]
    
    datacommon=np.genfromtxt(pos+"Commondata.txt")
    N=int(np.real(datacommon[3]))
    geometry=np.real(datacommon[4:(3*N+4)]).astype(int)
    d=np.real(datacommon[3*N+4])
    E_dir=np.real(datacommon[(3*N+5):(3*N+8)])
    k_dir=np.real(datacommon[(3*N+8):(3*N+11)])
    
    dec = 5

    cutnumber=13
    for filename in sorted(os.listdir(pos+"CoreStructure_verify"), key = lambda x: int(x[cutnumber:x.index(".txt")])):
        if filename.endswith(".txt"):
            
            nameit=int((filename[cutnumber:])[:-4])
            print(nameit)
            CoreStructure=np.genfromtxt(os.path.join(pos+"CoreStructure_verify","CoreStructure"+str(nameit)+".txt"),dtype=complex)
            Modelresults=np.genfromtxt(os.path.join(pos+"Model_output_verify","Model_results"+"it"+str(nameit)+".txt"),dtype=complex)
        
            diel=np.real(CoreStructure[(0):(3*N)])
            E_tot=(Modelresults[(0):(3*N)])
            zslice=5
            if(nameit >= int(it_start) and nameit <= int(it_end)):
                Shape(geometry, diel, d, iteration=nameit, position=pos+"ShapeSolid_verify/", decimal=dec, FullLattice=True)
                EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, iteration=nameit, Zslice=zslice,position=pos+"E-field_verify/")         #----------electric field intensity-----------
                #EField_slice_dirx(geometry, diel, d, k_dir, E_dir, E_tot, iteration=it, Zslice=zslice,position=pos+"E-field/")     #----------real Ex-----------------
                #EField_slice_diry(geometry, diel, d, k_dir, E_dir, E_tot, iteration=it, Zslice=zslice,position=pos+"E-field/")     #----------real Ey-----------------
                #EField_slice_dirz(geometry, diel, d, k_dir, E_dir, E_tot, iteration=it, Zslice=zslice,position=pos+"E-field/")     #----------real Ez-----------------
                #EField_slice_arrow(geometry, diel, d, k_dir, E_dir, E_tot, iteration=it, Zslice=zslice,position=pos+"E-field/")   #-----------Vector field (Ex.real, Ey.real, Ez.real) in a cross section----------
"""

