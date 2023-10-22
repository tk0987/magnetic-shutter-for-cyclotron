#++++++++++++++++++++++++++++++++++++++++++++++++++
# note the high inductance of such coil arrangement.
# here: all coils in paralell, 10 A in each one. This current creates B>150 mT inside of the system. More than enough
#
# such EM shutter was designed as an aid in flash measurements
#
# proper PSU is needed, but coil array making (PCB, and so on) is very simple

import magpylib as mg
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
import matplotlib.animation as animation
# import matplotlib
x_mm=0.0
y_mm=0.0
z_mm=-100.0

mag_no_1 = 6
mag_no_2=12
inner=[]
coil1 = mg.Collection()
angle_deg = 0.0

for zfgh in range(0,mag_no_1,1):
    for width in range(0,2000,1):
        x_mag = 55*np.cos(angle_deg)
        y_mag = 55*np.sin(angle_deg)
        winding = mg.current.Line(
            current=10/2000,
            vertices=[[0,width*5/2000+7.5,-100],[0,width*5/2000+7.5,100+width*5/2000],[0,-7.5-width*5/2000,100+width*5/2000],[0,-7.5-width*5/2000,-90-width*5/2000],[0,3-width*5/2000,-90-width*5/2000],[0,3-width*5/2000,-100]],
            position=(x_mag,y_mag,0),
            orientation= R.from_euler("xyz",(0,0,2*angle_deg)),
        )
        inner.append(winding)
    angle_deg+=2*np.pi/mag_no_1
    print(zfgh)
angle_deg=0.0

outer=[]
coil2 = mg.Collection()

for zgfd in range(0,mag_no_2,1):
    for width in range(0,2000,1):
        x_mag2 = 75*np.cos(angle_deg)
        y_mag2 = 75*np.sin(angle_deg)
        winding2 = mg.current.Line(
            current=10/2000,
            vertices=[[0,width*5/2000+7.5,-100],[0,width*5/2000+7.5,100+width*5/2000],[0,-7.5-width*5/2000,100+width*5/2000],[0,-7.5-width*5/2000,-90-width*5/2000],[0,3-width*5/2000,-90-width*5/2000],[0,3-width*5/2000,-100]],
            
            position=(x_mag2,y_mag2,0),
            orientation= R.from_euler("xyz",(0,0,2*angle_deg)),
        )
        outer.append(winding2)
        angle_deg+=2*np.pi/mag_no_2
    angle_deg+=2*np.pi/mag_no_2
    print(zgfd)

# mg.show(inner[:]+outer[:]) # uncomment if wanna see this beauty. it is horrible for cpu, though


FOI_size_mm = 40
step_mm = 2
FOI_x = np.linspace(-FOI_size_mm/2,FOI_size_mm/2,int(FOI_size_mm/step_mm))
FOI_y = np.linspace(-FOI_size_mm/2,FOI_size_mm/2,int(FOI_size_mm/step_mm))
FOI_b = list()

X,Y=np.meshgrid(FOI_x,FOI_y)



fieldz=np.ndarray(shape=(len(X),len(Y)))
fieldz.fill(0)
fieldX=np.ndarray(shape=(len(X),len(Y)))
fieldX.fill(0)
fieldY=np.ndarray(shape=(len(X),len(Y)))
fieldY.fill(0)
fieldZ=np.ndarray(shape=(len(X),len(Y)))
fieldZ.fill(0)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++contour fill loop
for column in range(0,len(X),1):
    for row in range(0,len(Y),1):
        sensor = mg.Sensor((X[column,row],Y[column,row],0))
        B=sensor.getB(outer[:]+inner[:],sumup=True)
        fieldz[column,row]=np.sqrt(B[0]**2+B[1]**2+B[2]**2)
        fieldX[column,row]=B[0]
        fieldY[column,row]=B[1]
        fieldZ[column,row]=B[2]
    print(100*column/len(X),"%")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++b total
fig, axs = plt.subplots(layout='constrained')
CS = axs.contourf(X, Y, fieldz, 85, cmap="gist_ncar")
axs.set_title('B_total, Z=0 (half coil)')
axs.set_xlabel('X [mm]')
axs.set_ylabel('Y [mm]')
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel('B_total [T]')
plt.show()
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++bx
fig, axs = plt.subplots(layout='constrained')
CS = axs.contourf(X, Y, fieldX, 85, cmap="gist_ncar")
axs.set_title('B_x, Z=0 (half coil)')
axs.set_xlabel('X [mm]')
axs.set_ylabel('Y [mm]')
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel('B_x [T]')
plt.show()
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++by
fig, axs = plt.subplots(layout='constrained')
CS = axs.contourf(X, Y, fieldY, 85, cmap="gist_ncar")
axs.set_title('B_y, Z=0 (half coil)')
axs.set_xlabel('X [mm]')
axs.set_ylabel('Y [mm]')
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel('B_y [T]')
plt.show()
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++bz
fig, axs = plt.subplots(layout='constrained')
CS = axs.contourf(X, Y, fieldZ, 85, cmap="gist_ncar")
axs.set_title('B_z, Z=0 (half coil)')
axs.set_xlabel('X [mm]')
axs.set_ylabel('Y [mm]')
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel('B_z [T]')
plt.show()



# ++++++++++++++++++++++++++++++++++++++++++++++++++ for trajectories:

# the energy&momentum conservation rules for velocities can be written as follows:

#     vz=(np.sqrt((vz_0**2)-(vx**2)-(vy**2))-az*time_step_s)

# where our particle moves with vz_0 at beginning, then - with acquiring vx and vy - vz reduces. in the term 'az*time_step_s' az is acceleration from lorentz force paralell to z axis. 
# time step is arbitrary, depending on what you want



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++linear plots in x=y=0 (only z)
fieldz=np.ndarray(shape=(100))
fieldX=np.ndarray(shape=(100))
fieldY=np.ndarray(shape=(100))
fieldZ=np.ndarray(shape=(100))
for zett in range(0,100,1):
    sensor = mg.Sensor((0,0,zett*400/100-200))
    B=sensor.getB(inner[:]+outer[:],sumup=True)
    fieldz[zett]=(np.sqrt(B[0]**2+B[1]**2+B[2]**2))
    fieldX[zett]=B[0]
    fieldY[zett]=B[1]
    fieldZ[zett]=B[2]
    if zett%5==0:
        print(zett/100+0.5)
Z=np.linspace(-200.,200.,100)

 



fig,ax=plt.subplots(layout='constrained')
plt.plot(Z,fieldz)
ax.set_title('B_total, X=0, Y=0 (Z axis)')
ax.set_xlabel('Z [mm]')
ax.set_ylabel('B [T]')
plt.show()
fig,ax=plt.subplots(layout='constrained')
plt.plot(Z,fieldX)
ax.set_title('B_X, X=0, Y=0 (Z axis)')
ax.set_xlabel('Z [mm]')
ax.set_ylabel('B [T]')
plt.show()
fig,ax=plt.subplots(layout='constrained')
plt.plot(Z,fieldY)
ax.set_title('B_Y, X=0, Y=0 (Z axis)')
ax.set_xlabel('Z [mm]')
ax.set_ylabel('B [T]')
plt.show()
fig,ax=plt.subplots(layout='constrained')
plt.plot(Z,fieldZ)
ax.set_title('B_Z, X=0, Y=0 (Z axis)')
ax.set_xlabel('Z [mm]')
ax.set_ylabel('B [T]')
plt.show()
