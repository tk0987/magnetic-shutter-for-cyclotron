"""
DC motor torque integration code

here - semi-halbach cylinder (circular array) rotor made of 6 magnets (5x5x5 mm, in 4 'rows'), intended to diminish field inside rotor/
/and make it stronger outside.

quite interesting geometry.

please note - i have just small experience in science (just MSc, not IT). this code may produce non-exact results, but the IDEA is clear.
Example for the IDEA making real is provided, hence...

there is not python/any_language module/library dedicated for such problems. so... just look inside, maybe you will find s'th usefull for you.

SIMULATION GEOMETRY: quadrupole rotor with (s*cking simple) 4 rectangle coils, creating also 4-pole magnetic field. 
There is no core inside the coils (it is far beyond my knowledge). The current is 1 [A], and the magnetisation of neodymium magnets is 1 [T] - in uT.

if the logic is correct, this code is easily scallable: multiply the results by (your_current)/1[A]. For the differences in rotor magnets, please use magpylib like here. 
"""

import magpylib as mg
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
import matplotlib.animation as animation

def euclid_dist(x0,x1,y0,y1,z0,z1):
    return np.sqrt((x0/1000.0-x1/1000.0)**2+(y0/1000.0-y1/1000.0)**2+(z0/1000.0-z1/1000.0)**2)

#integration variables
global integration_step_mm
integration_step_mm=2
global force
force=[]
# end of integration vars
radius=17 # radius of the coils in [mm]
temp = [] # all 'tempX' lists are for storing coords of the coils
temp2=[]
temp3=[]
temp4=[]
vertices=[[-15,-2.5,-20],[-15,-2.5,20],[15,-2.5,20],[15,-2.5,-20]] # basic vertices of the wires (coils), that will be used later.
vertices2=[[-2.5,-15,-20],[-2.5,-15,20],[-2.5,15,20],[-2.5,15,-20]] # here vertices orthogonal to the previous one(s)

noturns=100 # number of turns per coil/inductor

for turns in range(noturns): # yeah, building the coils
    shift=turns/noturns
    for i in range(len(vertices)):
        a,s,d=vertices[i][0],vertices[i][1],vertices[i][2]
        temp.append([a,s+shift,d])
    for i in range(len(vertices2)):
        a,s,d=vertices2[i][0],vertices2[i][1],vertices2[i][2]
        temp2.append([a+shift,s,d])

# placing properly basic coils/inductors
for ass in range(len(temp)):
    temp[ass][1]-=17.5
for ass in range(len(temp2)):
    temp2[ass][0]-=17.5
    
# copying basic coils/inductors into opposite side of (0,0,0)
for el in temp:
    temp3.append([el[0],el[1]+39,el[2]])
for el in temp2:
    temp4.append([el[0]+39,el[1],el[2]])
# numpy rulez
temp=np.asanyarray(temp)
temp2=np.asanyarray(temp2)
temp3=np.asanyarray(temp3)
temp4=np.asanyarray(temp4)

# plot the inductors
ax = plt.figure().add_subplot(projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.plot(temp[:,0],temp[:,1],temp[:,2],c='b')
ax.plot(temp2[:,0],temp2[:,1],temp2[:,2],c='y')
ax.plot(temp3[:,0],temp3[:,1],temp3[:,2],c='g')
ax.plot(temp4[:,0],temp4[:,1],temp4[:,2],c='r')
plt.show()

# function for making n points from just two. for now it doesn't have to work at any angle, as all wires are orthogonal in 3D
# improvement will be easy (euclidean), but i do not need it. time waste
def split_for_integral(array,integration_step_mm):
    
    array=np.asanyarray(array)
    dummy=[]
    for i in range(len(array)):
        if i<len(array)-1:
            if array[i,0]-array[i+1,0]!=0.0:
                sth=0.0
                increment=0
                while increment < abs(array[i,0]-array[i+1,0])/integration_step_mm:

                    dummy.append([array[i,0]+sth,array[i,1],array[i,2]])
                    increment+=1
                    if array[i,0]-array[i+1,0]>0.0:
                        sth-=integration_step_mm
                    else:
                        sth+=integration_step_mm                   
            if abs(array[i,1]-array[i+1,1])!=0.0:
                sth=0.0
                increment=0
                while increment < abs(array[i,1]-array[i+1,1])/integration_step_mm:

                    dummy.append([array[i,0],array[i,1]+sth,array[i,2]])
                    increment+=1
                    if array[i,1]-array[i+1,1]>0.0:
                        sth-=integration_step_mm
                    else:
                        sth+=integration_step_mm                    
            if abs(array[i,2]-array[i+1,2])!=0.0:
                sth=0.0
                increment=0
                while increment < abs(array[i,2]-array[i+1,2])/integration_step_mm:

                    dummy.append([array[i,0],array[i,1],array[i,2]+sth])
                    increment+=1
                    if array[i,2]-array[i+1,2]>0.0:
                        sth-=integration_step_mm
                    else:
                        sth+=integration_step_mm                    

        if i==len(array):
            if array[i,0]-array[-1,0]!=0.0:
                sth=0.0
                increment=0
                while increment < abs(array[i,0]-array[-1,0])/integration_step_mm:

                    dummy.append([array[i,0]+sth,array[i,1],array[i,2]])
                    increment+=1
                    if array[i,0]-array[-1,0]>0.0:
                        sth-=integration_step_mm
                    else:
                        sth+=integration_step_mm                    
            if abs(array[i,1]-array[-1,1])!=0.0:
                sth=0.0
                increment=0
                while increment < abs(array[i,1]-array[-1,1])/integration_step_mm:

                    dummy.append([array[i,0],array[i,1]+sth,array[i,2]])
                    increment+=1
                    if array[i,1]-array[-1,1]>0.0:
                        sth-=integration_step_mm
                    else:
                        sth+=integration_step_mm                  
            if abs(array[i,2]-array[-1,2])!=0.0:
                sth=0.0
                increment=0
                while increment < abs(array[i,2]-array[-1,2])/integration_step_mm:

                    dummy.append([array[i,0],array[i,1],array[i,2]+sth])
                    increment+=1
                    if array[i,2]-array[-1,2]>0.0:
                        sth-=integration_step_mm
                    else:
                        sth+=integration_step_mm
    return np.asanyarray(dummy)
# making coils usable for integration of lorentz forces 
temp=split_for_integral(temp,integration_step_mm)
temp2=split_for_integral(temp2,integration_step_mm)
temp3=split_for_integral(temp3,integration_step_mm)
temp4=split_for_integral(temp4,integration_step_mm)



length=40 # lenght in Z of coils [mm]
radius=17 # radius in XY of coils [mm]
mag_no_1 = 2 # number of coils facing each other (with center at X=0 or Y=0)




current=1 # for the integration of lorentz forces - very, very important

mag_no_2=6 # number of magnets in rotor. 6 magnets produce 4-pole field, that is weak inside of the rotor, but stronger OUTSIDE. plot it, see it. just like the halbach array

# how many rows of magnets:
halbach_mag_no_1 = list()
halbach_mag_no_2 = list()
halbach_mag_no_3 = list()
halbach_mag_no_4 = list()
# dimension of a single magnet
mag_size_mm = (5,5,5)

# for rotor magnets placing - there is a logical rule - both for placing them on XY plane, and also for their orientation in space. crucial variable
angle_deg=0.0

radius_mm_1 = 9 # spacing between the center of the magnet and [x,y]=[0,0]

angler=0.0 # this is rotation of rotor - an angle between x or y (doesn't matter) axis and main axis of motor.
for angler in range(0,110,1):
    angle_deg = 0.0
    z=-9
    # making magnets of a rotor 4x:
    for i in range(0,mag_no_2,1):
        x_mag = radius_mm_1*np.cos(angle_deg+angler/110.0*np.pi/2.0) # here you have angler variable. it is a quadrupole, so pi/2 is max rotation
        y_mag = radius_mm_1*np.sin(angle_deg+angler/110.0*np.pi/2.0)
        mag_pos = (x_mag,y_mag,z)
        magnetization = (1000,0,0)
        orientation = R.from_rotvec((0,0,5*angle_deg+angler/110.0*np.pi/2.0))
        magnet = mg.magnet.Cuboid(magnetization,mag_size_mm,mag_pos,orientation)
        halbach_mag_no_1.append(magnet)
        angle_deg+=2*np.pi/mag_no_2
    angle_deg=0.0
    z=-3
    for j in range(0,mag_no_2,1):
        x_mag = radius_mm_1*np.cos(angle_deg+angler/110.0*np.pi/2.0)
        y_mag = radius_mm_1*np.sin(angle_deg+angler/110.0*np.pi/2.0)
        mag_pos = (x_mag,y_mag,z)
        magnetization = (1000,0,0)
        orientation = R.from_rotvec((0,0,5*angle_deg+angler/110.0*np.pi/2.0))
        magnet = mg.magnet.Cuboid(magnetization,mag_size_mm,mag_pos,orientation)
        halbach_mag_no_2.append(magnet)
        angle_deg+=2*np.pi/mag_no_2
    angle_deg = 0.0
    z = 3
    for ind in range(0,mag_no_2,1):
        x_mag = radius_mm_1*np.cos(angle_deg+angler/110.0*np.pi/2.0)
        y_mag = radius_mm_1*np.sin(angle_deg+angler/110.0*np.pi/2.0)
        mag_pos = (x_mag,y_mag,z)
        magnetization = (1000,0,0)
        orientation = R.from_rotvec((0,0,5*angle_deg+angler/110.0*np.pi/2.0))
        magnet = mg.magnet.Cuboid(magnetization,mag_size_mm,mag_pos,orientation)
        halbach_mag_no_3.append(magnet)
        angle_deg+=2*np.pi/mag_no_2
    angle_deg=0.0
    z=9
    for jnd in range(0,mag_no_2,1):
        x_mag = radius_mm_1*np.cos(angle_deg+angler/110.0*np.pi/2.0)
        y_mag = radius_mm_1*np.sin(angle_deg+angler/110.0*np.pi/2.0)
        mag_pos = (x_mag,y_mag,z)
        magnetization = (1000,0,0)
        orientation = R.from_rotvec((0,0,5*angle_deg+angler/110.0*np.pi/2.0))
        magnet = mg.magnet.Cuboid(magnetization,mag_size_mm,mag_pos,orientation)
        halbach_mag_no_4.append(magnet)
        angle_deg+=2*np.pi/mag_no_2
    angle_deg = 0.0

    # mg.show(halbach_mag_no_1[:]+halbach_mag_no_1[:]+halbach_mag_no_3[:]+halbach_mag_no_4[:])
# here are the lists for acquiring 'dl' - as a np.ndarray for simplicity, and maybe speed?
    dl1=np.ndarray(shape=np.shape(temp))
    dl2=np.ndarray(shape=np.shape(temp2))
    dl3=np.ndarray(shape=np.shape(temp3))
    dl4=np.ndarray(shape=np.shape(temp4))
    # making dlX's:
    for i in range(len(temp)):
        dl1[i,0], dl1[i,1], dl1[i,2]=temp[i-1,0]-temp[i,0],temp[i-1,1]-temp[i,1],temp[i-1,2]-temp[i,2]
    dl1=np.asanyarray(dl1)

    for i in range(len(temp2)):
        dl2[i,0],dl2[i,1],dl2[i,2]=temp2[i-1,0]-temp2[i,0],temp2[i-1,1]-temp2[i,1],temp2[i-1,2]-temp2[i,2]
    dl2=np.asanyarray(dl2)

    for i in range(len(temp3)):
        dl3[i,0],dl3[i,1],dl3[i,2]=temp3[i-1,0]-temp3[i,0],temp3[i-1,1]-temp3[i,1],temp3[i-1,2]-temp3[i,2]
    dl3=np.asanyarray(dl3)

    for i in range(len(temp4)):
        dl4[i,0], dl4[i,1], dl4[i,2]=temp4[i-1,0]-temp4[i,0],temp4[i-1,1]-temp4[i,1],temp4[i-1,2]-temp4[i,2]
    dl4=np.asanyarray(dl4)

    """
    dTorque = current * SUM[(dl) *(dl x B)] +++++++++++++++ CRUCIAL LIKE THE GUNS AND SWORDS ++++++++++++++++++++++++++++++++
    
    where dl =[dx,dy,dz], and B=[Bx,By,Bz] 
    
    yeah - everywhere may be error. but: (dl x B) is just force, and when multiplied by dl - torque (it is summed still). at least - by dimension
    this is the simplest method of rectangles, and i believe there is no possibility to acquire better accuracy (if there are no errors) - when it comes to magnetic fields produced by magnets and coils

    there is all about lorentz force acting on wires of the coils. I've tried integrating forces 'between sources of magnetic field' - physically bullshit, numerically... bullshit too. But... i've tried.

    if there are wires only, with 1 A only, then it is just great motor... 
    i don't mind if i'm wrong or right, but i couldn't find better free solution for this problem.

    so: i'll make my brushless cnc spindle this way. should work. and bldc controller shouldn't be very complicated.
    """
    print(1) # start integrating torques from the 1st coil
    for j in range(len(dl1)):
        if dl1[j,0]>=0 and dl1[j,1]>=0 and dl1[j,2]>=0 and current>0:
            current=-1*current
        if dl1[j,0]<=0 and dl1[j,1]<=0 and dl1[j,2]<=0 and current<0:
            current=-1*current

        x,y,z=temp[j,0],temp[j,1],temp[j,2]
        sensor = mg.Sensor((x,y,z))
        B=sensor.getB(halbach_mag_no_1[:]+halbach_mag_no_2[:]+halbach_mag_no_3[:]+halbach_mag_no_4[:],sumup=True)
        force.append([temp[j],integration_step_mm/1000.0*np.cross([dl1[j,0]/1000.0,dl1[j,1]/1000.0,dl1[j,2]/1000.0],[B[0]/1000.0,B[1]/1000.0,B[2]/1000.0]),1,angler/110*np.pi/2])# there are [mm] and [mT], remember. it is dumb, but also reasonable. remember, or you will make scum!
    print(2)# start integrating torques from the 2nd coil
    for j in range(len(dl2)):
        if dl2[j,0]>=0 and dl2[j,1]>=0 and dl2[j,2]>=0 and current<0:
            current=-1*current
        if dl2[j,0]<=0 and dl2[j,1]<=0 and dl2[j,2]<=0 and current>0:
            current=-1*current

        x,y,z=temp2[j,0],temp2[j,1],temp2[j,2]
        sensor = mg.Sensor((x,y,z))
        B=sensor.getB(halbach_mag_no_1[:]+halbach_mag_no_2[:]+halbach_mag_no_3[:]+halbach_mag_no_4[:],sumup=True)
        force.append([temp2[j],integration_step_mm/1000.0*np.cross([dl2[j,0]/1000.0,dl2[j,1]/1000.0,dl2[j,2]/1000.0],[B[0]/1000.0,B[1]/1000.0,B[2]/1000.0]),2,angler/110*np.pi/2])
    print(3)# start integrating torques from the 3rd coil
    for j in range(len(dl3)):
        if dl3[j,0]>=0 and dl3[j,1]>=0 and dl3[j,2]>=0 and current>0:
            current=-1*current
        if dl3[j,0]<=0 and dl3[j,1]<=0 and dl3[j,2]<=0 and current<0:
            current=-1*current

        x,y,z=temp3[j,0],temp3[j,1],temp3[j,2]
        sensor = mg.Sensor((x,y,z))
        B=sensor.getB(halbach_mag_no_1[:]+halbach_mag_no_2[:]+halbach_mag_no_3[:]+halbach_mag_no_4[:],sumup=True)
        force.append([temp3[j],integration_step_mm/1000.0*np.cross([dl3[j,0]/1000.0,dl3[j,1]/1000.0,dl3[j,2]/1000.0],[B[0]/1000.0,B[1]/1000.0,B[2]/1000.0]),3,angler/110*np.pi/2])
    print(4)# start integrating torques from the 4th coil
    for j in range(len(dl4)):
        if dl4[j,0]>=0 and dl4[j,1]>=0 and dl4[j,2]>=0 and current<0:
            current=-1*current
        if dl4[j,0]<=0 and dl4[j,1]<=0 and dl4[j,2]<=0 and current>0:
            current=-1*current

        x,y,z=temp4[j,0],temp4[j,1],temp4[j,2]
        sensor = mg.Sensor((x,y,z))
        B=sensor.getB(halbach_mag_no_1[:]+halbach_mag_no_2[:]+halbach_mag_no_3[:]+halbach_mag_no_4[:],sumup=True)
        force.append([temp4[j],integration_step_mm/1000.0*np.cross([dl4[j,0]/1000.0,dl4[j,1]/1000.0,dl4[j,2]/1000.0],[B[0]/1000.0,B[1]/1000.0,B[2]/1000.0]),4,angler/110*np.pi/2])
    print(angler) # know the progress!
print('saving...') # write to file, as my approach can be wrong now and later. now circa 15 MB, quite low for numerical computing tasks...
f=open("rotor_negative_current.txt","w")
for line in force:
    f.writelines(str(line[0])+'\t'+str(line[1])+'\t'+str(line[2])+'\t'+str(line[3])+'\n')
f.close()
