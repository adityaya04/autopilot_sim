import numpy as np
from quad import Quad
from misc import *
import matplotlib.pyplot as plt
plt.style.use('seaborn-darkgrid')

Time = 300
num_time_slices = int(Time/dt) 
time_slices = np.linspace(0,Time,num_time_slices)

'''Creating a trajectory'''
size = int(num_time_slices/5)
x0 = np.zeros(size)
y0 = np.zeros(size)
SCALE = 50
vertex0 = np.array([0,0])
vertex1 = np.array([0.5, 0.5])*SCALE
vertex2 = np.array([-0.5, 0.5])*SCALE
vertex3 = np.array([-0.5, -0.5])*SCALE

x1 = np.linspace(vertex0[0], vertex1[0], size)
y1 = np.linspace(vertex0[1], vertex1[1], size)

x2 = np.linspace(vertex1[0], vertex2[0], size)
y2 = np.linspace(vertex1[1], vertex2[1], size)

x3 = np.linspace(vertex2[0], vertex3[0], size)
y3 = np.linspace(vertex2[1], vertex3[1], size)

x4 = np.linspace(vertex3[0], vertex0[0], size)
y4 = np.linspace(vertex3[1], vertex0[1], size)

x_coords = np.concatenate((x0,x1, x2, x3, x4))
y_coords = np.concatenate((y0,y1, y2, y3, y4))

'''Creating arrays filled with zeros to store
the data at each time step to be used for plotting'''
xF = np.zeros((12,num_time_slices))   # Each row contains the state of the quadrotor at a time step
U = np.zeros((4,num_time_slices))     # Each row contains the inputs to the quadrotor at a time step
Rot_D = np.zeros((3,num_time_slices)) # Each row contains the euler angles of the quadrotor body at a time step
VB = np.zeros((3,num_time_slices))    # Each row contains the body velocities of the quadrotor along the three axes at a time step
P_D = np.zeros((3,num_time_slices))   # Each row contains the desired position of the quadrotor at a time step

x0 = np.array([0,0,0,0,0,0,0,0,0,0,0,0]) # Initial state of the quadrotor

quad = Quad() # Initialising the quadrotor
quad.set_state(x0) # Setting the initial state of the quadrotor
xF[:,0] = x0

for i in range(0,num_time_slices):
    pd_i = np.array([x_coords[i],y_coords[i],-20.]) # Desired position of the quadrotor according to the triangle trajectory
    pd_b = quad.Ritob()@pd_i  # Desired position in body frame
    pb = quad.body_coords()
    wb = quad.body_angular_vel()
    vb = quad.body_vel()
    pd_b += -(sk(wb)@pd_b)*dt # Updating desired position
    del1 = pb - pd_b
    
    y1 = K2*del1 + vb
    y2 = vb
    y2 = np.array(y2)
    y1 = np.array(y1)
    
    u = -Sigma2(K2*y2 + Sigma1(K1*y1)) # Calculates the inputs to the inner loop
    
    thetad = -np.arcsin(u[0]/g)
    phid = np.arcsin(u[1]/(g*np.cos(thetad)))
    fd = mass*(g*np.cos(thetad)*np.cos(phid)-u[2]) 

    psid = 0  # Has to be zero for other equations to be correct
    n = quad.euler_anlges()
    nd = np.array([phid,thetad,psid])
    e_n = n-nd
    e_w = (quad.body_rate_to_euler_rate() @ wb) 
    Tbarb = -(Kpa @ e_n) - (Kda @ e_w)
    Winv = quad.euler_rate_body_rate()
    Tau = J @ (Winv@Tbarb + Jinv@sk(wb)@J@wb)
    [Tx,Ty,Tz] = Tau
    
    inp = np.array([fd,Tx,Ty,Tz]) # Thrust and torque inputs to the quadrotor
    
    quad.update_state(inp) # Updating the states of the quadrotor with the given input
    
    '''Updating the empty arrays with data at each time step'''
    P_D[:,i] = pd_i
    Rot_D[:,i] = [phid*RAD_TO_DEG,thetad*RAD_TO_DEG,psid*RAD_TO_DEG]
    xF[:,i] = quad.state
    VB[:,i] = quad.Rbtoi()@vb 
    U[:,i] = inp
    
    
   
    
    
thrust = U[0,:]
Taux = U[1,:]
Tauy = U[2,:]
Tauz = U[3,:]
PN = xF[0,:]
PE = xF[1,:]
alt = -xF[2,:]
Phi = xF[6,:]
Theta = xF[7,:]
Psi = xF[8,:]
Phi = RAD_TO_DEG*Phi
Theta = RAD_TO_DEG*Theta
Psi = RAD_TO_DEG*Psi
def plot_kinematics():
    figure, axis = plt.subplots(3, 3)
    axis[0,0].plot(time_slices,PN,'g', label = 'Current X')
    axis[0,0].plot(time_slices,P_D[0,:],'b', label = 'Desired')
    axis[0,0].legend()
    axis[0,1].plot(time_slices,PE,'g', label = 'Current Y')
    axis[0,1].plot(time_slices,P_D[1,:],'b', label = 'Desired')
    axis[0,1].legend()
    axis[0,2].plot(time_slices,alt,'g', label = 'Current altitude')
    axis[0,2].plot(time_slices,-P_D[2,:],'b', label = 'Desired')
    axis[0,2].legend()

    axis[1,0].plot(time_slices,Phi,'g', label = '$\phi$')
    axis[1,0].plot(time_slices,Rot_D[0,:],'b', label = 'Desired')
    axis[1,0].legend()
    axis[1,1].plot(time_slices,Theta,'g', label='\u03B8')
    axis[1,1].plot(time_slices,Rot_D[1,:],'b', label='Desired')
    axis[1,1].legend()
    axis[1,2].plot(time_slices,Psi,'g', label='$\psi$')
    axis[1,2].plot(time_slices,Rot_D[2,:],'b', label='Desired')
    axis[1,2].legend()


    axis[2,0].plot(time_slices,VB[0,:],'r')
    axis[2,1].plot(time_slices,VB[1,:],'r')
    axis[2,2].plot(time_slices,VB[2,:],'r')

    plt.setp(axis[2,:],xlabel= 'Time (s)')
    plt.setp(axis[0,0],ylabel = 'Inertial X (m)')
    plt.setp(axis[0,1],ylabel = 'Inertial Y (m)')
    plt.setp(axis[0,2],ylabel = 'Altitude (m)')

    plt.setp(axis[1,0],ylabel = '$\phi$ (deg)')
    plt.setp(axis[1,1],ylabel = '\u03B8 (deg)')
    plt.setp(axis[1,2],ylabel = '$\psi$ (deg)')

    plt.setp(axis[2,0],ylabel = '$v^i_x$ (m/s)')
    plt.setp(axis[2,1],ylabel = '$v^i_y$ (m/s)')
    plt.setp(axis[2,2],ylabel = '$v^i_z$ (m/s)')
def plot_dynamics():
    figure, axis = plt.subplots(2, 2)
    axis[0,0].plot(time_slices,thrust,'b')
    axis[0,1].plot(time_slices,Taux,'b')
    
    axis[1,0].plot(time_slices,Tauy,'b')
    axis[1,1].plot(time_slices,Tauz,'b')

    plt.setp(axis[1,:],xlabel= 'Time (s)')
    plt.setp(axis[0,0],ylabel = 'Thrust (N)')
    plt.setp(axis[0,1],ylabel = '\u03C4$_x$ (Nm)')

    plt.setp(axis[1,0],ylabel = '\u03C4$_y$ (Nm)')
    plt.setp(axis[1,1],ylabel = '\u03C4$_z$ (Nm)')

def plot_pos():
    figure, axis = plt.subplots(3, 1)
    axis[0].plot(time_slices,PN,'g', label = 'Current X')
    axis[0].plot(time_slices,P_D[0,:],'b', label = 'Desired')
    axis[0].legend()
    axis[1].plot(time_slices,PE,'g', label = 'Current Y')
    axis[1].plot(time_slices,P_D[1,:],'b', label = 'Desired')
    axis[1].legend()
    axis[2].plot(time_slices,alt,'g', label = 'Current altitude')
    axis[2].plot(time_slices,-P_D[2,:],'b', label = 'Desired')
    axis[2].legend()

    plt.setp(axis[2],xlabel= 'Time (s)')
    plt.setp(axis[0],ylabel = 'Inertial X (m)')
    plt.setp(axis[1],ylabel = 'Inertial Y (m)')
    plt.setp(axis[2],ylabel = 'Altitude (m)')
    
def plot_angles():
    figure, axis = plt.subplots(3, 1)
    axis[0].plot(time_slices,Phi,'g', label = '$\phi$')
    axis[0].plot(time_slices,Rot_D[0,:],'b', label = 'Desired')
    axis[0].legend()
    axis[1].plot(time_slices,Theta,'g', label='\u03B8')
    axis[1].plot(time_slices,Rot_D[1,:],'b', label='Desired')
    axis[1].legend()
    axis[2].plot(time_slices,Psi,'g', label='$\psi$')
    axis[2].plot(time_slices,Rot_D[2,:],'b', label='Desired')
    axis[2].legend()

    plt.setp(axis[2],xlabel= 'Time (s)')
    plt.setp(axis[0],ylabel = '$\phi$ (deg)')
    plt.setp(axis[1],ylabel = '\u03B8 (deg)')
    plt.setp(axis[2],ylabel = '$\psi$ (deg)')

def plot_vel():
    figure, axis = plt.subplots(3, 1)
    axis[0].plot(time_slices,VB[0,:],'r')
    axis[1].plot(time_slices,VB[1,:],'r')
    axis[2].plot(time_slices,VB[2,:],'r')
    
    plt.setp(axis[2],xlabel= 'Time (s)')
    plt.setp(axis[0],ylabel = '$v^i_x$ (m/s)')
    plt.setp(axis[1],ylabel = '$v^i_y$ (m/s)')
    plt.setp(axis[2],ylabel = '$v^i_z$ (m/s)')
def plot_trajectory():
    plt.plot(PN,PE,'g',label = 'Actual path')
    plt.plot(P_D[0,:],P_D[1,:],'b', label = 'Desired path')
    plt.xlabel("Inertial X (m)")
    plt.ylabel("Inertial Y (m)")
    plt.legend()
    
plot_kinematics()
plt.show()
