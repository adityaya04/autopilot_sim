import numpy as np
from quad import Quad
from misc import *
import matplotlib.pyplot as plt
plt.style.use('seaborn-darkgrid')

Time = 20
num_time_slices = int(Time/dt) 
time_slices = np.linspace(0,Time,num_time_slices)

'''Creating arrays filled with zeros to store
the data at each time step to be used for plotting'''
xF = np.zeros((12,num_time_slices))   # Each row contains the state of the quadrotor at a time step
U = np.zeros((4,num_time_slices))     # Each row contains the inputs to the quadrotor at a time step
Rot_D = np.zeros((3,num_time_slices)) # Each row contains the euler angles of the quadrotor body at a time step
VB = np.zeros((3,num_time_slices))    # Each row contains the body velocities of the quadrotor along the three axes at a time step
P_D = np.zeros((3,num_time_slices))   # Each row contains the desired position of the quadrotor at a time step

pd_i = np.array([0.,0.,-1.2]) # Desired position of the quadrotor

x0 = np.array([0,0,0,0,0,0,0,0,0,0,0,0]) # Initial state of the quadrotor

quad = Quad() # Initialising the quadrotor
quad.set_state(x0) # Setting the initial state of the quadrotor
xF[:,0] = x0

for i in range(0,num_time_slices):
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
figure, axis = plt.subplots(1, 2)

axis[0].plot(time_slices,alt,'g', label = 'Current altitude')
axis[0].plot(time_slices,-P_D[2,:],'b', label = 'Desired')
axis[0].set_title("Altitude")
axis[0].legend()
axis[1].plot(time_slices,thrust,'r')
axis[1].set_title("Thrust")

plt.setp(axis[0],xlabel='Time (s)',ylabel = 'Altitude (m)')
plt.setp(axis[1],xlabel='Time (s)',ylabel = 'Thrust (N)')

plt.show()
