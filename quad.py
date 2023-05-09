import numpy as np
from misc import *

class Quad():
    def __init__(self):
        self.state = np.zeros(12)
    
    def set_state(self,state):
        self.state = state
        #xb,yb,zb,u,v,w,phi,theta,psi,p,q,r
        
    def body_rate_to_euler_rate(self):
        phi = self.state[6]
        theta = self.state[7]
        A = np.array([[1 , np.sin(phi)*np.tan(theta) , np.cos(phi)*np.tan(theta)],
                      [0, np.cos(phi),-np.sin(phi)],
                      [0, np.sin(phi)/np.cos(theta), np.cos(phi)/np.cos(theta)]])
        return A
    def euler_rate_body_rate(self):
        A = self.body_rate_to_euler_rate()
        return np.linalg.inv(A)
    
    def Ritob(self):
        phi = self.state[6]
        theta = self.state[7]
        psi = self.state[8]
        R = Rot_phi(phi) @ Rot_theta(theta) @ Rot_psi(psi)
        return R
    def Rbtoi(self):
        phi = self.state[6]
        theta = self.state[7]
        psi = self.state[8]
        R = Rot_phi(phi) @ Rot_theta(theta) @ Rot_psi(psi)
        return np.transpose(R)
    def intertial_vel(self):
        return self.Rbtoi()@np.array([self.state[3],self.state[4],self.state[5]])
    def body_coords(self):
        return np.array([self.state[0],self.state[1],self.state[2]])
    def body_angular_vel(self):
        return np.array([self.state[9],self.state[10],self.state[11]])
    def body_vel(self):
        return np.array([self.state[3],self.state[4],self.state[5]])
    def euler_anlges(self):
        return np.array([self.state[6],self.state[7],self.state[8]])
    
    def update_state(self,u):
        dxk = np.zeros(12)
        pb = self.body_coords()
        wb = self.body_angular_vel()
        vb = self.body_vel()
        Tb = np.array([u[1],u[2],u[3]])
        [dxk[0],dxk[1],dxk[2]] = vb - sk(wb)@pb
        [dxk[3],dxk[4],dxk[5]] = -sk(wb)@vb + self.Ritob() @ np.array([0,0,g])-np.array([0,0,u[0]])/mass
        [dxk[6],dxk[7],dxk[8]] = self.body_rate_to_euler_rate() @ wb
        [dxk[9],dxk[10],dxk[11]]  = Jinv @ (Tb - sk(wb)@J@wb)
        self.state = self.state + dxk*dt
        self.state[6] = modulo2pi(self.state[6])
        self.state[7] = modulopi(self.state[7])
        self.state[8] = modulo2pi(self.state[8])
if __name__ == '__main__' :
    q = Quad()
    print(q.body_coords())