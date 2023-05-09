import numpy as np

'''Constants'''
dt = 0.01

RAD_TO_DEG = 180/np.pi
DEG_TO_RAD = np.pi/180

g = 9.81
mass = 0.486
Jxx = 820e-4
Jyy = 845e-4
Jzz = 1377e-4
Jxz = 0
J = np.diag([Jxx,Jyy,Jzz])
Jinv = np.linalg.inv(J)

Kpa = np.diag([7,7,7])
Kda = np.diag([2,2,2])

K1 = 1
K2 = 0.5
K3 = 3
M = np.array([[0.58,0.58,0.58],[3.8,3.8,3.8]])
L = 0.95*M
B = M-L
'''Elementary rotation matrices'''
def Rot_phi(angle):
        R = np.array([[1, 0 ,0],[0, np.cos(angle),np.sin(angle)],[0,-np.sin(angle),np.cos(angle)]])
        return R
def Rot_theta(angle):
    R = np.array([[np.cos(angle),0,-np.sin(angle)],[0, 1, 0],[np.sin(angle), 0, np.cos(angle)]])
    return R
def Rot_psi(angle):
    R = np.array([[np.cos(angle), np.sin(angle), 0],[-np.sin(angle), np.cos(angle), 0],[0, 0,1]])
    return R

'''Vector to skew-symmetric matrix'''
def sk(vec):
    A = np.array([[0,-vec[2],vec[1]],[vec[2],0,-vec[0]],[-vec[1],vec[0],0]])
    return(A)

'''Wrap an angle in radians to the range -pi to pi'''
def modulo2pi(angle):
    while angle > np.pi:
        angle -= 2*np.pi
    while angle < -np.pi:
        angle += 2*np.pi
    return angle

'''Wrap an angle in radians to the range -pi/2 to pi/2'''
def modulopi(angle):
    while angle > np.pi/2:
        angle -= np.pi
    while angle < -np.pi/2:
        angle += np.pi
    return angle

def sign(x) :
    if x > 0 :
        return 1
    elif x < 0 :
        return -1
    else :
        return 0

'''Saturation functions'''
def Sigma1(u):
    return np.array([sigma(u[0],0,0),sigma(u[1],0,1),sigma(u[2],0,2)])
def Sigma2(u):
    return np.array([sigma(u[0],1,0),sigma(u[1],1,1),sigma(u[2],1,2)])
def sigma(u,i,j):
    if abs(u) < L[i,j] :
        return u
    else :
        t = np.exp( -2*(abs(u) - L[i,j]) )
        return sign(u)*(L[i,j] + (B[i,j] - B[i,j]*t)/( 1 + (2*B[i,j]-1)*t ))
    

