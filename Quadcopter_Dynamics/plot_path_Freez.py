
import numpy as np
import sympy as smp
import matplotlib.pyplot as plt


# Read Data
data= np.load('Result_Path.npy')
xxf= data[-1][3:6]
xxi= data[ 0][3:6] 

ln=(xxf-xxi)
m_l= max(abs(ln))

lim_x_l, lim_x_r= min(xxi[0],xxf[0])-m_l, min(xxi[0],xxf[0])+ m_l
lim_y_l, lim_y_r= min(xxi[1],xxf[1])-m_l, min(xxi[1],xxf[1])+ m_l
lim_z_l, lim_z_r= min(xxi[2],xxf[2])-m_l, min(xxi[2],xxf[2])+ m_l


fig = plt.figure(0)

P_vec=np.empty((0,3), float)
X_vec=np.empty((0,3), float)
Y_vec=np.empty((0,3), float)
Z_vec=np.empty((0,3), float)
for i in range(len(data)):

    s1, s2, s3 = smp.sin(data[i][0]), smp.sin(data[i][1]), smp.sin(data[i][2])
    c1, c2, c3 = smp.cos(data[i][0]), smp.cos(data[i][1]), smp.cos(data[i][2])

    r11= c3*c2
    r12= c3*s2*s1-s3*c1
    r13= c3*s2*c1+s3*s1 
    r21= s3*c2
    r22= s3*s2*s1+c3*c1
    r23= s3*s2*c1-c3*s1
    r31= -s2
    r32= c2*s1
    r33= c2*c1
    R= smp.Matrix([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])
    xx= R@smp.Matrix([2,0,0])
    yy= R@smp.Matrix([0,2,0])
    zz= R@smp.Matrix([0,0,2])

    
    DOMP= np.array([data[i][3],data[i][4],data[i][5]])
    DOMX= np.array([xx[0],xx[1],xx[2]])
    DOMY= np.array([yy[0],yy[1],yy[2]])
    DOMZ= np.array([zz[0],zz[1],zz[2]])

    P_vec= np.append(P_vec, np.array([DOMP]), axis= 0)
    X_vec= np.append(X_vec, np.array([DOMX]), axis= 0)
    Y_vec= np.append(Y_vec, np.array([DOMY]), axis= 0)
    Z_vec= np.append(Z_vec, np.array([DOMZ]), axis= 0)


X_vec, Y_vec, Z_vec= 0.6*X_vec, 0.6*Y_vec, 0.6*Z_vec

ax = fig.add_subplot(111, projection='3d')

ax.quiver(P_vec[:,0], P_vec[:,1], P_vec[:,2], X_vec[:,0], X_vec[:,1], X_vec[:,2], color= 'blue')
ax.quiver(P_vec[:,0], P_vec[:,1], P_vec[:,2], Y_vec[:,0], Y_vec[:,1], Y_vec[:,2], color= 'black')
ax.quiver(P_vec[:,0], P_vec[:,1], P_vec[:,2], Z_vec[:,0], Z_vec[:,1], Z_vec[:,2], color= 'red')
ax.set_xlim3d(left=lim_x_l, right=lim_x_r)
ax.set_ylim3d(bottom=lim_y_l, top=lim_y_r)
ax.set_zlim3d(bottom=lim_z_l, top=lim_z_r)

ax.scatter3D(data[:,6], data[:,7],data[:,8], color = "red", s= 5)
plt.show()










