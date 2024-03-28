import sympy as smp
import numpy as np
from scipy.integrate import odeint
import math
import matplotlib.pyplot as plt

print("RUNNING .....")
dl_t= 0.1
g= 9.81
m= 0.468 
Ix= 4.856e-3
Iy= 4.856e-3
Iz= 8.801e-3
Ax= 0.25
Ay= 0.25
Az= 0.25

l= 0.225
k= 2.98e-6
b= 1.14e-7
Im= 3.357e-5
vw0,vw1, vw2= 10,1,10

# Defined Functions ##############################################################
# ODE Solver 
def dSdt(S,t,f3,tu1,tu2,tu3):

    t0,t1,t2,x0,x1,x2,t0d,t1d,t2d,x0d,x1d,x2d = S # it must be intial values 
    X = smp.Matrix([x0, x1, x2])
    T = smp.Matrix([t0, t1, t2])
    dX =smp.Matrix([x0d, x1d, x2d])
    dT =smp.Matrix([t0d, t1d, t2d])


    # making EQ
    s1, s2, s3 = smp.sin(t0), smp.sin(t1), smp.sin(t2)
    c1, c2, c3 = smp.cos(t0), smp.cos(t1), smp.cos(t2)

    # C Matrix
    c11= 0
    c12= (Iy-Iz)*(t1d*c1*s1+ t2d*s1*s1*c2)+(Iz-Iy)*(t2d*c1*c1*c2)- Ix*t2d*c2
    c13= (Iz-Iy)*t2d*c1*s1*c2*c2
    c21= (Iz-Iy)*(t1d*c1*s1+t2d*s1*c2)+(Iy-Iz)*t2d*c1*c1*c2+Ix*t2d*c2
    c22= (Iz-Iy)*t0d*c1*s1 
    c23= -Ix*t2d*s2*c2+Iy*t2d*s1*s1*s2*c2+Iz*t2d*c1*c1*s2*c2
    c31= (Iy-Iz)*t2d*c2*c2*s1*c1-Ix*t1d*c2
    c32= (Iz-Iy)*(t1d*c1*s1*s2+t0d*s1*s1*c2)+(Iy-Iz)*t0d*c1*c1*c2+Ix*t2d*s2*c2-Iy*t2d*s1*s1*s2*c2-Iz*t2d*c1*c1*s2*c2
    c33= (Iy-Iz)*t0d*c1*s1*c2*c2-Iy*t1d*s1*s1*c2*s2-Iz*t1d*c1*c1*c2*s2+Ix*t1d*c2*s2 
    C= smp.Matrix([[c11,c12,c13],[c21,c22,c23],[c31,c32,c33]])

    # J matrix -invers
    j11= Ix
    j12= 0
    j13= -Ix*s2
    j21= 0 
    j22= Iy*c1*c1+Iz*s1*s1
    j23= (Iy-Iz)*c1*s1*c2
    j31= -Ix*s2
    j32= (Iy-Iz)*c1*s1*c2
    j33= Ix*s2*s2+Iy*s1*s1*c2*c2+Iz*c1*c1*c2*c2
    J= smp.Matrix([[j11,j12,j13],[j21,j22,j23],[j31,j32,j33]]).inv()

    AA= smp.Matrix([[Ax,0,0],[0,Ay,0],[0,0,Az]])

    DOM00= J@(smp.Matrix([tu1,tu2,tu3])- (C@smp.Matrix([t0d,t1d,t2d])))
    DOM01= (smp.Matrix([0,0,-g]) +(f3/m)*smp.Matrix([c3*s2*c1+s3*s1, s3*s2*c1-c3*s1, c2*c1]) -(1/m)*(AA@smp.Matrix([x0d+vw0,x1d+vw1,x2d+vw2])))
    

    return [S[6],S[7],S[8],S[9],S[10],S[11],DOM00[0], DOM00[1],DOM00[2], DOM01[0],DOM01[1],DOM01[2]]


# Simulate Fligth based ob calculated torques and thrust 
def Dy_Run(f3,tu1,tu2,tu3,st):
    t= np.linspace(0,dl_t,11)  
    ans= odeint(dSdt, st, t=t, args= (f3,tu1,tu2,tu3))
    # print(ans)
    return(ans[-1]) # Return Final anguler Linear Pos/Velocity


# Algebric Trajectory : fit plolynomial (order 3)
def make_trj(xi, xf, vi, vf, t, dlt, A, B, tt): 
    global last1, last2, tf
    
    ti= tf-t
    # print(int(ti/dlt))
    B= (3.0/ti**2)*(xf- ti*vi -xi)-(1.0/ti)*(vf-vi)
    A= (ti**-2)*(vf-vi) -(2/ti**3)*(xf-xi-vi*ti)
    
    t=dlt
    domR= np.array([A*t**3 + B*t**2 + vi*t + xi])
    domV= np.array([3*A*t**2 + 2*B*t + vi])
    domAc= np.array([6*A*t+ 2*B])
    domF= m*domAc + m*g*np.array([0.0, 0.0, 1.0]) + (domV+[vw0,vw1,vw2])* np.matrix([[Ax,0,0],[0,Ay,0],[0,0,Az]])

    domtet= math.atan2(domF[0,0], domF[0,2])
    domphi= math.atan2(-math.cos(domtet)*domF[0,1], domF[0,2])  

    # take care of discontiuti ######
    if domtet < last1-math.pi: domtet += 2*math.pi
    if domtet > last1+math.pi: domtet -= 2*math.pi

    if domphi < last2-math.pi: domphi += 2*math.pi
    if domphi > last2+math.pi: domphi -= 2*math.pi
 
    domAg= np.array([domphi, domtet, 0.0])
    last1, last2= domtet, domphi 


    # Makin by Numeric diff for omega and alfa
    tg1= t+(dlt)
    domRg1= np.array([A*tg1**3 + B*tg1**2 + vi*tg1 + xi])
    domVg1= np.array([3*A*tg1**2 + 2*B*tg1 + vi])
    domAcg1= np.array([6*A*tg1+ 2*B])
    domFg1= m*domAc + m*g*np.array([0.0, 0.0, 1.0]) + (domV+[vw0,vw1,vw2]) * np.matrix([[Ax,0,0],[0,Ay,0],[0,0,Az]])
    domtetg1= math.atan2(domFg1[0,0],domFg1[0,2])
    domphig1= math.atan(-math.cos(domtetg1)*domFg1[0,1]/domFg1[0,2]) 


       
    tg2= t+2*(dlt)
    domRg2= np.array([A*tg2**3 + B*tg2**2 + vi*tg2 + xi])
    domVg2= np.array([3*A*tg2**2 + 2*B*tg2 + vi])
    domAcg2= np.array([6*A*tg2+ 2*B])
    domFg2= m*domAc + m*g*np.array([0.0, 0.0, 1.0]) + (domV+[vw0,vw1,vw2]) * np.matrix([[Ax,0,0],[0,Ay,0],[0,0,Az]])
    domtetg2= math.atan2(domFg2[0,0],domFg2[0,2])
    domphig2= math.atan(-math.cos(domtetg2)*domFg2[0,1]/domFg2[0,2])

  
    domOmg= np.array([(domphig1-domphi)/(tt[1]-tt[0]), (domtetg1-domtet)/(tt[1]-tt[0]), 0.0])
    domAlf= np.array([(domphig2- 2*domphig1+ domphi)/((tt[1]-tt[0])**2), (domtetg2- 2*domtetg1+ domtet)/((tt[1]-tt[0])**2), 0.0])


    return (domR, domV, domAc, domAg, domOmg, domAlf, A, B)






# Main call from here #######################################

def dynamics(xi,vi,xf,vf):
    global tf
    xxi, vvi= xi, vi
    xxf, vvf= xf, vf

    xf= xf+vf*dl_t

    tf=  5+ dl_t

    A= np.zeros(3)
    B= np.zeros(3)

    tt= np.linspace(0,tf,51)
   
    dlt= tt[1]-tt[0]

    R= np.empty((0,3), float)
    V , Ac, tet, tet_p, omg, alf= R, R, R, R, R, R
    cmd= np.empty((0,4), float) # [f3, tou]
    Rot= np.empty((0,4), float) # [w1,w2,w3,w4]
    save_path= np.empty((0,9), float)

    global last1, last2
    last1, last2 = 0.0, 0.0
    ode= np.zeros(12) # s= [t0,t1,t2,x0,x1,x2,t0d,t1d,t2d,x0d,x1d,x2d]
    cnt=0

    for t in tt:
        if cnt== len(tt)-1: break
        domR, domV, domAc, domAg, domOmg, domAlf, A, B= make_trj(xi, xf, vi, vf, t, dlt, A, B, tt) # aLL at time t 


    # Finding torque and Thrust 
        s1, s2, s3 = math.sin(domAg[0]), math.sin(domAg[1]), math.sin(domAg[2])
        c1, c2, c3 = math.cos(domAg[0]), math.cos(domAg[1]), math.cos(domAg[2])
        
        t0d, t1d, t2d= domOmg[0],domOmg[1],domOmg[2]


        # C Matrix (again global issue)
        c11= 0
        c12= (Iy-Iz)*(t1d*c1*s1+ t2d*s1*s1*c2)+(Iz-Iy)*(t2d*c1*c1*c2)- Ix*t2d*c2
        c13= (Iz-Iy)*t2d*c1*s1*c2*c2
        c21= (Iz-Iy)*(t1d*c1*s1+t2d*s1*c2)+(Iy-Iz)*t2d*c1*c1*c2+Ix*t2d*c2
        c22= (Iz-Iy)*t0d*c1*s1 
        c23= -Ix*t2d*s2*c2+Iy*t2d*s1*s1*s2*c2+Iz*t2d*c1*c1*s2*c2
        c31= (Iy-Iz)*t2d*c2*c2*s1*c1-Ix*t1d*c2
        c32= (Iz-Iy)*(t1d*c1*s1*s2+t0d*s1*s1*c2)+(Iy-Iz)*t0d*c1*c1*c2+Ix*t2d*s2*c2-Iy*t2d*s1*s1*s2*c2-Iz*t2d*c1*c1*s2*c2
        c33= (Iy-Iz)*t0d*c1*s1*c2*c2-Iy*t1d*s1*s1*c2*s2-Iz*t1d*c1*c1*c2*s2+Ix*t1d*c2*s2 
        C= smp.Matrix([[c11,c12,c13],[c21,c22,c23],[c31,c32,c33]])

        # J matrix -invers
        j11= Ix
        j12= 0
        j13= -Ix*s2
        j21= 0 
        j22= Iy*c1*c1+Iz*s1*s1
        j23= (Iy-Iz)*c1*s1*c2
        j31= -Ix*s2
        j32= (Iy-Iz)*c1*s1*c2
        j33= Ix*s2*s2+Iy*s1*s1*c2*c2+Iz*c1*c1*c2*c2
        J= smp.Matrix([[j11,j12,j13],[j21,j22,j23],[j31,j32,j33]])

    # Calculate Thrust 
        FF= np.array(m*domAc + m*g*np.array([0.0, 0.0, 1.0]) + (domV+[vw0,vw1,vw2]) * np.matrix([[Ax,0,0],[0,Ay,0],[0,0,Az]])).reshape(-1)
        if s1 != 0.0: 
            f3= -FF[1]/s1
        elif s2 != 0.0 :
            f3= FF[0]/s2  
        else: 
            f3= FF[2] 

        ttou= J@domAlf + C@domOmg
        domcmd= np.array([f3, ttou[0], ttou[1], ttou[2]])
        OK1, OK2, OK3, OK4= f3/(4.0*k), ttou[0]/(2.0*k*l), ttou[1]/(2.0*k*l),ttou[2]/(4.0*b) 
        try: 
            domRot= np.array([math.sqrt((OK1-OK3-OK4)),  math.sqrt((OK1-OK2+OK4)), math.sqrt((OK1+OK3-OK4)), math.sqrt((OK1+OK2+OK4))])
        except: 
            print("************   NOT FEASEBLE  ***********")


    # call fo ODE Solver

        DOMDOM= [ode[0],ode[1],ode[2], ode[3], ode[4], ode[5], domR[0][0], domR[0][1], domR[0][2]]
        save_path= np.append(save_path, np.array([DOMDOM]), axis= 0)

        f3,tu1,tu2,tu3= domcmd[0],domcmd[1],domcmd[2],domcmd[3]
        ode[0:3], ode[3:6], ode[6:9], ode[9:12]= domAg, domR, domOmg,domV
        ode= Dy_Run(f3,tu1,tu2,tu3,ode)
        


        R= np.append(R, np.array(domR), axis= 0)
        V= np.append(V, np.array(domV), axis= 0)
        Ac= np.append(Ac, np.array(domAc), axis= 0)
        tet= np.append(tet, np.array([domAg]), axis= 0)
        omg= np.append(omg, np.array([domOmg]), axis= 0)
        alf= np.append(alf, np.array([domAlf]), axis= 0)
        cmd= np.append(cmd, np.array([domcmd]), axis= 0)
        Rot= np.append(Rot, np.array([domRot]), axis= 0)

        xi= domR.reshape(xi.shape)
        vi= domV.reshape(vi.shape)
        cnt+=1
    ###############################################################
    

    # Save results 
    np.save("Result_Path", save_path)


    # plot Energy 
    Eng =np.array( [np.sum((Rot*Rot*Rot)[i]) for i in range(len(Rot))] ) 
    Tot_Eng= np.sum(Eng)
    print("Energy Consumption Estimation: ",Tot_Eng, "[J]")






    fig = plt.figure(figsize=(10, 7))
    ax = plt.axes(projection="3d")

    # Creating plot
    if np.dot(vvi, vvi) != 0 and np.dot(vvf, vvf) != 0:
        vvi = vvi / (math.sqrt(np.dot(vvi, vvi)))
        vvf = vvf / (math.sqrt(np.dot(vvf, vvf)))
    ln=(xxf-xxi)
    m_l= max(abs(ln))
    lim_x_l, lim_x_r= min(xxi[0],xxf[0])-m_l, min(xxi[0],xxf[0])+ m_l
    lim_y_l, lim_y_r= min(xxi[1],xxf[1])-m_l, min(xxi[1],xxf[1])+ m_l
    lim_z_l, lim_z_r= min(xxi[2],xxf[2])-m_l, min(xxi[2],xxf[2])+ m_l

    ax.scatter3D(R[:-2, 0], R[:-2, 1], R[:-2, 2], color="green")
    ax.quiver(xxf[0], xxf[1], xxf[2], vvf[0], vvf[1], vvf[2], color='black')
    ax.quiver(xxi[0], xxi[1], xxi[2], vvi[0], vvi[1], vvi[2], color='black')
    ax.set_xlim3d(left=lim_x_l, right=lim_x_r)
    ax.set_ylim3d(bottom=lim_y_l, top=lim_y_r)
    ax.set_zlim3d(bottom=lim_z_l, top=lim_z_r)
    plt.title("3D Trajectory (black arrows show initial and final velocity directions)")
    plt.show()






    plt.plot(tt[:-1],tet[0:len(tt),0]*57.3)
    plt.plot(tt[:-1],tet[0:len(tt),1]*57.3)
    plt.plot(tt[:-1],tet[0:len(tt),2]*57.3)
    plt.legend([" Roll", " Pitch", " Yaw"])
    plt.xlabel('Time [S]')
    plt.ylabel('Angel [Degree]')
    plt.show()


    # plt.plot(tt[:-1],Rot[0:len(tt),0]*9.55)
    # plt.plot(tt[:-1],Rot[0:len(tt),1]*9.55)
    # plt.plot(tt[:-1],Rot[0:len(tt),2]*9.55)
    # plt.plot(tt[:-1],Rot[0:len(tt),3]*9.55)
    # plt.title(" rpm vrs Time")
    # plt.show()


    plt.plot(tt[:-1],cmd[0:len(tt),0])
    plt.title(" Thrust vrs Time")
    plt.xlabel('Time [S]')
    plt.ylabel('Force [N]')
    plt.show()


    plt.plot(tt[:-1],cmd[0:len(tt),1])
    plt.plot(tt[:-1],cmd[0:len(tt),3])
    plt.plot(tt[:-1],cmd[0:len(tt),2])
    plt.title(" Torque vrs Time")
    plt.legend(["Causes Roll", " Causes Pitch", " Causes Yaw"])
    plt.xlabel('Time [S]')
    plt.ylabel('Torque [N.m]')
    plt.show()

    return 

