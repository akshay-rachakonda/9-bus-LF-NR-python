#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import cmath


# In[2]:


# bus 
Bus={"Bus":[0,1,2,3,4,5,6,7,8],"Type":[1,2,2,3,3,3,3,3,3],"Vm":[1.04,1.025 ,1.025 ,1,1,1,1,1,1],"vm_ang":[0,0,0,0,0,0,0,0,0],"Pgi":[0,163.0,85.0,0,0,0,0,0,0],"Qgi":[0,0,0,0,0,0,0,0,0],"Pli":[0,0,0,0,125,90,0,100,0],"Qli":[0,0,0,0,50,30,0,35,0],"Qmin":[0,-99,-99,0,0,0,0,0,0],"Qmax":[0,99,99,0,0,0,0,0]}


# In[3]:


# line data
Lines={"From Bus":[0,1,2,3,3,4,5,6,7],"To Bus":[3,6,8,4,5,6,8,7,8],"R":[0,0,0,0.01,0.017,0.032,0.039,0.0085,0.0119],"X":[0.0576,0.0625,0.0586,0.085,0.092,0.161,0.17,0.072,0.1008],"B/2":[0,0,0,0.1760/2,0.1580/2,0.3060/2,0.3580/2,0.1490/2,0.2090/2],"Tap":[1,1,1,1,1,1,1,1,1]}


# In[4]:


# Y bus formulation
def ybus():
    #from bus number to bus number index
    fb= Lines["From Bus"]
    tb= Lines["To Bus"]
    # resistance, reactnce, gound admittaance
    r=Lines["R"]
    x=Lines["X"]
    b=Lines["B/2"]

    #tap setting
    a=Lines["Tap"]

    # impedence
    z= r+ 1j*np.array(x)
    y=1/z
    b= 1j*np.array(b)

    #no. of buses / nodes
    nb= len(np.array(Bus["Bus"]))

    #no. of branches 
    nl= len(np.array(fb))


    # initialize Y bus
    Y= np.zeros([nb,nb],dtype=np.complex)

# python mein bus or line index 0 se shuru kiye hai taaki loop  mein index 0-8 tak jaye nahi to 9 par error aata hai

    # off diagonal elements of Y bus
    for i in range(nb):
        Y[fb[i]][tb[i]] = Y[fb[i]][tb[i]] - y[i]/a[i]
        Y[tb[i]][fb[i]] = Y[fb[i]][tb[i]]

    # diagonal elements

    for j in range(nb):
        for k in range(nl):
            if fb[k]==j:
                Y[j][j] = Y[j][j] + y[k]/(a[k]*a[k]) + b[k]
            elif tb[k]==j:
                Y[j][j] = Y[j][j] + y[k] + b[k]
    return (Y)


# In[5]:


def pol2rect(rho,theta):
#     print(rho*np.cos(theta))
    rect = rho*np.cos(theta)+ 1j*rho*np.sin(theta)
    return (rect)
#     print(rect)


# In[6]:


def loadflow(V,de,BMva):
    de=de
    V=V
    Y = ybus()
    Vm= pol2rect(V,de)
    fb= Lines["From Bus"]
    tb= Lines["To Bus"]
    b=Lines["B/2"]
    b= 1j*np.array(b)
    a= Lines["Tap"]
    nl= len(fb)            # no of branchs
    Pl= Bus["Pli"]
    Ql= Bus["Qli"]
    nb= len(Vm)
    Iij= np.zeros([nb,nb],dtype=np.complex)
    Sij= np.zeros([nb,nb],dtype=np.complex)
    
    
    I = np.matmul(Y,Vm)    
    Im= abs(I)
    ang_I= np.angle(I)

        # Line current flows
    for m in range(nb):
        for n in range(nl):
            if fb[n] == m:
                p = tb[n]
                Iij[m][p] = -(Vm[m] - Vm[p]*a[n])*Y[m][p]/a[n]**2 + (b[n]/a[n]**2)*Vm[m] # Y(m,n) = -y(m,n)..
                Iij[p][m] = -(Vm[p] - Vm[m]/a[n])*Y[p][m] + b[n]*Vm[p]
            elif tb[n] == m:
                p = fb[n]
                Iij[m][p] = -(Vm[m] - Vm[p]/a[n])*Y[p][m] + b[n]*Vm[m]
                Iij[p][m] = -(Vm[p] - Vm[m])*Y[m][p]/a[n]**2 + (b[n]/a[n]**2)*Vm[p]



    Iijr = np.real(Iij);        
    Iiji = np.imag(Iij);

        # Line power flows
    for m in range(nb):
        for n in range(nb):
            if m !=n:
                Sij[m][n] = Vm[m]*np.conj(Iij[m][n])*BMva
    Pij = np.real(Sij)
    Qij = np.imag(Sij)

        # Line Losses..

    Lij = np.zeros(nl,dtype=np.complex)
    for m in range(nl):
        p = fb[m]
        q = tb[m]
        Lij[m] = Sij[p][q] + Sij[q][p]

    Lpij = np.real(Lij)
    Lqij = np.imag(Lij)

        #Bus power injections

    Si= np.zeros(nb,dtype=np.complex)

    for i in range(nb):
        for k in range(nb):
            Si[i] = Si[i] + np.conj(Vm[i])* Vm[k]*Y[i][k]*BMva

    Pi = np.real(Si)
    Qi = -np.imag(Si)
    Pg = Pi+ Pl
    Qg = Qi+ Ql
    
    return(fb,tb,Pg,Qg)


# In[7]:


# newton raphson method
BMva = 100
Y=ybus()
bus= Bus["Bus"]
typ= Bus["Type"] # Type of Bus 1-Slack, 2-PV, 3-PQ
V= np.array(Bus["Vm"])
de= np.zeros(len(V))
Pg= np.array(Bus["Pgi"])/BMva 
Qg= np.array(Bus["Qgi"])/BMva
Pl = np.array(Bus["Pli"])/BMva       #PLi..
Ql = np.array(Bus["Qli"])/BMva       #QLi..
Qmin = np.array(Bus["Qmin"])/BMva    # Minimum Reactive Power Limit..
Qmax = np.array(Bus["Qmax"])/BMva 
nbus= max(bus)+1
P= Pg-Pl
Q= Qg-Ql
Psp=P
Qsp=Q
G= np.real(Y)
B= np.imag(Y)

pv=[] #index of PV buses
pq=[] #index of PQ buses
for i in range(len(typ)):
    if typ[i]==3:
        pq.append(i)
    else:
        pv.append(i)
pv= np.array(pv)
pq= np.array(pq)

npv = len(pv)                       # Number of PV buses..
npq = len(pq)                       # Number of PQ buses..

Tol=1
Iter=0

#matlab index     #python index
#1                0
#2                1     
#3                2
#4                3
#5                4
#6                5
#7                6
#8                7
#9                8  
while Tol>1e-8:
    
    P = np.zeros(nbus)
    Q = np.zeros(nbus)
    
    # Calculate P and Q
    for i in range(nbus): #  i,k goes from 0-8
        for k in range(nbus):
            P[i] = P[i] + V[i]* V[k]*(G[i][k]*np.cos(de[i]-de[k]) + B[i][k]*np.sin(de[i]-de[k]))
            Q[i] = Q[i] + V[i]* V[k]*(G[i][k]*np.sin(de[i]-de[k]) - B[i][k]*np.cos(de[i]-de[k]))
        
    # checking q limit violations
    if Iter>2 and Iter<=7:
        for n in range(1,nbus):
            if typ[n] ==2:
                QG = Q[n]+Ql[n]
                if QG < Qmin[n]:
                    V[n] = V[n] + 0.01
                elif QG > Qmax[n]:
                    V[n] = V[n] - 0.01
         
    
        
#     # Calculate change from specified value
    dPa = Psp-P
    dQa = Qsp-Q   
    
    dQ = np.zeros(npq)
    k=0
    for i in range(nbus): # dQ is a vector of 6 elements which refers to PQ bus
        if typ[i] == 3:
            dQ[k] = dQa[i]             
            k = k+1
    
#     # Change in active power of all buses except slack bus
    dP = dPa[1:nbus]
    
    # from change in active power is all bus except slack and reactive power change in PQ bus we 
    # create a mismatch vector
    
    #dP has 8 elements dQ has 6 elements M must have 14 elements
    
    M= np.zeros(len(dP)+len(dQ))
    for i in range(len(M)):
        if i<len(dP):
            M[i]=dP[i]
        else:
            M[i]=dQ[i-len(dP)]
            
    # Jacobian Matrix initialize
    J1 = np.zeros([nbus-1,nbus-1])
    J2 = np.zeros([nbus-1,npq])
    J3 = np.zeros([npq,nbus-1])
    J4 = np.zeros([npq,npq])
    
    #J1 - Derivative of Real Power Injections with Angles
    for i in range(len(J1)): # i 0-7, m 1-8, k 0-7,n1=1-8, n2= 0-8 
        m = i+1
        for k in range(len(J1)):
            n=k+1
            if n==m:
                for n in range(nbus): # diagonal
                    J1[i][k] = J1[i][k] + V[m]* V[n]*(-G[m][n]*np.sin(de[m]-de[n]) + B[m][n]*np.cos(de[m]-de[n]))
                J1[i][k] = J1[i][k] - V[m]**2*B[m][m]
            else:
                J1[i][k] = V[m]* V[n]*(G[m][n]*np.sin(de[m]-de[n]) - B[m][n]*np.cos(de[m]-de[n]))
                
                
    #J2 - Derivative of Real Power Injections with V..
    
    for i in range(len(J2)):
        m = i+1        
        for k in range(npq):
            n = pq[k]

            if n == m:
                for n in range(nbus):
                    J2[i][k] = J2[i][k] + V[n]*(G[m][n]*np.cos(de[m]-de[n]) + (B[m][n])*np.sin(de[m]-de[n]))
                J2[i][k] = J2[i][k] + V[m]*G[m][m]
            else:
                J2[i][k] = V[m]*(G[m][n]*np.cos(de[m]-de[n]) + B[m][n]*np.sin(de[m]-de[n]))
    
    #J3 - Derivative of Reactive Power Injections with Angles..
                
    for i in range(npq):
        m = pq[i]
        for k in range(nbus-1):
            n = k+1
            if n == m:
                for n in range(nbus):
                    J3[i][k] = J3[i][k] + V[m]* V[n]*(G[m][n]*np.cos(de[m]-de[n]) + B[m][n]*np.sin(de[m]-de[n]))
                
                J3[i][k] = J3[i][k] - V[m]**2*G[m][m]
            else:
                J3[i][k] = V[m]* V[n]*(-G[m][n]*np.cos(de[m]-de[n]) - B[m][n]*np.sin(de[m]-de[n]))     
        
    
        #J4 - Derivative of Reactive Power Injections with V..
    for i in range(npq):
        m = pq[i]
        for k in range(npq):
            n = pq[k]
            if n == m:
                for n in range(nbus):
                    J4[i][k] = J4[i][k] + V[n]*(G[m][n]*np.sin(de[m]-de[n]) - B[m][n]*np.cos(de[m]-de[n]))
                
                J4[i][k] = J4[i][k] - V[m]*B[m][m]
            else:
                J4[i][k]= V[m]*(G[m][n]*np.sin(de[m]-de[n]) - B[m][n]*np.cos(de[m]-de[n]))
                
                
    # Construct J matrix of size 14 by 14
    
    J = np.zeros([len(J1)+len(J3),len(J1)+len(J3)])
    
    for i in range(len(J)):
        if i< len(J1):
            for j in range(len(J)):
                if j < len(J1):
                    J[i][j] = J1[i][j]
                else:
                    J[i][j] = J2[i][j-len(J1)]
        else:
            for j in range(len(J)):
                if j< np.shape(J3)[1]:
                    J[i][j] = J3[i-len(J1)][j]
                else:
                    J[i][j] = J4[i-len(J1)][j-np.shape(J3)[1]]
                    
                    
    # J, M is correct
    
    #-----------------------------------------------------------------------------------------------------------
    #--------------------------------- ------------------------------------ --------------------- --------------
    #-----------------------------------------------------------------------------------------------------------
                    
    X = np.linalg.solve(J,M)
    
    dTh = X[0:nbus-1]
    dV = X[nbus-1:len(X)]
    
    # Update State Vectors (Voltage Angle & Magnitude)
    de[1:nbus] = dTh + de[1:nbus]
    k = 0
    for i in range(1,nbus): 
        if typ[i] ==3:
            V[i] = dV[k] + V[i]
            k = k+1
    
    Iter=Iter+1
    Tol= max(abs(M))

[fb,tb,Pij,Qij]= loadflow(V,de,BMva)
print(Pij)

