# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 16:10:05 2018

@author: allen
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
plt.style.use('ggplot')


class model(object):
    def __init__(self, k,D=0.1):
        self.k = k
        self.C = []
        self.Q = []
        self.W = []
        self.D = D
        self.A_c = []
        self.L = []
    def add_elements(self,c,Q,W,A_c=0,delx=20 ,L=100):
        self.C=np.array(c)
        self.Q=np.array(Q)
        self.W=np.array(W)
        self.L =np.array(L)        
        self.A_c=np.array(A_c)
        self.V_tot = self.A_c*self.L if len(self.A_c)==len(self.L) else self.A_c[:-1]*self.L
        self.delx= delx
        self.V = self.A_c*self.delx
    def solve(self,method='CMSS'):
        # build matrix
    
        if method=='CMSS':
            size = len(self.Q)
            C_new = np.zeros(size)
            for i in range(size):
                C_new[i] = (self.W[i]+(self.C[i]*self.Q[i]))/(self.Q[i]+self.k*self.V_tot[i])
            return C_new
        if method=='IMSS':
            size = len(self.C)-1
            C_new = []
            for i in range(size):
                n = int(self.L[i]/self.delx)
                D_aver = self.D*self.A_c[i]/self.delx
                low_diag = (-self.Q[i]-D_aver)*np.ones(n-1)
                up_diag = -D_aver*np.ones(n-1)
                main_diag = (self.Q[i]+D_aver+D_aver+self.k*self.A_c[i]*self.delx)*np.ones(n-2)
                main_diag = np.hstack((self.C[i]*self.Q[i]+self.Q[i]+self.k*self.A_c[0]*self.delx,main_diag, D_aver+self.k*self.A_c[i]*self.delx+self.C[i+1]*self.Q[i+1]))
                matrix = diags([up_diag,main_diag,low_diag],[1,0,-1])
                matrix += np.diag(low_diag,-1)
                matrix += np.diag(main_diag,0)
                rhs = np.hstack((self.W[i],np.zeros(n-1)))
                print(rhs.shape,matrix.shape)
                with open("results_IMSS.csv",'w') as f:
                    f.write(str(matrix))
                f.close()
                C_new.append(np.linalg.solve(matrix,rhs))
            return C_new
        if method == 'IMSS full':
            #try:
            #    self.Q.shape == self.W.shape+1
            #    self.C.shape == self.Q.shape
            #except:
            #    raise RuntimeError("The shape of Q or C should include boundary condition")
            # C is from the first to the second last
            #  the same for Q W,A
            size = len(self.C[:-1])
            n = [int(self.L[i]/self.delx) for i in range(size)]
            Q_1 = np.repeat(self.Q[:-1],n)
            A_1 = np.repeat(self.A_c,n)
            num = np.cumsum(n)-1
            num = np.hstack((0,num))
            W_1 = np.repeat([0,0,0,0,0],n)
            for i in range(size): 
                W_1[num[i]:num[i+1]:n[i]] =self.W[i]
            D_aver = self.D*A_1[1:]/self.delx
            low_diag = (-Q_1[1:]-D_aver)
            up_diag = -D_aver
            main_diag = Q_1[2:]+D_aver[1:]+D_aver[:-1]+self.k*A_1[1:-1]*self.delx
            #left BC is not correct
            main_diag = np.hstack((self.Q[1]+D_aver[0]+self.k*self.delx*self.A_c[0],main_diag, self.Q[-1]+D_aver[-1]+self.k*self.A_c[-1]*self.delx))
            matrix = diags([up_diag,main_diag,low_diag],[1,0,-1],format='csr')
            rhs_1 = W_1
            rhs_2 = np.zeros(len(W_1))
            rhs_2[0] = self.C[0]*self.Q[0]
            rhs = rhs_1+rhs_2
            C_new = spsolve(matrix,rhs)
            with open("results_full.csv",'w') as f:
                f.write(str(C_new))
            f.close()
            return C_new
    def draw(self):
        plt.figure()
        plt.plot(np.cumsum(self.L),self.solve())
        plt.show()
        
    def residence(self):
        # The resisdence time of a pollutant is the average time that a pollutant spends within a completely mixed system
        residence = np.zeros(len(self.Q))
        for i in range(len(self.Q)):
            residence[i] = self.V_tot[i]/(self.Q[i]+self.k*self.V_tot[i])
        return residence
    
    def resilience(self,t, to, loadtype='stepload'):
        # t is nd-array
        #resilience of a system is measured by how a system recovers after a pollution accident
        lab = np.zeros(len(self.Q))
        res = np.zeros((len(t),len(self.Q)))
        if loadtype == 'stepload':
            for i in range(len(self.Q)):
                for j in range(len(t)):
                    lab[i] = self.Q[i]/self.V_tot[i]+self.k
                    res[j,i] = self.W[i]/lab[i]/self.V_tot[i]*(1-np.exp(-lab[i]*(t[j]-to)))+self.C[i]
            return res
    
    def calibration(self, observed):
        obs = observed
        k = self.k
        def simulated(k):
            ml.k = k
            return ml.solve()
        def objective(k,o):
            #print((simulated(k)-o)**2).sum()
            #print(np.abs(simulated(k)-o)/o)
            return np.sqrt((np.abs(simulated(k)-o)**2).sum())
        bnds = (0, None)
        k = minimize(objective, x0=1e6, args=(obs) , method='Nelder-Mead' ).x
        print(objective(k,observed))
        #print((np.abs(simulated(k)-observed)/observed*100).sum())
        return k
'''
# Example of completely mixed system
ml = model(k=10)
C=[0.01,0.08,12]
Q = [250,155,350]
W = [0.01,0.26,0.36]
L = [0,100,200]
ml.add_elements(C,Q,W,L=L)
ml.solve(method='CMSS')

ml.draw()
ml.calibration(np.array([0.02,0.06,10]))
obs = np.array([0.2,0.6,10])
simulated(0.01)

#------------------------------------------------------------------------------
# Example of incompletely mixed system
Q = np.array([1000,1000])
D = np.array([2000])
C = np.array([1,1])
delx = 20
A_c = np.array([10])
L = np.array([100])
k=2
ml = model(k=k,D=D)
W = np.array([1000])
ml.add_elements(C,Q,W,A_c=A_c,delx=delx,L=L)
ml.solve(method='IMSS')
ml.residence(0)
t= np.arange(0,100,1)
ml.resilience(t,0,0)
'''




#------------------------------------------------------------------------------
#                             Read data
data = pd.read_excel("water quality.xlsx").iloc[:30,:]
# Split the first three year data as calibration data and the remaining as testing data
# unit for length is dm
# unit for time is second
# unit for weight is mg
C_COD = np.asarray(data['COD(mg/L)']).reshape(6,5)
C_NH = np.asarray(data['NH4(mg/L)']).reshape(6,5)
A = np.asarray(data['Cross-sectional Area(m^2)']).reshape(6,5)*10**(2)
C_COD_training = np.zeros((5,3))
C_COD_testing = np.zeros((5,2))
A_training = np.zeros((5,3))
C_NH_training = np.zeros((5,3))
C_NH_testing = np.zeros((5,2))
A_testing = np.zeros((5,2))
for i in range(5):
    C_COD_training[i,:] = C_COD[i,0:3]
    C_COD_testing[i,:] = C_COD[i,3:]
    C_NH_training[i,:] = C_NH[i,0:3]
    C_NH_testing[i,:] = C_NH[i,3:]
    A_training[i,:] = A[i,0:3]
    A_testing[i,:] = A[i,3:]
    
Q = np.asarray(data['Q(m^3/s)'])[::5]*10**(3)
L = np.asarray(data['Length(km)'])[::5][:-1]*1000*10
W_COD = np.array([4.632,4.849,4.537,5.362,5.221])*10**(6)
W_NH = np.array([0.786,0.863,0.842,1.112,0.965])*10**(6)
#-----------------------------------------------------------------------------
#                            Incompletely mixed system
c_out_training = []
c_out_testing = []
k=[]
Q = Q[:-1]
residence = []
#     COD
#    Training 
for i in range(len(C_COD_training[0,:])):
    ml = model(k=0.001)
    C = C_COD_training[:,i]
    W = W_COD
    A = A_training[:,i]
    ml.add_elements(C,Q,W,A_c=A,delx=100,L=L)
    #ml.solve()
    k.append(ml.calibration(C_COD[1:,i]))
    c_out_training.append(ml.solve(method='CMSS'))
k = np.asarray(k)
k_mean_COD = k.mean()
c_out_training = []
for i in range(len(C_COD_training[0,:])):
    ml = model(k=k_mean_COD)
    C = C_COD_training[:,i]
    W = W_COD
    A = A_training[:,i]
    ml.add_elements(C,Q,W,A_c=A,delx=100,L=L)
    c_out_training.append(ml.solve(method='CMSS'))
    residence.append(ml.residence())
#   Testing
for i in range(len(C_COD_testing[0,:])):
    ml = model(k=k_mean_COD)
    C = C_COD_testing[:,i]
    W = W_COD
    A = A_testing[:,i]
    ml.add_elements(C,Q,W,A_c=A,delx=100,L=L)
    residence.append(ml.residence())
    c_out_testing.append(ml.solve(method='CMSS'))
#       NH4
c_out_training = []
c_out_testing = []
k=[]
residence = []
#    Training 
C_NH_training[4,2] = 0.351987
C_NH[4,2] = 0.351987
for i in range(len(C_NH_training[0,:])):
    ml = model(k=1e-7)
    C = C_NH_training[:,i]
    W = W_NH
    A = A_training[:,i]
    ml.add_elements(C,Q,W,A_c=A,delx=100,L=L)
    ml.solve()
    k.append(ml.calibration(C_NH[1:,i]))
    ml = model(k=k[i])
    ml.add_elements(C,Q,W,A_c=A,delx=100,L=L)
    ml.solve()
    c_out_training.append(ml.solve(method='CMSS'))
k = np.asarray(k)
k_mean_NH = k[2]
for i in range(len(C_COD_training[0,:])):
    ml = model(k=k_mean_NH)
    C = C_NH_training[:,i]
    W = W_NH
    A = A_training[:,i]
    ml.add_elements(C,Q,W,A_c=A,delx=100,L=L)
    residence.append(ml.residence())
#   Testing
for i in range(len(C_NH_testing[0,:])):
    ml = model(k=k_mean_NH)
    C = C_NH_testing[:,i]
    W = W_NH
    A = A_testing[:,i]
    ml.add_elements(C,Q,W,A_c=A,delx=100,L=L)
    ml.solve()
    c_out_testing.append(ml.solve(method='CMSS'))
    residence.append(ml.residence())
#------------------------------------------------------------------------------
    #                      Plot
c_out = np.vstack((np.asarray(c_out_training),np.asarray(c_out_testing)))
plt.figure()
X = np.hstack((0,np.cumsum(L)))
plt.figure()
plt.subplot(231)
plt.title("year 2000")
plt.plot(X[1:],c_out[:,0],linestyle='',marker='o',label='modelled')
plt.plot(X[1:],C_NH[1:,0],linestyle='',marker='o',label='observed')
plt.legend()
plt.xlabel("distance    (m)")
plt.ylabel("concentration   (mg/l)")
plt.subplot(232)
plt.title("year 2001")
plt.plot(X[1:],c_out[:,1],linestyle='',marker='o',label='modelled')
plt.plot(X[1:],C_NH[1:,1],linestyle='',marker='o',label='observed')
plt.legend()
plt.xlabel("distance    (m)")
plt.ylabel("concentration   (mg/l)")
plt.subplot(233)
plt.title("year 2002")
plt.plot(X[1:],c_out[:,2],linestyle='',marker='o',label='modelled')
plt.plot(X[1:],C_NH[1:,2],linestyle='',marker='o',label='observed')
plt.legend()
plt.xlabel("distance    (m)")
plt.ylabel("concentration   (mg/l)")
plt.subplot(234)
plt.title("year 2003")
plt.plot(X[1:],c_out[:,3],linestyle='',marker='o',label='modelled')
plt.plot(X[1:],C_NH[1:,3],linestyle='',marker='o',label='observed')
plt.legend()
plt.xlabel("distance    (m)")
plt.ylabel("concentration   (mg/l)")
plt.subplot(235)
plt.title("year 2004")
plt.plot(X[1:],c_out[:,4],linestyle='',marker='o',label='modelled')
plt.plot(X[1:],C_NH[1:,4],linestyle='',marker='o',label='observed')
plt.legend()
plt.xlabel("distance    (m)")
plt.ylabel("concentration   (mg/l)")
plt.show()

'''
i=0
ml = model(k=0)
C = C_NH_training[:,i]
W = W_NH
A = A_training[:,i]
ml.add_elements(C,Q,W,A_c=A,delx=100,L=L)
(np.abs(ml.solve()-C_NH[1:,i])/ml.solve()*100).sum()
k.append(ml.calibration(C_NH[1:,i]))
ml = model(k=k[i])
ml.add_elements(C,Q,W,A_c=A,delx=100,L=L)
ml.solve()
'''
#------------------------------------------------------------------------------
#                              system recovery time
Q = Q[:-1]
A = A[:-1]
k_mean_COD = 0.0291/24/3600
k_mean_NH = 0.283/24/3600
t= np.arange(0,500,1)*3600    # modify here for COD it takes several seconds while for N it takes hours
C_COD_2002 = C_COD[:,2]
C_NH_2002 = C_NH[:,2]
ml = model(k_mean_COD)
ml.add_elements(C_COD_2002,Q,W_COD,A_c=A[:5,2],delx=100,L=L)
ml = model(k_mean_NH)
ml.add_elements(C_NH_2002,Q,W_NH,A_c=A[:5,2],delx=100,L=L)
res =ml.resilience(t,0)
plt.figure()
plt.subplot(231)
plt.title("Reactor 1")
plt.xlabel('time    (seconds)')
plt.ylabel('concentration NH (mg/l)')   #   modify pollutant here
plt.plot(t,res[:,0])
plt.subplot(232)
plt.plot(t,res[:,1])
plt.title("Reactor 2")
plt.xlabel('time    (seconds)')
plt.ylabel('concentration NH (mg/l)')
plt.subplot(233)
plt.plot(t,res[:,2])
plt.title("Reactor 3")
plt.xlabel('time    (seconds)')
plt.ylabel('concentration NH (mg/l)')
plt.subplot(234)
plt.plot(t,res[:,3])
plt.title("Reactor 4")
plt.xlabel('time    (seconds)')
plt.ylabel('concentration NH (mg/l)')
plt.subplot(235)
plt.plot(t,res[:,4])
plt.title("Reactor 5")
plt.xlabel('time    (seconds)')
plt.ylabel('concentration NH (mg/l)')
plt.show()
#------------------------------------------------------------------------------
#                           InCompletely mixed system full solution
delx=100
n = [int(L[i]/delx) for i in range(len(L))]
num = np.hstack((0,np.cumsum(n)))
tot_L = np.hstack((0,np.cumsum(L)))
X = np.repeat(L,n)
for i in range(len(n)):
    X[num[i]:num[i+1]:1]=np.linspace(tot_L[i],tot_L[i+1],n[i])
ml = model(k=0.0291/24/3600,D= 0.5*10**(2))
C = C_COD[:,4]
W = W_COD[:]
A = A[:-1,4]
Q=Q
ml.add_elements(C,Q,W,A_c=A,delx=100,L=L)
C_COD_IMSS = ml.solve(method='IMSS full')
plt.figure()
plt.plot(X,C_COD_IMSS)
plt.title("Response of COD in Yangtze River at Steady State")
plt.xlabel("distance    (m)")
plt.ylabel("concentration    (mg/l)")
plt.show()
ml = model(k=0.283/24/3600,D= 0.1*10**(2))
C = C_NH[:,4]
W = W_NH[:]
A = A[:-1,4]
Q=Q
ml.add_elements(C,Q,W,A_c=A,delx=100,L=L)
C_NH_IMSS = ml.solve(method='IMSS full')
plt.figure()
plt.plot(X,C_NH_IMSS)
plt.title("Response of Nitrogen in Yangtze River at Steady State")
plt.xlabel("distance    (m)")
plt.ylabel("concentration    (mg/l)")
plt.show()
#                        InCompletely mixed system
C = C_COD[:,0]
W = W_COD[:]


ml = model(k=1e-7)
ml.add_elements(C,Q,W,A_c=A,delx=1000,L=L)
C_COD_IMSS = np.asarray(ml.solve(method='IMSS'))
plt.plot(C_COD_IMSS[0])
#------------------------------------------------------------------------------

