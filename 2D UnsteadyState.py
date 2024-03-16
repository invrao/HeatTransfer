import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from openpyxl import Workbook
# Assigning values to variables
#=============================#
l = 10                  # length of billet
hL = 1                  # Left side convective coefficient1
hR = 1                  # Right side convective coefficient2
hT = 1                  # Top side convective coefficient2
hB = 1                  # Bottom side convective coefficient2
Ta = 30                 #Atmospheric temperature on left
Ts = 500                # Initial Surface temperature
p = 10                  # no. of nodes, the billet divided into, in X Direction
m = 10                  # no. of nodes, the billet divided into, in Y Direction
delx = l / (p - 1)      # change in distance X direction (distance between tow nodes) 
dely = l / (m - 1)      # change in distance Y Direction (distance between tow nodes)
delt = 0.02             # change in time (in sec)
rho = 1                 # density
total_time = 2.0        # total time
Tstep = int(total_time / delt)  # timestep/Iteration count
c=0.1;                 #carbon % for solidus and liquidus temps
#T_Solidus=1493;         #solidus temperature for carbon >0.1 and <0.18
#T_Liquidus=1534-80.39*c;#liquidus temperature for carbon >0.1 and <0.18

if (c<0.1) :
    T_Solidus = 1534-410 *c
elif (c<0.18 and c >=0.1) :
    T_Solidus = 1493
else:
    T_Solidus=1527-189.07*c

if (c<0.51) :
    T_Liquidus = 1534-80.39*c
else:
    T_Liquidus = 1540.82 -93.77*c


# Assigning initial values to arrays
#==================================#
Subdag = np.zeros(m*p) #Creating sub diagonal with size n
Dag = np.zeros(m*p)    #Creating diagonal with size n
Supdag = np.zeros(m*p) #Creating super diagonal with size n
D = np.zeros(m*p)      #Creating solution array D with size n
T = np.ones([m,p])*Ts  #Creating T(temperature) array for copying new values of T
k = np.zeros([m,p])       #thermal conductivity
cp=np.zeros([m,p]);      #specific heat
data = np.zeros((int(total_time / delt) + 1, (m*p) + 1))  # Add 1 for time column
data[0, 1:] = Ts  # Initial temperatures at 0 sec

#Creating a custom function "thomas algorithm"/ TRIDAG METHOD
#===========================================================#
#satrt of custom function
def thomas_algorithm(Subdag, Dag, Supdag, D):
    n = len(D)
    # Forward elimination ()
    #----------------------#
    for a in range(1, n):
        Subdag[a] = Subdag[a] / Dag[a - 1]
        Dag[a] = Dag[a] - (Subdag[a] * Supdag[a - 1])
        D[a] = D[a] - (Subdag[a] * D[a - 1])

    # Backward substitution
    #---------------------#
    X = np.zeros(n)
    X[-1] = D[-1] / Dag[-1]
    for i in range(n - 2, -1, -1):
        X[i] = (D[i] - Supdag[i] * X[i + 1]) / Dag[i]
    return X
#End of custom function

#Function for Calculating Thermal Conductiivty as fn of temperature
#==================================================================#
def Th_Conductivity():
    for i in range(m):
        for j in range(p):
            if T[i,j] < 801:
                k[i,j] = 59.4 - 0.0418 * T[i,j]
            elif 801 <= T[i,j] <= T_Solidus:
                k[i,j] = 18.4 + 0.0094 * T[i,j]
            else:
                THK = 18.4 + 0.0094 * T[i,j]
                k[i,j] = THK + (43 - THK) * (T[i,j] - T_Solidus) / (T_Liquidus - T_Solidus)
        
#End of Thermal conductivity custom function


#Function for Calculating Specific heat for each node
#====================================================#
def Specific_heat():
    for i in range(m):
        for j in range(p):
            if T[i,j] < 114.3:
                cp[i,j] = 0.49897
            elif 114.3 <= T[i,j] < 491.4:
                cp[i,j] = 0.456 + 2 * 0.000188 * T[i,j]
            elif 491.4 <= T[i,j] < 697.1:
                cp[i,j] = 0.268 + 2 * 0.000418 * T[i,j]
            elif 697.1 <= T[i,j] < 742.9:
                cp[i,j] = 1.431
            elif 742.9 <= T[i,j] < 868.6:
                cp[i,j] = 3.849 - 2 * 0.001883 * T[i,j]
            elif 868.6 <= T[i,j] < 1142.9:
                cp[i,j] = 0.648
            elif 1142.9 <= T[i,j] < T_Solidus:
                cp[i,j] = 0.268 + 2 * 0.000167 * T[i,j]
            elif T_Solidus <= T[i,j] < T_Liquidus:
                cp[i,j] = 0.268 + 2 * 0.000167 * T[i,j] + 272 / (T_Liquidus - T_Solidus)
            else:
                cp[i,j] = 0.787
        
#End of Specific Heat custom function

#=============================#
#For Horizontal Implicit Method
#=============================#
#Coner Nodes
def H_corner_nodes(delt):
    #Top left Node
    #define kR & kB
    i=0
    j=0
    kR = 2*(k[i,j]*k[i,j+1])/(k[i,j]+k[i,j+1])
    kB = 2*(k[i,j]*k[i+1,j])/(k[i,j]+k[i+1,j])
    
    Subdag[0] = 0
    Dag[0] = -(1)-((2*delt*hL)/(rho*cp[0,0]*delx)) - ((2*delt*kR)/(rho*cp[0,0]*delx*delx)) 
    Supdag[0] = (2*kR*delt) /(rho*cp[0,0]*delx*delx)
    D[0] =  - ((2*hT*delt)/(rho*cp[0,0]*dely))*Ta - ((2*hL*delt)/(rho*cp[0,0]*delx))*Ta - (((2*kB*delt))/(rho*cp[0,0]*dely*dely))*T[i+1,j] + ((2*delt*hT)/(rho*cp[0,0]*dely))*T[i,j] + ((2*delt*kB)/(rho*cp[0,0]*dely*dely))*T[i,j] - T[i,j]
    
    
    #Left Bottom Node
    n_bl = ((m-1)*p)
    kT = 2*(k[m-1,0]*k[m-2,0])/(k[m-1,0]+k[m-2,0])
    kR = 2*(k[m-1,0]*k[m-1,1])/(k[m-1,0]+k[m-1,1])

    Subdag[n_bl] = 0
    Dag[n_bl] = - ((2*hL*delt)/(rho*cp[m-1,0]*delx)) - ((2*kR*delt)/(rho*cp[m-1,0]*delx*delx))  -1
    Supdag[n_bl] =  (2*kR*delt)/(rho*cp[m-1,0]*delx*delx)
    D[n_bl] =   -(((2*hL*delt)/(rho*cp[m-1,0]*delx)) + (2*hB*delt)/(rho*cp[m-1,0]*dely))*Ta - (((2*kT*delt)/(rho*cp[m-1,0] *dely*dely))*T[m-2,0]) - (-  ((2*kT*delt)/(rho*cp[m-1,0]*dely*dely)) - ((2*hB*delt)/(rho*cp[m-1,0]*dely)))*T[m-1,0] -T[m-1,0]
    
    
    
    #Top Right Node
    i = 0
    j = p-1
    n_tr = p-1 #Top Rowlast node number, begining with 0
    kB = 2*(k[0,p-1]*k[1,p-1])/(k[0,p-1]+k[1,p-1])
    kL = 2*(k[0,p-1]*k[0,p-2])/(k[0,p-1]+k[0,p-2])
    
    Subdag[n_tr] =  (2*kL*delt)/(rho*cp[i,j]*delx*delx)
    Dag[n_tr] = -(1) - ((2*hR*delt)/(rho*cp[i,j]*delx)) - ((2*kL*delt)/(rho*cp[i,j]*delx*delx))  
    Supdag[n_tr] = 0
    D[n_tr] =   -(((2*hR*delt)/(rho*cp[i,j]*delx) + (2*hT*delt)/(rho*cp[i,j]*dely))*Ta + ((2*kB*delt)/(rho*cp[i,j] *dely*dely))*T[i+1,j]) +((2*kB*delt)/(rho*cp[i,j]*dely*dely))*T[i,j] +((2*hT*delt)/(rho*cp[i,j]*dely))*T[i,j] -( T[i,j])
    
    
    #Bottom Right Node
    i = m-1
    j = p-1
    n_br = (m*p)-1
    kT = 2*(k[m-1,p-1]*k[m-2,p-1])/(k[m-1,p-1]+k[m-2,p-1])
    kL = 2*(k[m-1,p-1]*k[m-1,p-2])/(k[m-1,p-1]+k[m-1,p-2])
    Subdag[n_br] = ((2*kL*delt) /(rho*cp[i,j]*delx*delx))
    Dag[n_br] = -((2 *hR* delt)/(rho*cp[i,j]*delx) ) -((2*kL*delt)/(rho*cp[i,j]*delx*delx))-1
    Supdag[n_br] = 0
    D[n_br] =   ((-2*hR*delt)/(rho*cp[i,j]*delx) + (-2*hB*delt)/(rho*cp[i,j]*dely))*Ta + ((-2*kT*delt)/(rho*cp[i,j] *dely*dely))*T[i-1,j] + (2*kT*delt/(rho*cp[i,j]*dely*dely))*T[i,j] + (2*hB*delt/(rho*cp[i,j]*dely))*T[i,j] - T[i,j]
    #D[n_br] = -((2*hB*Ta*delt)/(rho*cp[i,j]*dely)) + ((2*hR*Ta*delt)/(rho*cp[i,j]*delx)) + ((2*hB*delt)/(rho*cp[i,j]*dely))*T[i,j] + ((2*kT*delt)/(rho*cp[i,j]*dely*dely))*T[i,j] -((2*kT*delt)/(rho*cp[i,j]*delx*delx))*T[i,j] -T[i,j]
    #D[n_br] =  - (((2 *Ta* delt)/(rho*cp[m-1,p-1])) * ((hR/delx) + (hB/dely))) - ((2*kT*delt)/(rho*cp[m-1,p-1]*dely*dely))*T[m-2,p-1] + ((2 * delt*hB)/(delx*rho*cp[m-1,p-1]) )*T[m-1,p-1] + (2*kT * delt/(rho*cp[m-1,p-1])*dely*dely)*T[m-1,p-1] -T[m-1,p-1]
    
#Top Surface Node Function
def H_Top_Surface(delt):
    for j in range(1,p-1):
        #kR and kL
        i = 0
        kR = 2*(k[i,j]*k[i,j+1])/(k[i,j]+k[i,j+1])
        kL = 2*(k[i,j]*k[i,j-1])/(k[i,j]+k[i,j-1])
        kB = 2*(k[i,j]*k[i+1,j])/(k[i,j]+k[i+1,j])
        
        Subdag[j] = kL*delt / (rho*cp[i,j]*delx*delx)
        Dag[j] = -1- ((kL*delt) / (rho*cp[i,j]*delx*delx)) - ((kR*delt) / (rho*cp[i,j]*delx*delx)) 
        Supdag[j] = ((kR*delt) / (rho*cp[i,j]*delx*delx))
        D[j] =  -((Ta*(2*hT*delt / (rho*cp[i,j]*dely))) + ((2*kB*delt / (rho*cp[i,j]*dely*dely))*T[1,j])) -(1 -( (2*hT*delt) / (rho*cp[i,j]*dely) ) - ((2*kB*delt) / (rho*cp[i,j]*dely*dely)))*T[0,j]
       
def H_Bottom_Surface(delt):
    for j in range(1,(p-1)):
        i = m-1
        kR = 2*(k[i,j]*k[i,j+1])/(k[i,j]+k[i,j+1])
        kL = 2*(k[i,j]*k[i,j-1])/(k[i,j]+k[i,j-1])
        kT = 2*(k[i,j]*k[i-1,j])/(k[i,j]+k[i-1,j])
        n_bs = ((m-1)*p) + j #Node numbers for bottom surface nodes last row is m-1 index
        Subdag[n_bs] = kL*delt / (rho*cp[i,j]*delx*delx)
        Dag[n_bs] = -1- (kL*delt / (rho*cp[i,j]*delx*delx)) - (kR*delt / (rho*cp[i,j]*delx*delx)) 
        Supdag[n_bs] = kR*delt / (rho*cp[i,j]*delx*delx)
        D[n_bs] =  - (Ta*2*hB*delt / (rho*cp[i,j]*dely)) - (T[m-2,j]*2*kT*delt / (rho*cp[i,j]*dely*dely)) -(1-(( 2*hB*delt) / (rho*cp[i,j]*dely)+ (2*kT*delt / (rho*cp[i,j]*dely*dely)) ))* T[m-1,j]
       

def H_Left_Surface(delt):
    for i in range(1,m-1): 
        j=0
        n_ls= j+ (i*p)
        kB = (2*k[i,j]*k[i+1,j])/(k[i,j]+k[i+1,j])
        kR = (2*k[i,j]*k[i,j+1])/(k[i,j]+k[i,j+1])
        kT = (2*k[i,j]*k[i-1,j])/(k[i,j]+k[i-1,j])
        Subdag[n_ls] =0 
        Dag[n_ls]= -(1+((2*kR*delt)/(rho*cp[i,j]*delx*delx))+((2*hL*delt)/(rho*cp[i,j]*delx)))
        Supdag[n_ls] = ((2*kR*delt)/(rho*cp[i,j]*delx*delx))
        D[n_ls]= -(  ((kT*delt)/(rho*cp[i,j]*dely*dely))*T[i-1,j] +((2*hL*delt)/(rho*cp[i,j]*delx))*Ta + ((kB*delt)/(rho*cp[i,j]*dely*dely))*T[i+1,j]) + ((T[i,j]*kB*delt)/(rho*cp[i,j]*dely*dely)) +((T[i,j]*kT*delt)/(rho*cp[i,j]*dely*dely))-T[i,j]


def H_Right_Surface(delt):
    for i in range(1,m-1):
        j = p-1
        n_rs = ((i*p) + (p-1)) #index of last column is p-1
        
        kL = 2*(k[i,j]*k[i,j-1])/(k[i,j]+k[i,j-1])
        kB = 2*(k[i,j]*k[i+1,j])/(k[i,j]+k[i+1,j])
        kT = 2*(k[i,j]*k[i-1,j])/(k[i,j]+k[i-1,j])
        Supdag[n_rs] = 0
        Subdag[n_rs] = (2*kL*delt)/(rho*cp[i,j]*delx*delx)
        Dag[n_rs] = -1 -((2*hR*delt)/(rho*cp[i,j]*delx)) -(2*kL*delt/(rho*cp[i,j]*delx*delx))
        D[n_rs] =  - Ta*((2*hR*delt)/(rho*cp[i,j]*delx)) - T[i-1,p-1]*(kT*delt/(rho*cp[i,j]*dely*dely)) - (kB*delt/(rho*cp[i,j]*dely*dely))*T[i+1,p-1] +(kT*delt*T[i,p-1]/(rho*cp[i,j]*dely*dely)) + (kB*delt*T[i,p-1]/(rho*cp[i,j]*dely*dely))-T[i,p-1]
 
#Centr 
def H_Centre_Nodes(delt):
    for i in range(1,m-1):
        for j in range(1,p-1):
            n_cn= (i*p)+j
            
            kB = 2*(k[i,j]*k[i+1,j])/(k[i,j]+k[i+1,j])
            kT = 2*(k[i,j]*k[i-1,j])/(k[i,j]+k[i-1,j])
            kL = 2*(k[i,j]*k[i,j-1])/(k[i,j]+k[i,j-1])
            kR = 2*(k[i,j]*k[i,j+1])/(k[i,j]+k[i,j+1])
            Subdag[n_cn] = ((kL*delt)/(rho*cp[i,j]*delx*delx))
            Supdag[n_cn] = ((kR*delt)/(rho*cp[i,j]*delx*delx))
            Dag[n_cn] = -1 - ((kL*delt)/(rho*cp[i,j]*delx*delx)) - ((kR*delt)/(rho*cp[i,j]*delx*delx)) 
            D[n_cn] = - (T[i-1,j]*kT*delt/(rho*cp[i,j]*dely*dely)) - (T[i+1,j]*kB*delt/(rho*cp[i,j]*dely*dely)) + ((kT*delt/(rho*cp[i,j]*dely*dely))*T[i,j]) + ((kB*delt/(rho*cp[i,j]*dely*dely))*T[i,j]) -T[i,j] 
           #this centre node last T(i,j) sent to last and homogenity arrived, Why?

#Function to update the temperatures solved with TRIDAG Method for half DeltaT
def H_TUpdate():
    for K in range(m*p):
        q = int(K//p)
        r = int(K%p)
        T[q,r] = X[K]  # assign the value of X[K] to T[row,col]


#=============================#
#For Vertical Implicit Method
#=============================#
#Coner Nodes
def V_corner_nodes(delt):
    #Top left Node
    kR = 2*(k[0,0]*k[0,1])/(k[0,0]+k[0,1])
    kB = 2*(k[0,0]*k[1,0])/(k[0,0]+k[1,0])
    
    Subdag[0] = 0
    Dag[0] = -(1+ (2*kB*delt)/(rho*cp[0,0]*dely*dely)+(2*hT*delt)/(rho*cp[0,0]*dely))
    Supdag[0] = ((2*kB*delt)/(rho*cp[0,0] *dely*dely))
    D[0] =  -((2*hL*delt)/(rho*cp[0,0]*delx) + (2*hT*delt)/(rho*cp[0,0]*dely))*Ta - ((2*kR*delt)/(rho*cp[0,0]*delx*delx))*T[0,1] - T[0,0] + (((2*hL*delt)/(rho*cp[0,0]*delx))*T[0,0]) + ((2*kR*delt)/(rho*cp[0,0]*delx*delx))*T[0,0]

    #Left bottom Node
    n_bl = m-1
    kT = 2*(k[m-1,0]*k[m-2,0])/(k[m-1,0]+k[m-2,0])
    kR = 2*(k[m-1,0]*k[m-1,1])/(k[m-1,0]+k[m-1,1])
    Subdag[n_bl] = (2*kT*delt) /(rho*cp[0,0]*dely*dely)
    Dag[n_bl] = -1-(2 * delt/(rho*cp[0,0]) ) * (( (hB/dely)   + (kT/(dely*dely))))
    Supdag[n_bl] = 0
    D[n_bl] = -Ta*((2 * delt/(rho*cp[0,0]) ) * ((hL/delx) + (hB/dely))) - T[m-1,1]*(2*kR*delt/(rho*cp[0,0]*delx*delx)) +((2*delt*T[m-1,0]/(rho*cp[0,0]) ) * ((hL/delx) +(kR/(delx*delx)))) -T[m-1,0]

    #Top Right Node 
    n_tr = ((p-1)*m)
    kB = 2*(k[0,p-1]*k[1,p-1])/(k[0,p-1]+k[1,p-1])
    kL = 2*(k[0,p-1]*k[0,p-2])/(k[0,p-1]+k[0,p-2])
    Subdag[n_tr] = 0
    Dag[n_tr] = -1 -((2 * delt/(rho*cp[0,p-1]) ) * ( (hT/dely) + (kB/(dely*dely)) ))
    Supdag[n_tr] = (2*kB*delt) /(rho*cp[0,p-1]*dely*dely)
    D[n_tr] = -( Ta*((2 * delt/(rho*cp[0,p-1]) ) * ((hR/delx) + (hT/dely))) + T[0,p-2]*((2*kL*delt/(rho*cp[0,p-1]*delx*delx)))) +(((2 * delt*T[0,p-1])/(rho*cp[0,p-1]) ) * ((hR/delx) +(kL/(delx*delx)))) -T[0,p-1] 

    #Bottom Right Node 
    n_tr = (p*m)-1
    kT = 2*(k[m-1,p-1]*k[m-2,p-1])/(k[m-1,p-1]+k[m-2,p-1])
    kL = 2*(k[m-1,p-1]*k[m-1,p-2])/(k[m-1,p-1]+k[m-1,p-2])
    Subdag[n_tr] = (2*kT*delt) /(rho*cp[m-1,p-1]*dely*dely)
    Dag[n_tr] = -(1+((2 * delt/(rho*cp[m-1,p-1]) ) * ((hB/dely) + (kT/(dely*dely)))))
    Supdag[n_tr] = 0
    D[n_tr] = -( Ta*((2 * delt/(rho*cp[m-1,p-1]) ) * ((hR/delx) + (hB/dely))) + T[m-1,p-2]*((2*kL*delt/(rho*cp[m-1,p-1]*delx*delx)))) +((2 * delt*T[m-1,p-1]/(rho*cp[m-1,p-1]) ) * ((hR/delx)+(kL/(delx*delx)))) - T[m-1,p-1] 


#Top Surface Node Function
def V_Top_Surface(delt):
    for j in range(1,p-1):
        i = 0
        n_ts = ((j*m) + 0) #here row index is 0 for first row
        
        kR = 2*(k[i,j]*k[i,j+1])/(k[i,j]+k[i,j+1])
        kL = 2*(k[i,j]*k[i,j-1])/(k[i,j]+k[i,j-1])
        kB = 2*(k[i,j]*k[i+1,j])/(k[i,j]+k[i+1,j])
        
        Subdag[n_ts] = 0 
        Dag[n_ts] = -1 - (( 2*hT*delt) / (rho*cp[i,j]*dely) ) - ((2*kB*delt) / (rho*cp[i,j]*dely*dely))
        Supdag[n_ts] = ((2*kB*delt)/ (rho*cp[i,j]*dely*dely))
        D[n_ts] = -((Ta*(2*hT*delt )/ (rho*cp[i,j]*dely)) + ((T[0,j+1]*(kR*delt)) / (rho*cp[i,j]*delx*delx)) + (T[0,j-1]*(kL*delt )/ (rho*cp[i,j]*delx*delx))) + ((kL*delt*T[0,j]) / (rho*cp[i,j]*delx*delx)) + ((kR*delt*T[0,j]) / (rho*cp[i,j]*delx*delx)) -T[0,j]


def V_Bottom_Surface(delt):
    for j in range(1,(p-1)):
        i = m-1
        n_bs = ((j*m) + (m-1)) #here row index is m-1 for last row
        
        kR = 2*(k[i,j]*k[i,j+1])/(k[i,j]+k[i,j+1])
        kL = 2*(k[i,j]*k[i,j-1])/(k[i,j]+k[i,j-1])
        kT = 2*(k[i,j]*k[i-1,j])/(k[i,j]+k[i-1,j])
        Subdag[n_bs] = 2*kT*delt / (rho*cp[i,j]*dely*dely)
        Dag[n_bs] = -(1 + ( 2*hB*delt / (rho*cp[i,j]*dely) ) + (2*kT*delt / (rho*cp[i,j]*dely*dely)))
        Supdag[n_bs] = 0
        D[n_bs] = -((Ta*(2*hB*delt) / (rho*cp[i,j]*dely)) + T[i,j-1]*(kL*delt / (rho*cp[i,j]*delx*delx)) + T[i,j+1]*(kR*delt / (rho*cp[i,j]*delx*delx))) + (kL*delt*T[i,j] / (rho*cp[i,j]*delx*delx)) + (kR*delt*T[i,j] / (rho*cp[i,j]*delx*delx)) -(T[i,j])


def V_Left_Surface(delt):
    for i in range(1,m-1):
        j = 0
        
        kR = 2*(k[i,j]*k[i,j+1])/(k[i,j]+k[i,j+1])
        kB = 2*(k[i,j]*k[i+1,j])/(k[i,j]+k[i+1,j])
        kT = 2*(k[i,j]*k[i-1,j])/(k[i,j]+k[i-1,j])
        Subdag[i] = kT*delt/(rho*cp[i,j]*dely*dely)
        Dag[i] = -(1 + (kT*delt/(rho*cp[i,j]*dely*dely)) +(kB*delt/(rho*cp[i,j]*dely*dely)))
        Supdag[i] = kB*delt/(rho*cp[i,j]*dely*dely)
        D[i] = -(Ta*(2*hL*delt/(rho*cp[i,j]*delx)) +T[i,1]*(2*kR*delt/(rho*cp[i,j]*delx*delx))) + (2*hL*delt*T[i,0]/(rho*cp[i,j]*delx)) + (2*kR*delt*T[i,0]/(rho*cp[i,j]*delx*delx)) -T[i,0] 

def V_Right_Surface(delt):
    for i in range(1,m-1): 
        j=p-1
        n_rs= (p-1)*m +i
        #currently working
        kB = (2*k[i,j]*k[i+1,j])/(k[i,j]+k[i+1,j])
        kL= (2*k[i,j-1]*k[i,j])/(k[i,j-1]+k[i,j])
        kT = (2*k[i,j]*k[i-1,j])/(k[i,j]+k[i-1,j])
        Subdag[n_rs] = ((kT*delt)/(rho*cp[i,j]*dely*dely))
        Dag[n_rs]= -1- (kB*delt)/(rho*cp[i,j]*dely*dely)- (kT*delt)/(rho*cp[i,j]*dely*dely)
        Supdag[n_rs] = ((kB*delt)/(rho*cp[i,j]*dely*dely))
        D[n_rs]= -((2*kL*delt)/(rho*cp[i,j]*delx*delx))*T[i,j-1] -Ta*((2*hR*delt/(rho*cp[i,j]*delx))) + ((2*kL*delt*T[i,j])/(rho*cp[i,j]*delx*delx) )+((2*hR*delt*T[i,j])/(rho*cp[i,j]*delx))-T[i,j] 
 
def V_Centre_Nodes(delt):
    for j in range(1,p-1):
        for i in range(1,m-1):
            node_n= (j*m)+i
            
            kB = 2*(k[i,j]*k[i+1,j])/(k[i,j]+k[i+1,j])
            kT = 2*(k[i,j]*k[i-1,j])/(k[i,j]+k[i-1,j])
            kL = 2*(k[i,j]*k[i,j-1])/(k[i,j]+k[i,j-1])
            kR = 2*(k[i,j]*k[i,j+1])/(k[i,j]+k[i,j+1])
            
            Subdag[node_n] = kT*delt/(rho*cp[i,j]*dely*dely)
            Supdag[node_n] = kB*delt/(rho*cp[i,j]*dely*dely)
            Dag[node_n] = -1 - (kT*delt/(rho*cp[i,j]*dely*dely)) - (kB*delt/(rho*cp[i,j]*dely*dely))
            D[node_n] = - T[i,j-1]*(kL*delt/(rho*cp[i,j]*delx*delx)) -T[i,j+1]*(kR*delt/(rho*cp[i,j]*delx*delx)) + (kL*delt*T[i,j]/(rho*cp[i,j]*delx*delx)) + (kR*delt*T[i,j]/(rho*cp[i,j]*delx*delx))-T[i,j] 

#End of custom functions for nodes#
#=================================#

#Function to update the temperatures solved with TRIDAG Method for half DeltaT
def V_TUpdate():
    for K in range(m*p):
        q = int(K//m)
        r = int(K%m)
        T[r,q] = X[K]  # assign the value of X[K] to T[row,col]


#Main Claculations/Program Starts#
#================================#
for f in range (Tstep):
    Th_Conductivity()
    Specific_heat()
    
    #Calling functions for first DetlT/2 (Horizontal Implicit, Vertical Explicit)
    H_corner_nodes(delt/2)
    H_Top_Surface(delt/2)
    H_Bottom_Surface(delt/2)
    H_Left_Surface(delt/2)
    H_Right_Surface(delt/2)
    H_Centre_Nodes(delt/2)
    
    #Solving the problem using TriDiag
    ##################################
    # Thomas algorithm function for solving tridiagonal equation
    X = np.around(thomas_algorithm(Subdag, Dag, Supdag, D),2)
    
    H_TUpdate()
    
    #Calling functions for next DetlT/2 (Vertical Implicit, Horizontal Explicit)
    V_corner_nodes(delt/2)
    V_Top_Surface(delt/2)
    V_Bottom_Surface(delt/2)
    V_Left_Surface(delt/2)
    V_Right_Surface(delt/2)
    V_Centre_Nodes(delt/2)
    
    #Solving the problem using TriDiag
    ##################################
    # Thomas algorithm function for solving tridiagonal equation
    X = np.around(thomas_algorithm(Subdag, Dag, Supdag, D),2)
    
    V_TUpdate()

    V = T.flatten()
        
    #Saving the data to a variable data
    #---------------------------------#
    data[f + 1, 0] = (f + 1) * delt  # Store time values in first column
    
    data[f + 1, 1:] = V  # Store node temperatures in remaining columns
    
# Create a pandas DataFrame
#-------------------------#
df = pd.DataFrame(data, columns=["Time (s)"] + [f'T_{i}_{j}' for i in range(m) for j in range(p)]) 
#Exporting the data to external excel file
#----------------------------------------#
#columns = ['Time'] + [f'T_{i}_{j}' for i in range(m) for j in range(p)]

# Select rows (first two rows and remaining at cumulative 100th row) to write to Excel 
#(We can use custom name in place of "rows_to_write")
rows_to_write = pd.concat([df.head(2), df.iloc[10::10]])


# Specify the full path where you want to save the Excel file
#Here using 'r' as backslashes (\) in the file path being interpreted as escape characters in Python string literals.
# Or we can use double '\\'  2.....
file_path = r'C:\Users\invra\Desktop\2D_Unsteady_Modelling.xlsx'


# Specify the sheet name
sheet_name = 'ADI data'

# Export selected rows to an Excel file at the specified location with the specified sheet name
rows_to_write.to_excel(file_path, sheet_name=sheet_name, index=False)

#Printing the final temperature distribution
print(T)
#print(V)

# Final 2D array of temperatures plotting

plt.figure(figsize=(8, 6))
plt.imshow(T, cmap='hot', interpolation='nearest')
plt.colorbar(label='Temperature')
plt.title('Temperature Distribution at end of time sec')
plt.xlabel('X')
plt.ylabel('Y')
#plt.grid(True)
plt.show()

