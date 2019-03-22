from numpy import *
from numpy.linalg import norm, solve
from matplotlib.pyplot import *



# iptrack - interpolate track
#
# SYNTAX
# p=iptrack(filename)
#
# INPUT
# filename: data file containing exported tracking data on the standard
# Tracker export format
#
# mass_A
# t	x	y
# 0.0	-1.0686477620876644	42.80071293284619
# 0.04	-0.714777136706708	42.62727536827738
# ...
#
# OUTPUT
# p=iptrack(filename) returns the coefficients of a polynomial of degree 15
# that is the least square fit to the data y(x). Coefficients are given in
# descending powers.

import numpy as np

def iptrack(filename):
	data=np.loadtxt(filename,skiprows=2)
	return np.polyfit(data[:,1],data[:,2],15)






def trvalues(p,x):
	y=np.polyval(p,x)
	dp=np.polyder(p)
	dydx=np.polyval(dp,x)
	ddp=np.polyder(dp)
	d2ydx2=np.polyval(ddp,x)
	alpha=np.arctan(-dydx)
	R=(1.0+dydx**2)**1.5/d2ydx2
	return [y,dydx,d2ydx2,alpha,R]




# formel: (g*sin a(x) - (kv/m))/(1+ I_o/mr^2) der r er krumningsradius
p = iptrack("red1_data.bin")
def ode(p):
    g = 9.81
    m = 0.026
    k = 0.004
    h = 0.001
    n = 500
    rad_ball = 0.028
    v0=3.755013495E-1
    x0=1.671057728E-2
    vn=v0
    xn=x0
    x_array=[]
    v_array=[]
    x_array.append(xn)
    v_array.append(vn)    
    for i in range(n):
        trval = trvalues(p,xn)
        alpha = trval[3]
        f = (g*sin(alpha)- (k*vn/m))/(1.4)
        vn1=vn+h*vn*f #vn+1
        xn1=xn+h*vn*cos(alpha) #xn+1
        x_array.append(xn1)
        v_array.append(vn1)
        vn=vn1
        xn=xn1

    print("array x: \n {}".format(x_array))
    print("array v: \n {}".format(v_array))
        
ode(p)

def get_data(filename):
    data = np.genfromtxt(filename,skip_header = 2, usecols= (0, 2, 5))
    
    return data     

data = get_data("red3_alldata.bin")

def find_w(data):
    mass = 0.026

    v_1 = data[0][2]
    v_2 = data[9][2]
    y_1 = data[0][1]
    y_2 = data[9][1]
    k_1 = mass*(v_1**2)
    k_2 = mass*(v_2**2)
    potensial_diff = mass*9.81*(y_1-y_2)
    print("v_1: {} , v_2 : {} , y_1 : {}, y_2 {}".format(v_1,v_2,y_1,y_2))
    return k_1 + potensial_diff - k_2 

ydat=[3.75501350e-01,4.05722844e-01,4.89965153e-01,5.86403838e-01,6.60924023e-01,7.28988700e-01,7.81801743e-01,8.55190697e-01,9.54621625e-01,1.00488866e+00,1.04739765e+00]
xdat=[0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12]
print(np.trapz(ydat,x=xdat))
#print(get_data("red3_alldata.bin"))
