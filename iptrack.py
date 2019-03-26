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
    data = np.loadtxt(filename, skiprows=2)
    return np.polyfit(data[:, 1], data[:, 2], 15)


def trvalues(p, x):
    y = np.polyval(p, x)
    dp = np.polyder(p)
    dydx = np.polyval(dp, x)
    ddp = np.polyder(dp)
    d2ydx2 = np.polyval(ddp, x)
    alpha = np.arctan(-dydx)
    R = (1.0+dydx**2)**1.5/d2ydx2
    return [y, dydx, d2ydx2, alpha, R]


# formel: (g*sin a(x) - (kv/m))/(1+ I_o/mr^2) der r er krumningsradius
# p = iptrack("red1_data.bin")


def ode(p):
    g = 9.81
    m = 0.026
    k = 0.018
    h = 0.001
    n = 500
    rad_ball = 0.028
    v0 = 3.755013495E-1
    x0 = 1.671057728E-2
    vn = v0
    xn = x0
    x_array = []
    v_array = []
    x_array.append(xn)
    v_array.append(vn)
    for i in range(n):
        trval = trvalues(p, xn)
        alpha = trval[3]
        f = (g*sin(alpha) - (k*vn/m))/(1.4)
        vn1 = vn+h*vn*f  # vn+1
        xn1 = xn+h*vn*cos(alpha)  # xn+1
        x_array.append(xn1)
        v_array.append(vn1)
        vn = vn1
        xn = xn1

    print("array x: \n {}".format(x_array))
    print("array v: \n {}".format(v_array))


# ode(p)


def get_data(filename):
    data = np.genfromtxt(filename, skip_header=2, usecols=(0, 2, 5))

    return data


def find_w(data, start, end):
    mass = 0.026
    v_1 = data[start][2]
    v_2 = data[end][2]
    y_1 = data[start][1]
    y_2 = data[end][1]
    k_1 = 7/10*mass*(v_1**2)
    k_2 = 7/10*mass*(v_2**2)
    potensial_diff = mass*9.81*(y_1-y_2)
    # print("v_1: {} , v_2 : {} , y_1 : {}, y_2 {}".format(v_1, v_2, y_1, y_2))
    return k_1 + potensial_diff - k_2


def trapezoid(data, start, end):
    func = []
    xvalues = []
    for i in range(end+1):
        trueval = data[i+start][2]
        func.append(trueval**2)
        xvalues.append(data[i+start][0])
    return np.trapz(func, x=xvalues)


def find_k(data):
    return find_w(data, 0, 30)/trapezoid(data, 0, 30)


def find_avg(datalist):
    avg = 0
    for e in datalist:
        avg = avg + e
    return avg/3


def find_standard(datalist):
    standard = 0
    avg = find_avg(datalist)
    for e in datalist:
        standard = standard + (e-avg)**2
    standard = np.sqrt(1/2*standard)
    return standard


def find_usikkerhet(datalist):
    return find_standard(datalist)/sqrt(3)


data_red = get_data("red2_alldata.bin")
data_blue = get_data("blue1_alldata.bin")
data_green = get_data("green1_alldata.bin")

blue_k = find_k(data_blue)
green_k = find_k(data_green)
red_k = find_k(data_red)
k_list = [blue_k, green_k, red_k]
print("blue k: " + str(blue_k))
print("green k: " + str(green_k))
print("red k: " + str(red_k))
# print((blue_k+green_k+red_k)/3)
print("avg k: " + str(find_avg(k_list)))
print("usikkerhet k: " + str(find_usikkerhet(k_list)))
# avg_k = blue_k+green_k+red_k)/3
#standard_k = np.sqrt(1/2*((blue_k-avg_k)^2+(green_k-avg_k)^2+(red_k-avg_k)^2))
# usikkerhet
vred_list = [1.877434970, 1.709297180, 1.819931445]
print("avg v red: " + str(find_avg(vred_list)))
print("usikkerhet v red: " + str(find_usikkerhet(vred_list)))

vblue_list = [1.586149705E0, 1.609394185E0, 1.468681890E0]
print("avg v blue: " + str(find_avg(vblue_list)))
print("usikkerhet v blue: " + str(find_usikkerhet(vblue_list)))

vgreen_list = [1.423989411E0, 1.479436028E0, 1.401244473E0]
print("avg v green: " + str(find_avg(vgreen_list)))
print("usikkerhet v green: " + str(find_usikkerhet(vgreen_list)))
