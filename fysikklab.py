# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np
import math

c = 2/5  #
m = 29.6  # gram
g = 9.81

xmax = 0.98


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


def iptrack(filename):
    data = np.loadtxt(filename, skiprows=2, usecols=(0, 1, 2))
    return np.polyfit(data[:, 1], data[:, 2], 15)

# trvalues - track values
#
# SYNTAX
# [y,dydx,d2ydx2,alpha,R]=trvalues(p,x)
#
# INPUT
# p: the n+1 coefficients of a polynomial of degree n, given in descending
# order. (For instance the output from p=iptrack(filename).)
# x: ordinate value at which the polynomial is evaluated.
#
# OUTPUT
# [y,dydx,d2ydx2,alpha,R]=trvalues(p,x) returns the value y of the
# polynomial at x, the derivative dydx and the second derivative d2ydx2 in
# that point, as well as the slope alpha(x) and the radius of the
# osculating circle.
# The slope angle alpha is positive for a curve with a negative derivative.
# The sign of the radius of the osculating circle is the same as that of
# the second derivative.


def trvalues(p, x):
    y = np.polyval(p, x)
    dp = np.polyder(p)
    dydx = np.polyval(dp, x)
    ddp = np.polyder(dp)
    d2ydx2 = np.polyval(ddp, x)
    alpha = np.arctan(-dydx)
    R = (1.0+dydx**2)**1.5/d2ydx2
    return [y, dydx, d2ydx2, alpha, R]


def euler(h, n, x0, v0, file):
    m = 0.026
    g = 9.81
    k = 0.0051
    pol = iptrack(file)
    val = trvalues(pol, x0)
    v = [v0]
    x = [x0]
    y = [val[0]]
    f = []
    N = []
    for i in range(1, n):
        val = trvalues(pol, x[i-1])
        a = val[3]
        lf = k*v[i-1]/m
        vn = v[i-1] + h*(g*math.sin(a) - lf)/(1+c)
        xn = x[i-1] + h*v[i-1]*math.cos(a)
        y.append(val[0])
        v.append(vn)
        x.append(xn)
        friksjon = c*m*v[i-1]
        Nn = m*(g*math.cos(a) + vn**2/val[4])
        N.append(Nn)
        #fn = m*(g*math.sin(a) - (g*math.sin(a))/(1+c))
        f.append(friksjon)

    return x, v, y, f, N


def bane(file):
    p = iptrack(file)
    h = 0.01
    n = int(xmax/h)
    x = [h*i for i in range(n)]
    y = []
    for i in range(len(x)):
        y.append(trvalues(p, x[i])[0])
    return x, y


def energi(y, v):
    e = []
    for i in range(len(y)):
        e.append(m*(g*y[i] + 0.5*(1+c)*m*v[i]**2))

    return e


def snitt(liste):
    s = 0
    for i in range(len(liste)):
        s += liste[i]
    return s/len(liste)


def plottap(data):
    x = [i for i in range(len(data))]
    plt.plot(x, data, label='Mekanisk Energi')
    plt.xlabel('Periode')
    plt.ylabel('Energi [J]')
    plt.savefig('Energi.png', dpi=500)
    plt.show()


def plotrate(data):
    f = open(data)
    d = f.readlines()
    dat = []
    for j in d:
        j = j.strip()
        j = j.split('\t')
        dat.append(j)
    f.close()
    for i in range(len(dat)):
        for j in range(len(dat[i])):
            dat[i][j] = float(dat[i][j])

    energi = [[] for i in range(len(dat[0]))]
    for i in range(len(dat)):
        for j in range(len(dat[i])):
            energi[j].append(dat[i][j])
    snitten = [snitt(energi[i]) for i in range(len(energi))]

    rate = [100*(snitten[i]-snitten[i+1])/snitten[i]
            for i in range(len(snitten)-1)]
    y = [snitt(rate) for i in range(len(rate))]
    # plt.plot(x,snitten, label = 'Energi')
    x = [i for i in range(len(energi)-1)]

    font = {'fontname': 'Arial', 'size': '16'}

    plt.plot(x, rate, label='Tapsrate')
    plt.plot(x, y, linestyle='dashed', label='Gjennomsnitt')
    plt.xlabel('Periode', **font)
    plt.ylabel('Prosentvis tapt energi', **font)
    plt.legend(loc=4, fontsize=14)
    plt.tick_params(labelsize=12)
    plt.savefig('Energitap.png', dpi=500)
    plt.show()

    x = [i for i in range(len(energi))]
    plt.plot(x, snitten, label='Mekanisk Energi')
    plt.xlabel('Periode', **font)
    plt.ylabel('Energi [mJ]', **font)
    plt.legend(loc=1, fontsize=14)
    plt.tick_params(labelsize=12)
    plt.savefig('Energi.png', dpi=500)
    plt.show()


def getfile(file):
    f = open(file)
    d = f.readlines()
    dat = []
    for j in d:
        j = j.strip()
        j = j.split('\t')
        dat.append(j)
    f.close()
    for i in range(len(dat)):
        for j in range(len(dat[i])):
            dat[i][j] = float(dat[i][j])

    energi = [[] for i in range(len(dat[0]))]
    for i in range(len(dat)):
        for j in range(len(dat[i])):
            energi[j].append(dat[i][j])
    return energi[0], energi[2]


def plotfart(file):
    f = open(file)
    lines = f.readlines()
    lines = lines[3:]
    vx = []
    vy = []
    t = []
    x = []
    y = []

    for i in range(len(lines)):
        lines[i] = lines[i].strip()
        lines[i] = lines[i].split('\t')
        try:
            x.append(float(lines[i][1]))
            y.append(float(lines[i][2]))
            t.append(float(lines[i][0]))
            vx.append(float(lines[i][3]))
            vy.append(float(lines[i][4]))
        except:
            a = 1
    f.close()
    return t, vx, vy, x, y


def main():
    red = "red1_new"
    green = "ny_green1"
    blue = "blue1_allnew"
    file1 = 'one_nude.txt'  # lang
    file2 = 'two_nude.txt'  # fart data
    file3 = 'nude_data.txt'  # lang
    current = red

    n = 100000
    h = 0.00001
    t = [i*h for i in range(n)]


#    values = trvalues(p,xfart)

    font = {'fontname': 'Arial', 'size': '16'}
    lsize = 14
    tsize = 12

    x0, v0 = 0.05, 0
    x, v, y, f, N = euler(h, n, x0, v0, current)
    x1, y1 = bane(current)

    vt, vx, vy, xm, ym = plotfart(current)

    plt.plot(x1, y1, label='simulert')
    plt.plot(xm, ym, label='målt')
    plt.xlabel('x-posisjon [m]', **font)
    plt.ylabel('y-posisjon [m]', **font)
    plt.legend(loc=1, fontsize=lsize)
    plt.tick_params(labelsize=tsize)
    plt.savefig('Banegraf.png', dpi=500)
    plt.show()

    plt.plot(xm,vx, label ="simulert x-hastighet")
    plt.xlabel('x-posisjon [m]', **font)
    plt.ylabel('fart i x-retning [m/s]', **font)
    plt.legend(loc=1, fontsize=lsize)
    plt.tick_params(labelsize=tsize)
    plt.savefig('fartoverx.png', dpi=500)
    plt.show()

# t1, y1 = getfile(file3)
    #plt.plot(t, y, label='Simulert')
    #plt.plot(t1, y1, label='Målt')
    #plt.xlabel('Tid [s]', **font)
    #plt.ylabel('y-posisjon [m]', **font)
    #plt.legend(loc=1, fontsize=lsize)
    # plt.tick_params(labelsize=tsize)
    #plt.savefig('ygraf.png', dpi=500)
    # plt.show()

    plt.plot(x[:len(f)], f, label='Friksjon')
    plt.xlabel('x-posisjon [m]', **font)
    plt.ylabel('Friksjonskraft [N]', **font)
    plt.legend(loc=1, fontsize=lsize)
    plt.tick_params(labelsize=tsize)
    plt.savefig('friksjon.png', dpi=500)
    plt.show()

    plt.plot(x[:len(f)], N, label='Normalkraft')
    plt.axhline(linewidth=1, color='black')
    plt.xlabel('x-posisjon [m]', **font)
    plt.ylabel('Normalkraft [N]', **font)
    plt.ylim([-1, 2])
    plt.xlim([0, xmax])
    plt.legend(loc=8, fontsize=lsize)
    plt.tick_params(labelsize=tsize)
    plt.savefig('normalkraft.png', dpi=500)
    plt.show()

    # plotrate('Data.txt')

    plt.plot(t, v, label='Simulert')
    plt.plot(vt, vx, label="Målt x fart")
    # plt.plot(vt,vy, label = "y fart")
    plt.xlabel('Tid [s]', **font)
    plt.ylabel('Fart i x-retning [m/s]', **font)
    plt.legend(loc=1, fontsize=lsize)
    plt.tick_params(labelsize=tsize)
    plt.title("Ball høyde h3")
    plt.savefig('h3_fartsgraf.png', dpi=500)
    plt.show()




main()
