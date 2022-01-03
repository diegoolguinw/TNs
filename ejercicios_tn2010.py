# -*- coding: utf-8 -*-
"""ejercicios_tn2010.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1XSjBjO5R7gQ05H4BTnxjuMIvp2YQellv
"""

# Librerías a ocupar

import numpy as np

# Función de lado derecho
def fd(kM, p, t):
    return kM*p/(np.exp(p*t)-1)

def yMa(kM, p, c, t):
    return kM*(np.log(np.exp(p*t)-1)-p*t) + c

def omori_RK4(kM,p,c,t0,T,h,N,f):
    t = np.linspace(t0,T,N)
    y_RK4 = np.zeros(N)
    y_RK4[0] = yMa(t0)
    for i in range(N-1):
        g1 = f(t[i])
        g2 = f(t[i]+h/2)
        g3 = f(t[i]+h/2)
        g4 = f(t[i+1])
        y_RK4[i+1] = y_RK4[i] + h/6*(g1+2*g2+2*g3+g4)
    return y_RK4
