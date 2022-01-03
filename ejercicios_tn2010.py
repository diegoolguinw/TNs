# -*- coding: utf-8 -*-

# Librerías a ocupar

import numpy as np

# Función de lado derecho
def fd(kM, p, t):
    return kM*p/(np.exp(p*t)-1)

def yMa(kM, p, c, t):
    return kM*(np.log(np.exp(p*t)-1)-p*t) + c

def omori_RK4(kM,p,c,t0,T,h,f):
    N = int((T-t0)/h)
    t = np.linspace(t0,T,N)
    y_RK4 = np.zeros(N)
    y_RK4[0] = yMa(kM, p, c, t0)
    for i in range(N-1):
        g1 = f(kM, p, t[i])
        g2 = f(kM, p,t[i]+h/2)
        g3 = f(kM, p,t[i]+h/2)
        g4 = f(kM, p,t[i+1])
        y_RK4[i+1] = y_RK4[i] + h/6*(g1+2*g2+2*g3+g4)
    return y_RK4

def hint_1():
    print("Recuerde que debe imponer la misma condición inicial que las partes anteriores.")
    print("El método de RK4 se basa en calcular g1, g2, g3 y g4, luego calcular el siguiente elemento de la solución")
    print("con la recurrencia conocida.")
