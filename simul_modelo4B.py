#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LBT
Martín & José

Script para generar la curva IV simulando el circuito del modelo 4B (dos diodos
Schottky con resistencias en serie, en paralelo con una resistencia y un 
elemento SCLC)

Contiene funciones para generar dicha curva IV seleccionando todos los 
parámetros relevantes del modelo de circuito. 
Estos parámetros se guardan en la tupla guess_param, que después puede ser 
usada en el script fit_modelo4B.

Si se corre este script por sí mismo, se genera la curva IV para los parámetros
elegidos, y se calcula también su curva gamma. Luego, se grafican ambas en 
comparación con los datos de la muestra amorfa a 300K.

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from funciones_gamma import gamma_sin

from gammas_amorfa_sinplot import Vp_bajada300, Ip_bajada300, Vs_bajada300_mas, gamma_filt_bajada300_mas

muestra = 'Circuito Modelo 4B'

#%%

def generar_I(V, A, n, R, I_0, b, R_b, start):
    """
    Funcion de ajuste para el modelo de circuito propeusto. 
    Los parametros son:
    
    A:    factor del elemento SCLC
    n:    exponente del elemento SCLC
    R:    resistencia en paralelo a todos
    I_sb: corriente de saturacion del diodo Schottky (bottom)
    b:    inversa de V_termico por factor de idealidad del diodo Schottky
    R_b:  resistencia en serie al diodo Schottky (bottom)
    
    Para cada valor de V, determina el I que viene dado implicitamente por la
    ecuacion del circuito. El parámetro 'start' es el punto de partida para la
    resolución de la ecuación implícita.
    
    """

    
    I_sclc = A*V**n
    
    I_r = V/R
    
    V_d = lambda I: np.abs(V - (I - I_sclc - I_r)*R_b)

    I_schB = lambda I: I_0 * ( np.exp(b * np.sqrt(V_d(I))) - 1)
    
    
    implicit = lambda I: I - I_sclc - I_r - I_schB(I)
    
    I = fsolve(implicit, start)

    return I


def generar_IV(V, A, R, I_0, b, R_b, n=3): # guarda que el n está fijo
    """
    Idem a 'generar_I', pero recibe un array de voltajes y devuelve uno de 
    corrientes. 
    Aquí dejamos fijo el parámetro n del SCLC.
    Para cada valor de V, el valor de 'start' para la solución de la ecuación
    que determina I es el valor de I inmediatamente anterior.
    """
    
    I = np.zeros(len(V))
    I[0] = generar_I(V[0], A, n, R, I_0, b, R_b, 0)

    for i in range(len(V)):
    
        I[i] = generar_I(V[i], A, n, R, I_0, b, R_b, I[i-1])
        
    return I
    

#%% Generamos voltajes y corrientes
    
A = 0.3e0
n = 3     # este queda fijo, solo lo pongo acá para que aparezca en los labels
R = 4.3e0
I_0 = 1.77e0
b = 0.01e-1
R_b = 5e0

#guess_param = (A, n, R, I_0, b, R_b)
guess_param = (A, R, I_0, b, R_b)

if __name__ == '__main__':
    
    V_bot = 0.0
    V_top = 1.2
    
    V = np.linspace(V_bot, V_top, 1000)
    I = generar_IV(V, *guess_param)
        
        
    #%% Cálculo del gamma
    
    gamma = gamma_sin(V, I)
    
    exponente = 0.5
    Vs = np.abs(V)**exponente
    
    
    #%% Ploteo de la curva IV y del gamma
    
    
    params = ' $A = {0}$\n $R = {1}$\n $I_0 = {2:.2E}$\n $b = {3}$\n $R_b = {4}$\n $n = {5}$'.format(*guess_param, n)
    
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(V, I, label = 'simulacion\n'+ params)
    plt.plot(Vp_bajada300, Ip_bajada300, 'c-.', label = 'datos a 300K')
    plt.xlabel('$V$', fontsize = 15)
    plt.ylabel('$I$', fontsize = 15)
    plt.title(muestra + ': curva IV', fontsize = 15)
    plt.legend()
    plt.grid(True)
    
    plt.subplot(1,2,2)
    plt.plot(Vs, gamma, color = 'r', label = 'simulacion\n'+params)
    plt.plot(Vs_bajada300_mas, gamma_filt_bajada300_mas, 'c-.', label = 'datos a 300K')
    plt.xlabel('$|V|^{%.1f}$'%exponente, fontsize = 15)
    plt.ylabel('$\\gamma$', fontsize = 15)
    plt.title(muestra + ': curva $\\gamma$', fontsize = 15)
    plt.legend()
    plt.grid(True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    


