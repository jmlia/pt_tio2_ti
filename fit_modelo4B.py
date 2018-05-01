#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LBT
Martín & José

Script para el ajuste de IV simulando el circuito del modelo 4B (dos diodos
Schottky con resistencias en serie, en paralelo con una resistencia y un 
elemento SCLC)

Este programa hace un ajuste de los datos de la muestra amorfa a 300K, 
utilizando como función de ajuste 'generar_IV', definida en el script 
simul_modelo4B. Como parámetros iniciales, se utilizan los guess_param, 
definidos en ese mismo script. 


"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from funciones_gamma import gamma_sin
from simul_modelo4B import generar_IV, guess_param, n

from gammas_amorfa_sinplot import Vp_bajada300, Ip_bajada300

from gammas_amorfa_sinplot import Vs_bajada300_mas, gamma_filt_bajada300_mas

muestra = 'Circuito Modelo 4B'

    
#%%

V_bot = 0.0
V_top = 1.2

mask = Vp_bajada300['300'] > 0 

x_data = np.asarray(Vp_bajada300[mask]).ravel()
y_data = np.asarray(Ip_bajada300[mask]).ravel()


popt, pcov = curve_fit(generar_IV, x_data, y_data, p0 = guess_param)

V = np.linspace(V_bot, V_top, 1000)
I = generar_IV(V, *popt)
    
#%% Cálculo del gamma

gamma = gamma_sin(V, I)

exponente = 0.5
Vs = np.abs(V)**exponente

#%%

params = ' $A = {0}$\n $R = {1}$\n $I_0 = {2:.2E}$\n $b = {3}$\n $R_b = {4}$\n $n = {5}$'.format(*popt, n)

plt.figure()
plt.subplot(1,2,1)
plt.plot(V, I, label = 'ajuste\n'+ params)
plt.plot(Vp_bajada300, Ip_bajada300, 'c-.', label = 'datos a 300K')
plt.xlabel('$V$', fontsize = 15)
plt.ylabel('$I$', fontsize = 15)
plt.title(muestra + ': curva IV', fontsize = 15)
plt.legend()
plt.grid(True)

plt.subplot(1,2,2)
plt.plot(Vs, gamma, color = 'r', label = 'ajuste\n'+params)
plt.plot(Vs_bajada300_mas, gamma_filt_bajada300_mas, 'c-.', label = 'datos a 300K')
plt.xlabel('$|V|^{%.1f}$'%exponente, fontsize = 15)
plt.ylabel('$\\gamma$', fontsize = 15)
plt.title(muestra + ': curva $\\gamma$', fontsize = 15)
plt.legend()
plt.grid(True)

