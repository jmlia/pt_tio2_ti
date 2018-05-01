#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LBT
Martín & José

Script para el ajuste de IV simulando el circuito del modelo 4A (dos diodos
Schottky en paralelo con una rama con una resistencia y un elemento SCLC)

Este programa hace un ajuste de los datos de la muestra amorfa a 300K, 
utilizando como función de ajuste 'generar_IV', definida en el script 
simul_modelo4A. Como parámetros iniciales, se utilizan los guess_param, 
definidos en ese mismo script. 


"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from funciones_gamma import gamma_sin
from simul_modelo4A import generar_IV, guess_param

from gammas_amorfa_sinplot import Vp_bajada300, Ip_bajada300

from gammas_amorfa_sinplot import Vs_bajada300_mas, gamma_filt_bajada300_mas

muestra = 'Circuito Modelo 4A'

    
#%%

V_bot = 0.0
V_top = 1.2

mask = (Vp_bajada300['300'] > V_bot) * (Vp_bajada300['300'] < V_top) 

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

params = ' $A = {0}$\n $n = {1}$\n $R = {2:.2E}$\n $I_0 = {3:.2E}$\n $b = {4}$\n $'.format(*popt)


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

