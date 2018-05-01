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


import lmfit as lmf

#%%

V_bot = 0.0
V_top = 1.2

mask = Vp_bajada300['300'] > 0 

x_data = np.asarray(Vp_bajada300[mask]).ravel()
y_data = np.asarray(Ip_bajada300[mask]).ravel()


# En lugar de usar curve_fit, intentemos usar lmfit
#popt, pcov = curve_fit(generar_IV, x_data, y_data, p0 = guess_param)

# ver: https://lmfit.github.io/lmfit-py/intro.html

params = lmf.Parameters()
params.add("A", value = 0.41)
params.add("n", min = 0, value = 3, vary = False)
params.add("R", min = 0.0, value = 5.0)
params.add("I0", value = 0.1)
params.add("b", value = 0.0001)
params.add("Rb", min = 0.0, value = 100)


def Residual(params, x, data):
    A = params["A"]
    n = params["n"]
    R = params["R"]
    I_0 = params["I0"]
    b = params["b"]
    R_b = params["Rb"]

    model = generar_IV(x, A.value, R.value, I_0.value, b.value, R_b.value, n.value)
    return (data - model)

minner = lmf.Minimizer(Residual, params, fcn_args = (x_data, y_data), nan_policy = 'omit')
result = minner.minimize()

V = np.linspace(V_bot, V_top, 1000)
I = generar_IV(V, result.params["A"], result.params["R"], result.params["I0"],
               result.params["b"], result.params["Rb"], result.params["n"])
    
#%% Cálculo del gamma

gamma = gamma_sin(V, I)

exponente = 0.5
Vs = np.abs(V)**exponente

#%%

params = ''
# params = ' $A = {0}$\n $R = {1}$\n $I_0 = {2:.2E}$\n $b = {3}$\n $R_b = {4}$\n $n = {5}$'.format(result.params["A"].value, result.params["R"].value, result.params["I0"].value, result.params["b"].value, result.params["Rb"].value, result.params["n"].value)

fig, (ax_1, ax_2) = plt.subplots(1,2)

ax_1.plot(V, I, label = 'ajuste\n'+ params)
ax_1.plot(Vp_bajada300, Ip_bajada300, 'c-.', label = 'datos a 300K')
ax_1.set_xlabel('$V$', fontsize = 15)
ax_1.set_ylabel('$I$', fontsize = 15)
#ax_1.title(muestra + ': curva IV', fontsize = 15)
ax_1.legend()

ax_2.plot(Vs, gamma, color = 'r', label = 'ajuste\n'+params)
ax_2.plot(Vs_bajada300_mas, gamma_filt_bajada300_mas, 'c-.', label = 'datos a 300K')
ax_2.set_xlabel('$|V|^{%.1f}$'%exponente, fontsize = 15)
ax_2.set_ylabel('$\\gamma$', fontsize = 15)
#plt.title(muestra + ': curva $\\gamma$', fontsize = 15)
ax_2.legend()

fig.savefig('file.pdf', format = 'pdf')

