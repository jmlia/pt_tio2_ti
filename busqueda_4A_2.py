#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from funciones_gamma import gamma_sin
from simul_modelo4A import generar_IV
from gammas_amorfa_sinplot import Vp_bajada300, Ip_bajada300
from gammas_amorfa_sinplot import Vs_bajada300_mas, gamma_filt_bajada300_mas

output = 'Iteración {0:d}\n\tA = {1:.3E}' \
         '\n\tI0 = {2:.3E}\n\tb = {3:.3E}\n\tR = {4:.3E}\n\tResiduales: {5:.7f}\n\n'

V_bot = 0.0
V_top = 1.2

mask = (Vp_bajada300['300'] > V_bot) & (Vp_bajada300['300'] < V_top)
x_data = np.asarray(Vp_bajada300[mask]).ravel()
y_data = np.asarray(Ip_bajada300[mask]).ravel()

# Necesitamos el par (V_min, I(V_min)) para fijar R en cada paso...
#

n = np.argmin(x_data)
V_min = x_data[n]
I_min = y_data[n]

divisions = 5
arr_A = np.linspace(start = 0.1, stop = 10, num = divisions, endpoint = True)
arr_I0 = np.logspace(start = 1.0, stop = -5, num = divisions)
arr_b = np.linspace(start = 0.001, stop = 1, num = divisions)

# Mínimos.
minres = np.finfo('d').max
min_A = -1
min_I0 = -1
min_b = -1
i = 0
for A in arr_A:
    for I0 in arr_I0:
        for b in arr_b:

            # Dados estos valores, para cada iteración, hay que calcular R:
            # De acuerdo a la hoja:

            # R^{-1} = I_exp/V - A * V**(n-1) - I_s * b/sqrt(V).
            # Con (V,I_exp) el par que corresponde a V_min.

            R = 1.0/(I_min/V_min - A * V_min - (I0 * b)/np.sqrt(V_min))
            res = np.sum(np.abs(generar_IV(x_data, A, 2, R, I0, b) - y_data))

            if res < minres:
                minres = res
                min_A = A
                min_I0 = I0
                min_b = b
                
                print(output.format(i, A, I0, b, R, res))
                i += 1

# Calcula el R_min.
R = 1.0/(I_min/V_min - min_A * V_min - (min_I0 * min_b)/np.sqrt(V_min))

print('Mínimo alcanzado de suma |y_i - model_i| = {0:.7f}\n\tA = {1:.5f}\n\tR = {2:.5f}\n\tI0 = {3:.5f}\n\tb = {4:.5f}\n'.format(minres, min_A, R, min_I0, min_b))

# Cáclculo de gamma
fit = generar_IV(x_data, min_A, min_n, min_R, min_I0, min_b)
gamma_fit = gamma_sin(x_data, fit)

exponente = 0.5
Vs = np.abs(x_data)**exponente

#Esto debería graficar...

plt.figure()
plt.subplot(1,2,1)
plt.plot(x_data, y_data, 'o', label = 'Data')
plt.plot(x_data, fit, '-', label = 'Ajuste')
plt.grid(True)

plt.subplot(1,2,2)
plt.plot(Vs_bajada300_mas, gamma_filt_bajada300_mas, 'o', label = 'Data')
plt.plot(Vs, gamma_fit, '-', label = 'Ajuste')
plt.grid(True)
plt.show()
