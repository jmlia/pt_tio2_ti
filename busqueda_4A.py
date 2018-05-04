#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from funciones_gamma import gamma_sin
from simul_modelo4A import generar_IV, guess_param, n
from gammas_amorfa_sinplot import Vp_bajada300, Ip_bajada300
from gammas_amorfa_sinplot import Vs_bajada300_mas, gamma_filt_bajada300_mas

V_bot = 0.0
V_top = 1.2

mask = (Vp_bajada300['300'] > V_bot) * (Vp_bajada300['300'] < V_top)

x_data = np.asarray(Vp_bajada300[mask]).ravel()
y_data = np.asarray(Ip_bajada300[mask]).ravel()

# Acá se toman los parámetros que en el script
# `simul_modelo4A.py' están determinados como A, n, R, I_0, b.

# Algunas propiedades:
# A > 0
# 1 <= n <= 5, supongamos; y además entero.
# R > 0; y supongamos mayor a 1k, pero todo está en mA
#        por lo que no es necesario el factor 1000.
# I0;    debería andar en uA, que son 10^-3 mA.
# b > 0 (?)

arr_A = np.linspace(start = 0.1, stop = 10, num = 1000, endpoint = True)
arr_n = np.arange(start = 1, stop = 6)
arr_R = np.linspace(start = 1, stop = 1000, num = 1000)
arr_I0 = np.linspace(start = 1.0, stop = 10, num = 1000)
arr_b = np.linspace(start = 0.001, stop = 1, num = 1000)

i = 0
minres = np.finfo('d').max

# Mínimos.
min_A = -1
min_n = -1
min_R = -1
min_I0 = -1
min_b = -1

output = 'Iteración {0:d}\n\tA = {1:.3E}\n\tn = {2:d}\n\tR = {3:.3E}' \
         '\n\tI0 = {4:.3E}\n\tb = {5:.3E}\nResiduales: {6:.7f}\n\n'

for A in arr_A:
    for n in arr_n:
        for R in arr_R:
            for I0 in arr_I0:
                for b in arr_b:

                    res = np.sum(np.abs(generar_IV(x_data, A, n, R, I0, b) - y_data))

                    if res < minres:
                        minres = res
                        min_A = A
                        min_n = n
                        min_R = R
                        min_I0 = I0
                        min_b = b
          
                    print(output.format(i, A, n, R, I0, b, res))
                    i += 1

print('Mínimo alcanzado de suma |y_i - model_i| = {0:.7f}\n\tA = {1:.5f}\n\tn = {2:d}\n\tR = {3:.5f}\n\tI0 = {4:.5f}\n\tb = {5:.5f}\n'.format(minres, min_A, min_n, min_R, min_I0, min_b))

# Esto debería graficar...

# plt.figure()
# plt.subplots(1,1,1)
# plt.plot(x_data, y_data, 'o', label = 'Data')
# plt.plot(x_data, generar_IV(x_data, min_A, min_n, min_R, min_I0, min_b), '-', label = 'Ajuste')
# plt.grid(True)
# plt.show()
