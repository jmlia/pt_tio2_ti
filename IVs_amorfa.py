# -*- coding: utf-8 -*-
"""
Labo 7 
Martín & José
Script de análisis de las curvas IV, producidas al completar (o forzar la 
completición) el ciclo de estabilización en el LabVIEW. Se sigue la convención
de columnas y nombres de ese ciclo. 
Este en particular es para la muestra de Amorfa 9B, en las IVs que se 
produjeron en la bajada y subida de temperatura. 

Este script grafica las curvas IV para cada temperatura

Hecho usando pandas y DataFrames
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from funciones_gamma import central

#%%
# Importamos todos los datos para trabajar. Se hace una DataFrame para cada 
# temperatura con todas las mediciones de la curva IV, y luego una DataFrame
# para las tensiones, que contenga sus valores para cada temperatura (idem corrientes)

muestra = 'Amorfa 9B'

temps_bajada = [300, 260, 220, 180, 140, 100]
temps_subida = [100, 130, 160, 190, 210, 240]

path = '../datos/muestra/amorfa/' 
archivo_bajada = 'Pt_TiO2_Ti_9B_14_{:d}.00K_IV_08_11_20_42_54.txt'
archivo_subida = 'Pt_TiO2_Ti_9B_14_{:d}.00K_IV_09_11_13_26_50.txt'

# Estos son diccionarios de DataFrames, cuya key es la temperatura
datos_bajada = {}
# Estos son directamente DataFrames, y el nombre de las columnas corresponde
# a la temperatura
Vp_bajada = pd.DataFrame({})
Ip_bajada = pd.DataFrame({})

for temp in temps_bajada:
    key = '{:d}'.format(temp)  # para simplificar los slices posteriores
    datos_bajada[key] = pd.read_csv(path + archivo_bajada.format(temp))
    Vp_bajada[key] = datos_bajada[key]['V Pulso A (V)']
    Ip_bajada[key] = datos_bajada[key]['I Pulso (A)']*1000  # paso a mA
    
    
datos_subida = {}
Vp_subida = pd.DataFrame({})
Ip_subida = pd.DataFrame({})
for temp in temps_subida:
    key = '{:d}'.format(temp)
    datos_subida['{:d}'.format(temp)] = pd.read_csv(path + archivo_subida.format(temp))
    Vp_subida[key] = datos_subida[key]['V Pulso A (V)'] 
    Ip_subida[key] = datos_subida[key]['I Pulso (A)']*1000  # paso a mA
    
# Ahora una Dataframe con los datos de todas las temperaturas juntos

temperaturas_todas = sorted(temps_bajada + temps_subida, reverse = True)

# Los de 100K tienen un switcheo medio violento, los descartamos
del temperaturas_todas[-2:]

Vp_todas = pd.DataFrame({})
Ip_todas = pd.DataFrame({})

for temp in temperaturas_todas:
    
    key = '{:d}'.format(temp)
    
    try:
        Vp_todas[key] = Vp_bajada[key]
        Ip_todas[key] = Ip_bajada[key]
    
    except KeyError:
        Vp_todas[key] = Vp_subida[key]
        Ip_todas[key] = Ip_subida[key]

# Elimino las repeticiones sobre las curvas IV (cuando repite corrientes)
        
for temp in temperaturas_todas:
    
    key = '{:d}'.format(temp)
    
    imax = np.argmax(Ip_todas[key])
    
    Ip_todas[key] = Ip_todas[key][imax:]
    Vp_todas[key] = Vp_todas[key][imax:]

# Para la de 300K tengo que depurar un cacho más porque la I tiene una discontinuidad

index = 1496

Ip_todas['300'] = Ip_todas['300'][index:]
Vp_todas['300'] = Vp_todas['300'][index:]


#%% Graficamos las IVs
        
plt.figure()
for temp in temperaturas_todas:
    
    key = '{:d}'.format(temp)
    plt.plot(Vp_todas[key], Ip_todas[key], label = '{:d} K'.format(temp))

plt.xlabel('Voltaje (V)', fontsize = 15)
plt.ylabel('Corriente (mA)', fontsize = 15)
plt.title(muestra, fontsize = 15)
plt.legend()
plt.grid(True)


plt.figure()
for temp in temperaturas_todas:
    
    key = '{:d}'.format(temp)
    mask = Vp_todas[key] > 0
    
    plt.semilogy(Vp_todas[key][mask], np.abs(Ip_todas[key][mask]), label = '{:d} K'.format(temp))

plt.xlabel('Voltaje (V)', fontsize = 15)
plt.ylabel('Log Corriente (mA)', fontsize = 15)
plt.title(muestra, fontsize = 15)
plt.legend()
plt.grid(True)


#%% Ahora hago la extrapolación lineal al final de cada curva



def tangent(x, y, i):
    """
    Devuelve la recta tangente a la curva y(x) correspondiente a la posición
    dada por el indice i.
    La fórmula es la siguiente:
    
    r_x0 (x) = y'(x_0)(x - x_0) + y(x_0)
    
    Devuelve la función que a cada x le hace corresponder el valor de la recta 
    tangente al punto en cuestión.
    """
    
    m = central(x, y)[i]
    
    x_0 = x[i]
    y_0 = y[i]
    
    tangente = lambda x: m*(x - x_0) + y_0

    return tangente

def extrapole(x, y, reach):
    """
    Devuelve la extrapolación lineal de la función y(x).
    
    Toma un intervalo de valores desde x_max hacia atrás, dado por los índices
    desde el que corresponda a x_max (i_max) hasta i_max + reach. (Teneindo
    en cuenta que el array de tensiones (x) está ordenado de forma decreciente)
    
    Hace y devuelve el ajuste lineal que corresponde a los datos de y(x) en 
    ese intervalo.
    
    """
    
    i_max = np.argmax(x)
    
    x_data = x[i_max : imax + reach]
    y_data = y[i_max : imax + reach]
    
    linear = lambda x, m, b: m*x + b

    m, b = np.polyfit(x_data, y_data, deg = 1)
    extrap = lambda x: linear(x,m,b)

    return extrap

    












    
    
    