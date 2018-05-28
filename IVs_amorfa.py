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

logI = np.log(np.abs(Ip_todas))

#%% Graficamos las IVs
        
plt.figure()
for temp in temperaturas_todas:
    key = '{:d}'.format(temp)
    plt.semilogy(Vp_todas[key], Ip_todas[key], label = '{:d} K'.format(temp))

plt.xlabel('Voltaje (V)', fontsize = 15)
plt.ylabel('Corriente (mA)', fontsize = 15)
plt.title(muestra + ' bajada de temperatura', fontsize = 15)
plt.legend()
plt.grid(True)



















    
    
    