# -*- coding: utf-8 -*-
"""
Labo 7 
Martín & José
Script de análisis de las curvas IV, producidas al completar (o forzar la 
completición) el ciclo de estabilización en el LabVIEW. Se sigue la convención
de columnas y nombres de ese ciclo. 
Este en particular es para la muestra de Amorfa 9B, en las IVs que se 
produjeron en la bajada y subida de temperatura. 

Este script unicamente calcula las curvas gamma para cada temperatura

Este script aplica un filtro distinto según que region de la curva gamma sea.
Se eligen dos ventanas, y un rango de valores de Vs. Para valores bajos de Vs
se usa una ventana (chica) con un primer filtro, mientras que a valores altos
se usa otra ventana (grande) con un segundo filtro. Estos pueden ser savgol
o medfilt

Hecho usando pandas y DataFrames
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from funciones_gamma import *

#%%
# Importamos todos los datos para trabajar. Se hace una DataFrame para cada 
# temperatura con todas las mediciones de la curva IV, y luego una DataFrame
# para las tensiones, que contenga sus valores para cada temperatura (idem corrientes)

muestra = 'Amorfa 9B'
# le saque el 100K a todos porque en general molesta
temps_bajada = [260, 220, 180, 140] # hago 300 por separado porque tiene mas valores
#temps_subida = [100, 130, 160, 190, 210, 240]
temps_subida = [240, 210, 190, 160]
temps_subida2 = [130] # idem estos que tienen menos valores

#path = 'C:\\Users\\usuario\\Documents\\Usuarios\\Martin_Jose\\Labo 7\\datos\\muestra\\amorfa\\' 
#path = 'D:\\martin\\Documents\\FCEN\\Labo7\\datos\\muestra\\amorfa\\' 
path = '../datos/muestra/amorfa/'
archivo_bajada = 'Pt_TiO2_Ti_9B_14_{:d}.00K_IV_08_11_20_42_54.txt'
archivo_subida = 'Pt_TiO2_Ti_9B_14_{:d}.00K_IV_09_11_13_26_50.txt'

# Estos son diccionarios de DataFrames, cuya key es la temperatura
datos_bajada = {}
# Estos son directamente DataFrames, y el nombre de las columnas corresponde
# a la temperatura
Vp_bajada = pd.DataFrame({})
Ip_bajada = pd.DataFrame({})

# al de 300 lo agrego por separado porque sino falla for some reason
Vp_bajada300 = pd.DataFrame({})
Ip_bajada300 = pd.DataFrame({})
key = '{:d}'.format(300)  
datos_bajada[key] = pd.read_csv(path + archivo_bajada.format(300))
Vp_bajada300[key] = datos_bajada[key]['V Pulso A (V)']
Ip_bajada300[key] = datos_bajada[key]['I Pulso (A)']*1000  # paso a mA
    
for temp in temps_bajada:
    key = '{:d}'.format(temp)  # para simplificar los slices posteriores
    datos_bajada[key] = pd.read_csv(path + archivo_bajada.format(temp))
    Vp_bajada[key] = datos_bajada[key]['V Pulso A (V)']
    Ip_bajada[key] = datos_bajada[key]['I Pulso (A)']*1000  # paso a mA
    
    
datos_subida = {}
Vp_subida = pd.DataFrame({})
Ip_subida = pd.DataFrame({})

Vp_subida2 = pd.DataFrame({})
Ip_subida2 = pd.DataFrame({})
for temp in temps_subida:
    key = '{:d}'.format(temp)
    datos_subida['{:d}'.format(temp)] = pd.read_csv(path + archivo_subida.format(temp))
    Vp_subida[key] = datos_subida[key]['V Pulso A (V)'] 
    Ip_subida[key] = datos_subida[key]['I Pulso (A)']*1000  # paso a mA
    
for temp in temps_subida2:
    key = '{:d}'.format(temp)
    datos_subida['{:d}'.format(temp)] = pd.read_csv(path + archivo_subida.format(temp))
    Vp_subida2[key] = datos_subida[key]['V Pulso A (V)'] 
    Ip_subida2[key] = datos_subida[key]['I Pulso (A)']*1000  # paso a mA


#%% Cálculo de las gamma a cada temperatura sin filtrar
    
gamma_sin_bajada300 = pd.DataFrame({})
gamma_sin_bajada = pd.DataFrame({})

gamma_sin_bajada300['300'] = gamma_sin(Vp_bajada300['300'], Ip_bajada300['300'])

for temp in temps_bajada:
    key = '{:d}'.format(temp)   
    gamma_sin_bajada[key] = gamma_sin(Vp_bajada[key], Ip_bajada[key])

    
gamma_sin_subida = pd.DataFrame({})

for temp in temps_subida:
    key = '{:d}'.format(temp)   
    gamma_sin_subida[key] = gamma_sin(Vp_subida[key], Ip_subida[key])
    
gamma_sin_subida2 = pd.DataFrame({})

for temp in temps_subida2:
    key = '{:d}'.format(temp)   
    gamma_sin_subida2[key] = gamma_sin(Vp_subida2[key], Ip_subida2[key])


#%% Cálculo de las gamma a cada temperatura con medfilt 'a dos filtros'

exponente = 0.5    
Vs_bajada300 = np.abs(Vp_bajada300)**exponente
Vs_bajada = np.abs(Vp_bajada)**exponente

Vs_subida = np.abs(Vp_subida)**exponente
Vs_subida2 = np.abs(Vp_subida2)**exponente
    
window1 = 3     # donde se aplica el primer filtro,valores bajos de V con poca ventana
window2 = 45    # donde se aplica el segundo filtro, el resto de los valores de V
thresh_volt = 0.6  # el valor de tension que separa donde va cada filtro (en la escala de V**exponente)

gamma_filt_bajada300 = pd.DataFrame({})
gamma_filt_bajada = pd.DataFrame({})


mask_b = Vs_bajada300['300'] > thresh_volt # define donde se aplica el segundo filtro
gamma_filt_bajada300['300'] = gamma_medfilt(Vp_bajada300['300'], Ip_bajada300['300'], window1)
gamma_filt_bajada300['300'][mask_b] = gamma_medfilt(Vp_bajada300['300'], Ip_bajada300['300'], window2)


for temp in temps_bajada:
    key = '{:d}'.format(temp)
    
    mask_b = Vs_bajada[key] > thresh_volt # define donde se aplica el segundo filtro
    
    gamma_filt_bajada[key] = gamma_medfilt(Vp_bajada[key], Ip_bajada[key], window1)
    # primero hace el primer filtro en todos, pero después el segundo
    
    gamma_filt_bajada[key][mask_b] = gamma_medfilt(Vp_bajada[key], Ip_bajada[key], window2)

    
gamma_filt_subida = pd.DataFrame({})

for temp in temps_subida:
    key = '{:d}'.format(temp)   
    mask_s = Vs_subida[key] > thresh_volt
    gamma_filt_subida[key] = gamma_medfilt(Vp_subida[key], Ip_subida[key], window1)
    gamma_filt_subida[key][mask_s] = gamma_medfilt(Vp_subida[key], Ip_subida[key], window2)
    
gamma_filt_subida2 = pd.DataFrame({})

for temp in temps_subida2:
    key = '{:d}'.format(temp)   
    mask_s = Vs_subida2[key] > thresh_volt
    gamma_filt_subida2[key] = gamma_medfilt(Vp_subida2[key], Ip_subida2[key], window1)
    gamma_filt_subida2[key][mask_s] = gamma_medfilt(Vp_subida2[key], Ip_subida2[key], window2)

    
#%% separo los gammas segun la polaridad

mask_bajada_mas300 = Vp_bajada300 > 0
mask_bajada_menos300 = Vp_bajada300 < 0

mask_bajada_mas = Vp_bajada > 0
mask_bajada_menos = Vp_bajada < 0

mask_subida_mas = Vp_subida > 0
mask_subida_menos = Vp_subida < 0

mask_subida2_mas = Vp_subida2 > 0
mask_subida2_menos = Vp_subida2 < 0


# =============================================================================
# =============================================================================

Vs_bajada300_mas = Vs_bajada300[mask_bajada_mas300]
Vs_bajada_mas = Vs_bajada[mask_bajada_mas]

Vs_bajada300_menos = Vs_bajada300[mask_bajada_menos300]
Vs_bajada_menos = Vs_bajada[mask_bajada_menos]


Vs_subida_mas = Vs_subida[mask_subida_mas]
Vs_subida2_mas = Vs_subida2[mask_subida2_mas]    
 
Vs_subida_menos = Vs_subida[mask_subida_menos]
Vs_subida2_menos = Vs_subida2[mask_subida2_menos]       
    

gamma_sin_bajada300_mas = gamma_sin_bajada300[mask_bajada_mas300]
gamma_sin_bajada_mas = gamma_sin_bajada[mask_bajada_mas]

gamma_sin_bajada300_menos = gamma_sin_bajada300[mask_bajada_menos300]
gamma_sin_bajada_menos = gamma_sin_bajada[mask_bajada_menos]


gamma_sin_subida_mas = gamma_sin_subida[mask_subida_mas]
gamma_sin_subida2_mas = gamma_sin_subida2[mask_subida2_mas]

gamma_sin_subida_menos = gamma_sin_subida[mask_subida_menos]
gamma_sin_subida2_menos = gamma_sin_subida2[mask_subida2_menos]


gamma_filt_bajada300_mas = gamma_filt_bajada300[mask_bajada_mas300]
gamma_filt_bajada_mas = gamma_filt_bajada[mask_bajada_mas]

gamma_filt_bajada300_menos = gamma_filt_bajada300[mask_bajada_menos300]
gamma_filt_bajada_menos = gamma_filt_bajada[mask_bajada_menos]


gamma_filt_subida_mas = gamma_filt_subida[mask_subida_mas]
gamma_filt_subida2_mas = gamma_filt_subida2[mask_subida2_mas]

gamma_filt_subida_menos = gamma_filt_subida[mask_subida_menos]
gamma_filt_subida2_menos = gamma_filt_subida2[mask_subida2_menos]


