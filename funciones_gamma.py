# -*- coding: utf-8 -*-
"""
Labo 7 
Martín & José

Módulo de funciones para calcular el factor gamma
en todos los scripts importar como:
    
import funciones_gamma as gm 

ó

from funciones_gamma import *
"""

import numpy as np
from scipy.signal import savgol_filter
from scipy.signal import medfilt
import pandas as pd

#%% Funciones para calcular el gamma
    
def central(x, y, win = 1):
    '''
    Calcula mediante diferencia finita central, la derivada
    dy/dx. Para los valores de los bordes, repite los más cercanos.
    '''
    dydx = np.zeros(len(x))
    for i in range(1,len(x)-1):
        dydx[i] = (y[i+1] - y[i-1])/(x[i+1] - x[i-1])
    dydx[0] = dydx[1] # para mantener continuidad aprox
    dydx[-1] = dydx[-2] # idem
    
    return dydx 

def smooth(data, win):
    
    # Acá hay que agregar el parseo para evitar problemas...
    hw = win//2  # half-window
    N = data.shape[0]
    sdata = np.zeros(N)
    
    for i in range(hw, N - hw):        
        sdata[i] = np.mean(data[i - hw:i + hw])    

    return sdata

def gamma_smooth(V, I, win):
    '''
    Calcula la gamma smootheando con José los datos y la derivada. Usa el 
    filtro en una ventana de largo 'win'.
    '''
    
    V = np.asarray(V)
    I = np.asarray(I)
    
    logV = np.log(np.abs(V))
    logI = np.log(np.abs(I))
    
    logV_sm = savgol_filter(logV, win, 3)
    logI_sm = savgol_filter(logI, win, 3)
    
    gamma = central(logV_sm, logI_sm)
    gamma_sm = savgol_filter(gamma, win, 3)

    return gamma_sm


def gamma_savgol1(V, I, win):
    '''
    Calcula la gamma smootheando con savgol sólo la derivada. Usa el 
    filtro en una ventana de largo 'win', que debe ser impar
    '''
    assert (win % 2) == 1, 'La ventana no es impar'
    
    V = np.asarray(V)
    I = np.asarray(I)
    
    Va = np.abs(V)
    Ia = np.abs(I)
    
    gamma = Va / Ia / central(Ia, Va)
    gamma_sm = savgol_filter(gamma, win, 3)
    
    return gamma_sm   
    
def gamma_savgol2(V, I, win):
    '''
    Calcula la gamma smootheando con savgol los datos y la derivada. Usa el 
    filtro en una ventana de largo 'win', que debe ser impar
    '''
    assert (win % 2) == 1, 'La ventana no es impar'
    
    V = np.asarray(V)
    I = np.asarray(I)
    V_sm = np.abs(savgol_filter(V, win, 3))
    I_sm = np.abs(savgol_filter(I, win, 3))
    
    logV = np.log(V_sm)
    logI = np.log(I_sm)
    
    logV_sm = savgol_filter(logV, win, 3)
    logI_sm = savgol_filter(logI, win, 3)
    
    gamma = central(logV_sm, logI_sm)
    gamma_sm = savgol_filter(gamma, win, 3)
    
    return gamma_sm   
    
def gamma_savgol3(V, I, win):
    '''
    Calcula la gamma smootheando con savgol los datos y la derivada. Usa el 
    filtro en una ventana de largo 'win', que debe ser impar. No calcula 
    logaritmos sino directamente V/I*dI/dV
    '''
    assert (win % 2) == 1, 'La ventana no es impar'
    
    V = np.asarray(V)
    I = np.asarray(I)
    V_sm = np.abs(savgol_filter(V, win, 3))
    I_sm = np.abs(savgol_filter(I, win, 3))
    
    #gamma = central(V_sm, I_sm)*V_sm/I_sm
    # tomo dI/dV = (dV/dI)^-1
    gamma = V_sm / I_sm / central(I_sm, V_sm)
    gamma_sm = savgol_filter(gamma, win, 3)
    
    return gamma_sm


def gamma_medfilt(V, I, win):
    '''
    Calcula la gamma smootheando con median  los datos y la derivada. No calcula 
    logaritmos sino directamente V/I*dI/dV
    '''
    assert (win % 2) == 1, 'La ventana no es impar'
    
    V = np.asarray(V)
    I = np.asarray(I)
    V_sm = np.abs(medfilt(V, win))
    I_sm = np.abs(medfilt(I, win))
    
    #gamma = central(V_sm, I_sm)*V_sm/I_sm
    # tomo dI/dV = (dV/dI)^-1
    gamma = V_sm / I_sm / central(I_sm, V_sm)
    gamma_sm = medfilt(gamma, win)
    
    return gamma_sm


def gamma_sin(V, I):
    '''
    Calcula la gamma sin smoothear los datos y la derivada. No calcula 
    logaritmos sino directamente V/I*dI/dV
    '''
    
    V = np.asarray(V)
    I = np.asarray(I)
   
    #gamma = central(V_sm, I_sm)*V_sm/I_sm
    # tomo dI/dV = (dV/dI)^-1
    gamma = V / I / central(I, V)
    
    return gamma