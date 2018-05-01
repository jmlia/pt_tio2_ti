#!/usr/bin/env python
from lmfit import *
import numpy as np
#########################################
#        Aqui cargamos los archivos de la data medida
##########################################

#IV = np.loadtxt('IVs_Diego2.txt')
IV1 = np.loadtxt('IV_anatasa_16A_300K_22_11_positivos.txt')
IV = np.loadtxt('IV_anatasa_16A_300K_22_11_negativos.txt')

############################################################
#parte positiva de las medidas experimentales
############################################################

V1 = IV1[:,0]

I1 = IV1[:,1]

###############################################
#para la parte negativa donde tomamos valores absoluto
#################################################

V = np.abs(IV[:,0])

I = np.abs(IV[:,1])

###########################################

print ("------------------------------")
print ("------- Calculando ... -------")
print ("------------------------------")

###################################################
# Se define la funcion objeto y retorna el vector con los parametros optimizados
####################################################


#######################################################
#              La funcion para la parte positiva.
######################################################


def fcn2min(params, V, I):
	w = params.valuesdict()
	#model = ((V -  (I* w['R2']))/w['R1'])*(1.0 + (w['R1']/w['RPF'])*np.exp(w['C']*np.sqrt(V - (I* 
	modelplus = 1.0/w['RPF1']*V1*np.exp(w['a1']*np.sqrt(V1))+\
	((1.0/w['RT1'])*( V1 - w['RB1']*(I1 - (1.0/w['RPF1'])*V1*np.exp(w['a1']*np.sqrt(V1))) ))+\
	w['Is11']*(np.exp(w['b1']*(V1 - (w['RB1']* (I1 - (1.0/w['RPF1'])*V1*np.exp(w['a1']*np.sqrt(V1))))))-1.0)
	return modelplus - I1
	
params = Parameters()
########################################
#Estos son los parametros para la ecuacion model que esta comentado primero
#########################################

#params.add('R2',   value = 4.999999,  min = 0.09999, max = 10000)

#params.add('R1', value= 4.0, min= 0.1, max = 3000000000)

#params.add('RPF', value= 1000.89, min= 100.10995, max=10000000000)

#params.add('C', value= 8.99, min = 2. , max = 20)

######################################################################

params.add('RPF1', min= 10, max=8.0e8)

params.add('a1',   min = 1.0 , max = 1000)

params.add('Is11', min = 0.1e-5 , max = 5.0e-9)

params.add('b1',   min = 1.0 , max = 1000)

params.add('RB1',  min = 1.0, max = 100000)

params.add('RT1',  min = 1.0, max = 100000)

minner = Minimizer(fcn2min, params,fcn_args=(V, I), nan_policy='omit')


result = minner.minimize()

# calculate final result

final1 = I1 + result.residual

#derivada para calcular Gamma

dI_ex1 = np.zeros(I1.shape,np.float)

dI_ex1[0:-1] = np.diff(np.log(I1))/np.diff(np.log(V1))

dI_ex1[-1] = (np.log(I1[-1] - I1[-2]))/(np.log(V1[-1] - V1[-2]))

dI_teor1 = np.zeros(final1.shape,np.float)

dI_teor1[0:-1] = np.diff(np.log(final1))/np.diff(np.log(V1))

dI_teor1[-1] = (np.log(final1[-1] - final1[-2]))/(np.log(V1[-1] - V1[-2]))

# write error report

report_fit(result)

try:
	import matplotlib.pyplot as plt
	plt.plot(V1, I1, color='black', linewidth=2.5, linestyle='-', label='experimental_300+')
	plt.plot(V1, final1, color='red', linewidth=2.5, linestyle='--', label='Ajuste_Teorico_300+')
	plt.xlabel(r"$V (volt)$", fontsize = 15, color = (1,0,0))
	plt.ylabel(r"$I (mA)$", fontsize = 15, color = 'blue')
	#plt.xlabel('V (volt)')
	#plt.ylabel('I (mA)') 
	plt.legend()
	plt.show()
except ImportError:
	pass

try:
	import matplotlib.pyplot as plt
	plt.plot(V1**0.5, dI_ex1, color='black', linewidth=2.5, linestyle='-', label='Gamma_experimental_300+')
	plt.plot(V1**0.5, dI_teor1, color='red', linewidth=2.5, linestyle='--', label='Gamma_Ajuste_Teorico_300+')
	plt.xlabel(r"$V^{1/2} (volt)$", fontsize = 15, color = (1,0,0))
	plt.ylabel(r"$\gamma$", fontsize = 24, color = 'blue')
	#plt.xlabel('V\^0.5 (volt)')
	#plt.ylabel('Gamma') 
	plt.legend()
	plt.show()
except ImportError:
	pass
	
	
##############################################################################
#                             funcion parte negativa.
###############################################################################	


def fcn2min(params, V, I):
	w = params.valuesdict()
	#model = ((V -  (I* w['R2']))/w['R1'])*(1.0 + (w['R1']/w['RPF'])*np.exp(w['C']*np.sqrt(V - (I* 
	#modelplus = 1.0/w['RPF']*V*np.exp(w['a']*np.sqrt(V))+\
	#((1.0/w['RT'])*( V - w['RB']*(I - (1.0/w['RPF'])*V*np.exp(w['a']*np.sqrt(V))) ))+\
	 # w['Is1']*(np.exp(w['b']*(V - (w['RB']* (I - (1.0/w['RPF'])*V*np.exp(w['a']*np.sqrt(V))))))-1.0)
	  
	modelminus = 1.0/w['RPF']*V*np.exp(w['a']*np.sqrt(V))+\
	((1.0/w['RB'])*( V - w['RT']*(I - (1.0/w['RPF'])*V*np.exp(w['a']*np.sqrt(V))) ))+\
	 w['Is1']*(np.exp(w['b']*(V - (w['RT']* (I - (1.0/w['RPF'])*V*np.exp(w['a']*np.sqrt(V))))))-1.0)  
	#model = ((V)/w['R1'])*(1.0 + (w['R1']/w['RPF'])*np.exp(w['C']*np.sqrt(V)))
	return modelminus - I
	
	
####################################################################3	
# Creamos el conjunto de parametros a hacer calculado.
#####################################################################

params = Parameters()

########################################
#Estos son los parametros para la ecuacion model que esta comentado primero
#########################################

#params.add('R2',   value = 4.999999,  min = 0.09999, max = 10000)

#params.add('R1', value= 4.0, min= 0.1, max = 3000000000)

#params.add('RPF', value= 1000.89, min= 100.10995, max=10000000000)

#params.add('C', value= 8.99, min = 2. , max = 20)

######################################################################

params.add('RPF',  min= 1.0, max=2001000)

params.add('a',   min = 1.0 , max = 10)

params.add('Is1', min = 0.1e-5 , max = 5.0e-9)

params.add('b',   min = 1.0 , max = 100)

params.add('RB',  min = 1.0, max = 100000)

params.add('RT',  min = 1.0, max = 9219000)

minner = Minimizer(fcn2min, params,fcn_args=(V, I), nan_policy='omit')

result = minner.minimize()

# calculate final result

final = I + result.residual

#derivada para calcular Gamma

dI_ex = np.zeros(I.shape,np.float)

dI_ex[0:-1] = np.diff(np.log(I))/np.diff(np.log(V))

dI_ex[-1] = (np.log(I[-1] - I[-2]))/(np.log(V[-1] - V[-2]))

dI_teor = np.zeros(final.shape,np.float)

dI_teor[0:-1] = np.diff(np.log(final))/np.diff(np.log(V))

dI_teor[-1] = (np.log(final[-1] - final[-2]))/(np.log(V[-1] - V[-2]))

# write error report

report_fit(result)

# try to plot results

try:
	import matplotlib.pyplot as plt
	plt.plot(V, I, color='black', linewidth=2.5, linestyle='-', label='experimental_300-')
	plt.plot(V, final, color='red', linewidth=2.5, linestyle='--', label='Ajuste_Teorico_300-')
	plt.xlabel(r"$V (volt)$", fontsize = 15, color = (1,0,0))
	plt.ylabel(r"$I (mA)$", fontsize = 15, color = 'blue')
	#plt.xlabel('V (volt)')
	#plt.ylabel('I (mA)') 
	plt.legend()
	plt.show()
except ImportError:
	pass

try:
	import matplotlib.pyplot as plt
	plt.plot(V**0.5, dI_ex, color='black', linewidth=2.5, linestyle='-', label='Gamma_experimental_300-')
	plt.plot(V**0.5, dI_teor, color='red', linewidth=2.5, linestyle='--', label='Gamma_Ajuste_Teorico_300-')
	plt.xlabel(r"$V^{1/2} (volt)$", fontsize = 15, color = (1,0,0))
	plt.ylabel(r"$\gamma$", fontsize = 24, color = 'blue')
	#plt.xlabel('V\^0.5 (volt)')
	#plt.ylabel('Gamma') 
	plt.legend()
	plt.show()
except ImportError:
	pass
