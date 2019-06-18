#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:02:16 2019

@author: aripekka
"""

from __future__ import print_function, division

import sys
import xraylib as xrl

from scipy.optimize import root_scalar
from math import exp, pi, sqrt


_SHELL_MACROS = {'K': xrl.K_SHELL, 'L1' : xrl.L1_SHELL, 'L2' : xrl.L2_SHELL, 
                'L3' : xrl.L3_SHELL, 'M1' : xrl.M1_SHELL, 'M2' : xrl.M2_SHELL, 
                'M3' : xrl.M3_SHELL, 'M4' : xrl.M4_SHELL, 'M5' : xrl.M5_SHELL}

_VERSION = '0.1'

#python 2 compatibility 
if sys.version_info[0] < 3:
    input = raw_input

print('----------------------------------')
print('**  XAS SAMPLE OPTIMIZER v.'+ _VERSION +'  **')
print('----------------------------------')

#Ask for user input
input_str = {}
input_str['compound'] = input('Compound: ')
input_str['absorption_edge'] = input('Absorption edge (e.g. Cu-K, U-L3): ')
photons_per_second = float(input('Incident photon flux (background subtracted) (ph/s): '))

compound = xrl.CompoundParser(input_str['compound'])

#convert the absorption edge string to energy
edge_element, edge_shell = input_str['absorption_edge'].split('-')

edge_energy = xrl.EdgeEnergy(xrl.SymbolToAtomicNumber(edge_element), 
                             _SHELL_MACROS[edge_shell])


#calculate the mass attenuation coefficient below and above edge
murho = [0,0]

for i in range(compound['nElements']):
    murho[0] = murho[0] + compound['massFractions'][i]*xrl.CS_Total(compound['Elements'][i], edge_energy-0.001)
    murho[1] = murho[1] + compound['massFractions'][i]*xrl.CS_Total(compound['Elements'][i], edge_energy+0.001)
    
#calculate filler attenuation using an "average" nitrogen filter
filler_mux = xrl.CS_Total(7, edge_energy)*1.5*0.1

#Measurement specific values for background count rates, estimates are used
beta = 1.1
gamma = 0.03*exp(filler_mux)


#optimization condition f = 0
def f(x,mu0,mu,beta=1, gamma=0,filler_mux=0):
    aux1 = sqrt(8*beta*exp(-filler_mux))*sqrt(exp(mu0*x) + exp(mu*x) + gamma*(exp(2*mu0*x) + exp(2*mu*x)) )
    aux2 = 2*gamma*((1-mu0*x)*exp(2*mu0*x) + (1-mu*x)*exp(2*mu*x))
    aux3 = (2-mu0*x)*exp(mu0*x) + (2-mu*x)*exp(mu*x) 

    return aux1 + aux2 + aux3

#calculate the sample mass to cross sectional area of the sample container (in g/cm^2)
sample_mass2area = root_scalar(f,args=(murho[0],murho[1],beta,gamma,filler_mux),bracket=[0,5/max([murho[0],murho[1]])]).root

optimal_mux0 = murho[0]*sample_mass2area
optimal_mux = murho[1]*sample_mass2area

#note that in presence of filler this is an effective time ratio T'/T0
time_ratio = sqrt(exp(optimal_mux0) + exp(optimal_mux) + gamma*(exp(2*optimal_mux0) + exp(2*optimal_mux)))/sqrt(2*beta*exp(-filler_mux))

#estimate measurement time, rel_delta_mux = sought after statistical accuracy
rel_delta_mux = 0.01
T0 = 2*beta/photons_per_second*(1 + time_ratio)/rel_delta_mux**2/(optimal_mux - optimal_mux0)**2

print('')
print(u'SAMPLE INFO:')

print(u'    Energy of ' + edge_element + ' ' + edge_shell + ' edge: ' + str(edge_energy) + ' keV')
print(u'    Sample mass attenuation (cm²/g; below and above edge): %.2f, %.2f' % (murho[0],murho[1]))
print(u'    Assuming filler attenuation coefficient µx = %.2f (for 1 mm of Z=7, ρ = 1.5 g/cm³)' % filler_mux)

print('')
print(u'OPTIMAL SAMPLE:')

#print(u'    Optimal sample µx w/o filler (below edge, above edge):  %.3f, %.3f ' % (optimal_mux0, optimal_mux))
print(u'    Optimal sample µx w/ filler (below edge, above edge):    %.3f, %.3f ' % (optimal_mux0+filler_mux, optimal_mux+filler_mux))
print(u'    Optimal sample edge step Δµx:                            %.3f' % (optimal_mux-optimal_mux0))
print(u'    Transmission below and above edge w/ filler (%):         ' + '%.2f, %.2f' % (exp(-filler_mux-optimal_mux0)*100, exp(-filler_mux-optimal_mux)*100))
print(u'    Transmitted fluxes below and above edge (ph/s):          %i, %i ' % (exp(-filler_mux-optimal_mux0)*photons_per_second, exp(-filler_mux-optimal_mux)*photons_per_second))
print(u'    Optimal transmitted/direct beam measurement ratio        %.2f' % (time_ratio))
print(u'    Acquisition times per point to measure Δµx with %.2f ' % (rel_delta_mux*100) + '% accuracy: ')
print(u'        Direct beam:       %.2f s ' % (T0))
print(u'        Transmitted beam:  %.2f s ' % (T0*time_ratio))
print(u'        Total:             %.2f s ' % (T0 + T0*time_ratio))
print(u'    Mass per sample holder cross-sectional area:             %.2f mg/cm²' % (sample_mass2area*1000)) 
print(u'    Mass for standard ISO 7089 M8 washer:                    %.2f mg' % (sample_mass2area*1000*pi*(0.84/2)**2)) 
