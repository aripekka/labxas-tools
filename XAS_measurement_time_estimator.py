#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 11:52:57 2019

@author: aripekka
"""

from __future__ import print_function, division

import sys
import xraylib as xrl
from math import sqrt, log, pi

#python 2 compatibility 
if sys.version_info[0] < 3:
    input = raw_input

_SHELL_MACROS = {'K': xrl.K_SHELL, 'L1' : xrl.L1_SHELL, 'L2' : xrl.L2_SHELL, 
                'L3' : xrl.L3_SHELL, 'M1' : xrl.M1_SHELL, 'M2' : xrl.M2_SHELL, 
                'M3' : xrl.M3_SHELL, 'M4' : xrl.M4_SHELL, 'M5' : xrl.M5_SHELL}
    
_VERSION = '0.1'

print('--------------------------------------------')
print('**  XAS MEASUREMENT TIME ESTIMATOR v.'+ _VERSION +'  **')
print('--------------------------------------------')

#Ask for user input
inputs = {}
inputs['edge_string'] = input('Absorption edge (eg. Cu-K, U-L3): ')
inputs['direct_counts'] = float(input('Direct beam countrate (ph/s): '))
inputs['direct_bg_counts'] = float(input('Background rate for direct beam (ph/s): '))
inputs['transmitted_counts_below'] = float(input('Transmitted count rate BELOW the edge (ph/s): '))
inputs['transmitted_counts_above'] = float(input('Transmitted count rate ABOVE the edge (ph/s): '))
inputs['transmitted_bg_counts'] = float(input('Background rate for the transmitted beam (ph/s): '))

#calculate the optimal ratio between direct and transmitted beam measurements
temp1 = inputs['transmitted_counts_below']/(inputs['transmitted_counts_below'] - inputs['transmitted_bg_counts'])**2
temp2 = inputs['transmitted_counts_above']/(inputs['transmitted_counts_above'] - inputs['transmitted_bg_counts'])**2
temp3 = 2*inputs['direct_counts']/(inputs['direct_counts'] - inputs['direct_bg_counts'])**2

time_ratio = sqrt((temp1 + temp2)/temp3)

def estimate_measurement_time(precision):
    temp1 = 1/(log(inputs['transmitted_counts_below'] - inputs['transmitted_bg_counts']) - log(inputs['transmitted_counts_above'] - inputs['transmitted_bg_counts']))
    temp2 = 1/(inputs['direct_counts'] - inputs['direct_bg_counts'])
    temp3 = sqrt(2*inputs['direct_counts']) * sqrt(1 + time_ratio)

    T0 = (temp1*temp2*temp3)**2/precision**2   
    T = time_ratio*T0    
    T_tot = T0 + T
    
    return T, T0, T_tot

#estimate attenuation coefficients from the countrates
mux_below = -log(inputs['transmitted_counts_below'] - inputs['transmitted_bg_counts']) + log(inputs['direct_counts'] - inputs['direct_bg_counts'])
mux_above = -log(inputs['transmitted_counts_above'] - inputs['transmitted_bg_counts']) + log(inputs['direct_counts'] - inputs['direct_bg_counts'])
deltamux = mux_above - mux_below

time_estimates = {'10%' : estimate_measurement_time(0.10), '5%' : estimate_measurement_time(0.05), '1%' : estimate_measurement_time(0.01), 
                  '0.5%' : estimate_measurement_time(0.005), '0.1%' : estimate_measurement_time(0.001)}

#estimate the amount of element of interest in the sample

#convert the absorption edge string to energy
edge_element, edge_shell = inputs['edge_string'].split('-')
Z = xrl.SymbolToAtomicNumber(edge_element)
edge_energy = xrl.EdgeEnergy(Z,_SHELL_MACROS[edge_shell])

murho = [xrl.CS_Total(Z, edge_energy-0.001),
         xrl.CS_Total(Z, edge_energy+0.001)]

mass_surface_density = deltamux / (murho[1]-murho[0])

print('')
print('OUTPUT:')
print(u'    Estimated sample µx (below edge, above edge, edge step):    %.3f, %.3f, %.3f ' % (mux_below, mux_above, deltamux))
print(u'    Estimated sample µx (below edge, above edge, edge step):    %.3f, %.3f, %.3f ' % (mux_below, mux_above, deltamux))
print(u'    Optimal transmitted-to-direct beam acquisition time ratio:  %.2f' % time_ratio)
print(u'')
print(u'    Estimated measurement times per point for selected statistical accuracies:')
print(u'')
print(u'                               10%    5%     1%     0.5%     0.1%')
print(u'        Direct beam (s):       %.2f, %.2f, %.2f, %.2f, %.2f' %(time_estimates['10%'][1], time_estimates['5%'][1], time_estimates['1%'][1], time_estimates['0.5%'][1], time_estimates['0.1%'][1]))
print(u'        Transmitted beam (s):  %.2f, %.2f, %.2f, %.2f, %.2f' %(time_estimates['10%'][0], time_estimates['5%'][0], time_estimates['1%'][0], time_estimates['0.5%'][0], time_estimates['0.1%'][0]))
print(u'        Total (s):             %.2f, %.2f, %.2f, %.2f, %.2f' %(time_estimates['10%'][2], time_estimates['5%'][2], time_estimates['1%'][2], time_estimates['0.5%'][2], time_estimates['0.1%'][2]))
print(u'')
print(u'    Estimated mass per sample holder cross-sectional area:      %.2f mg/cm²' % (mass_surface_density*1000))
print(u'    Mass (if the holder is a standard ISO 7089 M8 washer):      %.2f mg' % (mass_surface_density*1000*pi*(0.84/2)**2))
    
