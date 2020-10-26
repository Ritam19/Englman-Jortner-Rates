######################################################################################################
#Englman-Jortner Rates Calculation                                                                   #
#(released 20/09/2020)                                                                               #
#                                                                                                    #
#How to run it?                                                                                      #
#                                                                                                    #
#Input:                                                                                              #
#                                                                                                    #
#All the input file are given, only the values should be modified:                                   #
#                                                                                                    #
#1-All the normal modes should be writing in "frequencies.py" file                                   #
#                                                                                                    #
#2-The electronic non-adiabatic coupling strength should be added in "ENCS.py" file                  #
#                                                                                                    #
#3-Run the program by: python EJ_rates.py                                                            #
#                                                                                                    #
#                                                                                                    #
#Output:                                                                                             #
#                                                                                                    #
#1-After execution , an “Output.dat” file is created containing 3 columns:                           #
#                                                                                                    #
#     dE(eV)          W_sc(1/s)          W_wc(1/s)                                                   #
#Where dE: is the range of energies from 0.1 to 2 eV                                                 #  
#      W_sc: is the strong coupling rate (1/s)                                                       #
#      W_wc: is the weak coupling rate (1/s)                                                         #
#                                                                                                    #
#                                                                                                    #
#2-To plot:                                                                                          #                                                                                 
#  -The output graph is including the plotting of dE vs W_sc and W_wc                                #
#                                                                                                    #
#   a)install "matplotlib" library                                                                   #
#   b)use a graphical porgram as gnuplot or Excel                                                    #
#                                                                                                    #
#3-After every execution,change the output file name since it is replaced for every new execution    #
#                                                                                                    #
#                                                                                                    #
######################################################################################################





#!/usr/bin/env python3

import matplotlib.pyplot as plt
import math
import sys
import warnings
import numpy as np

warnings.simplefilter(action='ignore', category=FutureWarning)


# Units
evtocm1=8065.54
autoev = 27.21
autotime=2.4188843265857*(10**-17)

# Printing options for floating point numbers
np.set_printoptions(precision=5, suppress=True)

#-----------#
# Variables #
#-----------#
freq = np.loadtxt('frequencies.py', usecols=(0), unpack=True)
w=(freq/evtocm1)/autoev

Coupling = float(open('ENCS.py').read())
C=(Coupling/evtocm1)/autoev

planck=1
boltz=(8.61733*(10**-5))/autoev

pi=np.pi

delta_default1=0.3
delta_default2=0.2
delta_default3=0.0

#---------------#
#Function to Em #
#---------------#
# Calculation of delta
def is_in_range(x, f1, f2):

    '''
    This function is to calculate delta which is the difference between the minima of the two states
    '''

    if x >= f1 and x < f2:
        return True
    else:
        return False
    
def calc_delta():

    """
    Write here what this function does.
    """

    with open('frequencies.py', 'r') as infile:
        delta = []
        for line in infile:
            line = float(line)
            if is_in_range(line, 1400, 15000):
                delta_i = delta_default1
                delta.append(delta_i)
            elif is_in_range(line, 800, 1400):
                delta_i = delta_default2
                delta.append(delta_i)
            else:
                delta_i = delta_default3
                delta.append(delta_i)
    return delta
delta=calc_delta()


# Convert from list to numpy array
delta = np.array(delta)

# Square directly using broadcasting
delta_sq = delta**2

# Function of Em
def calc_Em():
    Em = 0.5*np.sum(planck*delta_sq*w)
    return Em
Em=calc_Em()

#-------------------#
# Calculation of Ea #
#-------------------#
# Array of energies
dE_ev=np.arange(0.1,2.1,0.1)
dE=dE_ev/autoev

# Ea equation
Ea=((dE-Em)**2)/(4*Em)

#-------------------#
# Calculation of T* #
#-------------------#
# Mean over all the frequencies
w_mean=np.mean(w)

# Effective Temperature
a=1./np.tanh((planck*w_mean)/(2*293*boltz))
b=((1/(2*boltz))*planck*w_mean)
T_eff=a*b

#---------------------------------------#
# Function to calculate strong coupling #
#---------------------------------------#
def calc_W_sc():

    '''
    a: coupling**2/planck
    b: 2* pi /(Energy*boltzman constant* the effective tempreature)
    c: exponential term of -(energy/ the mean of all frequencies* planck constant)
    d: the strong coupling rate in 1/sec unit
    '''

    a=(C**2)/planck
    b=((2*pi)/(Em*boltz*T_eff))**0.5
    c=np.exp(-(Ea/(boltz*T_eff)))
    W_sc=(a*b*c)/autotime
    return W_sc
W_sc=calc_W_sc()

#-------------------------------------#
# Function to calculate weak coupling #
#-------------------------------------#
def is_in_range(x, f1):
    '''
    This function is to calculate delta (the difference between the two minima of two states) but over the high frequencies only (>=3000 cm-1)
    '''
    if x >= f1:
        return True
    else:
        return False
def calc_delta_m():
    with open('frequencies.py', 'r') as infile:
        delta_m=[]
        for line in infile:
            line = float(line)
            if is_in_range(line, 3000):
                delta_x = delta_default1
                delta_m.append(delta_x)
            else:
                pass
    return delta_m
delta_m = calc_delta_m()

# Convert from list to numpy array
delta_m = np.array(delta_m)

# Square directly using broadcasting
delta_m_sq = delta_m**2
            
#Calculation of the highest frequency
freq_index=np.where(freq >= 3000.0)
freq_m=(freq[freq_index]/evtocm1)/autoev
freq_max=np.max(freq_m)

#Calculation of DE_m
def calc_dE_m():
    dE_m = 0.5*np.sum(planck*delta_m_sq*freq_m)
    return dE_m
dE_m=calc_dE_m()

#Calculation of gamma DE
gamma = np.log((dE/dE_m)-1)

# Function of the weak coupling rates
def calc_W_wc():
    '''
    a= coupling*2/planck constant
    b=2*pi/(planck*the highest frequency* the range of energies(0.1 to 2.0 eV)
    c= exponential term of -(Energy/the mean of all frequencies* planck constant)
    d=the range of energies*planck constant*the highest frequency
    f= exponentional term of gamma in function of the energy
    W_wc=calculation of the weak coupling rate in 1/sec unit
    '''
    a=(C**2)/planck
    b=((2*pi)/(planck*freq_max*dE))**0.5
    c= np.exp(-(Em/(planck*w_mean)))
    d=dE*planck*freq_max
    e=-(gamma*d)
    f=np.exp(e)
    W_wc=(a*b*c*f)/autotime
    return W_wc
W_wc=calc_W_wc()


#------------------------------#
# Saving the Outputs in a file #
#------------------------------#
Output=np.column_stack((dE,W_sc,W_wc))
header= "dE(eV)                         W_sc(1/s)                       W_wc(1/s)\n"
results=np.savetxt('Output.dat', Output, fmt='%2.10e', delimiter= "\t\t", header=header)

#-----------#
# plotting  #
#-----------#
plt.plot(dE_ev,W_sc)
plt.plot(dE_ev,W_wc)
plt.yscale('log')
plt.xlabel('dE')
plt.ylabel('W_sc')
plt.show()

