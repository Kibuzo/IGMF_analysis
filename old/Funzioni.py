#coding=utf-8
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from gammapy.spectrum.models import Absorption
import astropy.units as u
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from scipy import stats
#---------------------Funzioni ausiliarie-----------------------------#
def Proiettaerrore(errorix,erroriy,datix,datiy): #Unused, l'effetto è sempre picolo
    a=[]
    dfdx=(datiy[1]-datiy[0])/(datix[1]-datix[0])
    a.append(erroriy[0]+(dfdx**2)*errorix[0]**2)
    for j in range (1,len(errorix)):
        dfdx=(datiy[j]-datiy[j-1])/(datix[j]-datix[j-1])
        print (dfdx**2)*errorix[j]**2
        a.append(erroriy[j]+(dfdx**2)*errorix[j]**2)
    return(a)
        
def deassorbi(energy, flux, model):
    a=[]
    for j in range (0,len(flux)):
        a.append(flux[j]*np.exp(model.spectral_index(energy[j]*u.GeV,epsilon=1e-1)))
    return(a)

def assorbi (energy,flux,model):
    a=[]
    for j in range (0,len(energy)):
        a.append(flux[j]*np.exp(-model.spectral_index(energy[j]*u.GeV,epsilon=1e-1)))
    return(a)

def PowerFermi(Energy,Gamma,A):
    Ecut=0.1
    return A*((Energy/Ecut)**(-Gamma))
    
def NewPowerLaw(Energy,Gamma,A):
    Ecut=1e5
    return A*(Energy**(-Gamma)*np.exp(-Energy/Ecut))
    
def FindBreak(Energy, Gamma, Ebreak):
    A=3.478585240408327e-08  #output from PowerFermi
    return A*(Energy**(-Gamma))*np.exp(-Energy/Ebreak)

def PowerFit(Energy,Prefactor,Index,Scale):
    return (Prefactor*(Energy/Scale)**Index)

def BrokenPowerFit(Energy,Prefactor,Index1,Break,Index2): #Thanks stackexchange
    y1=(Prefactor*(Energy/Break)**Index1) [Energy<Break]
    y2=(Prefactor*(Energy/Break)**Index2) [Energy>=Break]
    y=np.concatenate((y1,y2))
    return(y)

    
def propagaspettro(spettro,energia,dspettro,denergia):
    return (np.sqrt((energia**4)*(dspettro**2)+4*(spettro**2)*(energia**2)*(denergia**2)+4*(dspettro*denergia*spettro*(energia)**3)))
    