#ToDo: Powerlaw intrinseca compatibile con fermi, in particolare se aggiungo hess e veritas deve venire una power law, quando poi passi a broken power law aumenterà l'errore ma il fit verrà meglio. Cerca se qualcuno ha provato a stimare il redshift coi dati ottici. 
#ToDo2: Fai una curva luce annuale e semestrale; estrapola il fit alla power law di fermi
#gammapy fluxpointfitter gammapy.spectrum.model e via, per fare i fit combinati
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
from gammapy.spectrum import FluxPoints
from astropy.table import Table
from gammapy.spectrum.models import PowerLaw

savedir='/home/kibuzo/TesiMagistrale/TesiMia/plot/'
beta=6 #Activate last mod
SED=0  #Plot SED (if 1) or Flux (if 0)
ABSORBED=0 #Plot Absorbed model and pointa if 1, deabsorbed points and model if 0
###########################################################################################################
#------------------------------Blocco 0: Funzioni ausiliarie----------------------------------------------#
###########################################################################################################
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
        #a.append(flux[j]*np.exp(model.spectral_index(energy[j]*u.GeV,epsilon=1e-05)))
        a.append(flux[j]/model.evaluate(energy[j]*u.GeV,1))
    return(a)

def assorbi (energy,flux,model):
    a=[]
    for j in range (0,len(energy)):
        #a.append(flux[j]*np.exp(-model.spectral_index(energy[j]*u.GeV,epsilon=1e-5)+0.))
        a.append(flux[j]*model.evaluate(energy[j]*u.GeV,1))
    return(a)

def PowerFermi(Energy,Gamma,A):
    E0=0.1
    #Gamma=1.76
    return A*((Energy/E0)**(-Gamma))
    
def NewPowerLaw(Energy,Gamma,A):
    Ecut=1e5
    return A*(Energy**(-Gamma)*np.exp(-Energy/Ecut))
    
def FindBreak(Energy, Gamma, Ebreak):
    A=3.379e-08  #output from PowerFermi
    return A*(Energy**(-Gamma))*(np.exp(-Energy/Ebreak))
    
def FindBreakgamma(Energy, Gamma, Ebreak):
    A=1.5233713048782858e-08#3.478585240408327e-08  #output from PowerFermi
    return A*(Energy**(-Gamma))*(np.exp(-Energy/Ebreak))

def PowerFit(Energy,Prefactor,Index,Scale):
    return (Prefactor*(Energy/Scale)**Index)

def BrokenPowerFit(Energy,Prefactor,Index1,Break,Index2): #Thanks stackexchange
    y1=(Prefactor*(Energy/Break)**Index1) [Energy<Break]
    y2=(Prefactor*(Energy/Break)**Index2) [Energy>=Break]
    y=np.concatenate((y1,y2))
    return(y)

    
def propagaspettro(spettro,energia,dspettro,denergia):
    return (np.sqrt((energia**4)*(dspettro**2)+4*(spettro**2)*(energia**2)*(denergia**2)+4*(dspettro*denergia*spettro*(energia)**3)))
    
def proiettax(f,N0,min,max):
    Gamma=-1.76
    E0=15794.2
    inte=(max-min)
    f=f/inte #not so sure.
    return (E0*(f/N0)**(1/Gamma))

#A quanro pre possono tornare comode anche queste
def ErgToGeV(e):
    return(e*624.151)

def GeVToErg(e):
    return (e/624.151)
    
    
################################################################################################################
#----------------------------------------Blocco 1: Definizione dei vettori-------------------------------------#
################################################################################################################

if (beta==6):
    #$FGL psc_v18 Extended sources P8R3IRF
    #nphot=[2.03e-11,1.36e-11,9.09e-12,6.09e-12,4.077e-12,2.73e-12,1.82e-12]
    nphot=[2.07e-11,5.69e-12,2.11e-12,7.3e-13,2.688e-13,2.345e-14,9.9e-15]
    Fermix=np.array([6.94,13.41,25.9,50,96.54,186.38,359.84]) #in GeV
    Fermixmin=np.array([5,9.65,18.6,35.98,69.47,134.13,258.97])
    Fermixmax=np.array([9.65,18.6,35.98,69.47,134.13,258.97,500])
    #dNphot=np.array([5.075e-12,3.62e-12,2.54e-12,1.93e-12,1.51e-12,2.63e-12,1.86e-12])
    dNphot=np.array([4.5e-12,1.52e-12,5.9e-13,2.31e-13,9.95e-14,2.26e-14/2,1.01e-14/2])
    Fermidx=(Fermixmax-Fermixmin).tolist()
    fermiuplims=np.array([0,0,0,0,0,1,1])
    Fermits=[39.8573225013606,36.8792876649022,46.0275764309645,46.3142510222815,37.0156520052533,4.5468426120151,5.7454650115202]
    #nphot=(nphot/np.array(Fermidx)).tolist()
    #dNphot=dNphot/np.array(Fermidx)

if (beta==5):
    #catalogo 3FGL
    nphot=[1.29e-10,7.26e-11,4.8e-11,2.68e-11]
    Fermix=np.array([7.706, 18.094, 42.446, 99.745]) #in GeV
    Fermixmin=np.array([5,11.7,27.5,64.5])
    Fermixmax=np.array([11.7,27.5,64.5,151.514])
    dNphot=2*np.array([1.7e-11, 1.26e-11, 1.03e-11, 0.72e-11])
    Fermidx=(Fermixmax-Fermixmin).tolist()
    fermiuplims=np.array([0,0,0,0])
    nphot=(nphot/np.array(Fermidx)).tolist()
    dNphot=dNphot/np.array(Fermidx)

if (beta==4):
    #Solo la madonna sa perché gli errori sono così grandi, sarei per buttarlo via
    nphot=[2.77255134855615E-11, 1.85600030994195E-11, 1.24249584303576E-11, 8.31810116897444E-12, 5.56830024294583E-12,3.72779394149773E-12, 2.49545845721579E-12]
    Fermix=np.array([4.168, 8.048, 15.53, 30, 57.92, 118.27, 215.9]) #in GeV
    Fermixmin=np.array([3, 5.79, 11.18, 21.59, 41.68, 80.48,155.38])
    Fermixmax=np.array([5.79, 11.18, 21.59, 41.68, 80.48, 155.38,300])
    dNphot=2*np.array([1.74e-11, 1.38e-11, 1.03e-11, 0.72e-11, 0.66e-11, 0.62e-11,0.62e-11])
    Fermidx=(Fermixmax-Fermixmin).tolist()
    fermiuplims=np.array([0,0,0,0,0,0,0])
    nphot=(nphot/np.array(Fermidx)).tolist()
    dNphot=dNphot/np.array(Fermidx)
    
if (beta==3):
    #Dati fermi più bin
    nphot=[6.639E-11, 4.73191E-11, 3.59308e-11, 2.09613e-11, 1.6799e-11,1.59752e-11]
    Fermix=np.array([7.239, 12.450, 22.701, 38.490, 60.382, 139.266]) #in GeV
    Fermixmin=np.array([5.792, 9., 17., 30., 49., 74.])
    Fermixmax=np.array([9., 17., 30., 49., 74., 250.])
    dNphot=2*np.array([1.74e-11, 1.38e-11, 1.03e-11, 0.72e-11, 0.66e-11, 0.62e-11])
    Fermidx=(Fermixmax-Fermixmin).tolist()
    fermiuplims=np.array([0,0,0,0,0,0])
    nphot=(nphot/np.array(Fermidx)).tolist()
    dNphot=dNphot/np.array(Fermidx)

if (beta==2):
    #Dati fermi più bin
    nphot=[6.7704E-11,4.59135E-11,3.58885e-11,2.13419e-11,1.68467e-11,1.59015e-11]
    Fermix=np.array([7.234,12.617,22.701,38.490,60.381,139.266]) #in GeV
    Fermixmin=np.array([5.792,9.,17.,30.,49.,74.])
    Fermixmax=np.array([9.,17.,30.,49.,74.,250.])
    dNphot=np.array([1.7e-11,1.4e-11,0.75e-11,0.5e-11,0.6e-11,0.62e-11])
    Fermidx=(Fermixmax-Fermixmin).tolist()
    fermiuplims=np.array([0,0,0,0,0,0])
    nphot=(nphot/np.array(Fermidx)).tolist()
    dNphot=dNphot/np.array(Fermidx)

if (beta==1):
    #Dati fermi più bin
    nphot=[2.216391E-11,6.347918E-12,2.796748E-12,7.835567E-13,7.183973E-14]
    Fermix=np.array([5.688348,11.6045,23.67372,48.2955,98.52509])
    dNphot=np.array([5.767844E-12,1.811284E-12,7.00665E-13,2.487619E-13,5.773711E-14])
    Fermidx=[4,9,18,36,75]
    fermiuplims=np.array([0,0,0,0,0])

if (beta==0):
    #Meno bin ma più significance
    nphot=[3.082728E-11,8.744711E-12,1.742678E-12,4.677651E-13,6.171239E-14]
    dNphot=np.array([0,1.886901E-12,4.574961E-13,1.317449E-13,0])
    Fermix=np.array([4.403179,11.06028,27.78218,69.78568,175.3145]) #In GeV
    Fermidx=[]
    for j in range (1,len(Fermix)):
        Fermidx.append((Fermix[j]-Fermix[j-1])/2)
    Fermidx.append(np.max(Fermix)/2)
    fermiuplims=np.array([1,0,0,0,1])


#--------------Vettori di fermi-------------------------#
Fermiy=[9.577083E-13,1.714129E-12,2.155338E-12,3.650278E-12,3.038579E-12]
Nphot=np.array(nphot) #in fotoni per cm^2 per GeV
dnPhotsloppy=np.array([1.886901E-12,4.574961E-13,1.317449E-13])
Fermixsloppy=np.array([11.06028,27.78218,69.78568])
FermiNphotsloppy=[8.744711E-12,1.742678E-12,4.677651E-13]
Fermidy=np.array([1.309419163218391e-12,3.698684E-13,5.658294E-13,1.028092E-12,1.3563648604797577e-13])
Fermidx=np.array(Fermidx)


#--------------Vettori di veritas-----------------------#
veritasy=[5.626247237386407e-14,2.4498548272157202e-14,7.083053855910684e-15,3.185536918894964e-15,1.432665297776149e-15,3.486116852371258e-16,1.9873216216790803e-16]
Veritasy=np.array(veritasy)
Veritasdy=np.array([9.644669209754113e-12,3.2923184187082306e-12,3.2923184187082306e-12,1.0101356224329272e-12,2.51984282827237e-13,1.1241647549139845e-13,6.404098398059653e-14])*10**-3
Veritasx=np.array([0.21040591909975698,0.297501746151355,0.41945757279190776,0.5930886396919132,0.8385928812545467,1.18235989899921,1.671788208414575])*10**3
Veritasdx=[]
for j in range (1,len(Veritasx)):
    Veritasdx.append((Veritasx[j]-Veritasx[j-1])/2)
Veritasdx.append(np.max(Veritasx)/2)
Veritasdx=np.array(Veritasdx)


#--------------Vettori di Hess--------------------------#
hessy=[4.0e-15,1.5e-15,8.2e-16,2.5e-16,1.1e-16,6.0e-17,1.3e-17,2.2e-17,2.9e-18,2.4e-19,2.4e-18,4.8e-19]
Hessy=np.array(hessy)
Hessx=np.array([0.53,0.69,0.93,1.2,1.7,2.3,3.0,4.0,5.4,7.3,9.8,13])*10**3
Hessdy=np.array([1.95e-12,0.5e-12,0.185e-12,1.15e-13,0.475e-13,2.85e-14,1.245e-14,1.05e-14,0,0,0,0])*10**-3
Hessdx=[]
for j in range (1,len(Hessx)):
    Hessdx.append((Hessx[j]-Hessx[j-1])/2)
Hessdx.append(np.max(Hessdx)/2)
Hessdx=np.array(Hessdx)
Hessuplims=np.array([0,0,0,0,0,0,0,0,1,1,1,1])

#--------------Propagazione degli errori per SED-------#
deltafermi=propagaspettro(Nphot,Fermix,dNphot/2,Fermidx/2) #Corretto
deltaveritas=propagaspettro(Veritasy,Veritasx,Veritasdy/2,Veritasdx/2)
deltahess=propagaspettro(Hessy,Hessx,Hessdy/2,Hessdx/2)


#--------------Modelli di EBL--------------------------#
redshift=0.2
dominguez = Absorption.read_builtin('dominguez').table_model(redshift)
franceschini = Absorption.read_builtin('franceschini').table_model(redshift)
finke = Absorption.read_builtin('finke').table_model(redshift)

#--------------Deassorbimento vettori------------------#
flussotot=nphot+veritasy+hessy
enertot=Fermix.tolist()+Veritasx.tolist()+Hessx.tolist()
Enertot=np.array(enertot)
fra=deassorbi(enertot,flussotot,franceschini)
Fra=np.array(fra)
dom=deassorbi(enertot,flussotot,dominguez)
Dom=np.array(dom)
fin=deassorbi(enertot,flussotot,finke)
Fin=np.array(fin)
Fermideass=deassorbi(Fermixsloppy,FermiNphotsloppy,franceschini)

#-------------Vettore Totale----------------------------#
XFITBREAK=np.array(Veritasx.tolist()+Hessx.tolist())
XTOT=np.array(Fermix.tolist()+XFITBREAK.tolist())
#YTOT=np.array(FermiNphotsloppy+Veritasy.tolist()+Hessy.tolist())
YFITBREAK=deassorbi(XFITBREAK,np.array(Veritasy.tolist()+Hessy.tolist()),dominguez) #Vettore deassorbito dei cherenkov
YTOT=deassorbi(XTOT,np.array(nphot+Veritasy.tolist()+Hessy.tolist()),dominguez)#deassorbito totale
DYBREAK=np.array(Veritasdy.tolist()+Hessdy.tolist())
DYTOT=np.array(dNphot.tolist()+DYBREAK.tolist())

#####################################################################################################################
#---------------------------------------------Blocco 2: Fit---------------------------------------------------------#
#####################################################################################################################

# #-----Simple Power Law-----#
# initial=[3.79285579423e-12,  -1.606,   15.1]
# popt,pcov=curve_fit(PowerFit,Fermixsloppy,Fermideass,sigma=dnPhotsloppy,absolute_sigma= True,p0=None)
# xfermifit=np.logspace(np.log10(np.min(Fermix)), np.log10(np.max(enertot)),1000)
# yfermifit=PowerFit(xfermifit,popt[0],popt[1],popt[2])
# residuals = dnPhotsloppy - PowerFit(Fermixsloppy,popt[0],popt[1],popt[2])
# chi=np.sum(residuals**2/dnPhotsloppy**2)
# p=stats.chi2.pdf(chi , len(residuals)-1)
# print '\n'
# print 'Simple Power Law Fit results:\n'
# print 'Prefactor, Index, Scale:', popt[0], popt[1], popt[2]
# print 'Covariance matrix: \n',pcov
# print 'Chi squared: ', chi
# print 'p-value: ',p

#-------------------------------------------------Find Fermi normalization---------------------------------------#
initial=(1.5,5e-8)
poptnorm,pcov=curve_fit(PowerFermi,Fermix,Nphot,sigma=dNphot,absolute_sigma= True,p0=initial)
#poptnorm,pcov=curve_fit(PowerFermi,Fermix,deassorbi(Fermix,Nphot,dominguez),sigma=dNphot,absolute_sigma= True,p0=initial)
xnormfit=np.logspace(np.log10(np.min(Fermix)), np.log10(np.max(XTOT)),1000)
ynormfit=PowerFermi(xnormfit,poptnorm[0],poptnorm[1])
residuals = Nphot - PowerFermi(Fermix,poptnorm[0],poptnorm[1])
chinorm=np.sum((residuals/dNphot)**2)
p=stats.chi2.pdf(chinorm , len(residuals)-1)
print '\n'
print '\033[1m'+'Fermi band normalization fit result:\n'+'\033[0m'
print 'Index, Normalization', poptnorm[0], poptnorm[1]
print  'Chi squared/ndof: ', chinorm/4
print 'p-value: ',p
print 'Covariance matrix: \n',pcov

# #-----Find break energy------#
# initial=(2,1e4)
# popt,pcov=curve_fit(FindBreak,XFITBREAK[:-4],YFITBREAK[:-4],sigma=DYBREAK[:-4],absolute_sigma= True,p0=initial)
# xbreakfit=np.logspace(np.log10(np.min(Fermix)), np.log10(np.max(XTOT)),1000)
# ybreakfit=NewPowerLaw(xbreakfit,popt[0],popt[1])
# residuals = YFITBREAK[:-4] - FindBreak(XFITBREAK[:-4],popt[0],popt[1])
# chibreak=np.sum((residuals**2)/((DYBREAK[:-4])**2))
# p=stats.chi2.pdf(chibreak , len(residuals)-1)
# print '\n'
# print 'Cherenkov band break energy fit resiults:\n'
# print 'Index, Break energy:', popt[0], popt[1]
# print 'Covariance matrix: \n',pcov
# print  'Chi squared/ndof: ', chibreak/14
# print 'p-value: ',1-p

#-----------------------------------Fit power law ai dati cherenkov deassorbiti--------------------------------#
initial=(1.6,1e-10)
poptpl,pcovpl=curve_fit(PowerFit,XFITBREAK[:-4],YFITBREAK[:-4],sigma=DYBREAK[:-4],absolute_sigma= True,p0=None)
xcherfit=np.logspace(np.log10(np.min(Fermix)), np.log10(np.max(XTOT)),1000)
ycherfit=PowerFit(xcherfit,poptpl[0],poptpl[1],poptpl[2])
residuals = YFITBREAK[:-4] - PowerFit(XFITBREAK[:-4],poptpl[0],poptpl[1],poptpl[2])
chisq=chisquare(YFITBREAK[:-4], PowerFit(XFITBREAK[:-4],poptpl[0],poptpl[1],poptpl[2]))
chi=np.sum(residuals**2/DYBREAK[:-4]**2)
p=stats.chi2.pdf(chi , len(residuals)-1)
print '\n'
print '\033[1m''Power law fit for deabsorbed cherenkov data:\n'+'\033[0m'
print 'Prefactor, Index, Scale:', poptpl[0], poptpl[1], poptpl[2]
print  'Chi squared/ndof: ', chi/(len(YFITBREAK[:-4])-1)
print 'p-value: ',p
print 'Covariance matrix: \n',pcovpl

#-----------------------------------Fit broken power law ai dati cherenkov deassorbiti-------------------------#
initial=(-poptpl[1],1e4)
poptbpl,pcovbpl=curve_fit(FindBreak,XFITBREAK[:-4],YFITBREAK[:-4],sigma=DYBREAK[:-4],absolute_sigma= True,p0=initial)
xcherfitbpl=np.logspace(np.log10(np.min(Fermix)), np.log10(np.max(XTOT)),1000)
ycherfitbpl=FindBreak(xcherfitbpl,poptbpl[0],poptbpl[1])
residuals = YFITBREAK[:-4] - FindBreak(XFITBREAK[:-4],poptbpl[0],poptbpl[1])
chisq=chisquare(YFITBREAK[:-4], FindBreak(XFITBREAK[:-4],poptbpl[0],poptbpl[1]))
chi=np.sum(residuals**2/DYBREAK[:-4]**2)
p=stats.chi2.pdf(chi , len(residuals)-1)
print '\n'
print '\033[1m'+'Broken power law fit for deabsorbed cherenkov data:\n'+'\033[0m'
print 'Gamma,Ebreak:', poptbpl[0], poptbpl[1]
print  'Chi squared/ndof: ', chi/(len(YFITBREAK[:-4])-1)
print 'p-value: ',p
print 'Covariance matrix: \n',pcovbpl

#-----------------------------------Fit broken power law ai dati cherenkov non deassorbiti-------------------------#
initial=(-poptpl[1],1e4)
poptbpl,pcovbpl=curve_fit(FindBreak,XFITBREAK[:-4],assorbi(XFITBREAK[:-4],YFITBREAK[:-4],franceschini),sigma=DYBREAK[:-4],absolute_sigma= True,p0=initial)
xcherfitbpl2=np.logspace(np.log10(np.min(Fermix)), np.log10(np.max(XTOT)),1000)
ycherfitbpl2=FindBreak(xcherfitbpl,poptbpl[0],poptbpl[1])
residuals = YFITBREAK[:-4] - FindBreak(XFITBREAK[:-4],poptbpl[0],poptbpl[1])
chisq=chisquare(YFITBREAK[:-4], FindBreak(XFITBREAK[:-4],poptbpl[0],poptbpl[1]))
chi=np.sum(residuals**2/DYBREAK[:-4]**2)
p=stats.chi2.pdf(chi , len(residuals)-1)
print '\n'
print '\033[1m'+'Broken power law fit for cherenkov data:\n'+'\033[0m'
print 'Gamma,Ebreak:', poptbpl[0], poptbpl[1]
print  'Chi squared/ndof: ', chi/(len(YFITBREAK[:-4])-1)
print 'p-value: ',p
print 'Covariance matrix: \n',pcovbpl

#-------------------------------------------------NewPowerLaw-------------------------------------------#YTOT
initial=(1.6,1e-10)
popt,pcov=curve_fit(NewPowerLaw,XTOT[:-4],YTOT[:-4],sigma=DYTOT[:-4]/np.array((assorbi(XTOT[:-4],YTOT[:-4],dominguez))/np.array(YTOT[:-4])),absolute_sigma= True,p0=initial)
xfermifit=np.logspace(np.log10(np.min(Fermix)), np.log10(np.max(XTOT)),1000)
yfermifit=NewPowerLaw(xfermifit,popt[0],popt[1])
residuals = DYTOT[:-4] - NewPowerLaw(XTOT[:-4],popt[0],popt[1])
chisq=chisquare(YTOT[:-4], NewPowerLaw(XTOT[:-4],popt[0],popt[1]))
chi=np.sum(residuals**2/DYTOT[:-4]**2)
p=stats.chi2.pdf(chi , len(residuals)-1)
print '\n'
print '\033[1m'+'New Power Law Fit results:\n'+'\033[0m'
print 'Index, Normalization:', popt[0], popt[1]
print  'Chi squared/ndof: ', chi/19
print 'p-value: ',p
print 'Covariance matrix: \n',pcov


# # #------------------------------------------Broken Power Law-------------------------------------------#
# initial=(3.79285579423e-12,   -1.606,   100, -2.606)
# poptbreak,pcovbreak=curve_fit(BrokenPowerFit,Fermix,Nphot,sigma=dNphot,p0=initial)
# yfermibroken=BrokenPowerFit(xfermifit,poptbreak[0],poptbreak[1],poptbreak[2],poptbreak[3])
# residualsbreak = dNphot - PowerFit(Fermix,popt[0],popt[1],popt[2])
# chisqbreak = np.sum(residualsbreak**2/BrokenPowerFit(Fermix,poptbreak[0],poptbreak[1],poptbreak[2],poptbreak[3]))
# Chisqbreak=chisqbreak/(len(Fermix)-1)
# print'\n'
# print 'Broken Power Law Fit results'
# print 'Prefactor,Index1,Break,Index2:', poptbreak[0], poptbreak[1], poptbreak[2], poptbreak[3]
# print 'Covariance matrix: \n',pcovbreak
# print 'Chi squared: ',chisqbreak


######################################################################################################################
#-------------------------------------------------Blocco 3: Plot-----------------------------------------------------#
######################################################################################################################


#------------------------------------Modelli assorbiti da EBL-----------------------------------------#
dominguezfermi=assorbi(xfermifit,yfermifit,dominguez)
franceschinifermi=assorbi(xfermifit,yfermifit,franceschini)
finkefermi=assorbi(xfermifit,yfermifit,finke)

if (SED==1 and ABSORBED==1):
    #----------------------------------Plot SED dei modelli fittati assorbiti-----------------------------#
    #plt.plot(xfermifit,franceschinifermi*xfermifit**2,label='Fermi fit + Franceschini EBL')
    #plt.plot(xfermifit,dominguezfermi*xfermifit**2,label='Fermi fit + Dominguez EBL')
    #plt.plot(xfermifit,finkefermi*xfermifit**2,label='Fermi fit + Finke EBL')
    #plt.plot(xbreakfit,assorbi(xbreakfit,ybreakfit,dominguez)*xbreakfit**2,label='break fit')
    #plt.plot(xnormfit,assorbi(xnormfit,ynormfit,dominguez)*xnormfit**2,label='Fermi Normalization fit')
    plt.plot(xcherfit,assorbi(xcherfit,ycherfit,dominguez)*xcherfit**2,label='Cherenkov power law fit (Dominguez absorbed)')
    #plt.plot(xcherfitbpl,assorbi(xcherfitbpl,ycherfitbpl,dominguez)*xcherfitbpl**2,label='Cherenkov broken power law fit (Dominguez absorbed)')
    plt.plot(xfermifit,assorbi(xnormfit,ynormfit,dominguez)*xnormfit**2,label='Fermi band power law fit absorbed')
    # 



    #-----------------------------------Plot SED dei punti sperimentali---------------------------------#
    plt.errorbar(Fermix,nphot*Fermix**2,yerr=deltafermi,marker='s', markersize=2,uplims=fermiuplims,linestyle='',label='Fermi')
    plt.errorbar(Veritasx,Veritasy*Veritasx**2,yerr=deltaveritas,marker='^',linestyle='',label='Veritas',markersize=4)
    plt.errorbar(Hessx,Hessy*Hessx**2,yerr=deltahess,marker='.',uplims=Hessuplims,linestyle='',label='Hess',markersize=4)
    plt.ylabel('E$^2$dN/dE [GeVcm$^{-2}$s$^{-1}$]')
    plt.xscale("log")
    plt.yscale("log")

    plt.legend()
    plt.savefig(savedir + 'SEDabsorbed.png')
    plt.show()

if (SED==0 and ABSORBED==1):
    #-------Plot dei flussi dei punti sperimentali-------#
    plt.errorbar(Fermix,nphot,yerr=dNphot,marker='s', markersize=2,uplims=fermiuplims,linestyle='',label='Fermi')
    plt.errorbar(Veritasx,Veritasy,yerr=Veritasdy,marker='^',linestyle='',label='Veritas')
    plt.errorbar(Hessx,Hessy,yerr=Hessdy,marker='.',uplims=Hessuplims,linestyle='',label='Hess')
    # 

    #-------Plot dei flussi dei modelli assorbiti---------#
    # plt.plot(xfermifit,franceschinifermi,label='Fermi fit + Franceschini EBL')
    # plt.plot(xfermifit,dominguezfermi,label='Fermi fit + Dominguez EBL')
    # plt.plot(xfermifit,finkefermi,label='Fermi fit + Finke EBL')
    #plt.plot(xbreakfit,assorbi(xbreakfit,ybreakfit,dominguez)*xbreakfit**2,label='break fit')
    #plt.plot(xnormfit,assorbi(xnormfit,ynormfit,dominguez)*xnormfit**2,label='Fermi Normalization fit')
    plt.plot(xcherfit,assorbi(xcherfit,ycherfit,dominguez),label='Cherenkov power law fit (Dominguez absorbed)')
    #plt.plot(xcherfitbpl,assorbi(xcherfitbpl,ycherfitbpl,dominguez),label='Cherenkov broken power law fit (Dominguez absorbed)')
    plt.plot(xfermifit,assorbi(xnormfit,ynormfit,dominguez),label='Fermi band power law fit absorbed')
    plt.ylabel('dN/dE [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
    plt.xscale("log")
    plt.yscale("log")

    plt.legend()
    plt.savefig(savedir + 'fluxabsorbed.png')
    plt.show()
    
if (SED==0 and ABSORBED==0):

    #-------Plot dei flussi dei punti dessorbiti-------#
    plt.errorbar(Fermix,deassorbi(Fermix,nphot,dominguez),yerr=dNphot,marker='.', markersize=2,uplims=fermiuplims,linestyle=''
    ,label='Fermi dominguez deabsorbed')
    plt.errorbar(Veritasx,deassorbi(Veritasx,Veritasy,dominguez),yerr=Veritasdy*(deassorbi(Veritasx,Veritasy,dominguez)/Veritasy),marker='v',markersize=4,linestyle='',label='Veritas dominguez deabsorbed')
    plt.errorbar(Hessx,deassorbi(Hessx,Hessy,dominguez),yerr=Hessdy*(deassorbi(Hessx,Hessy,dominguez)/Hessy),marker='v', markersize=4,uplims=Hessuplims,linestyle='',label='Hess dominguez deabsorbed')
    # plt.errorbar(Fermix,Nphot,yerr=dNphot,marker='.',uplims=fermiuplims,linestyle=''
    # ,label='Fermi dominguez deabsorbed')
    # plt.errorbar(Veritasx,Veritasy,yerr=Veritasdy,marker='v',markersize=4,linestyle='',label='Veritas dominguez deabsorbed')
    # plt.errorbar(Hessx,Hessy,yerr=Hessdy,marker='v', markersize=4,uplims=Hessuplims,linestyle='',label='Hess dominguez deabsorbed')
    #plt.plot(xcherfit,ycherfit,label='Cherenkov power law fit')
    #plt.plot(xcherfitbpl,ycherfitbpl,label='Cherenkov broken power law fit')
    #plt.plot(xnormfit,ynormfit,label='Fermi band power law fit')
    plt.plot(xfermifit,yfermifit,label='Power law fit to deabsorbed MWL data'+'\n(Gamma='+str(-1.62)+')')
    plt.ylabel('dN/dE [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
    plt.xlabel('E [GeV]')
    plt.xscale("log")
    plt.yscale("log")

    plt.legend(loc='upper right')
    plt.savefig(savedir + 'fluxintrinsic.png')
    plt.show()
    
    
    #-------Plot dei flussi dei modelli intrinseci---------#
    # plt.plot(xfermifit,franceschinifermi,label='Fermi fit + Franceschini EBL')
    # plt.plot(xfermifit,dominguezfermi,label='Fermi fit + Dominguez EBL')
    # plt.plot(xfermifit,finkefermi,label='Fermi fit + Finke EBL')
    #plt.plot(xbreakfit,assorbi(xbreakfit,ybreakfit,dominguez)*xbreakfit**2,label='break fit')
    #plt.plot(xnormfit,assorbi(xnormfit,ynormfit,dominguez)*xnormfit**2,label='Fermi Normalization fit')
    # plt.plot(xcherfit,ycherfit,label='Cherenkov power law fit')
    # plt.plot(xcherfitbpl,ycherfitbpl2,label='Cherenkov broken power law fit')
    # plt.ylabel('dN/dE [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
    #plt.xscale("log")
    #plt.yscale("log")

    # plt.legend()
    # plt.show()
    

    #plt.savefig(savedir + 'FermiFit.png')

    



