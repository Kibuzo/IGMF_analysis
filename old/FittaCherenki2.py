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
from gammapy.utils.fitting import Fit
from iminuit import Minuit

#ToDo: estrapola la power law di fermi ad alta energia, assorbila e guarda se ti esce fuori la roba che osservi in veritas/hess (Fatto, no.), poi fallo con la logparabola (Fatto, meglio ma ancora no). Per CRPropa però non va bene perché la banda di fermi è influenzata dai secondari, quindi dovrai piuttosto usare i cherenkov deassorbiti come primario per alimentare la simulazione 
#ToDo: Sfuttare correttamente gli upper limit per fare i fit, la keyword dovrebbe essere t-statistic (Forse ci siamo)
#ToDo gravissimo: sembra che le tabelle di dominguez che hai te non funzionino con quelle del sito del tutorial https://docs.gammapy.org/dev/api/gammapy.spectrum.models.Absorption.html: prova a fare lo stesso plot

savedir='/home/kibuzo/TesiMagistrale/Python/Plots/'

redshift=0.2
dominguez = Absorption.read_builtin('dominguez').table_model(redshift)
franceschini = Absorption.read_builtin('franceschini').table_model(redshift)
finke = Absorption.read_builtin('finke').table_model(redshift)

def Proiettaerrore(errorix,erroriy,datix,datiy): #Unused, l'effetto è sempre picolo
    a=[]
    dfdx=(datiy[1]-datiy[0])/(datix[1]-datix[0])
    a.append(erroriy[0]+(dfdx**2)*errorix[0]**2)
    for j in range (1,len(errorix)):
        dfdx=(datiy[j]-datiy[j-1])/(datix[j]-datix[j-1])
        #print (dfdx**2)*errorix[j]**2
        a.append(erroriy[j]+(dfdx**2)*errorix[j]**2)
    return(a)

def LogParabola(Energy,N0,alpha,beta):
    Scale=14.731
    E=Energy/Scale
    return(N0*(E)**(-(alpha+beta*np.log10(E))))

def PowerFit(Energy,Prefactor,Index):
    Scale=14.731
    return (Prefactor*(Energy/Scale)**Index)
    
def deassorbi(energy, flux, model):
    a=[]
    for j in range (0,len(flux)):
        a.append(flux[j]*np.exp(model.spectral_index(energy[j]*u.GeV,epsilon=1e-1)))
    return(a)

def assorbi (energy,flux,model):
    a=[]
    for j in range (0,len(energy)):
        a.append(flux[j]*np.exp(-model.spectral_index(energy[j]*u.GeV,epsilon=1e-10)))
    return(a)

#$FGL psc_v18 Extended sources P8R3IRF
#nphot=[2.03e-11,1.36e-11,9.09e-12,6.09e-12,4.077e-12,2.73e-12,1.82e-12]
nphot=[2.07e-11,5.69e-12,2.11e-12,7.3e-13,2.688e-13,2.345e-14,9.9e-15]
Fermix=np.array([6.94,13.41,25.9,50,96.54,186.38,359.84]) #in GeV
Fermixmin=np.array([5,9.65,18.6,35.98,69.47,134.13,258.97])
Fermixmax=np.array([9.65,18.6,35.98,69.47,134.13,258.97,500])
#dNphot=np.array([5.075e-12,3.62e-12,2.54e-12,1.93e-12,1.51e-12,2.63e-12,1.86e-12])
dNphot=np.array([4.5e-12,1.52e-12,5.9e-13,2.31e-13,9.95e-14,2.26e-14/2,1.01e-14/2])/2
Fermidx=(Fermixmax-Fermixmin).tolist()
fermiuplims=np.array([0,0,0,0,0,1,1])
Fermits=[39.8573225013606,36.8792876649022,46.0275764309645,46.3142510222815,37.0156520052533,4.5468426120151,5.7454650115202]
#nphot=(nphot/np.array(Fermidx)).tolist()
#dNphot=dNphot/np.array(Fermidx)
#Provo a correggere l'errore dei punti con ts bassa con la seguente procedura: assumo che la TS dia il p-value, dal p value ricavo il chi quadro corrispondente che nel mio caso è circa 0.5 per entrambi poi uso la formula per il chi quadro con n gradi di libertà chi**2=(n-1)s**2sigma**-2 dove s è la varianza del sample, quindi l'incertezza che mi dà su nphot e sigma è la varianza vera. Poi sommo tutto in quadratura pregando il signore; non so se va bene e comunque non cambia abbastanza i risultati però è un inizio.
dNphot[5]=np.sqrt(dNphot[5]**2+(2.26e-14)**2)
dNphot[6]=np.sqrt(dNphot[6]**2+(1.01e-14)**2)

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

#--------------Deassorbimento vettori e vettori totali------------------#
flussotot=veritasy+hessy
enertot=Veritasx.tolist()+Hessx.tolist()
Enertot=np.array(enertot)
fra=deassorbi(enertot,flussotot,franceschini)
Fra=np.array(fra)
dom=deassorbi(enertot,flussotot,dominguez)
Dom=np.array(dom)
fin=deassorbi(enertot,flussotot,finke)
Fin=np.array(fin)
XTOT=np.array(Veritasx.tolist()+Hessx.tolist())
DYTOT=np.array(Veritasdy.tolist()+Hessdy.tolist())
XCOMBINED=np.array(Fermix.tolist()+Veritasx.tolist()+Hessx[:-4].tolist())
YCOMBINED=np.array(nphot+Veritasy.tolist()+Hessy[:-4].tolist())
#YCOMBINED=deassorbi(XCOMBINED,np.array(nphot+Veritasy.tolist()+Hessy[:-4].tolist()),dominguez)
DYCOMBINED=np.array(dNphot.tolist()+Veritasdy.tolist()+Hessdy[:-4].tolist())


#-----------------------------------Fit power law ai dati cherenkov deassorbiti--------------------------------#
#Fit Hess
initial=(1.0682018666081963e-16,-2.652093565978135)
poptpl,pcovpl=curve_fit(PowerFit,Hessx[:-4],deassorbi(Hessx[:-4],hessy[:-4],dominguez),sigma=Hessdy[:-4],absolute_sigma= True,p0=initial, maxfev=200000)
xcherfithess=np.linspace(np.min(Veritasx), np.max(Hessx),1000)
ycherfithess=PowerFit(xcherfithess,poptpl[0],poptpl[1])
residuals = hessy[:-4] - PowerFit(Hessx[:-4],poptpl[0],poptpl[1])
chisq=chisquare(hessy[:-4], PowerFit(Hessx[:-4],poptpl[0],poptpl[1]))
chi=np.sum(residuals**2/Hessdy[:-4]**2)
p=stats.chi2.pdf(chi , len(residuals)-1)
print '\n'
print '\033[1m''Power law fit for deabsorbed Hess data:\n'+'\033[0m'
print 'Prefactor, Index, Scale:', poptpl[0], poptpl[1]
print  'Chi squared/ndof: ', chi/(len(hessy[:-4])-1)
print 'p-value: ',p
print 'Covariance matrix: \n',pcovpl

#Fit Veritas
initial=(1.0682018666081963e-09,-2.652093565978135)
poptplV,pcovpl=curve_fit(PowerFit,Veritasx,deassorbi(Veritasx,veritasy,dominguez),sigma=Veritasdy,absolute_sigma= True,p0=initial)
xcherfitveritas=np.linspace(np.min(Veritasx), np.max(Hessx),1000)
ycherfitveritas=PowerFit(xcherfitveritas,poptplV[0],poptplV[1])
residuals = veritasy - PowerFit(Veritasx,poptplV[0],poptplV[1])
chisq=chisquare(veritasy, PowerFit(Veritasx,poptplV[0],poptplV[1]))
chi=np.sum(residuals**2/Veritasdy**2)
p=stats.chi2.pdf(chi , len(residuals)-1)
print '\n'
print '\033[1m''Power law fit for deabsorbed Veritas data:\n'+'\033[0m'
print 'Prefactor, Index, Scale:', poptplV[0], poptplV[1]
print  'Chi squared/ndof: ', chi/(len(veritasy)-1)
print 'p-value: ',p
print 'Covariance matrix: \n',pcovpl

# # #Fit complessivo
# initial=(1.0682018666081963e-09,-2.652093565978135)
# poptpl,pcovpl=curve_fit(PowerFit,XTOT[:-4],deassorbi(XTOT[:-4],flussotot[:-4],dominguez),sigma=DYTOT[:-4],absolute_sigma= True,p0=initial)
# xcherfit=np.linspace(np.min(XTOT), np.max(XTOT),3)
# ycherfit=PowerFit(xcherfit,poptpl[0],poptpl[1])
# residuals = flussotot[:-4] - PowerFit(XTOT[:-4],poptpl[0],poptpl[1])
# chisq=chisquare(flussotot[:-4], PowerFit(XTOT[:-4],poptpl[0],poptpl[1]))
# chi=np.sum(residuals**2/DYTOT[:-4]**2)
# p=stats.chi2.pdf(chi , len(residuals)-1)
# print '\n'
# print '\033[1m''Power law fit for deabsorbed cherenkov data:\n'+'\033[0m'
# print 'Prefactor, Index, Scale:', poptpl[0], poptpl[1]
# print  'Chi squared/ndof: ', chi/(len(flussotot[:-4])-1)
# print 'p-value: ',p
# print 'Covariance matrix: \n',pcovpl

print '\n\033[1m------Running minuit/migrad optimizer for HESS power law------\n'+'\033[0m'
def minim_hess(a,b):
    #return sum(((deassorbi(Hessx[:-4],Hessy[:-4],dominguez)-PowerFit(Hessx[:-4],a,b))**2)/(Hessdy[:-4])**2)
    return sum(((Hessy[:-4]-PowerFit(Hessx[:-4],a,b))**2)/(Hessdy[:-4])**2)
initial=(3.5597952256e-15,-2.68506943882)
m = Minuit(minim_hess, a=initial[0],b=initial[1],error_a=initial[0]/100.,error_b=initial[1]/100.,errordef=1)
m.migrad(precision=1e-09)  # run optimiser
print(m.values)
print(m.errors)
print(m.get_fmin())
print 'Reduced chi squared :'
print m.fval / (len(Hessy[:-4]) - 2)

print '\n\033[1m------Running minuit/migrad optimizer for VERITAS power law------\n'+'\033[0m'
def minim_veritas(a,b):
#    return sum(((deassorbi(Veritasx,Veritasy,dominguez)-PowerFit(Veritasx,a,b))**2)/((Veritasdy)**2))
    return sum(((Veritasy-PowerFit(Veritasx,a,b))**2)/((Veritasdy)**2))
initial=(1.0723380661e-15,-2.95147080071)
mv = Minuit(minim_veritas, a=initial[0],b=initial[1],error_a=initial[0]/100.,error_b=initial[1]/100.,errordef=1)
mv.migrad(precision=1e-09)  # run optimiser
print(mv.values)
print(mv.errors)
print(mv.get_fmin())
print 'Reduced chi squared :'
print mv.fval / (len(Veritasy) - 2)

print '\n\033[1m------Running minuit/migrad optimizer for combined cherenkov spectrum power law------\n'+'\033[0m'
def minim_tot(a,b):
    #return sum(((deassorbi(XTOT,flussotot[:-4],dominguez)-PowerFit(XTOT[:-4],a,b))**2)/((DYTOT[:-4])**2))
    return sum(((flussotot[:-4]-PowerFit(XTOT[:-4],a,b))**2)/((DYTOT[:-4])**2))
initial=(1.0723380661e-15,-2.95147080071)
mt = Minuit(minim_tot, a=initial[0],b=initial[1],error_a=initial[0]/100.,error_b=initial[1]/100.,errordef=1)
mt.migrad(precision=1e-09)  # run optimiser
print(mt.values)
print(mt.errors)
print(mt.get_fmin())
print 'Reduced chi squared :'
print mt.fval / (len(DYTOT[:-4]) - 2)

print '\n\033[1m------Running minuit/migrad optimizer for Fermi power law------\n'+'\033[0m'
def minim_fermi(a,b):
#    return sum(((deassorbi(Veritasx,Veritasy,dominguez)-PowerFit(Veritasx,a,b))**2)/((Veritasdy)**2))
    return sum(((nphot[:-2]-PowerFit(Fermix[:-2],a,b))**2)/((dNphot[:-2])**2))
initial=(1.0723380661e-15,-1.65147080071)
mf = Minuit(minim_fermi, a=initial[0],b=initial[1],error_a=initial[0]/100.,error_b=initial[1]/100.,errordef=1)
mf.migrad(precision=1e-09)  # run optimiser
print(mf.values)
print(mf.errors)
print(mf.get_fmin())
print 'Reduced chi squared :'
print mf.fval / (len(nphot) - 2)

print '\n\033[1m------Running minuit/migrad optimizer for Fermi logparabola------\n'+'\033[0m'
def minim_fermilp(a,b,c):
#    return sum(((deassorbi(Veritasx,Veritasy,dominguez)-PowerFit(Veritasx,a,b))**2)/((Veritasdy)**2))
    return sum(((nphot-LogParabola(Fermix,a,b,c))**2)/((dNphot)**2))
initial=(0.6223380661e-11,1.288,0.344)
mflp = Minuit(minim_fermilp, a=initial[0],b=initial[1],c=initial[2],error_a=initial[0]/100.,error_b=initial[1]/100.,error_c=initial[2]/100.,errordef=1)
mflp.migrad(precision=1e-20)  # run optimiser
print(mflp.values)
print(mflp.errors)
print(mflp.get_fmin())
print 'Reduced chi squared :'
print mflp.fval / (len(nphot) - 3)


print '\n\033[1m------Chi squared of combined fits: from fermi logparabola------'+'\033[0m'
chi2lp=np.sum(((YCOMBINED-LogParabola(XCOMBINED,mflp.values[0],mflp.values[1],mflp.values[2]) )**2)/(19*DYCOMBINED**2))
#chi2lp=np.sum(((YCOMBINED-assorbi(XCOMBINED,LogParabola(XCOMBINED,mflp.values[0],mflp.values[1],mflp.values[2]), dominguez) )**2)/(19*DYCOMBINED**2))
print chi2lp
print '\n\033[1m------Chi squared of combined fits: from fermi power law------'+'\033[0m'
chi2pl=np.sum((YCOMBINED-PowerFit(XCOMBINED,mf.values[0],mf.values[1]))**2/(19*DYCOMBINED**2))
print chi2pl

#plt.errorbar(XTOT, deassorbi(XTOT,flussotot,dominguez), yerr=DYTOT,linestyle='',marker='v', label='Dati cherenkov')
#plt.errorbar(Hessx, deassorbi(Hessx,Hessy,dominguez), yerr=Hessdy,linestyle='',marker='v', label='Dati HESS',markersize=5)
plt.errorbar(Hessx, Hessy, yerr=Hessdy,linestyle='',marker='v', label='Dati HESS',markersize=5)
#plt.errorbar(Veritasx, deassorbi(Veritasx,Veritasy,dominguez), yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITAS',markersize=5)
plt.errorbar(Veritasx, Veritasy, yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITAS',markersize=5)
plt.errorbar(Fermix,nphot,yerr=dNphot,marker='.', markersize=2,uplims=fermiuplims,linestyle=''
    ,label='Fermi')
#------------------Fit scipy--------------#
#plt.plot(xcherfit,ycherfit,label='fit combinato')
#plt.plot(xcherfithess,ycherfithess,label='fit hess')
#plt.plot(xcherfitveritas,ycherfitveritas,label='fit veritas')

#------------------Fit Minuit-------------#
plt.plot(Hessx,PowerFit(Hessx,m.values[0],m.values[1]),label='fit HESSMinuit')
plt.plot(Veritasx,PowerFit(Veritasx,mv.values[0],mv.values[1]),label='fit VERITASMinuit')
#plt.plot(Fermix,PowerFit(Fermix,mf.values[0],mf.values[1]),label='fit FermiMinuit power law')
XTOT=np.linspace(np.min(Fermix),np.max(Hessx))#np.array(Fermix[:-2].tolist()+Veritasx[:-4].tolist()+Hessx.tolist())
plt.plot(XTOT,LogParabola(XTOT,mflp.values[0],mflp.values[1],mflp.values[2]),label='Fermi logparabola estrapolato', linestyle='--')
plt.plot(XTOT,PowerFit(XTOT,mf.values[0],mf.values[1]), label='Fermi power law estrapolato', linestyle='--')
plt.xscale("log")
plt.yscale('log')
plt.ylim(1e-19,1e-10)
plt.xlabel('E[GeV]')
plt.ylabel('dN/dE [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
plt.legend()
plt.savefig(savedir+'noEBL.png')
plt.close()

plt.errorbar(Hessx, deassorbi(Hessx,Hessy,dominguez), yerr=Hessdy,linestyle='',marker='v', label='Dati HESS + dominguez EBL',markersize=5)
#plt.errorbar(Veritasx, deassorbi(Veritasx,Veritasy,dominguez), yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITAS',markersize=5)
plt.errorbar(Veritasx, deassorbi(Veritasx,Veritasy,dominguez), yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITAS + dominguez EBL',markersize=5)
plt.errorbar(Fermix,deassorbi(Fermix,nphot, dominguez),yerr=dNphot,marker='.', markersize=2,uplims=fermiuplims,linestyle=''
    ,label='Fermi analysis + dominguez ebl')
plt.plot(Hessx,PowerFit(Hessx,m.values[0],m.values[1]),label='fit HESSMinuit')
plt.plot(Veritasx,PowerFit(Veritasx,mv.values[0],mv.values[1]),label='fit VERITASMinuit')
#plt.plot(Fermix,PowerFit(Fermix,mf.values[0],mf.values[1]),label='fit FermiMinuit power law')
XTOT=np.linspace(np.min(Fermix),np.max(Hessx))#np.array(Fermix[:-2].tolist()+Veritasx[:-4].tolist()+Hessx.tolist())
plt.plot(XTOT,LogParabola(XTOT,mflp.values[0],mflp.values[1],mflp.values[2]),label='Fermi logparabola estrapolato', linestyle='--')
plt.plot(XTOT,PowerFit(XTOT,mf.values[0],mf.values[1]), label='Fermi power law estrapolato', linestyle='--')
plt.xscale("log")
plt.yscale('log')
plt.ylim(1e-19,1e-10)
plt.xlabel('E[GeV]')
plt.ylabel('dN/dE [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
plt.legend()
plt.savefig(savedir+'dominguez.png')
plt.close()


#--------------Plot di Logparabole estrapolate con e senza assorbimento---------------#
plt.errorbar(Hessx, Hessy, yerr=Hessdy,linestyle='',marker='v', label='Dati HESS',markersize=5)
plt.errorbar(Veritasx, Veritasy, yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITA',markersize=5)
plt.errorbar(Fermix,nphot,yerr=dNphot,marker='.', markersize=2,uplims=fermiuplims,linestyle=''
    ,label='Fermi analysis')
XTOT=np.linspace(np.min(Fermix),np.max(Hessx))#np.array(Fermix[:-2].tolist()+Veritasx[:-4].tolist()+Hessx.tolist())
plt.plot(XTOT,assorbi(XTOT,LogParabola(XTOT,mflp.values[0],mflp.values[1],mflp.values[2]),dominguez),label='Fermi logparabola estrapolato assorbito', linestyle='--')
plt.xscale("log")
plt.yscale('log')
plt.ylim(1e-19,1e-10)
plt.xlabel('E[GeV]')
plt.ylabel('dN/dE [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
plt.legend()
plt.savefig(savedir+'Logparabolassorbita.png')
plt.close()

plt.errorbar(Hessx, Hessy, yerr=Hessdy,linestyle='',marker='v', label='Dati HESS',markersize=5)
plt.errorbar(Veritasx, Veritasy, yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITA',markersize=5)
plt.errorbar(Fermix,nphot,yerr=dNphot,marker='.', markersize=2,uplims=fermiuplims,linestyle=''
    ,label='Fermi analysis')
XTOT=np.linspace(np.min(Fermix),np.max(Hessx))#np.array(Fermix[:-2].tolist()+Veritasx[:-4].tolist()+Hessx.tolist())
plt.plot(XTOT,LogParabola(XTOT,mflp.values[0],mflp.values[1],mflp.values[2]),label='Fermi logparabola estrapolato', linestyle='--')
plt.xscale("log")
plt.yscale('log')
plt.ylim(1e-19,1e-10)
plt.xlabel('E[GeV]')
plt.ylabel('dN/dE [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
plt.legend()
plt.savefig(savedir+'Logparabola.png')
plt.close()