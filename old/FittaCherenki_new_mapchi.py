#ACHTUNG: quando calcoli i chi quadri li devi calcolare sullo spettro assorbito, perché x-mu in generale cambia
#ToDo: Sono fortemente convinto che le due forme funzionali per il fit siano entrambe corrette in linea di principio. With that in mind invece che selezionare la fetta col chi quadro minimo calcola la cascade power in tutte le fette, calcola gli intervalli di confidenza in tutte le fette e seleziona la fetta che contiene il minimo della cascade power ammessa dal modello; in questo modo, se pur con larghezze diverse e fette diverse dovresti comunque trovare lo stesso valore della cascade power.
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
from scipy import integrate
from gammapy.spectrum import FluxPoints
from astropy.table import Table
from gammapy.spectrum.models import PowerLaw
from gammapy.utils.fitting import Fit
from iminuit import Minuit
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
Ezero=300.

savedir='/home/kibuzo/TesiMagistrale/Python/Plots/'

redshift=0.2
dominguez = Absorption.read_builtin('dominguez').table_model(redshift)
franceschini = Absorption.read_builtin('franceschini').table_model(redshift)
finke = Absorption.read_builtin('finke').table_model(redshift)

def calcolachi2(o,e,s):
        return np.sum(((o-e)**2)/s**2)
    
#------------------------------New chi2 calculator--------------------------------#
def Mapchi2(E,dy,fmin,fmax,cutmin,cutmax,gammamin,gammamax,stepf,stepcut,stepg): #metti fmin e fmax e stepf
    a=[]
    f_idx=-1
    map=np.zeros((len(np.arange(fmin,fmax,stepf)),len(np.arange (gammamin,gammamax,stepg)),len(np.arange(cutmin,cutmax,stepcut))))
    mappower=np.zeros((len(np.arange(fmin,fmax,stepf)),len(np.arange (gammamin,gammamax,stepg)),len(np.arange(cutmin,cutmax,stepcut))))
    
    for f in np.arange(fmin,fmax,stepf):
        f_idx+=1
        print 'step', f_idx, 'di', np.int((fmax-fmin)/stepf)
        c_idx=-1
        #--------------Scorri il vettore del cutoff---------------------------------#
        for ecut in np.arange(cutmin,cutmax,stepcut):
            c_idx+=1
            g_idx=-1
            #------------Scorri il vettore dell'indice power law--------------------#
            for gamma in np.arange (gammamin,gammamax,stepg):
                g_idx+=1
                #------------------------Chi---------------------------#
                #DYASS=DYCHER*(deassorbi(XCHER,YCHER,franceschini)/YCHER)
                residuals = YCHER - assorbi(XCHER,Cutoff(f,XCHER,gamma,ecut),franceschini)
                chi=np.sum(residuals**2/(DYCHER**2)) #Perché 4DY? perché DY è definito come la semilarghezza della barra di errore cioè sigma/2.
                a.append([f,gamma,ecut,chi])
                map[f_idx,g_idx,c_idx]=chi
                mappower[f_idx,g_idx,c_idx]=cascpower(gamma,f,ecut)[0]
    return(a,map,mappower)
    
def distribgamma(gammamin,gammamax,stepg,N0,Ecut):
    a=[]
    b=[]
    for gamma in np.arange(gammamin,gammamax,stepg):
        residuals=YCHER-assorbi(XCHER,Cutoff(N0,XCHER,gamma,Ecut),franceschini)
        chi=np.sum(residuals**2/DYCHER**2)
        a.append(gamma)
        b.append(chi)
    return a,b
    
#-------------------------Questo fitta N0 (in costruzione)------------------------------------#
def fittabreak(E,dy,cutmin,cutmax,gammamin,gammamax,stepcut,stepg):
    #--------------Scorri il vettore della normalizzazione----------------#
    a=[]
    c_idx=-1
    for cut in np.arange(cutmin,cutmax,stepcut):
        c_idx+=1
        g_idx=-1
        print 'step', c_idx, 'di', np.int((cutmax-cutmin)/stepcut)
        #------------Scorri il vettore dell'indice power law--------------------#
        for gamma in np.arange (gammamin,gammamax,stepg):
            g_idx+=1
            ecut_idx=-1
        #------------------------Fit---------------------------#
            def cutoff(E,f):
                return ((f*(E**-gamma))*np.exp(-E/cut))
            Y=deassorbi(XCHER,YCHER,dominguez) #sto cercando un cutoff intrinseco
            popt,pcov=curve_fit(cutoff,XCHER,Y,sigma=DYCHER*deassorbi(XCHER,YCHER,dominguez)/YCOMBINED)
            residuals = YCHER - assorbi(XCHER,cutoff(XCHER,popt[0]),dominguez)
            chi=np.sum(residuals**2/DYCHER**2)
            a.append([popt[0],gamma,cut,chi/17.])
    return(a)
        

def Cutoff(F0,E,Gamma,Ecut):
    #Se metto e/300 F0 è il flusso a 300 GeV
    return (F0*((E/Ezero)**-Gamma)*np.exp(-E/Ecut))
    
#Cascade power calculation
def cascpower(Gamma,F0,Ecut):
    def cutoff (E,F0,Gamma,Ecut):
        return (1e9*E*(F0*((E/Ezero)**-Gamma)*np.exp(-E/Ecut))) # Questa moltiplicazione 1e9 serve a far tornare una potenza, le unità dell'integranda sono 1/cm²s e l'integrale E/cm²s che dà densità di potenza, quindi in pratica integro una SED nuf(nu)
    power= integrate.quad(cutoff, 100, np.inf,args=(F0,Gamma,Ecut),epsabs=1e-3)
    return (power)
        
def LogParabola(Energy,N0,alpha,beta):
    Scale=14.731
    E=Energy/Scale
    return(N0*(E)**(-(alpha+beta*np.log10(E))))

def PowerFit(Energy,Prefactor,Index):
    Scale=14.731
    return (Prefactor*(Energy/Scale)**Index)

def PowerFit2(Energy,Prefactor,Index):
    return(Prefactor*(Energy**-Index))

def PaoloFit(energy,Prefactor,Index):
    return (Prefactor*(energy/Ezero)**-Index)


def PowerFitAss(Energy,Prefactor,Index):
    Scale=14.731
    return assorbi(Energy,(Prefactor*(Energy/Scale)**Index),dominguez)
    
def deassorbi(energy, flux, model):
    a=[]
    for j in range (0,len(flux)):
        #a.append(flux[j]*np.exp(model.spectral_index(energy[j]*u.GeV,epsilon=1e-05)))
        a.append(flux[j]/model.evaluate(energy[j]*u.GeV,1))
    return(a)

def assorbi (energy,flux,model):
    a=[]
    for j in range (0,len(energy)):
        a.append(flux[j]*model.evaluate(energy[j]*u.GeV,1))
    return(a)
    
def plottahess():
    #plt.loglog(XCOMBINED[:-4],deassorbi(XCOMBINED[:-4],YCOMBINED[:-4],dominguez),linestyle='',marker='v')
    #plt.loglog(XCOMBINED[:-4],PowerFit(XCOMBINED[:-4],mtc.values[0],mtc.values[1]))
    plt.errorbar(XCOMBINED,deassorbi(XCOMBINED,YCOMBINED,dominguez),yerr=DYCOMBINED,linestyle='',marker='v',markersize=2)
    plt.errorbar(Hessx,deassorbi(Hessx,Hessy,dominguez),yerr=Hessdy,linestyle='',marker='.',markersize=2)
    plt.plot(XCOMBINED,PowerFit(XCOMBINED,mtc.values[0],mtc.values[1]))
    plt.yscale('log')
    plt.xscale('log')
    #plt.loglog(XCOMBINED[:-4],PowerFit(XCOMBINED[:-4],4.10372356374e-12,-1.53332641727))
    plt.show()

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

# #--------------Vettori di veritas-----------------------#
# veritasy=[5.626247237386407e-14,2.4498548272157202e-14,7.083053855910684e-15,3.185536918894964e-15,1.432665297776149e-15,3.486116852371258e-16,1.9873216216790803e-16]
# Veritasy=np.array(veritasy)
# Veritasdy=np.array([9.644669209754113e-12,3.2923184187082306e-12,3.2923184187082306e-12,1.0101356224329272e-12,2.51984282827237e-13,1.1241647549139845e-13,6.404098398059653e-14])*10**-3#Non toccare, ho verificato col plot del paper e torna
# Veritasx=np.array([0.21040591909975698,0.297501746151355,0.41945757279190776,0.5930886396919132,0.8385928812545467,1.18235989899921,1.671788208414575])*10**3
# Veritasdx=[]
# for j in range (1,len(Veritasx)):
#     Veritasdx.append((Veritasx[j]-Veritasx[j-1])/2)
# Veritasdx.append(np.max(Veritasx)/2)
# Veritasdx=np.array(Veritasdx)
# 
# 
# #--------------Vettori di Hess--------------------------#
# hessy=[3.92e-15,1.47e-15,7.99e-16,2.48e-16,1.07e-16,5.88e-17,1.27e-17,2.17e-17,2.9e-18,2.4e-19,2.4e-18,4.8e-19]
# Hessy=np.array(hessy)
# Hessx=np.array([0.53, 0.69, 0.93, 1.2, 1.7, 2.3, 3.0, 4.0, 5.4, 7.3, 9.8, 13])*10**3
# Hessdy=np.array([1.95e-12,0.5e-12,0.2e-12,0.9e-13,0.43e-13,2.8e-14,1.245e-14,1.05e-14,0,0,0,0])*10**-3 #Non toccare, ho verificato col plot del paper e torna
# Hessdx=[]
# for j in range (1,len(Hessx)):
#     Hessdx.append((Hessx[j]-Hessx[j-1])/2)
# Hessdx.append(np.max(Hessdx)/2)
# Hessdx=np.array(Hessdx)
# Hessuplims=np.array([0,0,0,0,0,0,0,0,1,1,1,1])

#Vettori di Hess di Paolo
Hessx=np.array([0.538, 0.702, 0.934, 1.248, 1.67, 2.26, 3.018, 4.058, 5.4, 7.3, 9.8, 13])*10**3
Hessy=np.array([1.30311e-11, 7.44939e-12, 6.42965e-12, 2.95444e-12, 1.78543e-12, 1.35665e-12, 4.30136e-13, 1.32585e-12, 2.9e-14, 2.4e-15,2.4e-14,4.8e-15])*10**-3
Hessdy=np.array([2.64e-12, 2.52e-12, 1.46e-12, 1.09e-12, 8.09e-13, 7.42e-13, 5.06e-13, 6.96e-13,0,0,0,0])*10**-3
Hessdy=Hessdy*(assorbi(Hessx,Hessy,franceschini)/Hessy)
Hessy=np.array(assorbi(Hessx,Hessy,franceschini))
hessy=Hessy.tolist()

#Vettori di Veritas di Paolo
Veritasx=np.array([0.21, 0.299, 0.421, 0.592, 0.836, 1.18, 1.673])*10**3
Veritasy=np.array([7.22e-11, 4.12e-11, 1.76e-11, 1.26e-11, 9.87e-12, 4e-12, 3.49e-12])*10**-3
Veritasdy=np.array([1.25e-11, 5.74e-12, 2.42e-12, 2.67e-12, 2.06e-12, 1.27e-12, 1.17e-12])*10**-3
Veritasdy=Veritasdy*(assorbi(Veritasx,Veritasy,franceschini)/Veritasy)
Veritasy=np.array(assorbi(Veritasx,Veritasy,franceschini))
veritasy=Veritasy.tolist()

#--------------Deassorbimento vettori e vettori totali------------------#
flussotot=veritasy+hessy
Flussotot=np.array(flussotot)
#flussotot=deassorbi(np.array(Veritasx.tolist + Hessx.tolist),(veritasy+hessy),dominguez)
enertot=sorted(Veritasx.tolist()+Hessx.tolist())
Enertot=np.array(enertot)
XTOT=np.array(sorted(Veritasx.tolist()+Hessx.tolist()))
DYTOT=np.array(Veritasdy.tolist()+Hessdy.tolist())#ACHTUNG: VA RISCRITTO
XCHER=np.array(Veritasx.tolist()+Hessx[:-4].tolist())
YCHER=np.array(Veritasy.tolist()+Hessy[:-4].tolist())
DYCHER=np.array(Veritasdy.tolist()+Hessdy[:-4].tolist())
XCOMBINED=np.array(Fermix.tolist()[:-2]+Veritasx.tolist()+Hessx[:-4].tolist())
YCOMBINED=np.array(nphot[:-2]+Veritasy.tolist()+Hessy[:-4].tolist()) #ACHTUNG: ho tolto due punti (ULIMS)
#YCOMBINED=deassorbi(XCOMBINED,np.array(nphot+Veritasy.tolist()+Hessy[:-4].tolist()),dominguez)
DYCOMBINED=np.array(dNphot.tolist()[:-2]+Veritasdy.tolist()+Hessdy[:-4].tolist())


#-----------------------------------Fit power law ai dati cherenkov deassorbiti--------------------------------#
#Fit Hess
initial=(1.0682018666081963e-16,-2.652093565978135)
poptpl,pcovpl=curve_fit(PowerFitAss,Hessx[:-4],hessy[:-4],sigma=Hessdy[:-4],absolute_sigma= True,p0=initial)
xcherfithess=np.linspace(np.min(Veritasx), np.max(Hessx),1000)
ycherfithess=PowerFit(xcherfithess,poptpl[0],poptpl[1])
residuals = hessy[:-4] - PowerFit(Hessx[:-4],poptpl[0],poptpl[1])
chisq=chisquare(hessy[:-4], PowerFit(Hessx[:-4],poptpl[0],poptpl[1]))
chi=np.sum(residuals**2/Hessdy[:-4]**2)
p=stats.chi2.pdf(chi , len(residuals)-1)
print('\n')
print('\033[1m''Power law fit for deabsorbed Hess data:\n'+'\033[0m')
print('Prefactor, Index, Scale:', poptpl[0], poptpl[1])
print('Chi squared/ndof: ', chi/(len(hessy[:-4])-1))
print('p-value: ',p)
print('Covariance matrix: \n',pcovpl)

#Fit Veritas
initial=(1.0682018666081963e-09,-2.652093565978135)
poptplV,pcovpl=curve_fit(PowerFit,Veritasx,deassorbi(Veritasx,veritasy,dominguez),sigma=Veritasdy,absolute_sigma= True,p0=initial)
xcherfitveritas=np.linspace(np.min(Veritasx), np.max(Hessx),1000)
ycherfitveritas=PowerFit(xcherfitveritas,poptplV[0],poptplV[1])
residuals = veritasy - PowerFit(Veritasx,poptplV[0],poptplV[1])
chisq=chisquare(veritasy, PowerFit(Veritasx,poptplV[0],poptplV[1]))
chi=np.sum(residuals**2/Veritasdy**2)
p=stats.chi2.pdf(chi , len(residuals)-1)
print('\n')
print('\033[1m''Power law fit for deabsorbed Veritas data:\n'+'\033[0m')
print('Prefactor, Index, Scale:', poptplV[0], poptplV[1])
print('Chi squared/ndof: ', chi/(len(veritasy)-1))
print('p-value: ',p)
print('Covariance matrix: \n',pcovpl)

# #Fit complessivo
initial=(1.0682018666081963e-1,-1.652093565978135)
poptpl,pcovpl=curve_fit(PowerFit,Enertot[:-4],deassorbi(Enertot[:-4],flussotot[:-4],dominguez),sigma=DYTOT[:-4],absolute_sigma= True,p0=initial)
xcherfit=np.linspace(np.min(Enertot), np.max(Enertot),3)
ycherfit=PowerFit(xcherfit,poptpl[0],poptpl[1])
residuals = flussotot[:-4] - PowerFit(Enertot[:-4],poptpl[0],poptpl[1])
chisq=chisquare(flussotot[:-4], PowerFit(Enertot[:-4],poptpl[0],poptpl[1]))
chi=np.sum(residuals**2/DYTOT[:-4]**2)
p=stats.chi2.pdf(chi , len(residuals)-1)
print ('\n')
print ('\033[1m''Power law fit for deabsorbed cherenkov data:\n'+'\033[0m')
print ('Prefactor, Index, Scale:', poptpl[0], poptpl[1])
print ('Chi squared/ndof: ', chi/(len(flussotot[:-4])-1))
print ('p-value: ',p)
print ('Covariance matrix: \n',pcovpl)

print('\n\033[1m------Running minuit/migrad optimizer for HESS power law------\n'+'\033[0m')
def minim_hess(a,b):
    #return sum(((deassorbi(Hessx[:-4],Hessy[:-4],dominguez)-PowerFit(Hessx[:-4],a,b))**2)/((Hessdy[:-4])**2)) fitta male ma per colpa sua, apparentemente
    return np.sum(((Hessy[:-4]-PowerFit(Hessx[:-4],a,b))**2)/((Hessdy[:-4])**2))
initial=(1.0723380661e-15,-2.95147080071)
mh = Minuit(minim_hess, a=initial[0],b=initial[1],error_a=initial[0]/100.,error_b=initial[1]/100.,errordef=1)
mh.migrad(precision=1e-15)  # run optimiser
print((mh.values))
print((mh.errors))
print((mh.get_fmin()))
print('Reduced chi squared :')
print(mh.fval / (len(Hessy[:-4]) - 2))

print('\n\033[1m------Running minuit/migrad optimizer for VERITAS power law------\n'+'\033[0m')
def minim_veritas(a,b):
#    return sum(((deassorbi(Veritasx,Veritasy,dominguez)-PowerFit(Veritasx,a,b))**2)/((Veritasdy)**2))
    return sum(((Veritasy-PowerFit(Veritasx,a,b))**2)/((Veritasdy)**2))
initial=(1.0723380661e-15,-2.95147080071)
mv = Minuit(minim_veritas, a=initial[0],b=initial[1],error_a=initial[0]/100.,error_b=initial[1]/100.,errordef=1)
mv.migrad(precision=1e-09)  # run optimiser
print((mv.values))
print((mv.errors))
print((mv.get_fmin()))
print('Reduced chi squared :')
print(mv.fval / (len(Veritasy) - 2))

print('\n\033[1m------Running minuit/migrad optimizer for combined cherenkov spectrum power law------\n'+'\033[0m')
def minim_tot(a,b):
    return sum(((deassorbi(XCOMBINED,YCOMBINED,dominguez)-PowerFit(XCOMBINED,a,b))**2)/((DYCOMBINED)**2))
initial=(1.0723380661e-12,-1.95147080071)
mtc = Minuit(minim_tot, a=initial[0],b=initial[1],error_a=initial[0]/100.,error_b=initial[1]/100.,errordef=1)
mtc.migrad(precision=1e-9)  # run optimiser
print((mtc.values))
print((mtc.errors))
print((mtc.get_fmin()))
print('Reduced chi squared :')
res=[]
for i in range (0,len(DYCOMBINED[:-4])):
    res.append(PowerFit(XCOMBINED[i],mtc.values[0],mtc.values[1])-YCOMBINED[i])
res=np.array(res)
print(sum(res**2/(14*DYCOMBINED[:-4]**2)))# / (len(DYCOMBINED[:-4]) - 2))
print(mtc.fval / (len(DYCOMBINED) - 2))

print('\n\033[1m------Running minuit/migrad optimizer for Fermi power law------\n'+'\033[0m')
def minim_fermi(a,b):
#    return sum(((deassorbi(Veritasx,Veritasy,dominguez)-PowerFit(Veritasx,a,b))**2)/((Veritasdy)**2))
    return sum(((nphot[:-2]-PowerFit(Fermix[:-2],a,b))**2)/((dNphot[:-2])**2))
initial=(1.0723380661e-15,-1.65147080071)
mf = Minuit(minim_fermi, a=initial[0],b=initial[1],error_a=initial[0]/100.,error_b=initial[1]/100.,errordef=1)
mf.migrad(precision=1e-09)  # run optimiser
print((mf.values))
print((mf.errors))
print((mf.get_fmin()))
print('Reduced chi squared :')
print(mf.fval / (len(nphot) - 2))

print('\n\033[1m------Running minuit/migrad optimizer for Fermi logparabola------\n'+'\033[0m')
def minim_fermilp(a,b,c):
#    return sum(((deassorbi(Veritasx,Veritasy,dominguez)-PowerFit(Veritasx,a,b))**2)/((Veritasdy)**2))
    return sum(((nphot-LogParabola(Fermix,a,b,c))**2)/((dNphot)**2))
initial=(0.6223380661e-11,1.288,0.344)
mflp = Minuit(minim_fermilp, a=initial[0],b=initial[1],c=initial[2],error_a=initial[0]/100.,error_b=initial[1]/100.,error_c=initial[2]/100.,errordef=1)
mflp.migrad(precision=1e-20)  # run optimiser
print((mflp.values))
print((mflp.errors))
print((mflp.get_fmin()))
print('Reduced chi squared :')
print(mflp.fval / (len(nphot) - 3))


print('\n\033[1m------Chi squared of combined fits: from fermi logparabola------'+'\033[0m')
#chi2lp=np.sum(((YCOMBINED-LogParabola(XCOMBINED,mflp.values[0],mflp.values[1],mflp.values[2]) )**2)/(19*DYCOMBINED**2))
chi2lp=np.sum(((YCOMBINED-assorbi(XCOMBINED,LogParabola(XCOMBINED,mflp.values[0],mflp.values[1],mflp.values[2]), dominguez) )**2)/(19*DYCOMBINED**2))
print(chi2lp)
print('\n\033[1m------Chi squared of combined fits: from fermi power law------'+'\033[0m')
chi2pl=np.sum((YCOMBINED-assorbi(XCOMBINED,PowerFit(XCOMBINED,mf.values[0],mf.values[1]),dominguez))**2/(19*DYCOMBINED**2))
print(chi2pl)

#plt.errorbar(XTOT, deassorbi(XTOT,flussotot,dominguez), yerr=DYTOT,linestyle='',marker='v', label='Dati cherenkov')
#plt.errorbar(Hessx, deassorbi(Hessx,Hessy,dominguez), yerr=Hessdy,linestyle='',marker='v', label='Dati HESS',markersize=5)
plt.errorbar(Hessx, Hessy, yerr=Hessdy,linestyle='',marker='v', label='Dati HESS',markersize=5)
#plt.errorbar(Veritasx, deassorbi(Veritasx,Veritasy,dominguez), yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITAS',markersize=5)
plt.errorbar(Veritasx, Veritasy, yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITAS',markersize=5)
plt.errorbar(Fermix,nphot,yerr=dNphot,marker='.', markersize=2,uplims=fermiuplims,linestyle=''
    ,label='Fermi')

#------------------Fit Minuit-------------#
plt.plot(Hessx,PowerFit(Hessx,mh.values[0],mh.values[1]),label='fit HESSMinuit')
plt.plot(Veritasx,PowerFit(Veritasx,mv.values[0],mv.values[1]),label='fit VERITASMinuit')
XTOT=np.linspace(np.min(Fermix),np.max(Hessx))
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
plt.errorbar(Veritasx, deassorbi(Veritasx,Veritasy,dominguez), yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITAS + dominguez EBL',markersize=5)
plt.errorbar(Fermix,deassorbi(Fermix,nphot, dominguez),yerr=dNphot,marker='.', markersize=2,uplims=fermiuplims,linestyle=''
    ,label='Fermi analysis + dominguez ebl')
plt.plot(Hessx,PowerFit(Hessx,mh.values[0],mh.values[1]),label='fit HESSMinuit')
plt.plot(Veritasx,PowerFit(Veritasx,mv.values[0],mv.values[1]),label='fit VERITASMinuit')
XTOT=np.linspace(np.min(Fermix),np.max(Hessx))
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

#--------------Plot di power law estrapolata da fermi e assorbita---------------------#
plt.errorbar(Hessx, Hessy, yerr=Hessdy,linestyle='',marker='v', label='Dati HESS + dominguez EBL',markersize=5)
plt.errorbar(Veritasx, Veritasy, yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITAS + dominguez EBL',markersize=5)
plt.errorbar(Fermix,nphot,yerr=dNphot,marker='.', markersize=2,uplims=fermiuplims,linestyle=''
    ,label='Fermi analysis + dominguez ebl')
XTOT=np.linspace(np.min(Fermix),np.max(Hessx))
plt.plot(XTOT,assorbi(XTOT,PowerFit(XTOT,mf.values[0],mf.values[1]),dominguez), label='Fermi power law estrapolato + dominguez EBL', linestyle='--')
plt.xscale("log")
plt.yscale('log')
plt.ylim(1e-19,1e-10)
plt.xlabel('E[GeV]')
plt.ylabel('dN/dE [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
plt.legend()
plt.savefig(savedir+'PowerlawAssorbita.png')
plt.close()

#--------------Plot di Logparabole estrapolate con e senza assorbimento---------------#
plt.errorbar(Hessx, Hessy, yerr=Hessdy,linestyle='',marker='v', label='Dati HESS',markersize=5)
plt.errorbar(Veritasx, Veritasy, yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITA',markersize=5)
plt.errorbar(Fermix,nphot,yerr=dNphot,marker='.', markersize=2,uplims=fermiuplims,linestyle=''
    ,label='Fermi analysis')
XTOT=np.linspace(np.min(Fermix),np.max(Hessx))
plt.plot(XTOT,assorbi(XTOT,LogParabola(XTOT,mflp.values[0],mflp.values[1],mflp.values[2]),dominguez),label='Fermi logparabola estrapolato assorbito', linestyle='--')
plt.xscale("log")
plt.yscale('log')
plt.ylim(1e-19,1e-10)
plt.xlabel('E[GeV]')
plt.ylabel('dN/dE [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
plt.legend()
plt.savefig(savedir+'Logparabolassorbita.png')
plt.close()
#senza
plt.errorbar(Hessx, Hessy, yerr=Hessdy,linestyle='',marker='v', label='Dati HESS',markersize=5)
plt.errorbar(Veritasx, Veritasy, yerr=Veritasdy,linestyle='',marker='v', label='Dati VERITA',markersize=5)
plt.errorbar(Fermix,nphot,yerr=dNphot,marker='.', markersize=2,uplims=fermiuplims,linestyle=''
    ,label='Fermi analysis')
XTOT=np.linspace(np.min(Fermix),np.max(Hessx))
plt.plot(XTOT,LogParabola(XTOT,mflp.values[0],mflp.values[1],mflp.values[2]),label='Fermi logparabola estrapolato', linestyle='--')
plt.xscale("log")
plt.yscale('log')
plt.ylim(1e-19,1e-10)
plt.xlabel('E[GeV]')
plt.ylabel('dN/dE [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
plt.legend()
plt.savefig(savedir+'Logparabola.png')
plt.close()

# #------------------Parametri della mappa del chi quadro---------------------#
#step=(2e-12,0.005,2e2)#(5e-11,5e-3,5e2)
fval=(3.6e-14,4.2e-14)#(2e-10,3e-10)
gammaval=(0.8,1.9)#(1.56,1.63)
cutoff=(1e2,5e3)#[1e3,20e3]
npixel=100
nslices=10
step=(np.abs(fval[1]-fval[0])/nslices,np.abs(gammaval[0]-gammaval[1])/npixel,np.abs(cutoff[1]-cutoff[0])/npixel)


Mchi2,mappa,power=Mapchi2(XCHER,DYCHER,fval[0],fval[1],cutoff[0],cutoff[1],gammaval[0],gammaval[1],step[0],step[2],step[1])
Mchi2=np.transpose(np.array(Mchi2))

#------remember  fittabreak(E,dy,fmin,fmax,gammamin,gammamax,stepf,stepg)----#
#b=fittabreak(XCOMBINED,DYCOMBINED,cutoff[0],cutoff[1],gammaval[0],gammaval[1],step[2],step[1])
#b=fittabreakold(XCOMBINED,DYCOMBINED,fval[0],fval[1],gammaval[0],gammaval[1],step[0],step[1])
#print 'min chi squared',np.min(np.array(b)[:,3])
#row=np.where(np.array(b)[:,3] == np.min(np.array(b)[:,3]))[0][0]
#print 'min chi squared parameters E0, gamma, break', b[row]

#------------Diagnostic plots-------------#
def plottapuntide():
    Y=deassorbi(XCOMBINED,YCOMBINED,dominguez)
    DY=DYCOMBINED*deassorbi(XCOMBINED,YCOMBINED,dominguez)/YCOMBINED
    plt.errorbar(XCOMBINED, Y, yerr=DY,linestyle='',marker='v',markersize=2)

def plottacutoff(f,gamma,ebreak):
    plt.plot(XCOMBINED,Cutoff(f,XCOMBINED,gamma,ebreak))

def plottares(X,Y,YERR,fit):
    plt.errorbar(X,(fit-Y)/YERR, yerr=YERR, linestyle='',marker='v', markersize=2)

def plottaresidui():
    a=b
    plottares(XCOMBINED,deassorbi(XCOMBINED,YCOMBINED,dominguez),DYCOMBINED*(deassorbi(XCOMBINED,YCOMBINED,dominguez)/YCOMBINED),Cutoff(a[row][0],XCOMBINED,a[row][1],a[row][2]))
    plt.show()

# plt.xscale("log")
# plt.yscale('log')
# 
# plottapuntide()
# plottacutoff(b[row][0],b[row][1],b[row][2])
# plt.xscale("log")
# plt.yscale('log')
# plt.show()

# #How to plot the 2d histogram
# Define the position on which to slice the array in f0. There must be an easier way to do it but this appears to work
f0=np.int(np.round((np.array(Mchi2)[0][np.where(np.array(Mchi2)[3] == np.min(np.array(Mchi2)[3]))[0][0]]-fval[0])/step[0]))

# Slice the array: here i'm taking the value for E0 which gives the minimum chi² in the array
plt.imshow(mappa[f0,:,:],cmap='viridis',norm=LogNorm()) 


# Here i'm reshaping the axis labels. Since the array is made of integer indices i have to rescale them to show the correct values of gamma, Ecut. This necessarily depends on how you made the grid, so change it accordingly
plt.yticks(plt.yticks()[0][1:],np.round(np.linspace(gammaval[0],gammaval[1],num=len(plt.yticks()[0][1:])),decimals=3))
plt.xticks(plt.xticks()[0][1:],np.round(np.linspace(cutoff[0],cutoff[1],num=len(plt.xticks()[0][1:])),decimals=0))

#mappa[f0,:,:][mappa[f0,:,:]>18.03]=-1
plt.imshow(mappa[f0,:,:],cmap='viridis',norm=LogNorm()) 
plt.colorbar(ticks=[11.8,15.27,18])
plt.xlabel('Ecut [GeV]')
plt.ylabel('Spectral index ($\Gamma$)')
plt.title('Chi squared of fits in Gamma/Ecut parameters space')
# Finally show the plot and possibly be happy
plt.savefig(savedir+'Chiquadri.png')
#Create contour plot for future reference
g=np.arange(gammaval[0],gammaval[1],step=step[1])
cut=np.arange(cutoff[0],cutoff[1],step=step[2])
x,y=np.meshgrid(cut,g)
plt.close()

# Same for the 2d histogram of cascade power
minchi=np.min(Mchi2[3])
f0=np.int(np.round((np.array(Mchi2)[0][np.where(np.array(Mchi2)[3] == np.min(np.array(Mchi2)[3]))[0][0]]-fval[0])/step[0]))

# Slice the array: here i'm taking the value for E0 which gives the minimum chi² in the array
plt.imshow(power[f0,:,:],cmap='viridis') 
##-------------Print the minimum of the cascade power---------------------##
#1- slice the vectors around minimum E0
powerslice=power[f0,:,:]
chislice=mappa[f0,:,:]
#Search for the minimum of cascade power within the accepted chi squared values 
minimo=np.min(powerslice[np.where(chislice<minchi+6.25)])
#localize the point in the power cascade matrix
dove=np.where(power==minimo)
print ('minimum cascade power accepted by data=',np.round(minimo, decimals=2),'in coordinates N0=',f0*step[0]+fval[0],'gamma=',dove[1][0]*step[1]+gammaval[0],'Ecut=',dove[2][0]*step[2]+cutoff[0])
#Plot vertical lines
plt.plot((dove[2][0],dove[2][0]),(0,dove[1][0]),color='k',linestyle='--')
plt.plot((0,dove[2][0]),(dove[1][0],dove[1][0]),color='k',linestyle='--')


# Here i'm reshaping the axis labels. Since the array is made of integer indices i have to rescale them to show the correct values of gamma, Ecut. This necessarily depends on how you made the grid, so change it accordingly
#plt.yticks(plt.yticks()[0][1:],np.round(2-np.round(np.linspace(gammaval[0],gammaval[1],num=len(plt.yticks()[0][1:])),decimals=3),decimals=3)) #Plot con le y scalate
plt.yticks(plt.yticks()[0][1:],np.round(np.round(np.linspace(gammaval[0],gammaval[1],num=len(plt.yticks()[0][1:])),decimals=3),decimals=3))
plt.xticks(plt.xticks()[0][1:],np.round(np.linspace(cutoff[0],cutoff[1],num=len(plt.xticks()[0][1:])),decimals=0))

plt.imshow(power[f0,:,:],cmap='viridis',norm=LogNorm()) 
cbar=plt.colorbar()
cbar.set_label('Cascade power [eV/cm$^{2}$s]',rotation=90)
cs=plt.contour(mappa[f0,:,:],levels=(minchi+3.5,minchi+6.25))
plt.clabel(CS=cs)
plt.xlabel('Ecut [GeV]')
plt.ylabel('Spectral index ($\Gamma$)')
plt.title('Cascade power in $\Gamma$/Ecut parameters space')
# Finally show the plot and possibly be happy
plt.savefig(savedir+'Powers.png')
plt.show()

# # Funzione usata per manda il plot a paolo
# DYASS=DYCHER*(deassorbi(XCHER,YCHER,franceschini)/YCHER)
# poptpl,pcovpl=curve_fit(PaoloFit,XCHER,deassorbi(XCHER,YCHER,franceschini),sigma=DYASS)
# xcherfit=np.linspace(np.min(XCHER), np.max(XCHER),3)
# ycherfit=PaoloFit(XCHER,poptpl[0],poptpl[1])
# residuals = deassorbi(XCHER,YCHER,dominguez) - PowerFit2(XCHER,poptpl[0],poptpl[1])
# chisq=chisquare(deassorbi(XCHER,YCHER,dominguez), PowerFit2(XCHER,poptpl[0],poptpl[1]))
# chi=np.sum(residuals**2/DYASS**2)
# p=stats.chi2.pdf(chi , len(residuals)-1)
# print ('\n')
# print ('\033[1m''Power law fit for deabsorbed cherenkov data:\n'+'\033[0m')
# print ('Prefactor, Index, Scale:', poptpl[0], poptpl[1])
# print ('Chi squared/ndof: ', chi/(len(flussotot[:-4])-1))
# print ('p-value: ',p)
# print ('Covariance matrix: \n',pcovpl)
# plt.errorbar(XCHER, deassorbi(XCHER,YCHER,franceschini), yerr=DYASS/2,linestyle='',marker='v', label='Dati HESS',markersize=5)
# plt.plot(XCHER,ycherfit)
# plt.xscale('log')
# plt.yscale('log')
# plt.show()


