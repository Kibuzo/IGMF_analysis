#ACHTUNG: quando calcoli i chi quadri li devi calcolare sullo spettro assorbito, perché x-mu in generale cambia
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

#ToDo Urgentissimo: Devi ricontrollare gli errori di fermi e dei cherenkov once and for all. Non ti piace che la power law ottenuta col fit complessivo abbia un chi quadro più piccola di quella di fermi assorbito.
#ToDo: Trovare il break, see Da Vela (mail del 24 giugno), dovrebbe essere il metodo Lampton76. Devo fittare lo spettro con power law e devo trovare un risultato decente, poi mettere il break e trovare il minimo Ecut che mi rimane compatibile con le osservazioni (ipotesi più conservativa)
#toDo: Fit un parametro alla volta, che poi credo che sia sempre lampton76: purtroppo per come funzionano i vari ottimizzatori forse ti conviene fare dei cicli dove ridefinisci le funzioni di fit sbloccando i vari parametri o tirandoli a caso. In totale ne devi fittare 3.

#Prova a scegliere due parametri a mano (in un ciclo) e fittare il terzo con la power law
#Trova la tripla giusta e una volta trovata escluti N0 e plotta eb vs gamma. La cascade power la puoi considerare come l'integrale tra 100 gev e infinito del flusso

savedir='/home/kibuzo/TesiMagistrale/Python/Plots/'

redshift=0.21
dominguez = Absorption.read_builtin('dominguez').table_model(redshift)
franceschini = Absorption.read_builtin('franceschini').table_model(redshift)
finke = Absorption.read_builtin('finke').table_model(redshift)

def calcolachi2(o,e,s):
        return np.sum(((o-e)**2)/s**2)

#Mappachi dovrebbe restituire una mappa di chi quadro per i vari parametri. La mappa è tridimensionale, ma a te interessa solo trovare l'Ecut minimo consistente col 68% di CL per avere l'ipotesi più conservativa sul flusso di gamma secondari a terra (che verrà eventualmente ridotto dal campo magnetico, in questo senso è più conservativo: meno ne prevedi più piccola è la perdita imputabile al campo magnetico.)

#-------------------Old chi2 calculator------------------------------------------#
def Mappachi(E,dy,fmin,fmax,gammamin,gammamax,ecutmin,ecutmax,stepf,stepg,stepec):
    E0=14.731 #Identico alla scale nel powerfit
    map=np.zeros((np.int((fmax-fmin)/stepf)+1,np.int((gammamax-gammamin)/stepg),np.int((ecutmax-ecutmin)/stepec)))
    print 'Fittone gigantico multiparametro è iniziato, mettiti comodo'
    print 'Dimensioni array:', np.int((fmax-fmin)/stepf), np.int((gammamax-gammamin)/stepg), np.int((ecutmax-ecutmin)/stepec)
    f_idx=-1
    deassorbito=deassorbi(XCOMBINED,YCOMBINED,dominguez)
    for f in np.arange(fmin,fmax,stepf):
        f_idx+=1
        g_idx=-1
        ecut_idx=-1
        print 'step', f_idx, 'di', np.int((fmax-fmin)/stepf)
        for gamma in np.arange (gammamin,gammamax-stepg,stepg):
            g_idx+=1
            ecut_idx=-1
            for ecut in np.arange (ecutmin,ecutmax,stepec):
                ecut_idx+=1
                map[f_idx,g_idx,ecut_idx]=calcolachi2(assorbi(XCOMBINED,Cutoff(f,E,gamma,ecut),dominguez),YCOMBINED,dy)
    return (map)
    
#------------------------------New chi2 calculator--------------------------------#
def Mapchi2(E,dy,f,cutmin,cutmax,gammamin,gammamax,stepcut,stepg):
    #--------------Scorri il vettore del cutoff---------------------------------#
    a=[]
    c_idx=-1
    map=np.zeros((np.int((fmax-fmin)/stepf),np.int((gammamax-gammamin)/stepg),np.int((ecutmax-ecutmin)/stepec)))
    for ecut in np.arange(cutmin,cutmax,stepcut):
        c_idx+=1
        g_idx=-1
        print 'step', c_idx, 'di', np.int((cutmax-cutmin)/stepcut)
        #------------Scorri il vettore dell'indice power law--------------------#
        for gamma in np.arange (gammamin,gammamax,stepg):
            g_idx+=1
            ecut_idx=-1
        #------------------------Chi---------------------------#
            Y=deassorbi(XCOMBINED,YCOMBINED,dominguez) #sto cercando un cutoff intrinseco
            residuals = YCOMBINED - assorbi(XCOMBINED,Cutoff(f,XCOMBINED,gamma,ecut),dominguez)
            #chisq=chisquare(Y, cutoff(XCOMBINED,popt[0]))
            chi=np.sum(residuals**2/DYCOMBINED**2)
            a.append([f,gamma,ecut,chi])
            map[f,gamma,ecut]=chi
    return(a,map)
    
#----------------------------Questo fittava ecut, invece devi fittare N0------------------------#
def fittabreakold(E,dy,fmin,fmax,gammamin,gammamax,stepf,stepg):
    #--------------Scorri il vettore della normalizzazione----------------#
    a=[]
    f_idx=-1
    for f in np.arange(fmin,fmax,stepf):
        f_idx+=1
        g_idx=-1
        print 'step', f_idx, 'di', np.int((fmax-fmin)/stepf)
        #------------Scorri il vettore dell'indice power law--------------------#
        for gamma in np.arange (gammamin,gammamax,stepg):
            g_idx+=1
            ecut_idx=-1
        #------------------------Fit---------------------------#
            def cutoff(E,Ecut):
                return ((f*(E**-gamma))*np.exp(-E/Ecut))
            Y=deassorbi(XCOMBINED,YCOMBINED,dominguez) #sto cercando un cutoff intrinseco
            popt,pcov=curve_fit(cutoff,XCOMBINED,Y,sigma=DYCOMBINED*deassorbi(XCOMBINED,YCOMBINED,dominguez)/YCOMBINED)
            residuals = YCOMBINED - assorbi(XCOMBINED,cutoff(XCOMBINED,popt[0]),dominguez)
            #chisq=chisquare(Y, cutoff(XCOMBINED,popt[0]))
            chi=np.sum(residuals**2/DYCOMBINED**2)
            a.append([f,gamma,popt[0],chi/17.])
    return(a)

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
            Y=deassorbi(XCOMBINED,YCOMBINED,dominguez) #sto cercando un cutoff intrinseco
            popt,pcov=curve_fit(cutoff,XCOMBINED,Y,sigma=DYCOMBINED*deassorbi(XCOMBINED,YCOMBINED,dominguez)/YCOMBINED)
            residuals = YCOMBINED - assorbi(XCOMBINED,cutoff(XCOMBINED,popt[0]),dominguez)
            #chisq=chisquare(Y, cutoff(XCOMBINED,popt[0]))
            chi=np.sum(residuals**2/DYCOMBINED**2)
            a.append([popt[0],gamma,cut,chi/17.])
    return(a)
        

def Cutoff(F0,E,Gamma,Ecut):
    # E0=14.731
    # e=E/E0
    return (F0*((E)**-Gamma)*np.exp(-E/Ecut))
        
def LogParabola(Energy,N0,alpha,beta):
    Scale=14.731
    E=Energy/Scale
    return(N0*(E)**(-(alpha+beta*np.log10(E))))

def PowerFit(Energy,Prefactor,Index):
    Scale=14.731
    return (Prefactor*(Energy/Scale)**Index)


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
        #a.append(flux[j]*np.exp(-model.spectral_index(energy[j]*u.GeV,epsilon=1e-5)+0.))
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

#--------------Vettori di veritas-----------------------#
veritasy=[5.626247237386407e-14,2.4498548272157202e-14,7.083053855910684e-15,3.185536918894964e-15,1.432665297776149e-15,3.486116852371258e-16,1.9873216216790803e-16]
Veritasy=np.array(veritasy)
Veritasdy=np.array([9.644669209754113e-12,3.2923184187082306e-12,3.2923184187082306e-12,1.0101356224329272e-12,2.51984282827237e-13,1.1241647549139845e-13,6.404098398059653e-14])*10**-3#Non toccare, ho verificato col plot del paper e torna
Veritasx=np.array([0.21040591909975698,0.297501746151355,0.41945757279190776,0.5930886396919132,0.8385928812545467,1.18235989899921,1.671788208414575])*10**3
Veritasdx=[]
for j in range (1,len(Veritasx)):
    Veritasdx.append((Veritasx[j]-Veritasx[j-1])/2)
Veritasdx.append(np.max(Veritasx)/2)
Veritasdx=np.array(Veritasdx)


#--------------Vettori di Hess--------------------------#
hessy=[3.92e-15,1.47e-15,7.99e-16,2.48e-16,1.07e-16,5.88e-17,1.27e-17,2.17e-17,2.9e-18,2.4e-19,2.4e-18,4.8e-19]
Hessy=np.array(hessy)
Hessx=np.array([0.53,0.69,0.93,1.2,1.7,2.3,3.0,4.0,5.4,7.3,9.8,13])*10**3
Hessdy=np.array([1.95e-12,0.5e-12,0.2e-12,0.9e-13,0.43e-13,2.8e-14,1.245e-14,1.05e-14,0,0,0,0])*10**-3 #Non toccare, ho verificato col plot del paper e torna
Hessdx=[]
for j in range (1,len(Hessx)):
    Hessdx.append((Hessx[j]-Hessx[j-1])/2)
Hessdx.append(np.max(Hessdx)/2)
Hessdx=np.array(Hessdx)
Hessuplims=np.array([0,0,0,0,0,0,0,0,1,1,1,1])

#--------------Deassorbimento vettori e vettori totali------------------#
flussotot=veritasy+hessy
Flussotot=np.array(flussotot)
#flussotot=deassorbi(np.array(Veritasx.tolist + Hessx.tolist),(veritasy+hessy),dominguez)
enertot=sorted(Veritasx.tolist()+Hessx.tolist())
Enertot=np.array(enertot)
XTOT=np.array(sorted(Veritasx.tolist()+Hessx.tolist()))
DYTOT=np.array(Veritasdy.tolist()+Hessdy.tolist())#ACHTUNG: VA RISCRITTO
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

# #------------------Mappa del chi quadro
step=(5e-11,1e-1,1e3)
fval=(3e-10,8*5e-10)
gammaval=(1.,2.)
cutoff=[1e3,60e3]
#Mchi2=Mappachi(XCOMBINED,DYCOMBINED,fval[0],fval[1],gammaval[0],gammaval[1],cutoff[0],cutoff[1],step[0],step[1],step[2])
# print 'Coordinate del minimo:', np.where(Mchi2 == Mchi2.min())
# minchival=(np.where(Mchi2 == Mchi2.min())[0]*step[0]+fval[0],np.where(Mchi2 == Mchi2.min())[1]*step[1]+gammaval[0],np.where(Mchi2 == Mchi2.min())[2]*step[2]+cutoff[0])
# print 'fval=',minchival[0]
# print 'gamma=',minchival[1]
# print 'Cutoff=',minchival[2]

Mchi2,mappa=Mapchi2(XCOMBINED,DYCOMBINED,3.9920808393980515e-10,cutoff[0],cutoff[1],gammaval[0],gammaval[1],step[2],step[1])
Mchi2=np.transpose(np.array(Mchi2))

#------remember  fittabreak(E,dy,fmin,fmax,gammamin,gammamax,stepf,stepg)----#
b=fittabreak(XCOMBINED,DYCOMBINED,cutoff[0],cutoff[1],gammaval[0],gammaval[1],step[2],step[1])
#b=fittabreakold(XCOMBINED,DYCOMBINED,fval[0],fval[1],gammaval[0],gammaval[1],step[0],step[1])
print 'min chi squared',np.min(np.array(b)[:,3])
row=np.where(np.array(b)[:,3] == np.min(np.array(b)[:,3]))[0][0]
print 'min chi squared parameters E0, gamma, break', b[row]

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
    plottares(XCOMBINED,deassorbi(XCOMBINED,YCOMBINED,dominguez),DYCOMBINED*(deassorbi(XCOMBINED,YCOMBINED,dominguez)/YCOMBINED),Cutoff(a[row][0],XCOMBINED,a[row][1],a[row][2]))
    plt.show()

plt.xscale("log")
plt.yscale('log')

plottapuntide()
plottacutoff(b[row][0],b[row][1],b[row][2])
plt.xscale("log")
plt.yscale('log')
plt.show()




