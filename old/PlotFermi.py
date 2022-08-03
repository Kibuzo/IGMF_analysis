import numpy as np
import matplotlib.pyplot as plt

def PlotPL():
    Prefactor=4.778403E-15
    Index=-1.793774
    Scale=15794.25
    DeltaPF=6.038888E-16
    DeltaIDX=-0.1147944
    #Passo a GeV
    Prefactor=Prefactor*1e3
    Scale=Scale/1e3
    DeltaPF=DeltaPF*1e3
    E=np.logspace(np.log10(np.min(Fermix-Loxerr)),np.log10(np.max(Fermix+Hixerr)),num=100)
    plt.plot(E,((Prefactor)*(E/(Scale))**Index), label='Best fit PL model',color='r')
    #Sigmaup=((((Prefactor+DeltaPF))*(E/(Scale))**(Index+DeltaIDX)))
    Sigmaup=((Prefactor)*(E/(Scale))**Index)*(1+np.sqrt((DeltaPF/Prefactor)**2+((DeltaIDX*np.log((E)/Scale)))**2))
    Sigmadown=((Prefactor)*(E/(Scale))**Index)*(1-np.sqrt((DeltaPF/Prefactor)**2+((DeltaIDX*np.log((E)/Scale)))**2))
    #Sigmadown=((((Prefactor-DeltaPF))*(E/(Scale))**(Index-DeltaIDX)))
    plt.fill_between(E,np.maximum(Sigmaup,Sigmadown),np.minimum(Sigmaup,Sigmadown),color='r', alpha=0.2, label='1$\sigma$ CL (PL)')
    
def PlotLP():
    Norm=5.269057E-15*1e3
    Alpha=1.612989
    Beta=0.1323156
    Eb=15794.25/1e3
    DeltaNorm=7.674917E-16*1e3
    DeltaAlpha=0.176109
    DeltaBeta=0.09997965
    Corr=-0.01192512
    E=np.logspace(np.log10(np.min(Fermix-Loxerr)),np.log10(np.max(Fermix+Hixerr)),num=100)
    F=Norm*(E/Eb)**(-(Alpha+Beta*np.log(E/Eb)))
    plt.plot(E,F, color='b',label='Best fit, LP model')
    Factora=-Norm*((np.log(E/Eb)))*(E/Eb)**(-Alpha-Beta*np.log(E/Eb))
    Factorb=-Norm*((np.log(E/Eb))**2)*(E/Eb)**(-Alpha-Beta*np.log(E/Eb))
    Sigma=np.sqrt((DeltaNorm*F/Norm)**2+(DeltaAlpha*Factora)**2+(DeltaBeta*Factorb)**2+2*Corr*Factora*Factorb)
    plt.fill_between(E,np.maximum(F+Sigma,F-Sigma),np.minimum(F+Sigma,F-Sigma),alpha=0.2, color='b', label='1$\sigma$ CL (LP)')
    #plt.errorbar(Fermix,nphotLP,yerr=FermidxLP,xerr=[Loxerr,Hixerr],linestyle='',marker='.',color='k')
    
#$FGL psc_v18 Extended sources P8R3IRF
#nphot=[2.03e-11,1.36e-11,9.09e-12,6.09e-12,4.077e-12,2.73e-12,1.82e-12]
nphot=[2.07e-11,5.69e-12,2.11e-12,7.3e-13,2.688e-13,2.345e-14,9.9e-15]
nphotLP=[2.05e-11,5.68e-12,2.06e-12,7.14e-13,2.66e-13,2.35e-14,9.47e-15]
FermidxLP=[4.5e-12,1.51e-12,5.79e-13,2.73e-13,9.83e-14,2.24e-14,1e-14]
Fermix=np.array([6.94,13.41,25.9,50,96.54,186.38,359.84]) #in GeV
Fermixmin=np.array([5,9.65,18.6,35.98,69.47,134.13,258.97])
Fermixmax=np.array([9.65,18.6,35.98,69.47,134.13,258.97,500])
#dNphot=np.array([5.075e-12,3.62e-12,2.54e-12,1.93e-12,1.51e-12,2.63e-12,1.86e-12])
dNphot=np.array([4.5e-12,1.52e-12,5.9e-13,2.31e-13,9.95e-14,2.26e-14,1.01e-14])#/2
Fermidx=(Fermixmax-Fermixmin).tolist()
fermiuplims=np.array([0,0,0,0,0,1,1])
Fermits=[39.8573225013606,36.8792876649022,46.0275764309645,46.3142510222815,37.0156520052533,4.5468426120151,5.7454650115202]
#nphot=(nphot/np.array(Fermidx)).tolist()
#dNphot=dNphot/np.array(Fermidx)
#Provo a correggere l'errore dei punti con ts bassa con la seguente procedura: assumo che la TS dia il p-value, dal p value ricavo il chi quadro corrispondente che nel mio caso è circa 0.5 per entrambi poi uso la formula per il chi quadro con n gradi di libertà chi**2=(n-1)s**2sigma**-2 dove s è la varianza del sample, quindi l'incertezza che mi dà su nphot e sigma è la varianza vera. Poi sommo tutto in quadratura pregando il signore; non so se va bene e comunque non cambia abbastanza i risultati però è un inizio.
dNphot[5]=np.sqrt(dNphot[5]**2+(2.26e-14)**2)
dNphot[6]=np.sqrt(dNphot[6]**2+(1.01e-14)**2)
Loxerr=[]
Hixerr=[]
for j in range (0,6):
    if (j==0):
        Loxerr.append(Fermix[0]-5)
    else:
        Loxerr.append((Fermix[j]-Fermix[j-1])/2)
Loxerr.append(Loxerr[-1])
for j in range (0,6):
    Hixerr.append((Fermix[j+1]-Fermix[j])/2)
Hixerr.append(Hixerr[-1])
PlotPL()
PlotLP()
#plt.errorbar(Fermix,nphot,yerr=dNphot,xerr=[Loxerr,Hixerr],linestyle='',marker='.',color='k')
plt.xlabel('Energy[GeV]')
plt.ylabel('dN/dE [cm$^{-2}$s$^{-1}$GeV$^{-1}$]')
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.show()