#ToDo: Capire perché se il lower limit dell'estensione della power law si abbassa allora il montecarlo rappresenta di meno la PDF indipendentemente dal numero di sample. Questo è un problema puramente filosofico perché non sono interessato a E<100GeV e non esiste che qualcuno si metta a controllare questo montecarlo. Però magari vuoi essere honestissimo
import numpy as np
import matplotlib.pyplot as plt
import os
Cut=2080 #In GeV
#Cut=20000
Nzero=3.83e-14
nsample=1e5                                             #Montecarlo samples
extent=[200,2e4]                                        #estremi in energia della PDF
gridpoints=int(extent[1]-extent[0])*10                  #numero di punti della griglia per la discretizzazione della CDF
estep=(extent[1]-extent[0])/gridpoints                  #passo della griglia
egrid=np.linspace(extent[0],extent[1],num=gridpoints)   #definizione della griglia
filesplit=100
savedir='./MCsample/'
if not os.path.exists(savedir):
    os.mkdir(savedir)

def salvadati(energie,inizio,fine,dove):
    inizio=int(inizio)
    fine=int(fine)
    np.savetxt(dove,energie[inizio:fine])

def GridToE(gridpoint):
    return (extent[0]+gridpoint*estep)
    
def EtoGrid(Energy):
    return ((Energy-extent[0])/estep)

def Cutoff(E,N0,gamma,Ecut):
    return ((N0*(E/300.)**-gamma)*np.exp(-E/Ecut))

#Segmenta con rettangoli la funzione e  restituisce la CDF grigliata grid
def rettangoli(N0,Gamma,Ecut):
    CDF=[0]
    PDF=[]
    for x in egrid:
        areasup=estep*Cutoff(x,Nzero,Gamma,Ecut)
        areainf=estep*Cutoff(x+estep,Nzero,Gamma,Ecut)
        CDF.append((areasup+areainf)/2+CDF[-1])
        PDF.append((areasup+areainf)/2)
    return(CDF,PDF)

#Gli do uno xi generato a caso tra emin ed emax e lui mi ricava l'energia corrispondente nella CDF
def eval(xi,CDF):
    j=np.min(np.where(CDF>=xi)) #Funziona ed è veloce, che tu ci creda o no
    return (GridToE(j))
    
# ---------------MC starts ------------------ #
CDF,PDF=rettangoli(Nzero,1.5,Cut) #CDF, PDF
CDF=CDF/max(CDF)                  #Normalized CDF
f=[]
for j in range (0,int(nsample)):
    xi=np.random.random()        
    f.append(eval(xi,CDF))


# --------------- MC plot starts ------------------#
nbins=50
plt.close()
logbins=np.logspace(np.log10(extent[0]),np.log10(extent[1]),num=50)
data,bins=np.histogram(f,bins=nbins) #I bin adesso sono in energia perché eval returna energie
bincenters = 0.5*(bins[1:]+bins[:-1])
e=np.linspace(min(extent),max(extent),num=nbins)
func=Cutoff(e,Nzero,1.5,Cut)
binspace=np.logspace(np.log(extent[0]),np.log(extent[1]),num=nbins)
logbincenters=0.5*(binspace[1:]+binspace[:-1])

stepping=(extent[1]-extent[0])/nbins
#binpos=[100+(1./1.5)*(np.log(1.5*stepping)-np.log(1-np.exp(-1.5*stepping)))]
delta0=bins[1]-bins[0]
gamma=1.5
binpos=[]
#binpos=[(bins[0]+(1./gamma)-(bins[1]+1./gamma)*np.exp(-gamma*delta0))/(1-np.exp(-gamma*(delta0)))]
for j in range (0,nbins):
    x1=bins[j]
    x2=bins[j+1]
    deltax=x2-x1
    gamma=1.5
    numeratore=( x1 + (1./gamma) - (x2+1./gamma)*np.exp(-gamma*deltax) )
    denominatore=1-np.exp(-gamma*deltax)
    binpos.append(numeratore/denominatore)
binpos=np.array(binpos)
plt.figure(figsize=(10,5))
#Confronto plot di PDF valutata al centro del bin e istogramma montecarlo

# #New way
# #func=Cutoff(binpos,Nzero,1.5,Cut)
# plt.loglog(e,func,label='PDF')
# funcnorm=Cutoff(binpos,Nzero,1.5,Cut)
# plt.loglog(binpos,data*(np.sum(funcnorm)/nsample),label='Monte Carlo set',linestyle='',marker='.')

# #Old style
#plt.loglog(e,func,label='PDF')
funcnorm=Cutoff(bincenters,Nzero,1.5,Cut)
#plt.plot(bincenters,data/(np.sum(data)/np.sum(funcnorm)),label='Monte Carlo set',linestyle='',marker='.')
plt.bar(bincenters,width=bincenters[1]-bincenters[0],height=data/(np.sum(data)/np.sum(funcnorm)),label='Monte Carlo set',alpha=0.5)
plt.plot(bincenters,funcnorm,label='PDF', color='orange')

# stepping=(extent[1]-extent[0])/nbins
# binpos=[100+(1./1.5)*(np.log(1.5*stepping)-np.log(1-np.exp(-1.5*stepping)))]
# for j in range (0,49):
#     binpos.append(binpos[-1]+100+(1./1.5)*(np.log(1.5*stepping)-np.log(1-np.exp(-1.5*stepping))))
# plt.plot(binpos,func*(np.sum(data)/np.sum(func)),label='PDF')
# plt.hist(f,bins=50)
# plt.xscale('log')
# plt.yscale('log')
plt.legend(loc='lower left')
plt.title('Monte carlo rescaled sample (n=10$^5$) vs analytical PDF (Ecut='+str(Cut/1e3)+'TeV)')
plt.ylabel('Flux[GeV$^{-1}$s$^{-1}$cm$^{-2}$]')
plt.xlabel('Energy[GeV]')
plt.savefig(savedir+'MC_sample_complete.png')
plt.xscale('log')
plt.yscale('log')
plt.show()

# ------------- Plot diagnostico: hai fatto bene sampling della PDF? ----- #
# ics=np.linspace(extent[0],extent[1],num=grid-1)
# plt.loglog(ics,PDF)
# plt.loglog(ics,Cutoff(GridToE(ics),Nzero,1.5,Cut))
# plt.loglog(e,Cutoff(e,Nzero,1.5,Cut))

# plt.show()
for j in range (0,filesplit-1):
    salvadati(f,j*nsample/filesplit,(j+1)*nsample/filesplit,savedir+'mc'+str(j)+'.txt')
np.savetxt(savedir+'mc_complete.txt',f)

plt.figure(figsize=(10,5))
plt.plot(egrid,CDF[1:], label='Normalized CDF')
plt.plot(egrid,PDF/np.max(PDF), label='Rescaled PDF')
plt.title('Distribution functions')
plt.xlabel('Energy [GeV]')
plt.ylabel('normalized CDF/ Rescaled PDF')
#plt.xlim([egrid[0],egrid[-1]])
plt.xscale('log')
plt.legend()
#plt.yscale('log')
plt.show()
#ToDo: testare il montecarlo. Un modo è mettere delle barre di errore ai bin che dipendano esclusivamente dal numero di oggetti generati nei bin e successivamente calcolare il chi quadro tra distribuzione teorica e campione montecarlo. Do it over and over.


    

    