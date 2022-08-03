#Piccolo script che dimostra che prendere il centro dell'intervallo energetico binnato nel caso di power law è pericolodo
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

N0=5
E0=15000
Gamma=-1.6

def powerlaw(E):
    N0=5
    E0=15000
    Gamma=-1.6
    return((N0*(E/E0)**Gamma))

#Funzione che integra lo spettro di power law tra due estremi su un'ascissa x
def conta (min,max):
    a=0
    #for j in x:
    #    if (j>min and j<=max):
    #        a+=(N0*(j/E0)**Gamma)
    return(integrate.quad(powerlaw,min,max))[0]    #flusso nel bin
    
    

#Genero punti ditribuiti a power law
x=np.logspace(1,5,num=100)
y=[]
for j in x:
    y.append(N0*(j/E0)**Gamma)
    

#Binno la distribuzione power law
biny=[]
for a in range (1,5):
    biny.append(conta(10**(a),10**(a+1)))
binxmid=[4.5*1e1,4.5*1e2,4.5*1e3,4.5*1e4] #x1+Dx/2
binxright=((np.array(biny)/(N0*np.array(binxmid)*2))**(1/Gamma))*E0
biny=biny/(np.array(binxmid)*2) #Porto da integrale a flusso

plt.loglog(x,y)
plt.loglog(binxmid,biny,marker='v',linestyle='') #Nel mezzo del bin
plt.loglog(binxright,biny,marker='^',linestyle='') #Nel punto giusto
plt.show()
