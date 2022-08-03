import numpy as np
import matplotlib.pyplot as plt

def linfit (x,m,q):
    return m*x+q

def delta(dati):
    return (1/np.sqrt(dati))

Bergamo=np.array([372,423,537,623,761,997,1245,1472,1815,2136])
Lodi=np.array([482,559,658,739,811,853,928,963,1035,1123])
Milano=np.array([93,145,197,267,361,406,506,592,925,1146])
Piacenza=np.array([256,319,387,426,479,528,602,633,664,679])
Padova=np.array([144,162,175,198,216,255,273,296,373,439])
Toscana=np.array([19,38,61,79,113,166,208,264,320,364])
Toscana2=np.array([2,2,7,10,12,12,18,37,60,78,112,165,206,260,314,352,455,614,763,841,1024,1291,1482])
Lombardia=np.array([1520,1850,2251,2612,3420,4149,5469,5791,7280,8725])
Veneto=np.array([307,360,407,488,543,670,744,856,1023,1384])
giorni=np.linspace(0,len(Toscana2)-1, num=len(Toscana2))
piugiorni=np.linspace(0,30,num=31)

#poptuscany,pcovtuscany=curve_fit(linfit, giorni , np.log(Toscana))
#plt.plot(np.e**linfit(piugiorni,poptuscany[0],poptuscany[1]), label='fit Toscana, indice='+str(poptuscany[0]))
poptuscany2,pcovtuscany2=curve_fit(linfit, giorni , np.log(Toscana2))
plt.plot(np.e**linfit(piugiorni,poptuscany2[0],poptuscany2[1]), label='fit Toscana, indice='+str(poptuscany2[0]))
#poptpadova,pcovpadova=curve_fit(linfit, giorni , np.log(Padova))
#plt.plot(e**linfit(piugiorni,poptpadova[0],poptpadova[1]), label='fit Padova, indice='+str(poptpadova[0]))
#poptlodi,pcovlodi=curve_fit(linfit, giorni , np.log(Lodi))
#plt.plot(e**linfit(piugiorni,poptlodi[0],poptlodi[1]), label='fit Lodi, indice='+str(poptlodi[0]))
#poptlombardia,pcovlombardia=curve_fit(linfit, giorni , np.log(Lombardia))
#plt.plot(e**linfit(piugiorni,poptlombardia[0],poptlombardia[1]), label='fit Lombardia, indice='+str(poptlombardia[0]))

#plt.semilogy(Bergamo, label='Bergamo', marker='.', linestyle='')
#plt.plot(Lodi, label='Lodi', marker='.', linestyle='')
#plt.semilogy(Milano, Label='Milano', marker='.', linestyle='')
#plt.semilogy(Piacenza, label='Piacenza', marker='.', linestyle='')
#plt.plot(Padova, label='Padova', marker='.', linestyle='')
plt.plot(Toscana2, label='Toscana', marker='.', linestyle='')
#plt.plot(Lombardia, label='Lombardia', marker='.', linestyle='')
#plt.semilogy(Veneto, label='Veneto', marker='.', linestyle='')
plt.legend()
plt.xlabel('Giorni dal 25 febbraio')
plt.ylabel ('Casi diagnosticati')
plt.yscale('log')
plt.show()
