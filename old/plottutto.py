from pylab import *
savedir='/home/kibuzo/TesiMagistrale/Python/Outfigs/'

def countinside(radius, a):
    return np.sum(np.sqrt(a[:,6]**2+a[:,5]**2)<radius)
    
for i in range (3,9):
    Brms=float('1e-'+str(i)) 
    #filename='photon_electron_output'+str(Brms)+'.txt'
    filename='photon_electron_output'+str(Brms)+str(i)+'.txt'
    a=loadtxt(filename)
    #figure(figsize=(8,8))
    vecs=[]
    vettore=np.linspace(0,60,num=60)
    for j in vettore:
        vecs.append(countinside(j,a))
    plot(vecs,label='B='+str(Brms)+' ng')
    plot(vecs,label='B='+str(Brms)+str(i)+' ng')
    print 'angolo massimo di emissione= ' +str(np.round(np.arcsin(np.max(np.sqrt(a[:,17]**2+a[:,18]**2)/np.sqrt(a[:,19]**2+a[:,17]**2+a[:,18]**2)))*180/np.pi,decimals=2))+ ' gradi'

    
xlabel('Radius(Mpc)')
ylabel('Counts')
title('Number counts of photons as a fuction of radius')
legend(loc="upper left")
savefig(savedir+'countsrad_all.png')
close()

a=[]

for i in range (3,9):
    Brms=float('1e-'+str(i)) 
    #filename='photon_electron_output'+str(Brms)+'.txt'
    filename='photon_electron_output'+str(Brms)+str(i)+'.txt'
    a=loadtxt(filename)
    idx = np.abs(a[:,3]) == 22 #Seleziona gli indici che contengono fotoni
    #figure(figsize=(8,8)
    plot(a[idx,5],a[idx,6],marker='.',markersize=1,linestyle='',label=str(Brms)+'ng')
    s=np.unique(a[idx,2])
    print len(a[idx,2])
    print len(s)
    xlabel('Radius(Mpc)')
    ylabel('Counts')
    title('Events scatter plot')
    legend(loc="lower right")
    #savefig(savedir+'PSF_'+str(Brms)+'ng'+'.png')
    savefig(savedir+'PSF_'+str(Brms)+str(i)+'ng'+'.png')
    close()
    hist2d(a[idx,5],a[idx,6],bins=10)
    xlabel('Radius(Mpc)')
    ylabel('Counts/bin')
    title('Binned events')
    #savefig(savedir+'PSF_'+str(Brms)+'ng'+'_binned.png')
    savefig(savedir+'PSF_'+str(Brms)+str(i)+'ng'+'_binned.png')
    close()
    
for i in range (3,9):
    Brms=float('1e-'+str(i)) 
    #filename='photon_electron_output'+str(Brms)'.txt'
    filename='photon_electron_output'+str(Brms)+str(i)+'.txt'
    a_Brms=loadtxt(filename)
    a=np.concatenate((a,a_Brms)) #appiccica insieme tutti i file
idx = np.abs(a[:,3]) == 22 #Seleziona gli indici che contengono fotoni
#figure(figsize=(8,8)
#plot(a[idx,5],a[idx,6],marker='.',markersize=1,linestyle='',label=str(Brms)+'ng')
hist2d(a[idx,5],a[idx,6],bins=10)
s=np.unique(a[idx,2])
print len(a[idx,2])
print len(s)
xlabel('Radius(Mpc)')
ylabel('Counts')
title('Events scatter plot')
#legend(loc="lower right")
savefig(savedir+'PSF_avg.png')
close()



