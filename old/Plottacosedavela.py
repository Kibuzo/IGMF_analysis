from pylab import *
#savedir='/home/kibuzo/TesiMagistrale/Python/Outfigs/'
savedir='./Data/'
filedir='/home/kibuzo/davela/Output_for_Stefano/B1e-15G/'

def countinside(radius, a):
    return np.sum(np.sqrt(a[:,6]**2+a[:,7]**2)<radius)


#for i in range (3,9):
for i in range (100,199):
    filename=filedir+'HESSJ1943_spe_photindex1_6_Ecut_2TeV_B1e-6_n1e3ev_con5deg_obsize37Mpc_everything_sem_nrun_'+str(i)+'.txt'
    #Brms=float('1e-'+str(i)) 
    Brms=1e-06
    #filename=savedir+'photon_electron_output'+str(Brms)+'.txt'       #Usare per BRMS diversi
    #filename=savedir+'photon_electron_output'+str(Brms)+str(i)+'.txt' #Usare per BRMS uguali
    a=loadtxt(filename)
    #figure(figsize=(8,8))
    vecs=[]
    vettore=np.linspace(0,180,num=90)
    for j in vettore:
        vecs.append(countinside(j,a))
    plot(vecs,label='B='+str(Brms)+' ng')                     #Usare per BRMS diversi
    #plot(vecs,label='B='+str(Brms)+str(i)+' ng')              #Usare per BRMS uguali
    print 'angolo massimo di emissione= ' +str(np.round(np.arcsin(np.max(np.sqrt(a[:,19]**2+a[:,18]**2)/np.sqrt(a[:,19]**2+a[:,17]**2+a[:,18]**2)))*180/np.pi,decimals=2))+ ' gradi'

    
xlabel('Radius(Mpc)')
ylabel('Counts')
title('Number counts of photons as a fuction of radius')
legend(loc="upper left")
savefig(savedir+'countsrad_all.png')
close()

a=[]

for i in range (100,199):
    filename=filedir+'HESSJ1943_spe_photindex1_6_Ecut_2TeV_B1e-6_n1e3ev_con5deg_obsize37Mpc_everything_sem_nrun_'+str(i)+'.txt'
    Brms=float('1e-'+str(i)) 
    Brms=1e-06
    #filename='photon_electron_output'+str(Brms)+'.txt'        #Usare per BRMS diversi
    #filename=savedir+'photon_electron_output'+str(Brms)+str(i)+'.txt'  #Usare per BRMS uguali
    a=loadtxt(filename)
    idx = np.abs(a[:,3]) == 22 #Seleziona gli indici che contengono fotoni
    #figure(figsize=(8,8)
    plot(a[idx,7],a[idx,6],marker='.',markersize=1,linestyle='',label=str(Brms)+'ng')
    s=np.unique(a[idx,2])
    print len(a[idx,2])
    print len(s)
    xlabel('Radius(Mpc)')
    ylabel('Counts')
    title('Events scatter plot')
    legend(loc="lower right")
    #savefig(savedir+'PSF_'+str(Brms)+'ng'+'.png')              #Usare per BRMS diversi
    savefig(savedir+'PSF_'+str(Brms)+str(i)+'ng'+'.png')        #Usare per BRMS uguali
    close()
    # theta=np.arcsin(a[idx,8]/np.sqrt(a[idx,8]**2+a[idx,9]**2+a[idx,10]**2))
    # phi=np.arcsin(a[idx,9]/np.sqrt(a[idx,8]**2+a[idx,9]**2+a[idx,10]**2))    
    # or
    theta=np.arcsin(a[idx,7]/np.sqrt(a[idx,7]**2+a[idx,6]**2+a[idx,7]**2))
    phi=np.arcsin(a[idx,6]/np.sqrt(a[idx,7]**2+a[idx,6]**2+a[idx,7]**2))
    theta=theta*180/np.pi
    phi=phi*180/np.pi
    plot(theta,phi,marker='.',markersize=5,linestyle='',label=str(Brms)+'ng')
    xlabel('Theta')
    ylabel('Phi')
    title('Events scatter plot')
    legend(loc="lower right")   
    #savefig(savedir+'PSF_'+str(Brms)+'ng_thetaphi'+'.png')              #Usare per BRMS diversi
    savefig(savedir+'PSF_'+str(Brms)+str(i)+'ng_thetaphi'+'.png')        #Usare per BRMS uguali
    close()
    hist2d(a[idx,5],a[idx,6],bins=10)
    colorbar()
    xlabel('Radius(Mpc)')
    ylabel('Counts/bin')
    title('Binned events')
    #savefig(savedir+'PSF_'+str(Brms)+'ng'+'_binned.png')       #Usare per BRMS diversi
    savefig(savedir+'PSF_'+str(Brms)+str(i)+'ng'+'_binned.png') #Usare per BRMS uguali
    close()
    
for i in range (100,199):
    Brms=float('1e-'+str(i)) 
    Brms=1e-06
    #filename='photon_electron_output'+str(Brms)'.txt' 5          #Usare per BRMS diversi
   # filename=savedir+'photon_electron_output'+str(Brms)+str(i)+'.txt'   #Usare per BRMS uguali
    filename=filedir+'HESSJ1943_spe_photindex1_6_Ecut_2TeV_B1e-6_n1e3ev_con5deg_obsize37Mpc_everything_sem_nrun_'+str(i)+'.txt'
    a_Brms=loadtxt(filename)
    a=np.concatenate((a,a_Brms)) #appiccica insieme tutti i file
idx = np.abs(a[:,3]) == 22 #Seleziona gli indici che contengono fotoni
#figure(figsize=(8,8)
#plot(a[idx,5],a[idx,6],marker='.',markersize=1,linestyle='',label=str(Brms)+'ng')
hist2d(a[idx,7],a[idx,6],bins=100)
colorbar()
s=np.unique(a[idx,2])
print len(a[idx,2])
print len(s)
xlabel('Radius(Mpc)')
ylabel('Counts')
title('Events scatter plot')
#legend(loc="lower right")
savefig(savedir+'PSF_avg.png')
close()

# theta=np.arcsin(a[idx,8]/np.sqrt(a[idx,8]**2+a[idx,9]**2+a[idx,10]**2))
# phi=np.arcsin(a[idx,9]/np.sqrt(a[idx,8]**2+a[idx,9]**2+a[idx,10]**2))    
# or
d=D/Mpc
theta=np.arcsin((a[idx,7])/np.sqrt((d-a[idx,5])**2+(d-a[idx,6])**2+(d-a[idx,7])**2))
phi=np.arcsin((a[idx,6])/np.sqrt((d-a[idx,5])**2+(d-a[idx,6])**2+(d-a[idx,7])**2))
theta=theta*360/np.pi
phi=phi*360/np.pi
#plot(theta,phi,marker='.',markersize=5,linestyle='',label=str(Brms)+'ng')
#hist2d(theta,phi,bins=10,norm=mpl.colors.LogNorm())
hist2d(theta,phi,bins=100)
colorbar()
xlabel('Theta')
ylabel('Phi')
title('Binned events')
#legend(loc="lower right")   
#savefig(savedir+'PSF_'+str(Brms)+'ng_thetaphi'+'.png')              #Usare per BRMS diversi
savefig(savedir+'PSF_'+str(Brms)+str(i)+'ng_thetaphi_all'+'.png')        #Usare per BRMS uguali
close()



