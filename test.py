import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#input file
loopdat=np.loadtxt('..\\looptimelist.txt',skiprows=1)
maxdat=np.loadtxt('..\\maxT.timelist',skiprows=1)
data2=np.loadtxt('..\\aveT.timelist')
data3=np.loadtxt('..\\Tdis.txt')
#FIG1 sub0 control
zone=12
#FIG2 sub[1,5]control
tlist=[0,50,11,12,13]

time=loopdat[:,0]
flowrate=loopdat[:,2]
coreTin=loopdat[:,3]
coreTout=loopdat[:,4]
IHXTin=loopdat[:,5]
IHXTout=loopdat[:,6]
Qs=loopdat[:,7]
Tsin=loopdat[:,8]
Tsout=loopdat[:,9]
shc=151.0
powinput=loopdat[:,1]
pripow=shc*flowrate*(IHXTin-IHXTout)
corepow=shc*flowrate*(coreTout-coreTin)
shcs=4660.0
secpow=shcs*Qs*(Tsout-Tsin)

times=maxdat[:,0]
maxTfuel=maxdat[:,1]
maxTcoolant=maxdat[:,2]
maxTinner=maxdat[:,3]
maxTouter=maxdat[:,4]

fig,ax=plt.subplots(3,2)
plt.subplots_adjust(hspace=0.38)

ax[0,0].plot(time,flowrate,label='flowrate')
ax[0,0].set_ylabel('flowrate kg/s')
ax[0,0].set_xlabel('time/s')
ax[0,0].set_title( 'flowrate' )
ax[0,0].legend()

ax[0,1].plot(time,powinput,label='pow')
ax[0,1].set_ylabel('Pow/W')
ax[0,1].set_xlabel('time/s')
ax[0,1].set_title( 'Powerinput' )
ax[0,1].legend()

ax[1,0].plot(time,coreTin,label='coreTin')
ax[1,0].plot(time,coreTout,label='coreTout')
ax[1,0].plot(time,IHXTin,linestyle='-.',label='IHXTin')
ax[1,0].plot(time,IHXTout,linestyle='-.',label='IHXTout')
ax[1,0].set_ylabel('Temperature/K')
ax[1,0].set_xlabel('time/s')
ax[1,0].set_title( 'Temperature' )
ax[1,0].legend()

ax[1,1].plot(time,secpow,label='secpow(IHXscmdt)')
ax[1,1].plot(time,pripow,label='pripow(IHXpcmdt)')
ax[1,1].plot(time,corepow,label='corepow(Corecmdt)')
ax[1,1].set_ylabel('Pow/W')
ax[1,1].set_xlabel('time/s')
ax[1,1].set_title('Pow removed by IHX')
ax[1,1].legend()

#ax[2,0].plot(time,IHXTin,label='priInlet')
#ax[2,0].plot(time,IHXTout,label='priOutlet')
ax[2,0].plot(time,Tsin,label='secInlet')
ax[2,0].plot(time,Tsout,label='secOutlet')
ax[2,0].set_ylabel('Temperature/K')
ax[2,0].set_xlabel('time/s')
ax[2,0].set_title('IHX Temperature')
ax[2,0].legend()

ax[2,1].plot(times,maxTfuel,label='maxTfuel')
ax[2,1].plot(times,maxTcoolant,label='maxTcoolant')
ax[2,1].plot(times,maxTinner,label='maxTinner')
ax[2,1].plot(times,maxTouter,label='maxTouter')
ax[2,1].set_ylabel('Temperature/K')
ax[2,1].set_xlabel('time/s')
ax[2,1].set_title('max Temperature')
ax[2,1].legend()

time=data2[:,0]
Tfave=data2[:,4*zone-3]
Tcave=data2[:,4*zone-2]
Tinlet=data2[:,4*zone-1]
Toutlet=data2[:,4*zone]
zz=data3[0,1:]

fig,ax=plt.subplots(3,2)
ax[0,0].plot(time,Tfave,label='fuelTave')
ax[0,0].plot(time,Tcave,label='coolantTave')
ax[0,0].plot(time,Tinlet,label='Tinlet')
ax[0,0].plot(time,Toutlet,label='Toutlet')
ax[0,0].set_ylabel('Temperature/K')
ax[0,0].set_xlabel('time/s')
ax[0,0].set_title('zone='+str(zone)+' Temperature dynamic curve')
ax[0,0].legend()

ax1=ax.flatten()
for idx,it in enumerate(tlist):
    Tfuel=data3[3*it+1,1:]
    Tcoolant=data3[3*it+2,1:]
    ax1[idx+1].plot(zz,Tfuel,label='Tfuel')
    ax1[idx+1].plot(zz,Tcoolant,label='Tcoolant')
    ax1[idx+1].set_ylabel('Temperature/K')
    ax1[idx+1].set_xlabel('zz/m')
    ax1[idx+1].set_title('zone='+str(zone)+',t='+str(it)+' Temperature')
    ax1[idx+1].legend()

plt.show()