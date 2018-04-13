import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#dir='.\\test1.25MW1000s30mpipe\\'
#dir='.\\0403\\'
#dir='E:\\documents\\doctors degree\\software\\tansistant\\system\\test\\vvtest1.0\\'
dir='..\\'
file='looptimelist.txt'

loopdat=np.loadtxt(dir+file,skiprows=1)

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

fig,ax=plt.subplots(2,2)
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

ax[1,0].plot(time,coreTin,linestyle='--',label='coreTin')
ax[1,0].plot(time,coreTout,linestyle='--',label='coreTout')
ax[1,0].plot(time,IHXTin,linestyle='-.',label='IHXTin')
ax[1,0].plot(time,IHXTout,linestyle='-.',label='IHXTout')
ax[1,0].set_ylabel('Temperature/K')
ax[1,0].set_xlabel('time/s')
ax[1,0].set_title( 'Temperature' )
ax[1,0].legend()

ax[1,1].plot(time,pripow,label='pripow(IHXpcmdt)')
ax[1,1].plot(time,corepow,label='corepow(Corecmdt)')
ax[1,1].plot(time,secpow,label='secpow(IHXscmdt)')
ax[1,1].set_ylabel('Pow/W')
ax[1,1].set_xlabel('time/s')
ax[1,1].set_title('Pow removed by IHX')
ax[1,1].legend()


'''
file2='output.timelist'
redat=np.loadtxt(dir+file2,skiprows=10)
power=redat[:,4]
Tcmax=redat[:,11]
Tfmax=redat[:,14]
Tmaxinner=redat[:,16]
Tmaxouter=redat[:,17]

ax[1,0].plot(time,power,label='power')
ax[1,0].set_ylabel('power')
ax[1,0].set_xlabel('time/s')
ax[1,0].set_title( 'Power' )
ax[1,0].legend()

ax[1,1].plot(time,Tcmax,label='Tcmax')
ax[1,1].plot(time,Tfmax,label='Tfmax')
ax[1,1].plot(time,Tmaxinner,label='Tmaxinner')
ax[1,1].plot(time,Tmaxouter,label='Tmaxouter')
ax[1,1].set_ylabel('Temperature/K')
ax[1,1].set_xlabel('time/s')
ax[1,1].set_title( 'Temperature' )
ax[1,1].legend()
'''
plt.show()