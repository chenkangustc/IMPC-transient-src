import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#input file
loopdat=np.loadtxt('..\\looptimelist.txt',skiprows=1)
maxdat=np.loadtxt('..\\maxT.timelist',skiprows=1)
data2=np.loadtxt('..\\aveT.timelist')
data3=np.loadtxt('..\\Tdis.txt')
#exp input
exdir='E:\\documents\\doctors degree\\software\\tansistant\\system\EBR-II\\SHRT-17ex data\\'
exname=['high pressure inlet.txt',\
        'low pressure inlet.txt',\
        'Z-pipe inlet temperature.txt',\
        'peak fuel temperature.txt',\
        'peak coolant temperature.txt',\
        'IHX primary inlet temperature.txt',\
        'pump2 mass flow.txt']
exanl=['Exp','Exp','ANL','ANL','ANL','Exp','Exp']
#generel control
is_flag=[True,False,False,True]
is_exp=True
#FIG1 sub0 control
zone=12
#FIG2 sub[1,5]control
tlist=[0,50,100,200,300]
#data
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
#FIG1
if is_flag[0]==True:
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
#FIG2
if is_flag[1]==True:
    fig,ax=plt.subplots(3,2)
    plt.subplots_adjust(hspace=0.40)
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
        #ax1[idx+1].set_ylim(480.0,1300.0)
#FIG3
#time flowrate powinput  Tsin Qin
if is_flag[2]==True:
    fig,ax=plt.subplots(2,2)
    plt.subplots_adjust(hspace=0.27)
    ax[0,0].plot(time,flowrate)
    ax[0,0].set_ylabel('flowrate/(Kg/s)')
    ax[0,0].set_xlabel('time/s')  
    ax[0,0].set_title('Primary flowrate')  
    ax[0,0].legend()
    
    ax[0,1].plot(time,powinput)
    ax[0,1].set_ylabel('Power/W')
    ax[0,1].set_xlabel('time/s')  
    ax[0,1].set_title('Power total')  
    ax[0,1].legend()
    
    ax[1,0].plot(time,Qs)
    ax[1,0].set_ylabel('flowrate/(kg/s)')
    ax[1,0].set_xlabel('time/s')  
    ax[1,0].set_title('IHX intermedia flowrate')  
    ax[1,0].legend()
    
    ax[1,1].plot(time,Tsin)
    ax[1,1].set_ylabel('Temperature/K')
    ax[1,1].set_xlabel('time/s')  
    ax[1,1].set_title('IHX intermedia Tin')  
    ax[1,1].legend()
#FIG4:EBR-II
if is_flag[3]==True:
    fig,ax=plt.subplots(2,4)
    plt.subplots_adjust(wspace=0.28,hspace=0.38)

    #High pressure plenum inlet temperature
    ax[0,0].plot(time,coreTin,label='IMPC')
    #ax[0,0].plot(exdata[:,0],exdata[:,1],label='Exp')
    ax[0,0].set_ylabel('Temperature/K')
    ax[0,0].set_xlabel('Time/sec')  
    ax[0,0].set_title('high pressure inlet temperature')  
    ax[0,0].legend()
    ax[0,0].set_ylim(610.0,635.0)
    #Low pressure plenum inlet temperature
    ax[0,1].plot(time,coreTin,label='IMPC')
    ax[0,1].set_ylabel('Temperature/K')
    ax[0,1].set_xlabel('Time/sec')  
    ax[0,1].set_title('low pressure inlet temperature')  
    ax[0,1].legend()
    ax[0,1].set_ylim(610.0,635.0)
    #Z-pipe inlet termperature
    ax[0,2].plot(time,coreTout,label='IMPC')
    ax[0,2].set_ylabel('Temperature/K')
    ax[0,2].set_xlabel('Time/sec')  
    ax[0,2].set_title('Z-pipe inlet termperature')  
    ax[0,2].legend()   
    # ax[0,2].set_ylim(650.0,800.0)    
    #Peak fuel temperature
    ax[0,3].plot(time,maxTfuel,label='IMPC')
    ax[0,3].set_ylabel('Temperature/K')
    ax[0,3].set_xlabel('Time/sec')  
    ax[0,3].set_title('peak fuel temperature')  
    ax[0,3].legend()
    # ax[0,3].set_ylim(650.0,950.0) 
    #Peak coolant temperature
    ax[1,0].plot(time,maxTcoolant,label='IMPC')
    ax[1,0].set_ylabel('Temperature/K')
    ax[1,0].set_xlabel('Time/sec')  
    ax[1,0].set_title('peak coolant temperature')  
    ax[1,0].legend()
    # ax[1,0].set_ylim(650.0,950.0) 
    #IHX primary inlet temperature
    ax[1,1].plot(time,IHXTin,label='IMPC')
    ax[1,1].set_ylabel('Temperature/K')
    ax[1,1].set_xlabel('Time/sec')  
    ax[1,1].set_title('IHX primary inlet temperature')  
    ax[1,1].legend()
    # ax[1,1].set_ylim(595.0,750.0) 
    #Pump2 mass flow rate
    ax[1,2].plot(time,flowrate/2.,label='IMPC')
    ax[1,2].set_ylabel('Flowrate/(kg/s)')
    ax[1,2].set_xlabel('Time/sec')  
    ax[1,2].set_title('pump#2 mass flow rate')   
    ax[1,2].legend()
    #exp plot
    if is_exp==True:
        ax1=ax.flatten()
        for idx,iname in enumerate(exname):
            exdatab=np.loadtxt(exdir+iname,skiprows=4)
            dtime=exdatab[:,0]
            dvalue=exdatab[:,1]
            ax1[idx].plot(dtime,dvalue,label=exanl[idx])
            ax1[idx].legend()
            ax1[idx].yaxis.grid(True)  
            ax1[idx].xaxis.grid(True)  
        # exdata=np.array(exdatab)
        # for idx,iname in enumerate(exname):
            # dtime=exdata[idx,:,0]
            # dvalue=exdata[idx,:,1]
            # ax1[idx].plot(dtime,dvalue,label='exp')
plt.savefig(exdir+'test.png', format='png', bbox_inches='tight', transparent=True, dpi=600)
plt.show()
