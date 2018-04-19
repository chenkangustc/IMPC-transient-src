module imp_single_channel
    use constants
	use imp_property
    use imp_assm_global
    use imp_single_kernel
    use imp_assembly_header
    implicit none
     private
     public::driving_imp_steady
     public::driving_imp_Transient
	 public::driving_imp_THsteady   !only energy formula is calculated
	 public::driving_imp_THtransient!no pv step ,not simple
	 public::driving_imp_flowAlloc
	 public::driving_loop_flowAlloc
    contains
	subroutine driving_imp_flowAlloc(assm,flowrate)
		type(sys_assembly),allocatable::assm(:)
		real(KREAL),allocatable::flowrate(:)
		!local
		real(KREAL)::fuelArea,density!flow area of pin
		integer zone,n_pin
		integer i
		zone=size(assm)
		do  i=1,zone,1
			n_pin=assm(i)%geom%n_pin
			fuelArea=assm(i)%hydrau%aflow
			density=get_density(assm(i)%th_boundary%T%inlet)
			assm(i)%th_boundary%u%inlet=flowrate(i)/(n_pin*fuelArea*density)
		enddo
	end subroutine driving_imp_flowAlloc
	
	subroutine driving_loop_flowAlloc(assm,Qin)
		type(sys_assembly),intent(in out)::assm(:)
		real(KREAL),intent(in)::Qin
		!local
		real(KREAL),allocatable::flowrate(:)
		real(KREAL)::fuelArea,density!flow area of pin
		integer zone,n_pin,layer,nr
		integer i,j
		zone=size(assm)
		layer=size(assm(1)%thermal%velocity)
		nr=size(assm(1)%thermal%temperature,2)
		allocate(flowrate(zone))
		
		do i=1,zone,1
			flowrate(i)=Qin/zone!average
		enddo
		
		do  i=1,zone,1			
				n_pin=assm(i)%geom%n_pin
				fuelArea=assm(i)%hydrau%aflow
				density=get_density(assm(i)%th_boundary%T%inlet)
				assm(i)%th_boundary%u%inlet=flowrate(i)/(n_pin*fuelArea*density)
				density=get_density(assm(i)%th_boundary%T%outlet)
				assm(i)%th_boundary%u%outlet=flowrate(i)/(n_pin*fuelArea*density)			
			do j=1,layer,1
				density=get_density(assm(i)%thermal%temperature(j,nr))
				assm(i)%thermal%velocity(j)=flowrate(i)/(n_pin*fuelArea*density)
			enddo
		enddo
		!velocity
		
	end subroutine driving_loop_flowAlloc
	
subroutine driving_imp_steady(assm,power,fq_core)
    type(sys_assembly),intent(in out)::assm
    real(KREAL),intent(in)::power(:,:)
    real(KREAL),intent(in)::fq_core(:,:)
    !local rhoi/ui/Ti/dt/ap/pmodify/
    integer M,N,Ny,i,j
    real(KREAL):: flag,dt
    real(KREAL):: btotal,drho!判断因子
    real(KREAL),allocatable::pmodify(:)
    real(KREAL),allocatable::ui(:),Ti(:,:)
    real(KREAL),allocatable::rhoi(:,:),rhofi(:)
    real(KREAL),allocatable::ap(:)
    !write(*,*)'start steady calculation:'
    flag=0.0
    dt=1.0!为保证方程求解dt不为0，无具体意义+
    Ny=assm%mesh%Ny
    M=Ny+1
    N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
    allocate(pmodify(Ny),ui(Ny-1),Ti(M-1,N),rhoi(0:M,0:N),rhofi(0:M),ap(Ny-1))
 
    pmodify=0.0
    ui=assm%thermal%velocity
    Ti=assm%thermal%temperature
    rhoi=assm%property%rho
    rhofi=assm%property%rho(:,N)
    ap=0.0
    j=0
    drho=1.0
    call assm%pow%set(power,fq_core)
    do while(drho>assm%confactor_%alpha)  
       i=0
       btotal=1.0
      do while(btotal>assm%confactor_%sigma)
	   i=i+1
       call solve_momentum(assm,flag,rhofi,ui,dt,ap)!需要输入当前迭代步的rho,uz 
       call solve_pressureCorrection(assm,flag,ap,rhofi,dt,pmodify,btotal)
       call modify_PV(assm,ap,pmodify)
       !print*,'pv step=',i,' btotal=',btotal
      end do
      !print*,'velocity=',assm%th_boundary%u%inlet,assm%thermal%velocity,assm%th_boundary%u%outlet
      !print*,'pressure=',assm%th_boundary%p%inlet,assm%thermal%pressure,assm%th_boundary%p%outlet

      j=j+1
      call solve_temperature_rhoi(assm,flag,Ti,rhoi,dt)
      call update_property(assm,drho)!物性更新
      !print*,'property step=',j,' drho=',drho
    end do
  end subroutine driving_imp_steady
     
subroutine driving_imp_Transient(assm,power,fq_core,ltime,ctime)
    type(sys_assembly),intent(in out)::assm
    real(KREAL),intent(in)::power(:,:)
    real(KREAL),intent(in)::fq_core(:,:)
    real(KREAL),intent(in)::ltime
    real(KREAL),intent(in)::ctime
    !local rhoi/ui/Ti/dt/ap/pmodify/
    integer M,N,Ny,i,j
    real(KREAL):: dt
    real(KREAL):: flag
    real(KREAL):: btotal,drho!判断因子
    real(KREAL),allocatable::pmodify(:)
    real(KREAL),allocatable::ui(:),Ti(:,:)
    real(KREAL),allocatable::rhoi(:,:),rhofi(:)
    real(KREAL),allocatable::ap(:)
    !write(*,*)'start transient calculation:'
    flag=1.0
    dt=ctime-ltime
    Ny=assm%mesh%Ny
    M=Ny+1
    N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
    allocate(pmodify(Ny),ui(Ny-1),Ti(M-1,N),rhoi(0:M,0:N),rhofi(0:M),ap(Ny-1))
    
    pmodify=0.0
    Ti=assm%thermal%temperature
    ui=assm%thermal%velocity
    rhoi=assm%property%rho
    rhofi=assm%property%rho(:,N)
    ap=0.0
    !call assm%th_boundary%update(ctime)
    call assm%pow%set(power,fq_core)
    j=0
    drho=1.0
    do while(drho>assm%confactor_%alpha)
       i=0
       btotal=1.0
      do while(btotal>assm%confactor_%sigma)
       i=i+1
       call solve_momentum(assm,flag,rhofi,ui,dt,ap)!需要输入当前迭代步的rho,uz 
       call solve_pressureCorrection(assm,flag,ap,rhofi,dt,pmodify,btotal)
       call modify_PV(assm,ap,pmodify)
       !print*,'pv step=',i,' btotal=',btotal
      end do
      j=j+1     
      call solve_temperature(assm,flag,Ti,rhoi,dt)
      call update_property(assm,drho)!物性更新
      !print*,'density step=',j,' drho=',drho
    end do
  end subroutine driving_imp_Transient
subroutine driving_imp_THsteady(assm,power,fq_core)
    type(sys_assembly),intent(in out)::assm
    real(KREAL),intent(in)::power(:,:)
    real(KREAL),intent(in)::fq_core(:,:)
    !local rhoi/ui/Ti/dt/ap/pmodify/
    integer M,N,Ny,i,j
    real(KREAL):: flag,dt
    real(KREAL):: btotal,drho!判断因子
    real(KREAL),allocatable::pmodify(:)
    real(KREAL),allocatable::ui(:),Ti(:,:)
    real(KREAL),allocatable::rhoi(:,:),rhofi(:)
    real(KREAL),allocatable::ap(:)
    !write(*,*)'start steady calculation:'
    flag=0.0
    dt=1.0!为保证方程求解dt不为0，无具体意义+
    Ny=assm%mesh%Ny
    M=Ny+1
    N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
    allocate(pmodify(Ny),ui(Ny-1),Ti(M-1,N),rhoi(0:M,0:N),rhofi(0:M),ap(Ny-1))
 
    pmodify=0.0
    ui=assm%thermal%velocity
    Ti=assm%thermal%temperature
    rhoi=assm%property%rho
    rhofi=assm%property%rho(:,N)
    ap=0.0
    j=0
    drho=1.0
    call assm%pow%set(power,fq_core)
    call solve_temperature_rhoi(assm,flag,Ti,rhoi,dt)
	!cal Tfave & Tcave
	call cal_Tave(assm)
end subroutine driving_imp_THsteady
subroutine driving_imp_THtransient(assm,power,fq_core,ltime,ctime)
    type(sys_assembly),intent(in out)::assm
    real(KREAL),intent(in)::power(:,:)
    real(KREAL),intent(in)::fq_core(:,:)
    real(KREAL),intent(in)::ltime
    real(KREAL),intent(in)::ctime
    !local rhoi/ui/Ti/dt/ap/pmodify/
    integer M,N,Ny,i,j
    real(KREAL):: dt
    real(KREAL):: flag
    real(KREAL):: btotal,drho!判断因子
    real(KREAL),allocatable::pmodify(:)
    real(KREAL),allocatable::ui(:),Ti(:,:)
    real(KREAL),allocatable::rhoi(:,:),rhofi(:)
    real(KREAL),allocatable::ap(:)
    !write(*,*)'start transient calculation:'
    flag=1.0
    dt=ctime-ltime
    Ny=assm%mesh%Ny
    M=Ny+1
    N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
    allocate(pmodify(Ny),ui(Ny-1),Ti(M-1,N),rhoi(0:M,0:N),rhofi(0:M),ap(Ny-1))
    
    pmodify=0.0
    Ti=assm%thermal%temperature
    ui=assm%thermal%velocity
    rhoi=assm%property%rho
    rhofi=assm%property%rho(:,N)
    ap=0.0
    !call assm%th_boundary%update(ctime)
    call assm%pow%set(power,fq_core)
    j=0
    drho=1.0   
    call solve_temperature_rhoi(assm,flag,Ti,rhoi,dt)
	call cal_Tave(assm)
end subroutine driving_imp_THtransient

subroutine cal_Tave(assm)
	implicit none
	type(sys_assembly),intent(in out)::assm
	integer::j,k,N
	real(KREAL)::volumn,TVtotal
	real(KREAL)::TLtotal,Ltotal
	real(KREAL)::dr
	!计算元件的芯块平均温度和冷却剂平均温度
	N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
	dr=assm%geom%pellet/assm%mesh%Nf
	!fuel
	do j=1,assm%mesh%Ny,1
		volumn=0.0
		TVtotal=0.0
		do k=1,assm%mesh%Nf,1 !rod average					
			if (k==1) then
				TVtotal=TVtotal+assm%thermal%temperature(j,k)*PI*dr**2*assm%geom%height(j)
				volumn=volumn+PI*dr**2*assm%geom%height(j)
			else
				TVtotal=TVtotal+assm%thermal%temperature(j,k)*PI*((k*dr)**2-((k-1)*dr)**2)*assm%geom%height(j)
				volumn=volumn+PI*((k*dr)**2-((k-1)*dr)**2)*assm%geom%height(j)
			endif
		enddo
	enddo
	assm%thermal%Tfave=TVtotal/volumn
	!coolant
	TLtotal=0.0
	Ltotal=0.0
	do j=1,assm%mesh%Ny,1
		TLtotal=TLtotal+assm%thermal%temperature(j,N)*assm%geom%height(j)
		Ltotal=Ltotal+assm%geom%height(j)
	enddo
	assm%thermal%Tcave=TLtotal/Ltotal
end subroutine cal_Tave

end module imp_single_channel