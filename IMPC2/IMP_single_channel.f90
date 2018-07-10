module imp_single_channel
    use constants
	use imp_property
    use imp_assm_global
    use imp_single_kernel
    use imp_assembly_header
    use imp_loop_global
    implicit none
     private
     public::driving_imp_steady
     public::driving_imp_Transient
	 public::driving_imp_THsteady   !only energy formula is calculated
	 public::driving_imp_THtransient!no pv step ,not simple
	 public::driving_imp_flowAlloc
     public::update_property_rhoi
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
			density=get_density(assm(i)%property%Mtl_coolant,assm(i)%th_boundary%T%inlet)
			assm(i)%th_boundary%u%inlet=flowrate(i)/(n_pin*fuelArea*density)
		enddo
	end subroutine driving_imp_flowAlloc
	
subroutine driving_imp_steady(assm,power,fq_core)
    type(sys_assembly),intent(in out)::assm
    real(KREAL),intent(in)::power(:,:)
    real(KREAL),intent(in)::fq_core(:,:)
    !local rhoi/ui/Ti/dt/ap/pmodify/
    integer M,N,Ny,i,j
    real(KREAL):: flag,dt
    real(KREAL):: btotal,drho!判断因子
    real(KREAL),allocatable::pmodify(:)
    real(KREAL),allocatable::ui(:),Ti(:,:),hotTi(:,:)
    real(KREAL),allocatable::rhoi(:,:),rhofi(:)
    real(KREAL),allocatable::ap(:)
    !write(*,*)'start steady calculation:'
    flag=0.0
    dt=1.0!为保证方程求解dt不为0，无具体意义+
    Ny=assm%mesh%Ny
    M=Ny+1
    N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
    allocate(pmodify(Ny),ui(Ny-1),Ti(M-1,N),hotTi(M-1,N),rhoi(0:M,0:N),rhofi(0:M),ap(Ny-1))
 
    pmodify=0.0
    ui=assm%thermal%velocity
    Ti=assm%thermal%temperature
    hotTi=assm%hotthermal%temperature
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
      call solve_temperature_rhoi(assm,flag,Ti,hotTi,rhoi,dt)
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
    real(KREAL),allocatable::ui(:),Ti(:,:),hotTi(:,:)
    real(KREAL),allocatable::rhoi(:,:),rhofi(:)
    real(KREAL),allocatable::ap(:)
    !write(*,*)'start steady calculation:'
    flag=0.0
    dt=1.0!为保证方程求解dt不为0，无具体意义+
    Ny=assm%mesh%Ny
    M=Ny+1
    N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
    allocate(pmodify(Ny),ui(Ny-1),Ti(M-1,N),hotTi(M-1,N),rhoi(0:M,0:N),rhofi(0:M),ap(Ny-1))
 
    pmodify=0.0
    ui=assm%thermal%velocity
    Ti=assm%thermal%temperature
    hotTi=assm%hotthermal%temperature
    rhoi=assm%property%rho
    rhofi=assm%property%rho(:,N)
    ap=0.0
    j=0
    drho=1.0
    call assm%pow%set(power,fq_core)
    call solve_temperature_rhoi(assm,flag,Ti,hotTi,rhoi,dt)
    call update_property_rhoi(assm)
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
    real(KREAL),allocatable::ui(:),Ti(:,:),hotTi(:,:)
    real(KREAL),allocatable::rhoi(:,:),rhofi(:)
    real(KREAL),allocatable::ap(:)
    !write(*,*)'start transient calculation:'
    flag=1.0
    dt=ctime-ltime
    Ny=assm%mesh%Ny
    M=Ny+1
    N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
    allocate(pmodify(Ny),ui(Ny-1),Ti(M-1,N),hotTi(M-1,N),rhoi(0:M,0:N),rhofi(0:M),ap(Ny-1))
    
    pmodify=0.0
    Ti=assm%thermal%temperature
    hotTi=assm%hotthermal%temperature
    ui=assm%thermal%velocity
    rhoi=assm%property%rho
    rhofi=assm%property%rho(:,N)
    ap=0.0
    !call assm%th_boundary%update(ctime)
    call assm%pow%set(power,fq_core)
    j=0
    drho=1.0   
    call solve_temperature_rhoi(assm,flag,Ti,hotTi,rhoi,dt)
	call cal_Tave(assm)
    call update_property_rhoi(assm)
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
	volumn=0.0
	TVtotal=0.0
	do j=1,assm%mesh%Ny,1
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
	!Tcave
	TLtotal=0.0
	Ltotal=0.0
	do j=1,assm%mesh%Ny,1
		TLtotal=TLtotal+assm%thermal%temperature(j,N)*assm%geom%height(j)
		Ltotal=Ltotal+assm%geom%height(j)
	enddo
	assm%thermal%Tcave=TLtotal/Ltotal
	!Tcoolant
	do j=1,assm%mesh%Ny,1
	   assm%thermal%Tcoolant(j)=assm%thermal%temperature(j,N)
	enddo
	!Tfuel
	do j=1,assm%mesh%Ny,1
		TVtotal=0.0
		volumn=0.0
		do k=1,assm%mesh%Nf,1 !rod average					
			if (k==1) then
				TVtotal=TVtotal+assm%thermal%temperature(j,k)*PI*dr**2*assm%geom%height(j)
				volumn=volumn+PI*dr**2*assm%geom%height(j)
			else
				TVtotal=TVtotal+assm%thermal%temperature(j,k)*PI*((k*dr)**2-((k-1)*dr)**2)*assm%geom%height(j)
				volumn=volumn+PI*((k*dr)**2-((k-1)*dr)**2)*assm%geom%height(j)
			endif
		enddo
		assm%thermal%Tfuel(j)=TVtotal/volumn
	enddo
end subroutine cal_Tave

subroutine update_property_rhoi(assm)
    type(sys_assembly),intent(in out)::assm
        ! real(KREAL),allocatable::rho(:,:)!热物性
        ! real(KREAL),allocatable::shc(:,:)
        ! real(KREAL),allocatable::ctc(:,:)
    !local 
    integer::Ny,Nr,i,j
    Ny=size(assm%thermal%temperature,dim=1)
    Nr=size(assm%thermal%temperature,dim=2)

    do j=1,Ny,1
        do i=1,Nr,1
            if(i<=assm%mesh%Nf) then!fuel
                assm%property%rho(j,i)=get_density(assm%property%Mtl_fuel,assm%thermal%temperature(j,i))
                assm%property%shc(j,i)=get_shc(assm%property%Mtl_fuel,assm%thermal%temperature(j,i))
                assm%property%ctc(j,i)=get_conductivity(assm%property%Mtl_fuel,assm%thermal%temperature(j,i))
            elseif(i<=assm%mesh%Nf+assm%mesh%Ng) then!gap
                assm%property%rho(j,i)=get_density(assm%property%Mtl_gas,assm%thermal%temperature(j,i))
                assm%property%shc(j,i)=get_shc(assm%property%Mtl_gas,assm%thermal%temperature(j,i))
                assm%property%ctc(j,i)=get_conductivity(assm%property%Mtl_gas,assm%thermal%temperature(j,i))  
            elseif(i<=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns) then!clad
                assm%property%rho(j,i)=get_density(assm%property%Mtl_shell,assm%thermal%temperature(j,i))
                assm%property%shc(j,i)=get_shc(assm%property%Mtl_shell,assm%thermal%temperature(j,i))
                assm%property%ctc(j,i)=get_conductivity(assm%property%Mtl_shell,assm%thermal%temperature(j,i))  
            else!fluid
                assm%property%rho(j,Nr)=get_density(assm%property%Mtl_coolant,assm%thermal%temperature(j,Nr))
                assm%property%shc(j,Nr)=get_shc(assm%property%Mtl_coolant,assm%thermal%temperature(j,Nr))
                assm%property%ctc(j,Nr)=get_conductivity(assm%property%Mtl_coolant,assm%thermal%temperature(j,Nr))
                assm%property%dvs(j,Nr)=get_vis(assm%property%Mtl_coolant,assm%thermal%temperature(j,Nr),assm%property%rho(j,Nr))
            endif
        enddo
    enddo
    !up and down boundary
    assm%property%rho(0,:)=assm%property%rho(1,:)
    assm%property%shc(0,:)=assm%property%shc(1,:)
    assm%property%ctc(0,:)=assm%property%ctc(1,:)
    
    assm%property%rho(Ny+1,:)=assm%property%rho(Ny,:)
    assm%property%shc(Ny+1,:)=assm%property%shc(Ny,:)
    assm%property%ctc(Ny+1,:)=assm%property%ctc(Ny,:)
end subroutine update_property_rhoi

end module imp_single_channel