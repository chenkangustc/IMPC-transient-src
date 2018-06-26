module Imp_coremodle_header
    use constants
    use imp_property
    use imp_mathkerel
    use lapack_interface
    use Imp_coremodle_pow
    use imp_assm_global
  
    type coremodle
        !useful
        integer::Nzone !总组件数
        integer::Nflow
        integer::Nflowsemi!缩略模型半组件
        integer::Nsplit!计算时堆芯分成的份数
        integer::Nbranch!支路数
        real(KREAL)::Qtotal
        real(KREAL)::sigmaPass
        real(KREAL)::vqtotal!稳态冷却剂带走功率的计算值
        real(KREAL)::fric!
        integer,allocatable::fzone(:)!zones which should be allocated flow
        integer,allocatable::SAtable(:)
        integer::Mtl_fuel,Mtl_shell,Mtl_coolant,Mtl_gas
        integer::Nubundle
	    !geom 
        real(KREAL)::Ltotal
		real(KREAL),allocatable::Length(:)
        real(KREAL)::Rtube
		real(KREAL)::thicks
        real(KREAL)::theta
		real(KREAL)::area
		!mesh 
		integer::Ny
		!hydraulic
		real(KREAL)::De
		real(KREAL)::wet
        real(KREAL)::K
		real(KREAL)::Q
		real(KREAL)::beta
		!thermal
		real(KREAL),allocatable::Tf(:),Ts(:)
		real(KREAL),allocatable::htc(:)
		!material
		real(KREAL),allocatable::rhof(:),shcf(:),kf(:),visf(:)
		real(KREAL)::rhos,shcs
		!initial
		real(KREAL)::Ti
		!boundary
		real(KREAL)::Tfin,Tfout
    contains
        procedure,public::init=>init_coremodle
        procedure,public::alloc=>alloc_coremodle
        procedure,public::free=>free_coremodle   
		procedure,public::updatep=>update_property
		procedure,public::updateb=>update_boundary
		procedure,public::kerels=>kerel_steady
		procedure,public::kerelt=>kerel_transient
		procedure,public::calhtc=>cal_htc
		procedure,public::thcals=>cal_thermal_steady
		procedure,public::thcalt=>cal_thermal_transient		
		procedure,public::cbuoy=>cal_buoy		
		procedure,public::cbeta=>cal_beta		
		procedure,public::calpha=>cal_alpha		
		procedure,public::cfric=>cal_fric		
    end type coremodle

    private::init_coremodle
    private::alloc_coremodle
    private::free_coremodle
	private::update_property
	private::update_boundary
	private::kerel_steady
	private::kerel_transient
	private::cal_htc
	private::cal_thermal_steady
	private::cal_thermal_transient
	private::cal_buoy
	private::cal_beta
	private::cal_alpha
	private::cal_fric
  contains
      subroutine init_coremodle(this)
      implicit none
      class(coremodle),intent(in out)::this
      !local
      integer::i,Ny 
      real(KREAL)::Aflow,rcoremodle
      Ny=this%Ny
	  Aflow=this%Q
	  rcoremodle=this%Rtube
	  this%Tfin=600.0!K
	  this%Tfout=this%Ti
      do i=1,Ny,1
          this%length(i)=this%ltotal/Ny
		  this%Tf(i)=this%Ti
		  this%Ts(i)=this%Ti
          this%rhof(i)=get_density(this%Mtl_coolant,this%Tf(i))
		  this%shcf(i)=get_shc(this%Mtl_coolant,this%Tf(i))
		  this%kf(i)=get_conductivity(this%Mtl_coolant,this%Tf(i))
		  this%visf(i)=get_vis(this%Mtl_coolant,this%Tf(i))
      end do
	  !area
	  this%area=PI*rcoremodle*rcoremodle
      !de
	  this%wet=2*PI*rcoremodle
      this%de=4.0*Aflow/this%wet
	  call this%calhtc()
	  !beta
	  !this%beta=0.5*(this%fric*this%ltotal/this%de+this%K)*1.0/(this%rhof*this%area**2)
    end subroutine init_coremodle
	
    subroutine alloc_coremodle(this)
      implicit none
      class(coremodle),intent(in out)::this
      !local
      integer::Ny,Nflow,Nsemi,Nzone
      Ny=this%Ny
      Nflow=this%Nflow
      Nsemi=this%Nflowsemi
      Nzone=this%Nzone
      !integer,intent(in)::N
      !check allocated first
      call Free_coremodle(this)
      allocate(this%length(Ny))
	  allocate(this%Tf(Ny))
	  allocate(this%Ts(Ny))
	  allocate(this%htc(Ny))
	  allocate(this%rhof(Ny))
	  allocate(this%shcf(Ny))
	  allocate(this%visf(Ny))
	  allocate(this%kf(Ny))
      allocate(this%fzone(Nflow+Nsemi))
      !allocate(this%rhof(Ny))
      !allocate(this%T(Ny))
    end subroutine alloc_coremodle
    
    subroutine free_coremodle(this)
      implicit none
      class(coremodle),intent(in out)::this
      if(allocated(this%length)) deallocate(this%length)
	  if(allocated(this%Tf)) deallocate(this%Tf)
	  if(allocated(this%Ts)) deallocate(this%Ts)
	  if(allocated(this%htc)) deallocate(this%htc)
	  if(allocated(this%rhof)) deallocate(this%rhof)
	  if(allocated(this%shcf)) deallocate(this%shcf)
	  if(allocated(this%kf)) deallocate(this%kf)
	  if(allocated(this%visf)) deallocate(this%visf)
      if(allocated(this%fzone)) deallocate(this%fzone)
      if(allocated(this%SAtable)) deallocate(this%SAtable)
    end subroutine free_coremodle
	
	subroutine cal_htc(this)
		implicit none
		class(coremodle),intent(in out)::this
		!local
		integer::i
		real(KREAL)::Nu
		associate(Ny=>this%Ny,    &
		          area=>this%area,&
				  Q=>this%Q,      &
				  wet=>this%wet,  &
				  De=>this%De,    &
				  rhof=>this%rhof,&
				  shcf=>this%shcf,&
				  kf=>this%kf,    &
				  visf=>this%visf)
		do i=1,Ny,1
			Nu=get_Nusselt_tube(this%Mtl_coolant,this%Nubundle,area,wet,De,rhof(i),Q,visf(i),shcf(i),kf(i))!secondary
			this%htc(i)=Nu*kf(i)/De
		end do
		end associate
	end subroutine cal_htc
	
	subroutine update_property(this)
		implicit none
		class(coremodle),intent(in out)::this
		!local
		integer::i,Ny
		Ny=this%Ny
		do i=1,Ny,1
			this%rhof(i)=get_density(this%Mtl_coolant,this%Tf(i))
			this%shcf(i)=get_shc(this%Mtl_coolant,this%Tf(i))
			this%kf(i)=get_conductivity(this%Mtl_coolant,this%Tf(i))
			this%visf(i)=get_vis(this%Mtl_coolant,this%Tf(i))
		enddo
	end subroutine update_property
	
	subroutine update_boundary(this,Q,Tfin)
		implicit none
		class(coremodle),intent(in out)::this
		real(KREAL),intent(in)::Q,Tfin
		this%Q=Q
		this%Tfin=Tfin
	end subroutine update_boundary
	
	subroutine kerel_steady(this)
		implicit none
		class(coremodle),intent(in out)::this
		!local
		integer::i,j,info,Ny
		real(KREAL)::aw,ap,ae,as,apz,S
		real(KREAL)::ms!m per tube per Nyy of shell
		real(KREAL),allocatable::A(:,:),B(:)
		real(KREAL),allocatable::Tfl(:),Tsl(:)

		associate( Rtube=>this%Rtube,&
				  thicks=>this%thicks,&
				  length=>this%length,&
				  rhos=>this%rhos,&
				  rhof=>this%rhof,&
				  shcs=>this%shcs,&
				  shcf=>this%shcf,&
				  htc=>this%htc,  &
				  Q=>this%Q,      &
				  Tfin=>this%Tfin)
        Ny=this%Ny
        allocate(A(2*Ny,2*Ny))
		allocate(B(2*Ny))
		allocate(Tfl(Ny),Tsl(Ny))
		!dt=current-last
		A=0.0
		B=0.0
		Tfl=this%Tf
		Tsl=this%Ts

		!calculate h0
		do j=1,2,1!different part
			do i=1,Ny,1!different height
			!rpq=(this%areap+this%Ntube*PI*(this%Rtube+this%thickt)**2)/PI!radiau**2 r2 of primary 
			ms=PI*((Rtube+thicks)**2-Rtube**2)*Length(i)*rhos
			select case(j)
			  case(1)!fluid			  
			  ae=htc(i)*PI*2*Rtube*Length(i)
			  !apz=areas*this%length(i)*this%rhos(i)*this%shcs(i)/dt
			  apz=0.0
			  ap=Q*shcf(i)+ae+apz
			  if(i==1)then
				as=0
				S=apz*Tfl(i)+Q*shcf(i)*Tfin+pow/Ny
			  else
				as=Q*shcf(i)
				S=apz*Tfl(i)+pow/Ny
                A(i,i-1)=-as
			  endif
			  A(i,i)=ap
			  A(i,i+Ny)=-ae!tube
			  B(i)=S			  
			  case(2)!solid
			  aw=htc(i)*PI*2*Rtube*Length(i)
			  !apz=mv*this%shcv/dt
			  apz=0.0
			  ap=aw+apz
			  S=apz*Tsl(i)
			  A(i+Ny,i+Ny)=ap
			  A(i+Ny,i)=-aw
			  B(i+Ny)=S
			end select
			end do
        enddo
        
		call self_DGESV(A,B,info)
		end associate
		do i=1,2*Ny,1
			if(i<=Ny)then
				this%Tf(i)=B(i)
			else
				this%Ts(i-Ny)=B(i)
            endif
        end do
        this%Tfout=this%Tf(Ny)	
	end subroutine kerel_steady
	
	subroutine kerel_transient(this,last,current)
		implicit none
		class(coremodle),intent(in out)::this
		real(KREAL),intent(in)::last,current
		!local
		integer::i,j,info,Ny
        real(KREAL)::dt
		real(KREAL)::aw,ap,ae,as,apz,S
		real(KREAL)::ms!m per tube per Nyy of shell
		real(KREAL),allocatable::A(:,:),B(:)
		real(KREAL),allocatable::Tfl(:),Tsl(:)

		associate( Rtube=>this%Rtube,&
				  thicks=>this%thicks,&
				  length=>this%length,&
                  area=>this%area,&
				  rhos=>this%rhos,&
				  rhof=>this%rhof,&
				  shcs=>this%shcs,&
				  shcf=>this%shcf,&
				  htc=>this%htc,  &
				  Q=>this%Q,      &
				  Tfin=>this%Tfin)
           
        Ny=this%Ny
		
        allocate(A(2*Ny,2*Ny))
		allocate(B(2*Ny))
		allocate(Tfl(Ny),Tsl(Ny))
		dt=current-last
		A=0.0
		B=0.0
		Tfl=this%Tf
		Tsl=this%Ts
		!calculate h0
		do j=1,2,1!different part
			do i=1,Ny,1!different height
			ms=PI*((Rtube+thicks)**2-Rtube**2)*Length(i)*rhos
			select case(j)
			  case(1)!fluid			  
			  ae=htc(i)*PI*2*Rtube*Length(i)
			  apz=area*length(i)*rhof(i)*shcf(i)/dt
			  !apz=0.0
			  ap=Q*shcf(i)+ae+apz
			  if(i==1)then
				as=0
				S=apz*Tfl(i)+Q*shcf(i)*Tfin+pow/Ny
			  else
				as=Q*shcf(i)
				S=apz*Tfl(i)+pow/Ny
                A(i,i-1)=-as
			  endif
			  A(i,i)=ap
			  A(i,i+Ny)=-ae!tube
			  B(i)=S			  
			  case(2)!solid
			  aw=htc(i)*PI*2*Rtube*Length(i)
			  apz=ms*shcs/dt
			  !apz=0.0
			  ap=aw+apz
			  S=apz*Tsl(i)
			  A(i+Ny,i+Ny)=ap
			  A(i+Ny,i)=-aw
			  B(i+Ny)=S
			end select
			end do
        enddo
        
		call self_DGESV(A,B,info)
		end associate
		do i=1,2*Ny,1
			if(i<=Ny)then
				this%Tf(i)=B(i)
			else
				this%Ts(i-Ny)=B(i)
            endif
        end do
        this%Tfout=this%Tf(Ny)	
	end subroutine kerel_transient
	
	subroutine cal_thermal_steady(this)
		implicit none
		class(coremodle),intent(in out)::this
		!local
		integer::i
		real(KREAL)::sigma,sums,sf,ss
		real(KREAL),allocatable::Tfi(:),Tsi(:)
		sigma=1.0

		!local
		associate(Ny=>this%Ny,&
				  Tf=>this%Tf,&
				  Ts=>this%Ts)
        allocate(Tfi(Ny),Tsi(Ny))
		do while(sigma>=0.001)
			Tfi=Tf
			Tsi=Ts
			call this%updatep()
			call this%calhtc()
			call this%kerels()
			sums=0.0
			do i=1,Ny,1
				sf=abs((Tf(i)-Tfi(i))/Tfi(i))
				ss=abs((Ts(i)-Tsi(i))/Tsi(i))
				sums=sums+sf+ss
			enddo
			sigma=sums
		end do
		end associate
	end subroutine cal_thermal_steady
	
	subroutine cal_thermal_transient(this,last,current)
		implicit none
		class(coremodle),intent(in out)::this
		!real(KREAL),intent(in)::Q,Tfin
		real(KREAL),intent(in)::last,current
		
		call this%updatep()!update property
		call this%calhtc()
		!call this%updateb(Q,Tfin)
		call this%kerelt(last,current)	
	end subroutine cal_thermal_transient
    
    function cal_buoy(this) result (buoy)
		implicit none
		class(coremodle),intent(in out)::this
        real(KREAL)::buoy,dbuoy,rhoa
        integer::i,j
        integer::Nfluid,Nzone
        buoy=0.0
        dbuoy=0.0
        associate(Ny=>assm1(1)%mesh%Ny,&
                  Nf=>assm1(1)%mesh%Nf,&  
                  Ng=>assm1(1)%mesh%Ng,&  
                  Ns=>assm1(1)%mesh%Ns,&  
                  Nzone=>this%Nzone,&  
                  Nbranch=>this%Nbranch,&  
                  rho=>assm1(1)%property%rho,&
                  length=>assm1(1)%geom%height)
            Nfluid=Nf+Ng+Ns+1
            do i=1,Ny,1
                rhoa=0.0
                do j=1,Nzone,1
                    rhoa=rhoa+assm1(j)%property%rho(i,Nfluid)/Nzone
                enddo
                dbuoy=-GRAVG*sin(90.0/180.0*PI)*rhoa*length(i)
                buoy=buoy+dbuoy
            enddo
            buoy=buoy/Nbranch
    end associate
    end function cal_buoy
    
    function cal_fric(this) result(frics)
		class(coremodle),intent(in out)::this
        real(KREAL)::frics
        real(KREAL)::ltotal
        real(KREAL)::Re,visa,De,area
        integer::i,j
        integer::Nzone,Nfluid
        ! integer::Ftype,Frtype
        associate(Ny=>assm1(1)%mesh%Ny,&
                  Nf=>assm1(1)%mesh%Nf,&
                  Ng=>assm1(1)%mesh%Ng,&          
                  Ns=>assm1(1)%mesh%Ns,&          
                  Npin=>assm1(1)%geom%N_pin,&  
                  height=>assm1(1)%geom%height,&  
                  flowrate=>this%Qtotal,&
                  areap=>assm1(1)%hydrau%aflow,&
                  fric=>assm1(1)%hydrau%fric,&
                  Ftype=>assm1(1)%property%Mtl_coolant,&
                  Frtype=>assm1(1)%thermal%Frtype,&
                  Dep=>assm1(1)%hydrau%De)
                  
            Nzone=size(assm1)
            Nfluid=Nf+Ng+Ns+1
            ltotal=sum(height)
            De=Dep*Npin*Nzone
            area=areap*Npin*Nzone
          
            visa=0.0
            do i=1,Nzone,1
                do j=1,Ny,1
                    visa=visa+assm1(i)%property%dvs(j,Nfluid)*height(j)/(ltotal*Nzone)
                enddo
            enddo
            Re=flowrate*De/(visa*area)
            frics=get_fric_pin(Ftype,Frtype,Re)
            this%fric=frics
            ! assm1(:)%hydrau%fric=frics
            ! if(Re>=2050.) then
                ! frics=0.1875/Re**2
            ! else
                ! frics=76.5/Re
            ! endif
        end associate
    end function cal_fric
    
    function cal_beta(this) result (beta)
        implicit none
		class(coremodle),intent(in out)::this
        real(KREAL)::beta,rhoa,De
        real(KREAL)::ltotal,area
        integer::i,j
        integer::Nzone,Nfluid
        beta=0.0
        rhoa=0.0
        associate(Ny=>assm1(1)%mesh%Ny,&
                  Nf=>assm1(1)%mesh%Nf,&
                  Ng=>assm1(1)%mesh%Ng,&          
                  Ns=>assm1(1)%mesh%Ns,&          
                  Npin=>assm1(1)%geom%N_pin,&          
                  ! rho=>assm1(1)%property%rho,&
                  height=>assm1(1)%geom%height,&
                  areap=>assm1(1)%hydrau%aflow,&
                  fric=>assm1(1)%hydrau%fric,&
                  Dep=>assm1(1)%hydrau%De,&
                  Nbranch=>this%Nbranch,&
                  K=>assm1(1)%hydrau%K)
            Nzone=size(assm1)
            Nfluid=Nf+Ng+Ns+1
            De=Dep*Npin*Nzone
            area=areap*Npin*Nzone
            ltotal=sum(height)
            rhoa=0.0
            do i=1,Nzone,1
                do j=1,Ny,1
                    rhoa=rhoa+assm1(i)%property%rho(j,Nfluid)*height(j)/(ltotal*Nzone)
                enddo
            enddo
            fric=this%cfric()
            beta=-Nbranch*fric*ltotal/(2*De*rhoa*area**2)-K/(2*rhoa*area**2)
    end associate
    end function cal_beta
    
    function cal_alpha(this) result(alpha)
        implicit none
		class(coremodle),intent(in out)::this
        real(KREAL)::alpha
        real(KREAL)::area,ltotal
        integer::Nzone,Npin
        associate(height=>assm1(1)%geom%height,&
                  Npin=>assm1(1)%geom%N_pin,&
                  areap=>assm1(1)%hydrau%aflow)
            Nzone=size(assm1)
            ltotal=sum(height)
            area=areap*Npin*Nzone
            alpha=ltotal/area
        end associate
    end function cal_alpha
end module Imp_coremodle_header