module Imp_coremodle_header
    use constants
    use imp_property
    use imp_mathkerel
    use lapack_interface
    use Imp_coremodle_pow
  
    type coremodle
        !useful
        integer::Nflow
        integer::Nflowsemi!缩略模型半组件
        integer,allocatable::fzone(:)!zones which should be allocated flow
        real(KREAL)::sigmaPass
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
        real(KREAL)::fric
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
          this%rhof(i)=get_density(this%Tf(i))
		  this%shcf(i)=get_shc_LBE(this%Tf(i))
		  this%kf(i)=get_conductivity_LBE(this%Tf(i))
		  this%visf(i)=get_vis_LBE()
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
      integer::Ny,Nflow,Nsemi
      Ny=this%Ny
      Nflow=this%Nflow
      Nsemi=this%Nflowsemi
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
			Nu=get_Nusselt_PIPE_tube(area,wet,De,rhof(i),Q,visf(i),shcf(i),kf(i))!secondary
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
			this%rhof(i)=get_density(this%Tf(i))
			this%shcf(i)=get_shc_LBE(this%Tf(i))
			this%kf(i)=get_conductivity_LBE(this%Tf(i))
			this%visf(i)=get_vis_LBE()
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
    

end module Imp_coremodle_header