module Imp_pipe_header
    use constants
    use imp_property
    use imp_mathkerel
    use lapack_interface
  
    type pipe
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
		real(KREAL)::rhos,shcs,ks
		!initial
		real(KREAL)::Ti
		!boundary
		real(KREAL)::Tfin,Tfout
        logical::is_Tb
        real(KREAL),allocatable::Tb(:),Bq(:)
    contains
        procedure,public::init=>init_pipe
        procedure,public::alloc=>alloc_pipe
        procedure,public::free=>free_pipe   
		procedure,public::updatep=>update_property
		procedure,public::updateb=>update_boundary
		procedure,public::kerels=>kerel_steady
		procedure,public::kerelt=>kerel_transient
		procedure,public::calhtc=>cal_htc
		procedure,public::thcals=>cal_thermal_steady
		procedure,public::thcalt=>cal_thermal_transient		
		procedure,public::cbuoy=>cal_buoy		
		procedure,public::cbeta=>cal_beta		
		procedure,public::cfric=>cal_fric		
    end type pipe

    private::init_pipe
    private::alloc_pipe
    private::free_pipe
	private::update_property
	private::update_boundary
	private::kerel_steady
	private::kerel_transient
	private::cal_htc
	private::cal_thermal_steady
	private::cal_thermal_transient
	private::cal_buoy
	private::cal_beta
	private::cal_fric
  contains
      subroutine init_pipe(this)
      implicit none
      class(pipe),intent(in out)::this
      !local
      integer::i,Ny 
      real(KREAL)::Aflow,rpipe
      Ny=this%Ny
	  Aflow=this%Q
	  rpipe=this%Rtube
	  this%Tfin=600.0!K
	  this%Tfout=this%Ti
      this%Bq=0.0
      !property
      this%rhos=get_density_304()
      this%shcs=get_shc_304()
      this%ks=get_conductivity_304()
      do i=1,Ny,1
          this%length(i)=this%ltotal/Ny
		  this%Tf(i)=this%Ti
		  this%Ts(i)=this%Ti
          this%rhof(i)=get_density_Na(this%Tf(i))
		  this%shcf(i)=get_shc_Na(this%Tf(i))
		  this%kf(i)=get_conductivity_Na(this%Tf(i))
		  this%visf(i)=get_vis_Na(this%Tf(i),this%rhof(i))
      end do
	  !area
	  this%area=PI*rpipe*rpipe
      !de
	  this%wet=2*PI*rpipe
      this%de=4.0*Aflow/this%wet
	  call this%calhtc()
	  !beta
	  !this%beta=0.5*(this%fric*this%ltotal/this%de+this%K)*1.0/(this%rhof*this%area**2)
    end subroutine init_pipe
	
    subroutine alloc_pipe(this)
      implicit none
      class(pipe),intent(in out)::this
      !local
      integer::Ny
      Ny=this%Ny
      !integer,intent(in)::N
      !check allocated first
      call Free_pipe(this)
      allocate(this%length(Ny))
	  allocate(this%Tf(Ny))
	  allocate(this%Ts(Ny))
	  allocate(this%htc(Ny))
	  allocate(this%rhof(Ny))
	  allocate(this%shcf(Ny))
	  allocate(this%visf(Ny))
	  allocate(this%kf(Ny))
	  allocate(this%Bq(Ny))
	  if(this%is_Tb) allocate(this%Tb(Ny))
      !allocate(this%rhof(Ny))
      !allocate(this%T(Ny))
    end subroutine alloc_pipe
    
    subroutine free_pipe(this)
      implicit none
      class(pipe),intent(in out)::this
      if(allocated(this%length)) deallocate(this%length)
	  if(allocated(this%Tf)) deallocate(this%Tf)
	  if(allocated(this%Ts)) deallocate(this%Ts)
	  if(allocated(this%htc)) deallocate(this%htc)
	  if(allocated(this%rhof)) deallocate(this%rhof)
	  if(allocated(this%shcf)) deallocate(this%shcf)
	  if(allocated(this%kf)) deallocate(this%kf)
	  if(allocated(this%visf)) deallocate(this%visf)
	  if(allocated(this%Tb)) deallocate(this%Tb)
	  if(allocated(this%Bq)) deallocate(this%Bq)
    end subroutine free_pipe
	
	subroutine cal_htc(this)
		implicit none
		class(pipe),intent(in out)::this
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
			Nu=get_Nusselt_Na_tube(area,wet,De,rhof(i),Q,visf(i),shcf(i),kf(i))!secondary
			this%htc(i)=Nu*kf(i)/De
		end do
		end associate
	end subroutine cal_htc
	
	subroutine update_property(this)
		implicit none
		class(pipe),intent(in out)::this
		!local
		integer::i,Ny
		Ny=this%Ny
		do i=1,Ny,1
			this%rhof(i)=get_density_Na(this%Tf(i))
			this%shcf(i)=get_shc_Na(this%Tf(i))
			this%kf(i)=get_conductivity_Na(this%Tf(i))
			this%visf(i)=get_vis_Na(this%Tf(i),this%rhof(i))
		enddo
	end subroutine update_property
	
	subroutine update_boundary(this,Q,Tfin)
		implicit none
		class(pipe),intent(in out)::this
		real(KREAL),intent(in)::Q,Tfin
		this%Q=Q
		this%Tfin=Tfin
	end subroutine update_boundary
	
	subroutine kerel_steady(this)
		implicit none
		class(pipe),intent(in out)::this
		!local
		integer::i,j,info,Ny
		real(KREAL)::aw,ap,ae,as,an,apz,S,ab
		real(KREAL)::ms!m per tube per Nyy of shell
		real(KREAL)::kfn,kfs
		real(KREAL),allocatable::A(:,:),B(:)
		real(KREAL),allocatable::Tfl(:),Tsl(:)

		associate( Rtube=>this%Rtube,&
				  thicks=>this%thicks,&
				  length=>this%length,&
				  rhos=>this%rhos,&
				  rhof=>this%rhof,&
				  shcs=>this%shcs,&
				  shcf=>this%shcf,&
                  ks=>this%ks,&
				  htc=>this%htc,  &
				  Bq=>this%Bq,  &
				  Q=>this%Q,      &
				  area=>this%area,&
				  kf=>this%kf,    &
				  Tb=>this%Tb,&
				  is_Tb=>this%is_Tb,&
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
			  if(i==1.or.i==Ny) then
				ae=htc(i)*PI*2*Rtube*Length(i)
				!apz=areas*this%length(i)*this%rhos(i)*this%shcs(i)/dt
				apz=0.0
				ap=Q*shcf(i)+ae+apz
				if(i==1)then
					as=0
					S=apz*Tfl(i)+Q*shcf(i)*Tfin
				else
					as=Q*shcf(i)
					S=apz*Tfl(i)
					A(i,i-1)=-as
				endif
			  else
				kfn=2*kf(i)*kf(i+1)/(kf(i)+kf(i+1))
				kfs=2*kf(i)*kf(i-1)/(kf(i)+kf(i-1))
				an=area*kfn/this%length(i)
			  	as=Q*shcf(i)+area*kfs/this%length(i)
				ae=htc(i)*PI*2*Rtube*Length(i)
				apz=0.0
				ap=an+as+ae+apz
				S=apz*Tfl(i)
				A(i,i+1)=-an
				A(i,i-1)=-as
			  endif
			  A(i,i)=ap
			  A(i,i+Ny)=-ae!tube
			  B(i)=S			  
			  case(2)!solid
			  aw=htc(i)*PI*2*Rtube*Length(i)
              if(is_Tb==.TRUE.) then
                ab=ks*PI*2*(Rtube+thicks)*Length(i)/(thicks/2.0)
			  else
                ab=0.0
              endif
              !apz=mv*this%shcv/dt
			  apz=0.0
			  ap=aw+apz+ab
              if(is_Tb==.True.) then
			     S=apz*Tsl(i)+ab*Tb(i)+Bq(i)
              else
                 S=apz*Tsl(i)+Bq(i)
              endif
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
		class(pipe),intent(in out)::this
		real(KREAL),intent(in)::last,current
		!local
		integer::i,j,info,Ny
        real(KREAL)::dt
		real(KREAL)::aw,ap,ae,as,an,apz,S,ab
		real(KREAL)::ms!m per tube per Nyy of shell
		real(KREAL)::kfn,kfs
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
				  Bq=>this%Bq,      &
				  Q=>this%Q,      &
				  kf=>this%kf,    &
				  ks=>this%ks,    &
				  Tb=>this%Tb,    &
				  is_Tb=>this%is_Tb,    &
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
			  if(i==1.or.i==Ny) then
				ae=htc(i)*PI*2*Rtube*Length(i)
				!apz=areas*this%length(i)*this%rhos(i)*this%shcs(i)/dt
				apz=0.0
				ap=Q*shcf(i)+ae+apz
				if(i==1)then
					as=0
					S=apz*Tfl(i)+Q*shcf(i)*Tfin
				else
					as=Q*shcf(i)
					S=apz*Tfl(i)
					A(i,i-1)=-as
				endif
			  else
				kfn=2*kf(i)*kf(i+1)/(kf(i)+kf(i+1))
				kfs=2*kf(i)*kf(i-1)/(kf(i)+kf(i-1))
				an=area*kfn/this%length(i)
			  	as=Q*shcf(i)+area*kfs/this%length(i)
				ae=htc(i)*PI*2*Rtube*Length(i)
				apz=area*length(i)*rhof(i)*shcf(i)/dt
				ap=an+as+ae+apz
				S=apz*Tfl(i)
				A(i,i+1)=-an
				A(i,i-1)=-as
			  endif
			  A(i,i)=ap
			  A(i,i+Ny)=-ae!tube
			  B(i)=S			  
			  case(2)!solid
			  aw=htc(i)*PI*2*Rtube*Length(i)
              if(is_Tb==.TRUE.) then
                ab=ks*PI*2*(Rtube+thicks)*Length(i)/(thicks/2.0)
              else
                ab=0.0
              endif
			  apz=ms*shcs/dt
			  !apz=0.0
			  ap=aw+apz
              if(is_Tb==.TRUE.) then
			      S=apz*Tsl(i)+ab*Tb(i)+Bq(i)
              else
                  S=apz*Tsl(i)+Bq(i)
              endif
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
		class(pipe),intent(in out)::this
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
        call this%updatep()
		end associate
	end subroutine cal_thermal_steady
	
	subroutine cal_thermal_transient(this,last,current)
		implicit none
		class(pipe),intent(in out)::this
		!real(KREAL),intent(in)::Q,Tfin
		real(KREAL),intent(in)::last,current
		call this%calhtc()
		!call this%updateb(Q,Tfin)
		call this%kerelt(last,current)	
        call this%updatep()!update property
	end subroutine cal_thermal_transient

    function cal_buoy(this) result (buoy)
		implicit none
		class(pipe),intent(in out)::this
        real(KREAL)::buoy,dbuoy
        integer::i
        buoy=0.0
        dbuoy=0.0
        associate(Ny=>this%Ny,&
                  rhof=>this%rhof,&
                  length=>this%Length,&
                  theta=>this%theta)
            do i=1,Ny,1
               dbuoy=-GRAVG*sin(theta/180.0*PI)*rhof(i)*length(i)
               buoy=buoy+dbuoy
            enddo
        end associate
    end function cal_buoy
        
    function cal_fric(this) result(fric)
        implicit none
		class(pipe),intent(in out)::this
        real(KREAL)::fric
        real(KREAL)::Re,visa,De
        integer::i
        associate(De=>this%De,&
                  Ny=>this%Ny,&
                  vis=>this%visf,&
                  flowrate=>this%Q,&
                  flowarea=>this%area)
            visa=0.0
            do i=1,Ny,1
                visa=visa+vis(i)/Ny
            enddo
            Re=flowrate*De/(visa*flowarea)
            if(Re>=1082) then
                fric=0.0055*(1.+(20000*1e-5/De+1e6/Re)**(1./3.))
            else
                fric=64./Re
            endif
        end associate
    end function cal_fric
    
    function cal_beta(this) result (beta)
        implicit none
		class(pipe),intent(in out)::this
        real(KREAL)::beta,rhoa
        integer::i
        beta=0.0
        rhoa=0.0
        associate(Ny=>this%Ny,&
                  fric=>this%fric,&
                  rho=>this%rhof,&
                  ltotal=>this%ltotal,&
                  De=>this%De,&
                  area=>this%area,&
                  K=>this%K)
            do i=1,Ny,1
                rhoa=rhoa+rho(i)/Ny
            enddo
            fric=this%cfric()
            beta=-fric*ltotal/(2*De*rhoa*area**2)-K/(2*rhoa*area**2)
        end associate
    end function cal_beta
end module Imp_pipe_header