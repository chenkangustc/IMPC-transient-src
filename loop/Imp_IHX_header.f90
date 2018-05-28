!type::IHX
!
!type subroutine
! 			private::init_IHX
! 			private::alloc_IHX
! 			private::free_IHX
! 			private::cal_thermal_steady
! 			private::cal_thermal_transient
!*******************************************************************************
module Imp_IHX_header
	use constants
	use lapack_interface
	use Imp_property
	use Imp_mathkerel

	implicit none
	
	type IHX
		!geom
		real(KREAL),allocatable::Length(:)
		real(KREAL)::Lsingle
		real(KREAL),allocatable::zz(:)
		real(KREAL)::Rtube,Rp,Rpa,thickt,thickv
		real(KREAL)::Plength!中心节距
		real(KREAL)::Areap
		integer::Ntube
		real(KREAL)::AreaTubeSingle,AreaTubeTotal
		!mesh
		integer::N
		!hydraulic
		real(KREAL)::Dep,Fricp,Kricp,Wetp,betap,Wets,Des!beta is used to calculate the flowrate
		real(KREAL)::Qp,Qs
		!material
		real(KREAL),allocatable::rhop(:),rhos(:)
		real(KREAL),allocatable::shcp(:),shcs(:)
        real(KREAL),allocatable::visp(:),viss(:)
		real(KREAL),allocatable::kp(:),ks(:)
		real(KREAL)::rhot,rhov
		real(KREAL)::shct,shcv
		!thermal
		real(KREAL),allocatable::htp(:),hvp(:),hts(:)
		real(KREAL),allocatable::Tp(:),Ts(:),Tt(:),Tv(:)
		!initial
		real(KREAL)::Ti
		!boundary
		real(KREAL)::Tpin,Tsin
        real(KREAL)::Tpout,Tsout
        integer::Nftime,NTtime
        real(KREAL),allocatable::flowtable(:,:),Tintable(:,:)
        !control
        logical::is_Tintable,is_flowtable
	  contains
		procedure,public::init=>init_IHX
        procedure,public::alloc=>alloc_IHX
        procedure,public::free=>free_IHX
		procedure,public::thcals=>cal_thermal_steady
		procedure,public::thcalt=>cal_thermal_transient
        procedure,public::update=>update_property
        procedure,public::CalHtc=>cal_htc
        procedure,public::kerel=>cal_thermal_kerel
		procedure,public::kerelTransient=>cal_thermal_kerel_transient
		procedure,public::updateBoundary=>update_Boundary
		procedure,public::cbuoy=>cal_buoy
		procedure,public::cbeta=>cal_beta
		procedure,public::cfric=>cal_fric
	end type IHX
	
	private::init_IHX
    private::alloc_IHX
    private::free_IHX
	private::cal_thermal_steady
	private::cal_thermal_transient
	private::cal_thermal_kerel
	private::update_property
	private::cal_htc
	private::cal_thermal_kerel_transient
	private::update_Boundary
	private::cal_buoy
	private::cal_beta
  contains
	subroutine init_IHX(this)
		implicit none
		class(IHX),intent(in out)::this
		!local
		integer::i,N,Ntube
		real(KREAL)::Lsingle,rtube,Qp,Kiouter,Rpa
        real(KREAL)::Ts,Tp
		!real(KREAL),allocatable::temperature(:)
		
		Lsingle=this%Lsingle
		Qp=this%Qp
		Kiouter=this%Kricp/this%N
		rtube=this%Rtube
		Ntube=this%Ntube
		N=size(this%Tp)
		!allocate(temperature(1:N))
		
		this%Tp(:)=this%Ti
		this%Ts(:)=this%Ti
		this%Tt(:)=this%Ti
		this%Tv(:)=this%Ti
		this%AreaTubeSingle=PI*rtube*rtube
		this%AreaTubeTotal=Ntube*this%AreaTubeSingle
		
		this%Tpin=663.15!K
        this%Tpout=this%Ti
        this%Tsout=this%Ti
		
		this%Rp=sqrt((this%areap+this%Ntube*PI*(this%Rtube+this%thickt)**2)/PI)!for mass_shell calculation
		this%Rpa=sqrt((this%Areap/this%Ntube+PI*(this%Rtube+this%thickt)**2)/PI)!single modle ridiau of primary
		this%Wetp=2*PI*(this%Rtube+this%thickt+this%Rpa)
	    this%Dep=4*this%Areap/(this%Wetp*this%Ntube)	
		this%betap=0.0	
		this%Wets=2*PI*(this%Rtube+this%thickt)
		this%Des=4*this%AreaTubeSingle/this%Wets
        !property
        this%rhot=get_density_304()
        this%shct=get_shc_304()
        this%rhov=get_density_304()
        this%shcv=get_shc_304()

		do i=1,N,1
			Tp=this%Tp(i)
            Ts=this%Ts(i)
			this%rhop(i)=get_density_Na(Tp)
			this%shcp(i)=get_shc_Na(Tp)
            this%visp(i)=get_vis_Na(Tp,this%rhop(i))
			this%kp(i)=get_conductivity_Na(Tp)		
			this%rhos(i)=get_density_Na(Ts)
			this%shcs(i)=get_shc_Na(Ts)
		    this%viss=get_vis_Na(Ts,this%rhos(i))
			this%ks(i)=get_conductivity_Na(Ts)
			this%Length(i)=Lsingle/N
			if(i==1)then
				this%zz(i)=this%Length(i)/2
			else
				this%zz(i)=this%zz(i-1)+this%Length(i)
			endif
			!betap
			!this%betap=this%betap+0.5*(this%Fricp*this%length(i)/this%Dep+Kiouter)*1/(this%rhop(i)*this%Areap**2)
		end do
		call this%CalHtc()
		!if(allocated(temperature)) deallocate(temperature)
	end subroutine init_IHX
	
	subroutine alloc_IHX(this)
      implicit none
      class(IHX),intent(in out)::this
      !local
      integer::N
      N=this%N
      !integer,intent(in)::N
      !check allocated first
      call Free_IHX(this)
      allocate(this%length(1:N))
	  allocate(this%zz(1:N))
      allocate(this%rhop(1:N))
	  allocate(this%shcp(1:N))
      allocate(this%visp(1:N))
	  allocate(this%kp(1:N))
	  allocate(this%rhos(1:N))
	  allocate(this%shcs(1:N))
      allocate(this%viss(1:N))
	  allocate(this%ks(1:N))
      allocate(this%Tp(1:N))
	  allocate(this%Ts(1:N))
	  allocate(this%Tt(1:N))
	  allocate(this%Tv(1:N))
	  
	  allocate(this%htp(1:N))
	  allocate(this%hvp(1:N))
	  allocate(this%hts(1:N))
      if(this%is_flowtable) allocate(this%flowtable(2,this%Nftime))
      if(this%is_Tintable)  allocate(this%Tintable(2,this%NTtime))

    end subroutine alloc_IHX
    
    subroutine free_IHX(this)
      implicit none
      class(IHX),intent(in out)::this
      if(allocated(this%length)) deallocate(this%length)
	  if(allocated(this%zz)) deallocate(this%zz)
      if(allocated(this%rhop)) deallocate(this%rhop)
	  if(allocated(this%shcp)) deallocate(this%shcp)
      if(allocated(this%visp)) deallocate(this%visp)
	  if(allocated(this%kp)) deallocate(this%kp)
	  if(allocated(this%rhos)) deallocate(this%rhos)
	  if(allocated(this%shcs)) deallocate(this%shcs)
      if(allocated(this%viss)) deallocate(this%viss)
	  if(allocated(this%ks)) deallocate(this%ks)
      if(allocated(this%Tp)) deallocate(this%Tp)
	  if(allocated(this%Ts)) deallocate(this%Ts)
	  if(allocated(this%Tt)) deallocate(this%Tt)
	  if(allocated(this%Tv)) deallocate(this%Tv)
	  if(allocated(this%htp)) deallocate(this%htp)
	  if(allocated(this%hvp)) deallocate(this%hvp)
	  if(allocated(this%hts)) deallocate(this%hts)
	  if(allocated(this%flowtable)) deallocate(this%flowtable)
	  if(allocated(this%Tintable)) deallocate(this%Tintable)
    end subroutine free_IHX
	
    subroutine update_property(this)
		implicit none
		class(IHX),intent(in out)::this
		!local
		integer::i,N
		N=this%N
		do i=1,N,1
			this%rhop(i)=get_density_Na(this%Tp(i))
			this%shcp(i)=get_shc_Na(this%Tp(i))
			this%kp(i)=get_conductivity_Na(this%Tp(i))
			this%visp(i)=get_vis_Na(this%Tp(i),this%rhop(i))
			
            this%rhos(i)=get_density_Na(this%Ts(i))
			this%shcs(i)=get_shc_Na(this%Ts(i))
			this%ks(i)=get_conductivity_Na(this%Ts(i))
			this%viss(i)=get_vis_Na(this%Ts(i),this%rhos(i))
		enddo
	end subroutine update_property
    
	subroutine cal_thermal_steady(this)
		implicit none
		class(IHX),intent(in out)::this
		!local
		integer::i,N
		real(KREAL)::sigma,sums,sp,ss,st,sv
		real(KREAL),allocatable::Tpi(:),Tsi(:),Tti(:),Tvi(:)
		N=this%N
		sigma=1.0
		allocate(Tpi(N),Tsi(N),Tti(N),Tvi(N))

		do while(sigma>=0.001)!该迭代为物性迭代
			Tpi=this%Tp
			Tsi=this%Ts
			Tti=this%Tt
			Tvi=this%Tv
			call this%update()
			call this%calhtc()
			call this%kerel()
			sums=0.0
			do i=1,N,1
				sp=abs((this%Tp(i)-Tpi(i))/Tpi(i))
				ss=abs((this%Ts(i)-Tsi(i))/Tsi(i))
				st=abs((this%Tt(i)-Tti(i))/Tti(i))
				sv=abs((this%Tv(i)-Tvi(i))/Tvi(i))
				sums=sums+sp+ss+st+sv
			enddo
			sigma=sums
		end do
		call this%update()
	end subroutine cal_thermal_steady
	
	subroutine cal_htc(this)
		implicit none
		class(IHX),intent(in out)::this
		!local
		integer::i,N
		real(KREAL)::Nust,Nupt,Nupv
		real(KREAL)::Areap_,Areas_,Qp_,Qs_
		real(KREAL)::Dtubeo
		N=this%N
		Areap_=this%Areap/this%Ntube
		Areas_=this%AreaTubeSingle
		Qp_=this%Qp/this%Ntube
		Qs_=this%Qs/this%Ntube
		Dtubeo=(this%Rtube+this%thickt)*2.0
		!get_Nusselt_IHX_shell(flowarea,wet,De,rho,flowrate,vis,shc,conductivity)
		do i=1,N,1
			Nust=get_Nusselt_Na_tube(Areas_,this%wets,this%Des,this%rhos(i),Qs_,this%viss(i),this%shcs(i),this%ks(i))!secondary
			Nupt=get_Nusselt_Na_bundle(this%Plength/Dtubeo,Areap_,this%wetp,this%Dep,this%rhop(i),Qp_,this%visp(i),this%shcp(i),this%kp(i))!primary, assume that the htc between p and t is as same as that between p and v
			Nupv=get_Nusselt_Na_bundle(this%Plength/Dtubeo,Areap_,this%wetp,this%Dep,this%rhop(i),Qp_,this%visp(i),this%shcp(i),this%kp(i))	
			this%hts(i)=Nust*this%ks(i)/this%Des
			this%htp(i)=Nupt*this%kp(i)/this%Dep
			this%hvp(i)=Nupv*this%kp(i)/this%Dep
		end do
	end subroutine cal_htc
	
	
	subroutine cal_thermal_kerel(this)
		implicit none
		class(IHX),intent(in out)::this
		!integer,intent(in)::key!if key==0 steady,else if key==1,transient
		!local
		integer::i,j,info
		integer::N
		real(KREAL)::an,aw,ap,ae,as,apz,S
		real(KREAL)::areap,areas
		real(KREAL)::mt,mv!m per tube per Ny of tube or shell
		real(KREAL)::kpn,kps,ksn,kss
		!real(KREAL)::dt
		real(KREAL)::Qsa,Qpa
		real(KREAL),allocatable::A(:,:),B(:)
		real(KREAL),allocatable::Tpl(:),Tsl(:),Ttl(:),Tvl(:)
		
		N=this%N
		allocate(A(4*N,4*N))
		allocate(B(4*N))
		allocate(Tpl(N),Tsl(N),Ttl(N),Tvl(N))
		
		Qsa=this%Qs/this%Ntube
		Qpa=this%Qp/this%Ntube
        areap=this%areap/this%Ntube!p 
		areas=this%AreaTubeSingle!s
		!dt=current-last
		A=0.0
		B=0.0
		 Tpl=this%Tp
		 Tsl=this%Ts
		 Ttl=this%Tt
		 Tvl=this%Tv

		!calculate h0
		do j=1,4,1!different part
			do i=1,N,1!different height
			!cmp=this%rho(i)!(C*M)
			mt=PI*((this%Rtube+this%thickt)**2-this%Rtube**2)*this%Length(i)*this%rhot
		    !assmue that the IHX is cylindrical
			!rpq=(this%areap+this%Ntube*PI*(this%Rtube+this%thickt)**2)/PI!radiau**2 r2 of primary 
			mv=PI*((this%Rp+this%thickv)**2-this%Rp**2)*this%Length(i)*this%rhov

			select case(j)
			  case(1)!primary
			  if(i==1.or.i==N) then!inlet and outlet
				aw=this%htp(i)*PI*2*(this%Rtube+this%thickt)*this%Length(i)
				ae=this%hvp(i)*PI*2*this%Rpa*this%Length(i)
				!apz=areap*this%length(i)*this%rhop(i)*this%shcp(i)/dt
				apz=0.0
				ap=Qpa*this%shcp(i)+aw+ae+apz
				if(i==N)then
					an=0.0
					S=apz*Tpl(i)+Qpa*this%shcp(i)*this%Tpin
				else
					an=Qpa*this%shcp(i)
					S=apz*Tpl(i)
					A(i,i+1)=-an
				endif
			  else!inner:i=[2,Ny-1]
				kpn=2*this%kp(i)*this%kp(i+1)/(this%kp(i)+this%kp(i+1))
				kps=2*this%kp(i)*this%kp(i-1)/(this%kp(i)+this%kp(i-1))
				an=Qpa*this%shcp(i)+areap/this%Length(i)*kpn
				as=areas/this%Length(i)*kps
			  	aw=this%htp(i)*PI*2*(this%Rtube+this%thickt)*this%Length(i)
				ae=this%hvp(i)*PI*2*this%Rpa*this%Length(i)
				apz=0.0
				ap=an+as+aw+ae+apz				
				S=apz*Tpl(i)
				A(i,i+1)=-an
				A(i,i-1)=-as
			  endif		
			  A(i,i)=ap
			  A(i,i+2*N)=-aw!tube
			  A(i,i+3*N)=-ae!shell
			  B(i)=S
			  case(2)!secondary
			  if(i==1.or.i==N) then
				ae=this%hts(i)*PI*2*this%Rtube*this%Length(i)
				!apz=areas*this%length(i)*this%rhos(i)*this%shcs(i)/dt
				apz=0.0
				ap=Qsa*this%shcs(i)+ae+apz
				if(i==1)then
					as=0
					S=apz*Tsl(i)+Qsa*this%shcs(i)*this%Tsin
				else
					as=Qsa*this%shcs(i)
					S=apz*Tsl(i)
					A(i+N,i-1+N)=-as
				endif
			  else
			  	ksn=2*this%ks(i)*this%ks(i+1)/(this%ks(i)+this%ks(i+1))
				kss=2*this%ks(i)*this%ks(i-1)/(this%ks(i)+this%ks(i-1))
			  	an=areas/this%length(i)*ksn
				as=Qsa*this%shcs(i)+areas/this%length(i)*kss
				ae=this%hts(i)*PI*2*this%Rtube*this%Length(i)
				apz=0.0
				ap=an+as+ae+apz
				S=apz*Tsl(i)
				A(i+N,i+1+N)=-an
				A(i+N,i-1+N)=-as
			  endif
              A(i+N,i+N)=ap			  
			  A(i+N,i+2*N)=-ae!tube
			  B(i+N)=S
			  
			  case(3)!tube
			  aw=this%hts(i)*PI*2*this%Rtube*this%Length(i)
			  ae=this%htp(i)*PI*2*(this%Rtube+this%thickt)*this%Length(i)
			  !apz=mt*this%shct/dt
			  apz=0.0
			  ap=aw+ae+apz
			  S=apz*Ttl(i)
			  A(i+2*N,i+2*N)=ap
			  A(i+2*N,i+N)=-aw
			  A(i+2*N,i)=-ae
		      B(i+2*N)=S
			  case(4)!shell
			  aw=this%hvp(i)*PI*2*this%Rpa*this%Length(i)
			  !apz=mv*this%shcv/dt
			  apz=0.0
			  ap=aw+apz
			  S=apz*Tvl(i)
			  A(i+3*N,i+3*N)=ap
			  A(i+3*N,i)=-aw
			  B(i+3*N)=S
			end select
			end do
        enddo
        
		call self_DGESV(A,B,info)
		
		do i=1,4*N,1
			if(i<=N)then
				this%Tp(i)=B(i)
			elseif(i<=2*N)then
				this%Ts(i-N)=B(i)
			elseif(i<=3*N)then
				this%Tt(i-2*N)=B(i)
			else
				this%Tv(i-3*N)=B(i)
            endif
        end do
        this%Tpout=this%Tp(1)
        this%Tsout=this%Ts(N)
	end subroutine cal_thermal_kerel
	subroutine cal_thermal_transient(this,last,current)
		implicit none
		class(IHX),intent(in out)::this
		!real(KREAL),intent(in)::Qp,Qs,Tpin,Tsin
		real(KREAL),intent(in)::last,current

		call this%calhtc()
		!call this%updateBoundary(Qp,Qs,Tpin,Tsin)
		call this%kerelTransient(last,current)	
        !update property 
		call this%update()       
	end subroutine cal_thermal_transient
	
	subroutine cal_thermal_kerel_transient(this,last,current)
		implicit none
		class(IHX),intent(in out)::this
		real(KREAL),intent(in)::last,current
		!local
		integer::i,j,info
		integer::N
		real(KREAL)::an,aw,ap,ae,as,apz,S
		real(KREAL)::areap,areas
		real(KREAL)::mt,mv!m per tube per Ny of tube or shell
        real(KREAL)::kpn,kps,ksn,kss
		real(KREAL)::dt
		real(KREAL)::Qsa,Qpa
		real(KREAL),allocatable::A(:,:),B(:)
		real(KREAL),allocatable::Tpl(:),Tsl(:),Ttl(:),Tvl(:)
		
		N=this%N
		allocate(A(4*N,4*N))
		allocate(B(4*N))
		allocate(Tpl(N),Tsl(N),Ttl(N),Tvl(N))
		
		Qsa=this%Qs/this%Ntube
		Qpa=this%Qp/this%Ntube
        areap=this%areap/this%Ntube
		areas=this%AreaTubeSingle
		dt=current-last
		A=0.0
		B=0.0
		Tpl=this%Tp
		Tsl=this%Ts
		Ttl=this%Tt
		Tvl=this%Tv
		!calculate h0
		do j=1,4,1!different part
			do i=1,N,1!different height
			!cmp=this%rho(i)!(C*M)
			mt=PI*((this%Rtube+this%thickt)**2-this%Rtube**2)*this%Length(i)*this%rhot
		    !assmue that the IHX is cylindrical
			!rpq=(this%areap+this%Ntube*PI*(this%Rtube+this%thickt)**2)/PI!radiau**2 r2 of primary 
			mv=PI*((this%Rp+this%thickv)**2-this%Rp**2)*this%Length(i)*this%rhov
			select case(j)
			  case(1)!primary
			  if(i==1.or.i==N) then!inlet and outlet
				aw=this%htp(i)*PI*2*(this%Rtube+this%thickt)*this%Length(i)
				ae=this%hvp(i)*PI*2*this%Rpa*this%Length(i)
				apz=areap*this%length(i)*this%rhop(i)*this%shcp(i)/dt
				ap=Qpa*this%shcp(i)+aw+ae+apz
				if(i==N)then
					an=0.0
					S=apz*Tpl(i)+Qpa*this%shcp(i)*this%Tpin
				else
					an=Qpa*this%shcp(i)
					S=apz*Tpl(i)
					A(i,i+1)=-an
				endif
			  else!inner:i=[2,Ny-1]
				kpn=2*this%kp(i)*this%kp(i+1)/(this%kp(i)+this%kp(i+1))
				kps=2*this%kp(i)*this%kp(i-1)/(this%kp(i)+this%kp(i-1))
				an=Qpa*this%shcp(i)+areap/this%Length(i)*kpn
				as=areas/this%Length(i)*kps
			  	aw=this%htp(i)*PI*2*(this%Rtube+this%thickt)*this%Length(i)
				ae=this%hvp(i)*PI*2*this%Rpa*this%Length(i)
				apz=areap*this%length(i)*this%rhop(i)*this%shcp(i)/dt
				ap=an+as+aw+ae+apz				
				S=apz*Tpl(i)
				A(i,i+1)=-an
				A(i,i-1)=-as
			  endif		
			  A(i,i)=ap
			  A(i,i+2*N)=-aw!tube
			  A(i,i+3*N)=-ae!shell
			  B(i)=S
			  case(2)!secondary
			  if(i==1.or.i==N) then
				ae=this%hts(i)*PI*2*this%Rtube*this%Length(i)
				apz=areas*this%length(i)*this%rhos(i)*this%shcs(i)/dt
				ap=Qsa*this%shcs(i)+ae+apz
				if(i==1)then
					as=0
					S=apz*Tsl(i)+Qsa*this%shcs(i)*this%Tsin
				else
					as=Qsa*this%shcs(i)
					S=apz*Tsl(i)
					A(i+N,i-1+N)=-as
				endif
			  else
			  	ksn=2*this%ks(i)*this%ks(i+1)/(this%ks(i)+this%ks(i+1))
				kss=2*this%ks(i)*this%ks(i-1)/(this%ks(i)+this%ks(i-1))
			  	an=areas/this%length(i)*ksn
				as=Qsa*this%shcs(i)+areas/this%length(i)*kss
				ae=this%hts(i)*PI*2*this%Rtube*this%Length(i)
				apz=areas*this%length(i)*this%rhos(i)*this%shcs(i)/dt
				ap=an+as+ae+apz
				S=apz*Tsl(i)
				A(i+N,i+1+N)=-an
				A(i+N,i-1+N)=-as
			  endif
              A(i+N,i+N)=ap			  
			  A(i+N,i+2*N)=-ae!tube
			  B(i+N)=S
			  
			  case(3)!tube
			  aw=this%hts(i)*PI*2*this%Rtube*this%Length(i)
			  ae=this%htp(i)*PI*2*(this%Rtube+this%thickt)*this%Length(i)
			  apz=mt*this%shct/dt
			  ap=aw+ae+apz
			  S=apz*Ttl(i)
			  A(i+2*N,i+2*N)=ap
			  A(i+2*N,i+N)=-aw
			  A(i+2*N,i)=-ae
		      B(i+2*N)=S
			  case(4)!shell
			  aw=this%hvp(i)*PI*2*this%Rpa*this%Length(i)
			  apz=mv*this%shcv/dt
			  ap=aw+apz
			  S=apz*Tvl(i)
			  A(i+3*N,i+3*N)=ap
			  A(i+3*N,i)=-aw
			  B(i+3*N)=S
			end select
			end do
        enddo
        
		call self_DGESV(A,B,info)
		
		do i=1,4*N,1
			if(i<=N)then
				this%Tp(i)=B(i)
			elseif(i<=2*N)then
				this%Ts(i-N)=B(i)
			elseif(i<=3*N)then
				this%Tt(i-2*N)=B(i)
			else
				this%Tv(i-3*N)=B(i)
			endif
        end do	
        this%Tpout=this%Tp(1)
        this%Tsout=this%Ts(N)
	end subroutine cal_thermal_kerel_transient
	
	subroutine update_Boundary(this,Qp,Qs,Tpin,Tsin)
		implicit none
		class(IHX),intent(in out)::this
		real(KREAL),intent(in)::Qp,Qs,Tpin,Tsin
		this%Qp=Qp
		this%Qs=Qs
		this%Tpin=Tpin
		this%Tsin=Tsin		
	end subroutine update_Boundary
    
    function cal_buoy(this) result (buoy)
        class(IHX),intent(in out)::this
        real(KREAL)::buoy,dbuoy
        integer::i
        buoy=0.0
        dbuoy=0.0
        associate(Ny=>this%N,&
                  rhop=>this%rhop,&
                  length=>this%Length)
            do i=1,Ny,1
               dbuoy=-GRAVG*sin(-90.0/180.0*PI)*rhop(i)*length(i)
               buoy=buoy+dbuoy
            enddo
        end associate
    end function cal_buoy
    
    function cal_fric(this) result(fric)
        implicit none
		class(IHX),intent(in out)::this
        real(KREAL)::fric
        real(KREAL)::Re,visa,De
        integer::i
        associate(Dep=>this%Dep,&
                  Ny=>this%N,&
                  Ntube=>this%Ntube,&
                  vis=>this%visp,&
                  flowrate=>this%Qp,&
                  flowarea=>this%areap)
            visa=0.0
            do i=1,Ny,1
                visa=visa+vis(i)/Ny
            enddo
            De=Dep*Ntube
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
		class(IHX),intent(in out)::this
        real(KREAL)::beta,rhoa,De,area
        integer::i
        beta=0.0
        rhoa=0.0
        associate(Ny=>this%N,&
                  Ntube=>this%Ntube,&
                  fric=>this%Fricp,&
                  rho=>this%rhop,&
                  ltotal=>this%lsingle,&
                  Des=>this%Des,&
                  areap=>this%areap,&
                  K=>this%Kricp)
            De=Des*Ntube
            area=areap
            fric=this%cfric()
            do i=1,Ny,1
                rhoa=rhoa+rho(i)/Ny
            enddo
            beta=-fric*ltotal/(2*De*rhoa*area**2)-K/(2*rhoa*area**2)
        end associate
    end function cal_beta
    

end module Imp_IHX_header