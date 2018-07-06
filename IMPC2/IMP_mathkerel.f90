module imp_mathkerel
    use constants
    implicit none
    private
    public::tdma
    public::get_convection
    public::get_hyconstant
    public::get_nusselt_tube
    public::get_nusselt_bundle
    public::get_fric_IHX
    public::get_fric_pin
    public::get_fric_pipe
    
contains
    subroutine tdma(N,A,B,u)
      real(KREAL):: A(N,N)
      real(KREAL):: B(N)
      real(KREAL):: u(N)
      integer i,N
      real(KREAL),dimension(N)::aa,bb,cc,dd,x,y

      do i=1,N,1
          if(i==1)then
          aa(1)=0.0
          bb(1)=A(1,1)
          cc(1)=A(1,2)
          dd(1)=B(1)
          elseif (i==N) then
           aa(N)=A(N,N-1)
           bb(N)=A(N,N)
           cc(N)=0.0
           dd(N)=B(N)
          else
           aa(i)=A(i,i-1)
           bb(i)=A(i,i)
           cc(i)=A(i,i+1)
           dd(i)=B(i)
          endif
      enddo
      
      x(1)=cc(1)/bb(1)
      y(1)=dd(1)/bb(1)
      
      do i=2,N,1
          x(i)=cc(i)/(bb(i)-x(i-1)*aa(i))
          y(i)=(dd(i)-y(i-1)*aa(i))/(bb(i)-x(i-1)*aa(i))
      enddo
      
      u(N)=y(N)
      
      do i=N-1,1,-1
          u(i)=y(i)-x(i)*u(i+1)
      enddo
    endsubroutine tdma    
    
    subroutine get_convection(lenth,flow_area,wetted_perimeter,density,velocity,viscosity,capacity,conductivity,convection) !lenth是特征长度,取水力学直径  Nu=(h*x)/k
     real(KREAL):: lenth
     real(KREAL):: flow_area
     real(KREAL):: wetted_perimeter
     real(KREAL):: density
     real(KREAL):: velocity
     real(KREAL):: De
     real(KREAL):: viscosity,capacity,convection,conductivity
     real(KREAL):: Pr,Re,Pe,Nu
     
     call get_Nusselt(flow_area,wetted_perimeter,density,velocity,viscosity,capacity,conductivity,Nu)  
     !Nu=get_Nusselt_Na_bundle(flow_area,wetted_perimeter,density,velocity,viscosity,capacity,conductivity)
     convection=Nu*conductivity/lenth
    end subroutine get_convection
    
    subroutine get_Nusselt(flow_area,wetted_perimeter,density,velocity,viscosity,capacity,conductivity,Nu)!液态重金属的努赛尔数
     real(KREAL):: flow_area
     real(KREAL):: wetted_perimeter    
     real(KREAL):: density
     real(KREAL):: velocity
     real(KREAL):: De
     real(KREAL):: Re 
     real(KREAL):: viscosity
     real(KREAL):: capacity
     real(KREAL):: conductivity
     real(KREAL):: Pr
     real(KREAL):: Pe
     real(KREAL):: Nu
     
     De=4*flow_area/wetted_perimeter
     Re=4*density*velocity*De/viscosity
     Pr=viscosity*capacity/conductivity
     Pe=Re*Pr
     Nu=5.0+2.5D-2*Pe**0.8
    end subroutine get_Nusselt
	!========================================================================
	! water
    !========================================================================
	function get_Nusselt_water_tube(Nutype,flowarea,wet,De,rho,flowrate,vis,shc,conductivity) result(Nu)
		integer,intent(in)::Nutype
		real(KREAL),intent(in)::flowarea,wet,De
		real(KREAL),intent(in)::rho,vis,shc,conductivity
		real(KREAL),intent(in)::flowrate
		!local
		real(KREAL)::Re,Pr,Pe,Nu
		Re=flowrate*De/(vis*flowarea)
		Pr=vis*shc/conductivity
		Pe=Re*Pr
        select case(Nutype)
        case(1)
		Nu=0.023*Re**0.8*Pr**0.333
        end select
	end function get_Nusselt_water_tube
	!========================================================================
	! LBE
    !========================================================================
	function get_Nusselt_LBE_bundle(Nutype,pd,flowarea,wet,De,rho,flowrate,vis,shc,conductivity) result(Nu)
		integer,intent(in)::Nutype
		real(KREAL),intent(in)::pd
		real(KREAL),intent(in)::flowarea,wet,De
		real(KREAL),intent(in)::rho,vis,shc,conductivity
		real(KREAL),intent(in)::flowrate
		!local
		real(KREAL)::Re,Pr,Pe,Nu
		Re=flowrate*De/(vis*flowarea)
		Pr=vis*shc/conductivity
		Pe=Re*Pr
        select case(Nutype)
        case(1)
		Nu=7.55*pd-20.*(pd)**(-13)+3.67/(90.*(pd)**2)*Pe**(0.56+0.19*pd)
        if (Nu<0) print*,'Nu is .lt. 0.0,the Nu equation cannot be used here'
        end select
    end function get_Nusselt_LBE_bundle

	function get_Nusselt_LBE_tube(Nutype,flowarea,wet,De,rho,flowrate,vis,shc,conductivity) result(Nu)
		integer,intent(in)::Nutype
		real(KREAL),intent(in)::flowarea,wet,De
		real(KREAL),intent(in)::rho,vis,shc,conductivity
		real(KREAL),intent(in)::flowrate
		!local
		real(KREAL)::Re,Pr,Pe,Nu
		Re=flowrate*De/(vis*flowarea)
		Pr=vis*shc/conductivity
		Pe=Re*Pr
        select case(Nutype)
        case(1)
		Nu=4.8+0.0156*Re**0.85*Pr**0.93
        end select
	end function get_Nusselt_LBE_tube
    
    function get_fric_LBE_IHX(Frtype,Re) result(fric)
        integer,intent(in)::Frtype
        real(KREAL),intent(in)::Re
        real(KREAL)::fric
        select case(Frtype)
        case(1)
        fric=0.3164/Re**0.25
        end select
    end function
    
    function get_fric_LBE_pin(Frtype,Re) result(fric)
        integer,intent(in)::Frtype
        real(KREAL),intent(in)::Re
        real(KREAL)::fric
        select case(Frtype)
        case(1)
        fric=70.399/Re+20.722/Re**0.133+0.01!吕P41
        case(2)
        fric=0.3164/Re**0.25!chen zhao
        case(3)
        fric=0.55/Re**0.25
        end select
    end function
    
    function get_fric_LBE_pipe(Frtype,Re) result(fric)
        integer,intent(in)::Frtype
        real(KREAL),intent(in)::Re
        real(KREAL)::fric
        select case(Frtype)
        case(1)
        fric=0.3164/Re**0.25
        end select  
	end function
    !========================================================================
	!       Na
    !========================================================================
    function get_Nusselt_Na_tube(Nutype,flowarea,wet,De,rho,flowrate,vis,shc,conductivity) result(Nu)
		integer,intent(in)::Nutype
		real(KREAL),intent(in)::flowarea,wet,De
		real(KREAL),intent(in)::rho,vis,shc,conductivity
		real(KREAL),intent(in)::flowrate
		!local
		real(KREAL)::Re,Pr,Pe,Nu
		Re=flowrate*De/(vis*flowarea)
		Pr=vis*shc/conductivity
		Pe=Re*Pr
        select case(Nutype)
        case(1)
		Nu=4.8+0.025*Pe**0.8!Argonne
        !Nu=4.5+0.018*Pe**0.8
        end select
	end function get_Nusselt_Na_tube
    
    function get_Nusselt_Na_bundle(Nutype,pd,flowarea,wet,De,rho,flowrate,vis,shc,conductivity) result(Nu)
		integer,intent(in)::Nutype
		real(KREAL),intent(in)::pd
		real(KREAL),intent(in)::flowarea,wet,De
		real(KREAL),intent(in)::rho,vis,shc,conductivity
		real(KREAL),intent(in)::flowrate
        !local
		real(KREAL)::Re,Pr,Pe,Nu
		Re=flowrate*De/(vis*flowarea)
		Pr=vis*shc/conductivity
		Pe=Re*Pr
        select case(Nutype)
        case(1)
        !Argonne
        Nu=5.0+0.025*Pe**0.8
        end select
        !West formula
        ! if(pd>=1.05.and.pd<=1.15) then
            ! if(Pe<=150.)then
                ! Nu=4.496*(-16.15+24.96*pd-8.55*pd**2)
            ! elseif(Pe>=150.0.and.Pe<=1000.)then
                ! Nu=(-16.15+24.96*pd-8.55*pd**2)*Pe**0.3
            ! else
                ! print*,'Pe is out of range of Nu'
            ! endif
        ! elseif(pd>=1.15.and.pd<=1.30)then
            ! Nu=4.0+0.16*pd**5.+0.33*pd**3.8*(Pe/100.)**0.86
        ! else
            ! print*,'PD is out of range of Nu'
        ! endif
	end function get_Nusselt_Na_bundle
    
    function get_fric_Na_IHX(Frtype,De,Re) result(fric)
        integer,intent(in)::Frtype
        real(KREAL),intent(in)::De,Re
        real(KREAL)::fric
        select case(Frtype)
        case(1)
        if(Re>=1082) then
            fric=0.0055*(1.+(20000*1e-5/De+1e6/Re)**(1./3.))
        else
            fric=64./Re
        endif
        end select
    end function
    
    function get_fric_Na_pin(Frtype,Re) result(fric)
        integer,intent(in)::Frtype
        real(KREAL),intent(in)::Re
        real(KREAL)::fric
        select case(Frtype)
        case(1)
        if(Re>=2050.) then
            fric=0.1875/Re**2
        else
            fric=76.5/Re
        endif
        end select
    end function
    
    function get_fric_Na_pipe(Frtype,De,Re) result(fric)
        integer,intent(in)::Frtype
        real(KREAL),intent(in)::De,Re
        real(KREAL)::fric
        select case(Frtype)
        case(1)
        if(Re>=1082) then
            fric=0.0055*(1.+(20000*1e-5/De+1e6/Re)**(1./3.))
        else
            fric=64./Re
        endif
        end select
	end function
    !========================================================================
	!       total
    !========================================================================
    function get_Nusselt_tube(Ftype,Nutype,flowarea,wet,De,rho,flowrate,vis,shc,conductivity) result(Nu)
		integer,intent(in)::Ftype,Nutype
		real(KREAL),intent(in)::flowarea,wet,De
		real(KREAL),intent(in)::rho,vis,shc,conductivity
		real(KREAL),intent(in)::flowrate
        !
        real(KREAL)::Nu
        select case(Ftype)
        case(101)
        Nu=get_Nusselt_LBE_tube(Nutype,flowarea,wet,De,rho,flowrate,vis,shc,conductivity) 
        case(102)!Na
        Nu=get_Nusselt_Na_tube(Nutype,flowarea,wet,De,rho,flowrate,vis,shc,conductivity) 
        case(103)
        Nu=get_Nusselt_water_tube(Nutype,flowarea,wet,De,rho,flowrate,vis,shc,conductivity)
        end select
    end function
    
    function get_Nusselt_bundle(Ftype,Nutype,pd,flowarea,wet,De,rho,flowrate,vis,shc,conductivity) result(Nu)
		integer,intent(in)::Ftype,Nutype
		real(KREAL),intent(in)::pd
		real(KREAL),intent(in)::flowarea,wet,De
		real(KREAL),intent(in)::rho,vis,shc,conductivity
		real(KREAL),intent(in)::flowrate
        !
        real(KREAL)::Nu
        select case(Ftype)
        case(101)
        Nu=get_Nusselt_LBE_bundle(Nutype,pd,flowarea,wet,De,rho,flowrate,vis,shc,conductivity)
        case(102)!Na
        Nu=get_Nusselt_Na_bundle(Nutype,pd,flowarea,wet,De,rho,flowrate,vis,shc,conductivity)
        end select
    end function
    
    function get_fric_IHX(Ftype,Frtype,De,Re) result(fric)
        integer,intent(in)::Ftype,Frtype
        real(KREAL),intent(in)::De,Re
        real(KREAL)::fric
        select case(Ftype)
        case(101)
        fric=get_fric_LBE_IHX(Frtype,Re)
        case(102)!Na
        fric=get_fric_Na_IHX(Frtype,De,Re)
        end select
    end function
    
    function get_fric_pin(Ftype,Frtype,Re) result(fric)
        integer,intent(in)::Ftype,Frtype
        real(KREAL),intent(in)::Re
        real(KREAL)::fric
        select case(Ftype)
        case(101)
        fric=get_fric_LBE_pin(Frtype,Re)
        case(102)!Na
        fric=get_fric_Na_pin(Frtype,Re)
        end select
    end function
    
    function get_fric_pipe(Ftype,Frtype,De,Re) result(fric)
        integer,intent(in)::Ftype,Frtype
        real(KREAL),intent(in)::De,Re
        real(KREAL)::fric
        select case(Ftype)
        case(101)!LBE
        fric=get_fric_LBE_pipe(Frtype,Re)
        case(102)!Na
        fric=get_fric_Na_pipe(Frtype,De,Re)
        end select
    end function
    !
    subroutine get_hyconstant(rc,pd,Aflow,wet,de)
       real(KREAL):: rc,p,pd !r是包壳外半径 p是对边距
       real(KREAL):: Aflow,Ashell,Atotal,wet,de
       
           p=pd*2*rc
           Ashell=3.14*rc*rc
           Atotal=0.5*sqrt(3.0)*p*p
           Aflow=Atotal-Ashell
           wet=2*3.14*rc+p/sqrt(3.0)*6
           de=4*Aflow/wet
    endsubroutine get_hyconstant
end module imp_mathkerel
