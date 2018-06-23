module imp_property
    use constants
    contains
    !=============================================================
    !LBE
    !=============================================================
    function get_density_LBE(tin) result(density)
     real(KREAL),intent(in)::tin
     real(KREAL)::density
     density=11096-1.3326*tin
    end function get_density_LBE
	
	function get_shc_LBE(Tin) result(shc)
	     real(KREAL),intent(in)::tin
		 real(KREAL)::shc
		 shc=159.0-0.0272*Tin+7.12*10.**(-6)*Tin**2!J/(kg*K)
	end function get_shc_LBE
	
	function get_conductivity_LBE(Tin) result(conductivity)
	     real(KREAL),intent(in)::tin
		 real(KREAL)::conductivity
		 conductivity=3.61+1.517*10.**(-2)*Tin-1.741*10.**(-6)*Tin**2
	end function get_conductivity_LBE
	
	function get_vis_LBE() result(vis)
        real(KREAL)::vis,Tin
		Tin=500.0!K
		vis=3.26*10.**(-3)-6.26*10.**(-6)*Tin+4.63*10.**(-9)*Tin**2
	end function get_vis_LBE
    !=============================================================
    !WATER
	!<传热与传质>P579
	!饱和水在500K时候的热物性
    !=============================================================
	function get_density_water(tin) result(density)
     real(KREAL),intent(in)::tin
     real(KREAL)::density
     density=1000.0!kg/m3
    end function get_density_water
	
	function get_shc_water(Tin) result(shc)
	     real(KREAL),intent(in)::tin
         real(KREAL)::shc
		 shc=4660.0!J/(kg*K)
	end function get_shc_water
	
	function get_conductivity_water(Tin) result(conductivity)
	     real(KREAL),intent(in)::tin
         real(KREAL)::conductivity
		 conductivity=0.642!W/(m*K)
	end function get_conductivity_water
	
	function get_vis_water() result(vis)
         real(KREAL)::vis
		 vis=118.0*10.**(-6)!Pa*s
	end function get_vis_water
    !=============================================================
    ! liquid Na
    !=============================================================
    function get_density_Na(Tin) result(rho)
	    real(KREAL),intent(in)::tin
		real(KREAL)::rho
        !local
        real(KREAL),allocatable::dtem(:)
        real(KREAL),allocatable::drho(:)       
        integer::Nrho
        Nrho=23
        allocate(dtem(Nrho))
        allocate(drho(Nrho))
        dtem=[400.,500.,600.,700.,800.,900.,1000.,1100.,1200.,1300.,1400.,1500.,1600.,1700.,1800.,1900.,2000.,2100.,2200.,2300.,2400.,2500.,2503.7]
            drho=[919.,897.,874.,852.,828.,805.,781.,756.,732.,706.,680.,653.,626.,597.,568.,537.,504.,469.,431.,387.,335.,239.,219.]
        ! if(Tin<dtem(1).or.Tin>dtem(Nrho)) write(*,fmt="('Tin is out of range of Na rho table')")
        if(Tin<dtem(1)) then
            rho=(drho(2)-drho(1))/(dtem(2)-dtem(1))*(Tin-dtem(1))+drho(1)
        elseif(Tin>dtem(Nrho)) then
            rho=(drho(Nrho)-drho(Nrho-1))/(dtem(Nrho)-dtem(Nrho-1))*(Tin-dtem(Nrho))+drho(Nrho)
        else
            do i=1,Nrho-1,1
                if(Tin>=dtem(i).and.Tin<=dtem(i+1)) then
                    rho=(drho(i+1)-drho(i))/(dtem(i+1)-dtem(i))*(Tin-dtem(i))+drho(i)
                    exit
                endif
            enddo  
        endif
        !Cui Manman
         ! rho=950.076-0.2298*Tin-1.4605*1e-5*Tin**2+5.5379*1e-9*Tin**3!kg/m3
	end function get_density_Na
    
    function get_shc_Na(Tin) result(shc)
	    real(KREAL),intent(in)::tin
		real(KREAL)::shc
        !local
        integer::Nar
        real(KREAL),allocatable::argonne(:,:)
        real(KREAL),allocatable::dtem(:)
        real(KREAL),allocatable::dar(:)
        Nar=24
        allocate(argonne(2,Nar))
        allocate(dtem(Nar))
        allocate(dar(Nar))
        argonne=reshape([371.,1383,400.,1372.,500.,1334.,600.,1301.,700.,1277.,800.,1260.,900.,1252.,&
                       &1000.,1252.,1100.,1261.,1200.,1279.,1300.,1305.,1400.,1340.,1500.,1384.,1600.,1437.,&
                       &1700.,1500.,1800.,1574.,1900.,1661.,2000.,1764.,2100.,1926.,2200.,2190.,2300.,2690.,&
                       &2400.,4012.,2469.,8274.,2500.,39279],[2,24])
        dtem=argonne(1,:)
        dar=argonne(2,:)
        ! if(Tin<dtem(1).or.Tin>dtem(Nar)) write(*,fmt="('Tin is out of range of Na shc table')")
        
        if(Tin<dtem(1)) then
            shc=(dar(2)-dar(1))/(dtem(2)-dtem(1))*(Tin-dtem(1))+dar(1)
        elseif(Tin>dtem(Nar)) then
            shc=(dar(Nar)-dar(Nar-1))/(dtem(Nar)-dtem(Nar-1))*(Tin-dtem(Nar))+dar(Nar)
        else
            do i=1,Nar-1,1
                if(Tin>=dtem(i).and.Tin<=dtem(i+1)) then
                    shc=(dar(i+1)-dar(i))/(dtem(i+1)-dtem(i))*(Tin-dtem(i))+dar(i)
                    exit
                endif
            enddo  
        endif
        ! do i=1,Nar-1,1
           ! if(Tin>=dtem(i).and.Tin<=dtem(i+1)) then
               ! shc=(dar(i+1)-dar(i))/(dtem(i+1)-dtem(i))*(Tin-dtem(i))+dar(i)
               ! exit
           ! endif
           ! if(Tin>dtem(Nar)) shc=dar(Nar)
        ! enddo
        !Cui Manman
		! shc=1436.05+(4.625*1e-5*Tin-0.5802)*Tin!J/(kg*K)
	end function get_shc_Na
    
    function get_conductivity_Na(Tin) result(conductivity)
	    real(KREAL),intent(in)::tin
		real(KREAL)::conductivity
		! conductivity=92.25-0.058*Tin+1.17*1e-5*Tin**2!W/(m*K)
		!Argonne
        conductivity=124.67-0.1138*Tin+5.5226*1e-5*Tin**2-1.1842*1e-8*Tin**3!W/(m*K)
	end function get_conductivity_Na
    
    function get_vis_Na(Tin,rho) result(vis)
	    real(KREAL),intent(in)::tin
        real(KREAL),intent(in)::rho
		real(KREAL)::vis
		if(Tin<=773) then
           vis=0.1235*1e-4*rho**(1./3.)*exp(0.697*rho/Tin)!Pa*s
        else
           vis=0.0851*1e-4*rho**(1./3.)*exp(1.04*rho/Tin)!Pa*s
        endif
	end function get_vis_Na
    !================================================================
    !U5Fs
    !================================================================
    function get_density_U5Fs(Tin) result(rho)
	    real(KREAL),intent(in)::Tin
		real(KREAL)::rho
		rho=18200.0!kg/m3
	end function get_density_U5Fs
    
    function get_shc_U5Fs(Tin) result(shc)
	    real(KREAL),intent(in)::Tin
		real(KREAL)::shc!J/(Kg*K)
        if(273<Tin<833) then
            shc=139.61-1978.5*1e-5*Tin+3156.6*1e-7*Tin**2
        elseif(833<Tin<913) then
            shc=463.64-2615.0*1e-4*Tin
        else
            shc=2401.26-3512.6*1e-3*Tin+1401.2*1e-6*Tin**2
        endif
	end function get_shc_U5Fs
    
    function get_conductivity_U5Fs(Tin) result(conductivity)
	    real(KREAL),intent(in)::Tin
		real(KREAL)::conductivity
		conductivity=14.1+2.98*1e-2*Tin+3.01*1e-6*Tin**2!W/(m*K)
	end function get_conductivity_U5Fs
    !================================================================
    !316L
    !================================================================
    function get_density_316L(Tin) result(rho)
		real(KREAL),intent(in)::Tin
        real(KREAL)::rho
		rho=7980.!kg/m3
	end function get_density_316L
    
    function get_shc_316L(Tin) result(shc)
		real(KREAL),intent(in)::Tin
        real(KREAL)::shc!J/(Kg*K)
        shc=502.
	end function get_shc_316L
    
    function get_conductivity_316L(Tin) result(conductivity)
		real(KREAL),intent(in)::Tin
        real(KREAL)::conductivity
		conductivity=20.9!W/(m*K)
	end function get_conductivity_316L
    !================================================================
    !304
    !================================================================
    function get_density_304(Tin) result(rho)
		real(KREAL),intent(in)::Tin
        real(KREAL)::rho
		rho=7930.!kg/m3
	end function get_density_304
    
    function get_shc_304(Tin) result(shc)
		real(KREAL),intent(in)::Tin
        real(KREAL)::shc!J/(Kg*K)
        shc=500.
	end function get_shc_304
    
    function get_conductivity_304(Tin) result(conductivity)
		real(KREAL),intent(in)::Tin
        real(KREAL)::conductivity
		conductivity=21.5!W/(m*K)
	end function get_conductivity_304
    !================================================================    
    !density
    function get_density(num,Tin) result(rho)
        integer,intent(in)::num
        real(KREAL),intent(in),optional::Tin
        real(KREAL)::rho
        select case(num)
            case(101)
            rho=get_density_LBE(Tin)
            case(102)
            rho=get_density_water(Tin)
            case(103)
            rho=get_density_Na(Tin)
            case(301)
            rho=get_density_U5Fs(Tin)
            case(202)
            rho=get_density_316L(Tin)
            case(203)
            rho=get_density_304(Tin)
        end select
    end function get_density
    !shc
    function get_shc(num,Tin) result(shc)
        integer,intent(in)::num
        real(KREAL),intent(in)::Tin
        real(KREAL)::shc
        select case(num)
            case(101)
            shc=get_shc_LBE(Tin)
            case(102)
            shc=get_shc_water(Tin)
            case(103)
            shc=get_shc_Na(Tin)
            case(301)
            shc=get_shc_U5Fs(Tin)
            case(202)
            shc=get_shc_316L(Tin)
        end select
    end function get_shc
    !conductivity
    function get_conductivity(num,Tin) result(conductivity)
        integer,intent(in)::num
        real(KREAL),intent(in)::Tin
        real(KREAL)::conductivity
        select case(num)
            case(101)
            conductivity=get_conductivity_LBE(Tin)
            case(102)
            conductivity=get_conductivity_water(Tin)
            case(103)
            conductivity=get_conductivity_Na(Tin)
            case(301)
            conductivity=get_conductivity_U5Fs(Tin)
            case(202)
            conductivity=get_conductivity_316L(Tin)
            case(203)
            conductivity=get_conductivity_304(Tin)
        end select      
    end function get_conductivity
    !vis
    function get_vis(num,Tin,rho) result(vis)
        integer,intent(in)::num
        real(KREAL),intent(in),optional::Tin
        real(KREAL),intent(in),optional::rho
        real(KREAL)::vis
        select case(num)
            case(101)
            vis=get_vis_LBE()
            case(102)
            vis=get_vis_water()
            case(103)
            if (PRESENT(Tin).and.PRESENT(rho)) then
                vis=get_vis_Na(Tin,rho)
            else
                print*,'vis_Na need Tin and rho'
            endif
        end select        
    end function get_vis
    !================================================================
end module imp_property
