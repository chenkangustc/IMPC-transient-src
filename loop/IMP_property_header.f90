module imp_property
    use constants
    contains
    !=============================================================
    !LBE
    !=============================================================
    function get_density(tin) result(density)
     real(KREAL),intent(in)::tin
     real(KREAL)::density
     density=11096-1.3326*tin
    end function get_density
	
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
		rho=950.076-0.2298*Tin-1.4605*1e-5*Tin**2+5.5379*1e-9*Tin**3!kg/m3
	end function get_density_Na
    
    function get_shc_Na(Tin) result(shc)
	    real(KREAL),intent(in)::tin
		real(KREAL)::shc
		shc=1436.05+(4.625*1e-5*Tin-0.5802)*Tin!J/(kg*K)
	end function get_shc_Na
    
    function get_conductivity_Na(Tin) result(conductivity)
	    real(KREAL),intent(in)::tin
		real(KREAL)::conductivity
		conductivity=92.25-0.058*Tin+1.17*1e-5*Tin**2!W/(m*K)
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
    function get_density_316L() result(rho)
		real(KREAL)::rho
		rho=7980.!kg/m3
	end function get_density_316L
    
    function get_shc_316L() result(shc)
		real(KREAL)::shc!J/(Kg*K)
        shc=502.
	end function get_shc_316L
    
    function get_conductivity_316L() result(conductivity)
		real(KREAL)::conductivity
		conductivity=20.9!W/(m*K)
	end function get_conductivity_316L
    !================================================================
    !304
    !================================================================
    function get_density_304() result(rho)
		real(KREAL)::rho
		rho=7930.!kg/m3
	end function get_density_304
    
    function get_shc_304() result(shc)
		real(KREAL)::shc!J/(Kg*K)
        shc=500.
	end function get_shc_304
    
    function get_conductivity_304() result(conductivity)
		real(KREAL)::conductivity
		conductivity=21.5!W/(m*K)
	end function get_conductivity_304       

end module imp_property
