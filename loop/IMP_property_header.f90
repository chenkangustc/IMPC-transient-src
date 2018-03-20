module imp_property
    use constants
    contains
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
	!<传热与传质>P579
	!饱和水在500K时候的热物性
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
end module imp_property
