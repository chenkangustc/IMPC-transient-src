	module Imp_driving_THcore
		use constants
		use imp_property
		use Imp_assm_global
		use imp_single_channel
		implicit none
		contains
		subroutine driving_THcore_steady(Qin,Tin,assembly,Tout)
			real(KREAL),intent(in)::assembly(:,:)!power(zone,layer) W
			real(KREAL),intent(in)::Qin,Tin
			real(KREAL),intent(out)::Tout	!Tave of core
			!local
			integer  :: i,j,k 
            integer  :: nr,na,M,N
			real(KREAL),allocatable::pow(:,:),fq_core(:,:)
			real(KREAL)::density,flowrate

			fq_core=1.0D0		
			nr = SIZE(assembly, dim=1)!zone                           
			na = SIZE(assembly, dim=2)!layer                          
			M=size(assm1(1)%thermal%temperature,dim=1)
			N=size(assm1(1)%thermal%temperature,dim=2)		
            allocate(pow(M,N),fq_core(M,N))
			pow=0.0
			fq_core=1.0
		
			do i=1,nr,1
				assm1(i)%th_boundary%T%inlet=Tin
			enddo

			call driving_loop_flowAlloc(assm1,Qin)

			do i=1,nr,1!zone start
				do j=1,assm1(i)%mesh%ny,1!dy
					do k=1,N,1
						if(k<=assm1(i)%mesh%Nf) pow(j,k)=assembly(i,j+assm1(i)%mesh%layer_bottom)/(assm1(i)%geom%N_fuelpin*assm1(i)%geom%height(j)*3.14159*assm1(i)%geom%pellet**2)
					enddo
				enddo
			enddo
			
			do i=1,nr,1
				if (assm1(i)%th_boundary%u%inlet==0.0) then
					assm1(i)%thermal%velocity=0.0
					assm1(i)%th_boundary%u%outlet=0.0
					assm1(i)%thermal%temperature=assm1(i)%th_boundary%T%inlet
					assm1(i)%th_boundary%T%outlet=assm1(i)%th_boundary%T%inlet
					assm1(i)%property%rho=get_density(assm1(i)%th_boundary%T%inlet)
				else
						call driving_imp_steady(assm1(i),pow,fq_core)
				endif	
			enddo !zone		
			!Tout volum ave
			Tout=0.0
			do i=1,nr,1
				density=get_density(assm1(i)%th_boundary%T%inlet)
				flowrate=assm1(i)%th_boundary%u%inlet*(assm1(i)%geom%n_pin*assm1(i)%hydrau%aflow*density)
				Tout=Tout+assm1(i)%th_boundary%T%outlet*flowrate/Qin	
			enddo
		end subroutine driving_THcore_steady
	end module Imp_driving_THcore