	module Imp_driving_THcore
		use constants
		use imp_property
		use Imp_assm_global
		use imp_single_channel
		use imp_loop_global
		implicit none
		contains
		subroutine driving_TH_core(transient_flag,Qin,Tin,assembly,Tout,last,current)
			logical,intent(in)::transient_flag
			real(KREAL),intent(in)::assembly(:,:)!power(zone,layer) W
			real(KREAL),intent(in)::Qin,Tin
			real(KREAL),intent(out)::Tout	!Tave of core
			real(KREAL),intent(in),optional::last,current
			!local
			integer  :: i,j,k 
            integer  :: nr,na,M,N
            integer  :: Nzone,izone
			real(KREAL),allocatable::pow(:,:),fq_core(:,:)
			real(KREAL)::density,flowrate
            
			fq_core=1.0D0
            Nzone=core%Nflow+core%Nflowsemi!需要计算的zone的数量
            
			nr = SIZE(assembly, dim=1)!zone                           
			na = SIZE(assembly, dim=2)!layer                          
			M=size(assm1(1)%thermal%temperature,dim=1)
			N=size(assm1(1)%thermal%temperature,dim=2)		
            allocate(pow(na,N),fq_core(na,N))
			pow=0.0
			fq_core=1.0
		
			do i=1,Nzone,1
				assm1(core%fzone(i))%th_boundary%T%inlet=Tin
			enddo
			
			call driving_loop_flowAlloc(assm1,Qin)

			do i=1,Nzone,1!zone start
                izone=core%fzone(i)
				do j=1,assm1(izone)%mesh%ny,1!dy
					do k=1,N,1
						if(k<=assm1(izone)%mesh%Nf) pow(j,k)=assembly(i,j+assm1(izone)%mesh%layer_bottom)/(assm1(izone)%geom%N_fuelpin*assm1(izone)%geom%height(j)*3.14159*assm1(izone)%geom%pellet**2)
					enddo
				enddo

				! if (assm1(izone)%th_boundary%u%inlet==0.0) then
					! assm1(izone)%thermal%velocity=0.0
					! assm1(izone)%th_boundary%u%outlet=0.0
					! assm1(izone)%thermal%temperature=assm1(izone)%th_boundary%T%inlet
					! assm1(izone)%th_boundary%T%outlet=assm1(izone)%th_boundary%T%inlet
					! assm1(izone)%property%rho=get_density(assm1(izone)%th_boundary%T%inlet)
				! else
                if (transient_flag)then
                    if(PRESENT(last).and.PRESENT(current)) call driving_imp_THtransient(assm1(izone),pow,fq_core,last,current)
                else
                    call driving_imp_THsteady(assm1(izone),pow,fq_core)					
                endif
				! endif	
			enddo !zone		
			!Tout volum ave
			Tout=0.0
			do i=1,Nzone,1
				density=get_density(assm1(izone)%th_boundary%T%inlet)
				flowrate=assm1(izone)%th_boundary%u%inlet*(assm1(izone)%geom%n_pin*assm1(izone)%hydrau%aflow*density)
				Tout=Tout+assm1(izone)%th_boundary%T%outlet*flowrate/Qin	
			enddo
			core%Tfout=Tout
		end subroutine driving_TH_core
	end module Imp_driving_THcore