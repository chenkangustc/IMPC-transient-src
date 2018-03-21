	module driving_THcore
		use contants
		use Imp_assm_global
		use imp_single_channel
		implicit none
		contains
		subroutine driving_THcore_steady(Qin,Tin,assembly)
			real(KREAL),intent(in)::assembly(:,:)!power(zone,layer) W
			real(KREAL),intent(in)::Qin,Tin
			!local
			integer  :: i,j,k 
			real(KREAL),allocatable::pow(:,:),fq_core(:,:)

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
						if(k<=assm1(i)%mesh%Nf) power(j,k)=assembly(i,j+assm1(i)%mesh%layer_bottom)/(assm1(i)%geom%N_fuelpin*assm1(i)%geom%height(j)*3.14159*assm1(i)%geom%pellet**2)
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
					if (transient_flag)  then
						call driving_imp_transient(assm1(i),power, fq_core,last_, current_)
					else
						call driving_imp_steady(assm1(i),power,fq_core)
					end if
				endif	
			enddo !zone		
		end subroutine driving_THcore_steady
	end module driving_THcore