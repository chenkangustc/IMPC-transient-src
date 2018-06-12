	module Imp_driving_THcore
		use constants
		use imp_property
		use Imp_assm_global
		use imp_single_channel
		use imp_loop_global
        use imp_re_input_global
		implicit none
        public::driving_TH_core
        public::driving_loop_flowAlloc
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
			real(KREAL)::density,flowrate,Qinpart
            real(KREAL),allocatable::powersingle(:)
            real(KREAL)::zonepower
            real(KREAL)::qtotal,dqi

            fq_core=1.0D0
            Nzone=core%Nzone!需要计算的zone的数量
            
			nr = SIZE(assembly, dim=1)!zone                           
			na = SIZE(assembly, dim=2)!layer                          
			M=size(assm1(1)%thermal%temperature,dim=1)
			N=size(assm1(1)%thermal%temperature,dim=2)		
            allocate(pow(na,N),fq_core(na,N))
            allocate(powersingle(Nzone))
			pow=0.0
			fq_core=1.0
		
			do i=1,Nzone,1
				assm1(i)%th_boundary%T%inlet=Tin
			enddo
			Qinpart=Qin/core%Nsplit
			call driving_loop_flowAlloc(assm1,Qinpart)
            !powertotal=0.0
			do izone=1,Nzone,1!zone start
				do j=1,assm1(izone)%mesh%ny,1!dy
					do k=1,N,1
						! if(k<=assm1(izone)%mesh%Nf) pow(j,k)=assembly(izone,j+assm1(izone)%mesh%layer_bottom)/(assm1(izone)%geom%N_fuelpin*assm1(izone)%geom%height(j)*PI*assm1(izone)%geom%pellet**2)!W/m3
						if(k<=assm1(izone)%mesh%Nf) pow(j,k)=assembly(izone,j)/(assm1(izone)%geom%N_fuelpin*assm1(izone)%geom%height(j)*PI*assm1(izone)%geom%pellet**2)!W/m3
                    enddo
                    !powertotal=powertotal+assembly(izone,j)
				enddo
                if (transient_flag)then
                    if(PRESENT(last).and.PRESENT(current)) call driving_imp_THtransient(assm1(izone),pow,fq_core,last,current)
                else
                    call driving_imp_THsteady(assm1(izone),pow,fq_core)					
                endif
                ! powersingle(izone)=151.0*assm1(izone)%hydrau%Qf*(assm1(izone)%th_boundary%T%outlet-assm1(izone)%th_boundary%T%inlet)
			enddo !zone		
			!Tout volum ave
			Tout=0.0
			do izone=1,Nzone,1
				density=get_density_Na(assm1(izone)%th_boundary%T%inlet)
				flowrate=assm1(izone)%hydrau%Qf*assm1(izone)%geom%N_pin
				Tout=Tout+assm1(izone)%th_boundary%T%outlet*flowrate/Qinpart
            enddo
            Tout=Tout+core%sigmaPass*Qinpart*Tin/Qinpart
			core%Tfout=Tout
            !v&v power
            
            qtotal=0.0
            do izone=1,Nzone,1
                associate(Qf=>assm1(izone)%hydrau%Qf,&
                          N_pin=>assm1(izone)%geom%N_pin,&
                          Ny=>assm1(izone)%mesh%ny,&
                          Nf=>assm1(izone)%mesh%nf,&
                          Ng=>assm1(izone)%mesh%ng,&
                          Ns=>assm1(izone)%mesh%ns,&
                          shc=>assm1(izone)%property%shc,&
                          temperature=>assm1(izone)%thermal%temperature,&
                          Tin=>assm1(izone)%th_boundary%T%inlet,&
                          Tout=>assm1(izone)%th_boundary%T%outlet)
                dqi=0.0
                do j=1,Ny,1
                    if(j==1) then
                        dqi=dqi+Qf*N_pin*shc(j,Nf+Ng+Ns+1)*(temperature(j,Nf+Ng+Ns+1)-Tin)
                    elseif(j==Ny) then
                        dqi=dqi+Qf*N_pin*shc(j,Nf+Ng+Ns+1)*(Tout-temperature(j,Nf+Ng+Ns+1))
                    else
                        dqi=dqi+Qf*N_pin*shc(j,Nf+Ng+Ns+1)*(temperature(j,Nf+Ng+Ns+1)-temperature(j-1,Nf+Ng+Ns+1))
                    endif
                enddo
                end associate
                qtotal=qtotal+dqi
            enddo
            core%vqtotal=qtotal
        end subroutine driving_TH_core       
        
        subroutine driving_loop_flowAlloc(assm,Qin)
            type(sys_assembly),intent(in out)::assm(:)
            real(KREAL),intent(in)::Qin
            !local
            integer zone
            integer i
            zone=core%Nzone
            !以zone为统计，而非SA
            !Qave=Qin*(1-core%sigmaPass)/(core%Nflow+core%Nflowsemi/2)!average
            do i=1,zone,1
                assm1(i)%saflag=core%SAtable(i)
                assm(i)%hydrau%Qf=Qin*reInputdata%sa(assm1(i)%saflag)%flowdis/assm(i)%geom%N_pin
            enddo
            ! do i=1,zone,1
                ! izone=core%fzone(i)
                ! if (i<=core%Nflow) then
                    ! flowrate=Qave
                ! else
                    ! flowrate=Qave/2.0
                ! endif
                ! assm(izone)%hydrau%Qf=flowrate/assm(izone)%geom%N_pin
            ! enddo
        end subroutine driving_loop_flowAlloc
    end module Imp_driving_THcore