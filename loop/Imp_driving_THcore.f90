	module Imp_driving_THcore
		use constants
		use imp_property
		use Imp_assm_global
		use imp_single_channel
		use imp_loop_global
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

            fq_core=1.0D0
            Nzone=core%Nflow+core%Nflowsemi!需要计算的zone的数量
            
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
			do i=1,Nzone,1!zone start
                izone=core%fzone(i)
				do j=1,assm1(izone)%mesh%ny,1!dy
					do k=1,N,1
						if(k<=assm1(izone)%mesh%Nf) pow(j,k)=assembly(izone,j+assm1(izone)%mesh%layer_bottom)/(assm1(izone)%geom%N_fuelpin*assm1(izone)%geom%height(j)*PI*assm1(izone)%geom%pellet**2)!W/m3
                    enddo
                    !powertotal=powertotal+assembly(izone,j)
				enddo
                zonepower=0.0
                zonepower=sum(assembly(izone,:))
				if (zonepower==0.0) then
					assm1(izone)%thermal%temperature=assm1(izone)%th_boundary%T%inlet
					assm1(izone)%th_boundary%T%outlet=assm1(izone)%th_boundary%T%inlet
				else
                    if (transient_flag)then
                        if(PRESENT(last).and.PRESENT(current)) call driving_imp_THtransient(assm1(izone),pow,fq_core,last,current)
                    else
                        call driving_imp_THsteady(assm1(izone),pow,fq_core)					
                    endif
				endif	
                powersingle(i)=151.0*assm1(izone)%hydrau%Qf*(assm1(izone)%th_boundary%T%outlet-assm1(izone)%th_boundary%T%inlet)
			enddo !zone		
			!Tout volum ave
			Tout=0.0
			do i=1,Nzone,1
                izone=core%fzone(i)
				density=get_density_Na(assm1(izone)%th_boundary%T%inlet)
				flowrate=assm1(izone)%hydrau%Qf*assm1(izone)%geom%N_pin
				Tout=Tout+assm1(izone)%th_boundary%T%outlet*flowrate/Qinpart
			enddo
            Tout=Tout+core%sigmaPass*Qinpart*Tin/Qinpart
			core%Tfout=Tout
		end subroutine driving_TH_core
        
        subroutine driving_loop_flowAlloc(assm,Qin)
            type(sys_assembly),intent(in out)::assm(:)
            real(KREAL),intent(in)::Qin
            !local
            real(KREAL)::flowrate
            real(KREAL)::fuelArea,density!flow area of pin
            real(KREAL)::Qave
            integer zone,n_pin,layer,nr,izone
            integer i,j
            zone=size(core%fzone)
            layer=size(assm(1)%thermal%velocity)
            nr=size(assm(1)%thermal%temperature,2)
            !以zone为统计，而非SA
            Qave=Qin*(1-core%sigmaPass)/(core%Nflow+core%Nflowsemi/2)!average
            do i=1,zone,1
                izone=core%fzone(i)
                if (i<=core%Nflow) then
                    flowrate=Qave
                else
                    flowrate=Qave/2.0
                endif
                assm(izone)%hydrau%Qf=flowrate/assm(izone)%geom%N_pin
            enddo
        end subroutine driving_loop_flowAlloc
    end module Imp_driving_THcore