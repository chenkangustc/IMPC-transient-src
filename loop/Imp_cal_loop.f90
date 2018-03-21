module Imp_cal_loop
	use constants
    use Imp_loop_global
	use Imp_driving_THcore
	implicit none
	contains
	subroutine driving_loop_steady(assembly)
		implicit none
		real(KREAL),intent(in)::assembly(:,:)! (zone, layer), in W, 各组件功率;
        !local
        real(KREAL)::sigma
		real(KREAL)::coreTin,coreTout,coreQin
        sigma=1.0
		coreTin=0.0
		coreTout=0.0
		coreQin=0.0
        do while(sigma>0.001)
            coreTin=PipePR%Tfout
			coreQin=PipePR%Q
            print *,'core cal'
			!	 driving_THcore_steady(Qin,Tin,assembly,Tout)
		    call driving_THcore_steady(coreQin,coreTin,assembly,coreTout)
		    pipeRI%Tfin=coreTout
            print*,'pipe cal'
		    call PipeRI%thCals()
		    IHX1%Tpin=PipeRI%Tfout
		    call IHX1%thCals()
		    PipeIP%Tfin=IHX1%Tpout
		    call PipeIP%thCals()
            PipePR%Tfin=PipeIP%Tfout
            call PipePR%thCals()
            sigma=abs((PipePR%Tfout-core%Tfin)/core%Tfin)
            print*,sigma
        enddo
	end subroutine driving_loop_steady
	
	subroutine driving_loop_transient(last,current)
        implicit none
        real(KREAL),intent(in)::last,current
        
            call cal_loop_hydraulic(current)    
		    call core%thCalt(last,current)
		    pipeRI%Tfin=core%Tfout
            print*,'pipe cal'
		    call PipeRI%thCalt(last,current)
		    IHX1%Tpin=PipeRI%Tfout
		    call IHX1%thCalt(last,current)
		    PipeIP%Tfin=IHX1%Tpout
		    call PipeIP%thCalt(last,current)
            PipePR%Tfin=PipeIP%Tfout
            call PipePR%thCalt(last,current)
            core%Tfin=PipePR%Tfout    
            
	end subroutine driving_loop_transient
	
	subroutine cal_loop_hydraulic(current)
		implicit none
		real(KREAL),INTENT(in)::current!time
		!local
		real(KREAL)::alpha,beta,LAsum
		real(KREAL)::flowrate
		integer i!i is the num of the device
		LAsum=PipePR%ltotal/PipePR%Area+core%Ltotal/core%Area+PipeRI%ltotal/PipeRI%Area+IHX1%Lsingle/IHX1%Areap+PipeIP%ltotal/PipeIP%Area
		alpha=LAsum+pump1%rho*pump1%yita*pump1%I*pump1%omegae**2/pump1%Qe**2
		call cal_beta(beta,0)
		flowrate=alpha/(beta*current+alpha/pump1%Qe)
        call set_flowrate(flowrate)
	end subroutine cal_loop_hydraulic
	
	subroutine cal_beta(beta,formula)
		implicit none
		real(KREAL),intent(in out)::beta
		integer,intent(in)::formula!if formula==1 then K/fric else he/we
		select case(formula)
        case(0)!rho keep constant
			beta=pump1%rho*9.80*pump1%He/pump1%Qe**2
		case(1)!rho changes	
            beta=PipePR%beta+core%beta+PipeRI%beta+IHX1%betap+PipeIP%beta
		end select
    end subroutine cal_beta
    
	subroutine test_loop_hydraulic()
		implicit none
		!local
		real(KREAL)::alpha,beta,LAsum
		real(KREAL)::flowrate,ffe
		integer i,current!i is the num of the device
		!LAsum=PipePR%Length/PipePR%Area+core%Length/core%Area+PipeRI%Length/PipeRI%Area+IHX1%Lsingle/IHX1%AreaOuter+PipeIP%Length/PipeIP%Area
		LAsum=959.0!m-1
		pump1%Qe=125.0
        pump1%He=41.9
		pump1%yita=0.756
		pump1%omegae=49*PI
		pump1%I=3.7
		pump1%rho=1000.0
		alpha=LAsum+pump1%rho*pump1%yita*pump1%I*pump1%omegae**2/pump1%Qe**2
		call cal_beta(beta,0)
		
		open(1,file='.\flowrate.txt')
		do i=0,25,1
			current=i
			flowrate=alpha/(beta*current+alpha/pump1%Qe)
            call set_flowrate(flowrate)
			print*, flowrate
            ffe=flowrate/pump1%Qe
			write(1,100) current,ffe
			100 Format(I6,' ',F10.2)
		enddo
		close(1)
    end subroutine test_loop_hydraulic
    
    subroutine set_flowrate(flowrate)
        implicit none
        real(KREAL),intent(in)::flowrate
            pump1%Q=flowrate
            PipePR%Q=flowrate
            core%Q=flowrate
            PipeRI%Q=flowrate
            IHX1%Qp=flowrate
            PipeIP%Q=flowrate
    end subroutine set_flowrate
end module Imp_cal_loop
