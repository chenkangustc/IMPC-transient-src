module Imp_cal_loop
	use constants
    use Imp_loop_global
	use Imp_driving_THcore
	use Imp_inputcard
	use Imp_driving_syspost
	implicit none
	contains
	subroutine driving_loop_steady(assembly)
		implicit none
		real(KREAL),intent(in)::assembly(:,:)! (zone, layer), in W, 各组件功率;
        !local
		real(KREAL)::current
        real(KREAL)::sigma
		real(KREAL)::coreTin,coreTout,coreQin
		real(KREAL)::powinput
        integer::num
		logical:: transient_flag
        is_natural=.FALSE.
		transient_flag=.FALSE.
		current=0.0
        sigma=1.0
		coreTin=0.0
		coreTout=0.0
		coreQin=0.0
		powinput=0.0
        num=0
        Qloop=pump1%Nbranch*Pump1%Qe
        core%Qtotal=Qloop
		!cal total input power for output
		call cal_total_inputpow(assembly,powinput)
		!write(*,*)'driving loop steady:'
        write(*,fmt="('------------------------------------------------------------------------------')") 
        write(*,fmt="('  ','num','   ','sigma','      ','coreTin','  ','coreTout','  ','IHXTpin','   ','IHXTpout','   ','IHXTsin','  ','IHXTsout')") 
        write(*,fmt="('------------------------------------------------------------------------------')") 
        do while(sigma.gt.1.0D-6)
            num=num+1
            coreTin=PipePR%Tfout
			coreQin=pump1%Nbranch*Pump1%Qe
            !print *,'core cal'
			!	 driving_THcore_steady(Qin,Tin,assembly,Tout)
		    call driving_TH_core(transient_flag,coreQin,coreTin,assembly,coreTout)
		    pipeRI%Tfin=coreTout
           ! print*,'pipe cal steady'
            if(is_Bq) call cal_Bq()
		    call PipeRI%thCals()
		    IHX1%Tpin=PipeRI%Tfout
		    call IHX1%thCals()
		    PipeIP%Tfin=IHX1%Tpout
		    call PipeIP%thCals()
            PipePR%Tfin=PipeIP%Tfout
            call PipePR%thCals()
            sigma=abs((PipePR%Tfout-coreTin)/coreTin)
            ! print*,'num=',num,'sigma=',sigma,'coreTin=',coreTin,'coreTout=',coreTout
            ! write(*,fmt="('num=',I4,'  sigma=',F10.7,'  coreTin=',F8.2,'  coreTout=',F8.2,'IHXTpin=',F8.2,'IHXTpout=',F8.2)") num,sigma,coreTin,coreTout,IHX1%Tpin,IHX1%Tpout
            write(*,fmt="(I4,'|',F10.7,6F10.2)")num,sigma,coreTin,coreTout,IHX1%Tpin,IHX1%Tpout,IHX1%Tsin,IHX1%Tsout
        enddo
		!call driving_output_steady()
		write(unit=file_t,fmt="(F6.1,' ',F10.1,8F8.2,2F15.1)") current,powinput,Pump1%Qe,coreTin,coreTout,IHX1%Tpin,IHX1%Tpout,IHX1%Qs,IHX1%Tsin,IHX1%Tsout,core%vqtotal,IHX1%vqtotal
	end subroutine driving_loop_steady
	
	subroutine driving_loop_transient(assembly,last,current)
        implicit none
		real(KREAL),intent(in)::assembly(:,:)!power(zone,layer)
        real(KREAL),intent(in)::last,current
        !local
		real(KREAL)::coreTin,coreTout,coreQin
		real(KREAL)::powinput
		integer::i,j,nr,na
		logical::transient_flag,is_table
		
		transient_flag=.TRUE.
        is_table=.TRUE.

		call cal_total_inputpow(assembly,powinput)
        call cal_loop_hydraulic(is_table,is_natural,last,current,Qloop)		
        call update_secflow(current)
        call update_Tsin(current)
        coreQin=pump1%Nbranch*Qloop
        core%Qtotal=coreQin
        coreTin=PipePR%Tfout
        core%Tfin=PipePR%Tfout
        ! print *,'core cal'
        !driving_THcore_steady(Qin,Tin,assembly,Tout)
        call driving_TH_core(transient_flag,coreQin,coreTin,assembly,coreTout,last,current)
        pipeRI%Tfin=coreTout
        core%Tfout=coreTout
        ! print*,'pipe cal transient'
        if(is_Bq) call cal_Bq()
        call PipeRI%thCalt(last,current)
        IHX1%Tpin=PipeRI%Tfout
        call IHX1%thCalt(last,current)
        PipeIP%Tfin=IHX1%Tpout
        call PipeIP%thCalt(last,current)
        PipePR%Tfin=PipeIP%Tfout
        call PipePR%thCalt(last,current) 
        write(unit=file_t,fmt="(F6.1,' ',F10.1,8F8.2)") current,powinput,Qloop,coreTin,coreTout,IHX1%Tpin,IHX1%Tpout,IHX1%Qs,IHX1%Tsin,IHX1%Tsout
        write(unit=file_hy,fmt="(F6.1,' ',F8.2,F8.4,2F15.1,F8.2,F10.6)") current,PipeRI%velocity,PipeRI%fric,PipeRI%Re,PipeRI%rhof(1),PipeRI%Q,PipeRI%area
          
	end subroutine driving_loop_transient
    
	subroutine cal_loop_hydraulic(is_table,is_natural,last,current,flowrate)
        logical,intent(in)::is_table
        logical,intent(in out)::is_natural
		real(KREAL),INTENT(in)::last!time
		real(KREAL),INTENT(in)::current!time
		real(KREAL),INTENT(out)::flowrate
        !local
        real(KREAL)::flowf,flows
        flowf=flowrate
        flows=flowrate
        if(is_natural==.FALSE.) then
            call cal_hydraulic_fir(is_table,current,flowf)
            call cal_hydraulic_sec(last,current,flows)
            if(flows<flowf) then
                flowrate=flowf
            else
                flowrate=flows
                is_natural=.TRUE.
            endif
        else
            call cal_hydraulic_sec(last,current,flows)
            flowrate=flows
        endif
    end subroutine cal_loop_hydraulic
    
	subroutine cal_hydraulic_fir(is_table,current,flowrate)
		implicit none
        logical,intent(in)::is_table
		real(KREAL),INTENT(in)::current!time
		real(KREAL),INTENT(out)::flowrate
		!local
		real(KREAL)::alpha,beta,LAsum
		integer i!i is the num of the device
        real(KREAL)::crotate
		if(is_table.eq..FALSE.)then
            !if(pump1%Q>=0.10*pump1%Qe) then
            !	LAsum=PipePR%ltotal/PipePR%Area+core%Ltotal/core%Area+PipeRI%ltotal/PipeRI%Area+IHX1%Lsingle/IHX1%Areap+PipeIP%ltotal/PipeIP%Area
            !	alpha=LAsum+pump1%rho*pump1%yita*pump1%I*pump1%omegae**2/pump1%Qe**2
            !	call cal_beta(beta,0)
            !	flowrate=alpha/(beta*current+alpha/pump1%Qe)
            !else
            !	flowrate=pump1%Q
            !endif
            flowrate=pump1%Qe
        else
            associate(dtime=>pump1%rotate(1,:),&
                      rotate=>pump1%rotate(2,:),&
                      Ntime=>pump1%Ntime,&
                      Qe=>pump1%Qe)
            do i=1,Ntime-1,1
                if(current>=dtime(i).and.current<=dtime(i+1)) then
                    crotate=(rotate(i+1)-rotate(i))/(dtime(i+1)-dtime(i))*(current-dtime(i))+rotate(i)
                    exit
                endif
                if(current>dtime(Ntime)) crotate=rotate(Ntime)
            enddo
            !if (crotate>=0.024*Qe) then
                flowrate=crotate/rotate(1)*Qe
            !else 
            !    flowrate=0.024*Qe
            !endif
            end associate
        endif
        call set_flowrate(flowrate)
	end subroutine cal_hydraulic_fir
    
	subroutine cal_hydraulic_sec(last,current,flowrate)
        implicit none
		real(KREAL),INTENT(in)::last!time
		real(KREAL),INTENT(in)::current!time
		real(KREAL),INTENT(out)::flowrate
        !local
        real(KREAL)::alpha,beta,buoy
        real(KREAL)::dt
        dt=current-last
        alpha=PipePR%ltotal/PipePR%Area+core%calpha()+PipeRI%ltotal/PipeRI%Area+IHX1%Lsingle/IHX1%Areap+PipeIP%ltotal/PipeIP%Area
        beta=PipePR%cbeta()+PipeRI%cbeta()+PipeIP%cbeta()+IHX1%cbeta()+core%cbeta()
        buoy=PipePR%cbuoy()+PipeRI%cbuoy()+PipeIP%cbuoy()+IHX1%cbuoy()+core%cbuoy()
        !this%betap=this%betap+0.5*(this%Fricp*this%length(i)/this%Dep+Kiouter)*1/(this%rhop(i)*this%Areap**2)
        !solve
        flowrate=(alpha*flowrate+buoy*dt)/(alpha-beta*flowrate*dt)
        call set_flowrate(flowrate)
   end subroutine cal_hydraulic_sec
    
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
            !core%Q=flowrate
            PipeRI%Q=flowrate
            IHX1%Qp=flowrate
            PipeIP%Q=flowrate
    end subroutine set_flowrate
	
	subroutine cal_total_inputpow(powinput,powtotal)
		implicit none
		real(KREAL),intent(in)::powinput(:,:)
		real(KREAL),intent(out)::powtotal
		!local
		integer::nr,na,i,j
		nr=size(powinput,dim=1)
		na=size(powinput,dim=2)
		powtotal=0.0
		do i=1,nr,1
			do j=1,na,1
				powtotal=powtotal+powinput(i,j)
			enddo
		enddo
	end subroutine cal_total_inputpow
    subroutine update_Tsin(current)
        real(KREAL),intent(in)::current
        !local
        integer::i
        associate(dtime=>IHX1%Tintable(1,:),&
                  Tsin=>IHX1%Tintable(2,:),&
                  Ntime=>IHX1%NTtime)
        if (IHX1%is_Tintable) then
            do i=1,Ntime-1,1
                if(current>=dtime(i).and.current<=dtime(i+1)) then
                    IHX1%Tsin=(Tsin(i+1)-Tsin(i))/(dtime(i+1)-dtime(i))*(current-dtime(i))+Tsin(i)
                    exit
                endif
                if(current>dtime(Ntime)) IHX1%Tsin=Tsin(Ntime)
            enddo            
        endif
        end associate
    end subroutine update_Tsin
    subroutine update_secflow(current)
        real(KREAL),intent(in)::current
        !local
        integer::i
        associate(dtime=>IHX1%flowtable(1,:),&
                  Qs=>IHX1%flowtable(2,:),&
                  Ntime=>IHX1%Nftime)
        if (IHX1%is_flowtable) then
            do i=1,Ntime-1,1
                if(current>=dtime(i).and.current<=dtime(i+1)) then
                    IHX1%Qs=(Qs(i+1)-Qs(i))/(dtime(i+1)-dtime(i))*(current-dtime(i))+Qs(i)
                    exit
                endif
                if(current>dtime(Ntime)) IHX1%Qs=Qs(Ntime)
            enddo          
        endif
        end associate
    end subroutine update_secflow
    subroutine cal_Bq()
        integer::i,Num
        real(KREAL)::Bqtotal
        Num=5!cool pool 前5个节点分配功率
        associate(Ny=>PipeRI%Ny,&
                  Bq=>PipeRI%Bq,&
                  Rtube=>PipeRI%Rtube,&
                  thicks=>PipeRI%thicks,&
                  Length=>PipeRI%Length,&
                  Ts=>PipeRI%Ts,&
                  ctc=>PipeRI%ks)
        do i=1,Ny,1
            Bq(i)=ctc*2*PI*(Rtube+thicks)*Length(i)*(IHX1%Tpout-Ts(i))/(thicks/2.0)
        enddo
        Bqtotal=sum(Bq)
        do i=1,Num,1
            PipeIP%Bq(i)=-Bqtotal/Num
            ! PipeRI%Bq(i)=Bqtotal/Num
        enddo
        end associate
    end subroutine cal_Bq
end module Imp_cal_loop
