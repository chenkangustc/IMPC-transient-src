module Imp_THonly
	use global_state
    use imp_assm_global
    use TH2NK_interface_loop
    use imp_loop_global
    use imp_re_input_global
	implicit none
    contains
    subroutine driving_THloop()

        real(KREAL)  :: power(core%Nzone,reInputdata%Ny)
		real(KREAL)  :: powerSteady(core%Nzone,reInputdata%Ny)
        real(KREAL)  :: fq(core%Nzone,reInputdata%Ny)
        
        logical  :: transient_flag = .FALSE.
        real(KREAL)  :: Tfuel(core%Nzone,reInputdata%Ny) 
        real(KREAL)  :: Tcoolant(core%Nzone,reInputdata%Ny) 
        real(KREAL)  :: Rhocoolant(core%Nzone,reInputdata%Ny) 
        real(KREAL)  :: toutlet
        real(KREAL)  :: max_Tfuel 
        real(KREAL)  :: max_Tcoolant 
        real(KREAL)  :: min_Rhocoolant 
        real(KREAL)  :: last 
        real(KREAL)  :: current 
        real(KREAL)  :: powtotal
		!local
		integer Nradial,i_zone,nTime,i,j
        integer::Nzone,izone
        real(KREAL) ::tTotal,dtime
        Nzone=core%Nzone
        Tfuel = 0.0; Tcoolant = 0.0; Rhocoolant = 0.0; 
        max_Tfuel = 0.0; max_Tcoolant = 0.0; min_Rhocoolant = 0.0; 
        last = 0.0; current = 0.0;
	    
        transient_flag=.TRUE.
        power=0.0
        do izone=1,Nzone,1
            assm1(izone)%saflag=core%SAtable(izone)
            power(izone,:)=tpower1%pow(2,1)*1.0e6*reInputdata%sa(assm1(izone)%saflag)%powdis/reInputdata%Ny
        enddo
        if(transient_flag==.FALSE.) then
            if (ns%feedback%is_loop)  then
	          call Perform_TH_loop(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)  
			endif  
        else       
			transient_flag=.FALSE.
			call Perform_TH_loop(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)  	   
			!power=0.8*power
			powerSteady=power
			transient_flag=.TRUE.
       		tTotal=900.0
			nTime=900
			dtime=tTotal/nTime
            write(*,fmt="('------------------------------------------------------------------------------')")
            write(*,fmt="(' ','time','   ','maxTfuel','  ','maxTcoolant',' ','coreTin',' ','coreTout','   ','IHXTpin','   ','IHXTpout')")
            write(*,fmt="('------------------------------------------------------------------------------')")
			do i=1,nTime,1
				current=current+dtime
                call set_pow(current,powtotal)
                do izone=1,Nzone,1
                    assm1(izone)%saflag=core%SAtable(izone)
                    power(izone,:)=powtotal*1.0e6*reInputdata%sa(assm1(izone)%saflag)%powdis/reInputdata%Ny
                enddo
				!call get_pow(current,power,powerSteady)
				call Perform_TH_loop(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)  
				write(*,fmt="(F5.1,'|',6F10.2)")  current,max_Tfuel,max_Tcoolant,core%Tfin,core%Tfout,IHX1%Tpin,IHX1%Tpout
                last=current
			enddo
        endif
            
    end subroutine driving_THloop
	
	subroutine get_pow(current,pow,powS)
		implicit none
		real(KREAL),intent(in)::current,powS(:,:)
		real(KREAL),intent(in out)::pow(:,:)
		!local
        integer::i,j,nr,na
        real(KREAL)::powtotal,powStotal
        powtotal=0.0
        powStotal=0.0
        nr=size(pow,dim=1)
        na=size(pow,dim=2)
		!if(current<=400.0) pow=-powS*(current-400.0)/400.0
		!pow=powS
		pow=0.0
        do i=1,nr,1
            do j=1,na,1
                powtotal=powtotal+pow(i,j)
                powStotal=powStotal+powS(i,j)
            enddo
        enddo
	end subroutine get_pow
    
    subroutine set_pow(current,powtotal)
        real(KREAL),intent(in)::current
        real(KREAL),intent(in out)::powtotal
        !local
        real(KREAL)::cpow
        integer::i
        associate(dtime=>tpower1%pow(1,:),&
                  Ntime=>tpower1%Ntime   ,&
                  dpow=>tpower1%pow(2,:))
        do i=1,Ntime-1,1
            if(current>=dtime(i).and.current<=dtime(i+1)) then
                cpow=(dpow(i+1)-dpow(i))/(dtime(i+1)-dtime(i))*(current-dtime(i))+dpow(i)
                exit
            endif
            if(current>dtime(Ntime)) cpow=dpow(Ntime)
        enddo
        powtotal=cpow
        end associate
        ! do i=1,Ntime-1,1
            ! if(current>=dtime(i).and.current<=dtime(i+1)) then
                ! crotate=(rotate(i+1)-rotate(i))/(dtime(i+1)-dtime(i))*(current-dtime(i))+rotate(i)
                ! exit
            ! endif
            ! if(current>dtime(Ntime)) crotate=rotate(Ntime)
        ! enddo
    end subroutine set_pow
    end module Imp_THonly