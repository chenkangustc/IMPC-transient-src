module testNK2loop
	use global_state
    use imp_assm_global
    use TH2NK_interface_loop
    use imp_loop_global
    use imp_timer_global
    use output_visit,               only : Print_vtk_files
	implicit none
    contains
    subroutine driving_testNK2loop()

        real(KREAL)  :: power(ns%state%zone, ns%state%layer)
		real(KREAL)  :: powerSteady(ns%state%zone, ns%state%layer)
        real(KREAL)  :: fq(ns%state%zone, ns%state%layer)
        
        logical  :: transient_flag = .FALSE.
        real(KREAL)  :: Tfuel(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: Tcoolant(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: Rhocoolant(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: toutlet
        real(KREAL)  :: max_Tfuel 
        real(KREAL)  :: max_Tcoolant 
        real(KREAL)  :: min_Rhocoolant 
        real(KREAL)  :: last 
        real(KREAL)  :: current 
		!local
		integer Nradial,i_zone,nTime,i
        real(KREAL) ::tTotal,dtime
        
        Tfuel = 0.0; Tcoolant = 0.0; Rhocoolant = 0.0; 
        max_Tfuel = 0.0; max_Tcoolant = 0.0; min_Rhocoolant = 0.0; 
        last = 0.0; current = 0.0;
	    
	   transient_flag=.TRUE.
	   open(unit=1,file='.\output\powDistribution.txt')
       read(1,100) power
	   100 Format(F15.5)
	   close(1)

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
       		  tTotal=timer1%ttotal
			  nTime=timer1%Nt
			  dtime=tTotal/nTime
              write(*,fmt="('------------------------------------------------------------------------------')")
              write(*,fmt="(' ','time','   ','maxTfuel','  ','maxTcoolant',' ','coreTin',' ','coreTout','   ','IHXTpin','   ','IHXTpout')")
              write(*,fmt="('------------------------------------------------------------------------------')")
			  do i=1,nTime,1
				  current=current+dtime
				  call get_pow(current,power,powerSteady)
				  call Perform_TH_loop(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)  
                  write(*,fmt="(F5.1,'|',6F10.2)")  current,max_Tfuel,max_Tcoolant,core%Tfin,core%Tfout,IHX1%Tpin,IHX1%Tpout
				  last=current
                  ! call Print_vtk_files (is_adjoint=.FALSE., is_transient=.TRUE., tidx=id_Total, ctime=ctime)
                  call Print_vtk_files (.FALSE., .TRUE., i, current)
              enddo
	   endif
	      
    !   i_zone=19
	   !Nradial=assm1(i_zone)%mesh%Nf+assm1(i_zone)%mesh%Ng+assm1(i_zone)%mesh%Ns+1
	   !print*,'maxTfuel=',max_Tfuel,'maxTcoolant=',max_Tcoolant
	   !print*,'Tcoolant=',Tfuel(i_zone,:)
	   !print*,'zone=',i_zone,'Tinlet=',assm1(i_zone)%th_boundary%T%inlet
	   !print*,'zone=',i_zone,'Temperature=',assm1(i_zone)%thermal%temperature(:,Nradial)
    !   print*,'zone=',i_zone,'uinlet=',assm1(i_zone)%th_boundary%u%inlet
	   !print*,'zone=',i_zone,'velocity=',assm1(i_zone)%thermal%velocity
	   !print*,'zone=',i_zone,'pressure=',assm1(i_zone)%thermal%pressure
    !   read(*,*)
    end subroutine driving_testNK2loop
	
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
		! pow=powS
        if(current>=1.and.current<2.) then
            pow=pows*(1.0-current*0.1)
        elseif(current>=2..and.current<4.) then
            pow=pows*(0.8-(current-2.0)*0.35)
        elseif(current>=4..and.current<=14.) then
            pow=pows*(0.1-(current-4.0)*0.005)
        else
            pow=pows*0.05
        endif
		!pow=0.0
        do i=1,nr,1
            do j=1,na,1
                powtotal=powtotal+pow(i,j)
                powStotal=powStotal+powS(i,j)
            enddo
        enddo
	end subroutine get_pow
end module testNK2loop