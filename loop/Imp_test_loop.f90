
!****************************************************************************
!
!  PROGRAM: system
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    module Imp_test_loop
    use Imp_driving_presys
	use Imp_cal_loop
    use Imp_timer_global
    use Imp_inputcard
    use Imp_driving_syspost
    use Imp_test_IHX_thermal
    use Imp_test_PIPE_thermal
    implicit none
    contains
    subroutine driving_loop_interface()
    !local
    integer::i,Nt
	real(KREAL)::dt
	!call test_loop_hydraulic()
	!call test_IHX_thermal()
    !call test_PIPE_thermal()
	associate(Nt=>timer1%Nt,ttotal=>timer1%ttotal,last=>timer1%last,current=>timer1%current)

    call driving_presys()
    dt=ttotal/Nt
	last=0.0
	current=0.0
    
	!call driving_loop_steady()
    
    write(unit=FILE_O,fmt="('  time',' flowrate',' coreTout',' IHXTin',' IHXTout')")
	do i=1,Nt,1
        current=last+dt
!		call driving_loop_transient(last,current)
        write(unit=FILE_O,fmt="(F6.1,' ',4F8.2)") current,core%Q,core%Tfout,IHX1%Tpin,IHX1%Tpout
        last=current
	enddo
	end associate
    call driving_postsys()
    end subroutine driving_loop_interface
    end module Imp_test_loop

