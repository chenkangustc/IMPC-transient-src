module Imp_driving_syspost
	use Imp_loop_global!common data
    use Imp_inputcard!IOfile
	use Imp_driving_presys, only:driving_free_loop
    implicit none
	contains
	subroutine driving_output_steady()
		implicit none
		integer    ::i
		real(KREAL)::LL
		!outputIHX
        write(unit=FILE_O,fmt="('  height','   Ts  ',' Ttube','  Tp ','  Tshell')")
        do i=0,IHX1%N+1,1
          if(i==0)then
              write(unit=FILE_O,fmt="('  0.0000',F8.1,'        ',F8.1)") IHX1%Tsin,IHX1%Tpout
          elseif(i==IHX1%N+1)then
              write(unit=FILE_O,fmt="(F8.4,F8.1,'        ',F8.1)") IHX1%Lsingle,IHX1%Tsout,IHX1%Tpin          
          else
              write(unit=FILE_O,fmt="(F8.4,4F8.1)") IHX1%zz(i),IHX1%Ts(i),IHX1%Tt(i),IHX1%Tp(i),IHX1%Tv(i)
          end if
        enddo
        !output core
        write(unit=FILE_O,fmt="('  height','  Tf  ','  Tshell')")
        LL=0.0
        do i=0,core%Ny+1,1
          if(i==0)then
              write(unit=FILE_O,fmt="('  0.0000',F8.1)") core%Tfin
          elseif(i==core%Ny+1)then
              write(unit=FILE_O,fmt="(F8.4,F8.1)") core%Ltotal,core%Tfout         
          else
			 if(i==1)then 
				LL=core%length(i)/2.0
			 else
				LL=LL+(core%length(i)+core%length(i-1))/2.0
             end if
              write(unit=FILE_O,fmt="(F8.4,2F8.1)") LL,core%Tf(i),core%Ts(i)
          end if
        enddo 


	end subroutine driving_output_steady
	
	subroutine driving_postsys()
        implicit none
        integer::i
        real(KREAL)::LL
        !outputIHX
        ! write(unit=FILE_O,fmt="('  height','   Ts  ',' Ttube','  Tp ','  Tshell')")
        ! do i=0,IHX1%N+1,1
          ! if(i==0)then
              ! write(unit=FILE_O,fmt="('  0.0000',F8.1,'        ',F8.1)") IHX1%Tsin,IHX1%Tpout
          ! elseif(i==IHX1%N+1)then
              ! write(unit=FILE_O,fmt="(F8.4,F8.1,'        ',F8.1)") IHX1%Lsingle,IHX1%Tsout,IHX1%Tpin          
          ! else
              ! write(unit=FILE_O,fmt="(F8.4,4F8.1)") IHX1%zz(i),IHX1%Ts(i),IHX1%Tt(i),IHX1%Tp(i),IHX1%Tv(i)
          ! end if
        ! enddo
        !output core
        ! write(unit=FILE_O,fmt="('  height','  Tf  ','  Tshell')")
        ! LL=0.0
        ! do i=0,core%Ny+1,1
          ! if(i==0)then
              ! write(unit=FILE_O,fmt="('  0.0000',F8.1)") core%Tfin
          ! elseif(i==core%Ny+1)then
              ! write(unit=FILE_O,fmt="(F8.4,F8.1)") core%Ltotal,core%Tfout         
          ! else
			 ! if(i==1)then 
				! LL=core%length(i)/2.0
			 ! else
				! LL=LL+(core%length(i)+core%length(i-1))/2.0
             ! end if
              ! write(unit=FILE_O,fmt="(F8.4,2F8.1)") LL,core%Tf(i),core%Ts(i)
          ! end if
        ! enddo 
        !free
        call driving_free_loop()
        !close
        close(FILE_I)
        close(FILE_O)
		close(FILE_T)
		close(FILE_maxT)
		close(FILE_aveT)
    end subroutine driving_postsys
end module Imp_driving_syspost