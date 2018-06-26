module Imp_driving_syspost
	use Imp_loop_global!common data
    use Imp_assm_global!core
	use Imp_inputcard!IOfile
	use Imp_driving_presys, only:driving_free_loop
    
	implicit none
	contains
	subroutine loop_output_steady()
		integer::ny,n_zone
		integer::i,j
		ny=assm1(1)%mesh%ny
		n_zone=assm1(1)%mesh%n_zone
		write(unit=File_o,fmt="(<Ny>F8.2)") (assm1(1)%mesh%z(i,1),i=1,ny)
		do i=1,n_zone,1
			write(unit=File_o,fmt="(<Ny>F8.2)") assm1(i)%thermal%Tfuel(:)
			write(unit=File_o,fmt="(<Ny>F8.2)") assm1(i)%thermal%Tcoolant(:)
		enddo

	end subroutine loop_output_steady
	
	subroutine loop_output_transient(current)
		real(KREAL),intent(in)::current
		!local
		integer::nr,nave,ny
		integer::i,izone
		real(KREAL),allocatable::aveT(:)
		nr=assm1(1)%mesh%n_zone
        ny=assm1(1)%mesh%ny
		nave=4*nr
		allocate(aveT(nave))
		do i=1,nr,1
			aveT(4*i-3)=assm1(i)%thermal%Tfave
			aveT(4*i-2)=assm1(i)%thermal%Tcave
			aveT(4*i-1)=assm1(i)%th_boundary%T%inlet
			aveT(4*i)=assm1(i)%th_boundary%T%outlet
		enddo
		write(unit=File_aveT,fmt="(F6.1,' ',<Nave>F8.2)") current,aveT!(aveT(i),i=1,Nave)
		!looptimelist.txt
		!Tdis.txt
		izone=izoneTdis!
		write(unit=file_disT,fmt="(F6.1,' ',<Ny>F8.2)") current,(assm1(1)%mesh%z(i,1),i=1,ny)
		write(unit=file_disT,fmt="(F6.1,' ',<Ny>F8.2)") current,assm1(izone)%thermal%Tfuel(:)
		write(unit=file_disT,fmt="(F6.1,' ',<Ny>F8.2)") current,assm1(izone)%thermal%Tcoolant(:)
        izone=izoneTdis2
        write(unit=file_disT,fmt="(F6.1,' ',<Ny>F8.2)") current,(assm1(1)%mesh%z(i,1),i=1,ny)
		write(unit=file_disT,fmt="(F6.1,' ',<Ny>F8.2)") current,assm1(izone)%thermal%Tfuel(:)
		write(unit=file_disT,fmt="(F6.1,' ',<Ny>F8.2)") current,assm1(izone)%thermal%Tcoolant(:)
		
	end subroutine loop_output_transient
	
	subroutine driving_postsys()
        implicit none
        integer::i
        real(KREAL)::LL
        !free
        call driving_free_loop()
        !close
        close(FILE_I)
        close(FILE_O)
		close(FILE_T)
		close(FILE_HY)
		close(FILE_maxT)
		close(FILE_aveT)
		close(FILE_disT)
    end subroutine driving_postsys
end module Imp_driving_syspost

	! subroutine driving_output_steady()
		! implicit none
		! integer    ::i
		! real(KREAL)::LL
		! !outputIHX
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
        ! !output core
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


	! end subroutine driving_output_steady