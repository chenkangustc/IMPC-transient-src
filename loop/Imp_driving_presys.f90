module Imp_driving_presys
	use Imp_loop_global
	use Imp_inputcard
    use imp_re_input_global
    use loop_vtk_global
    ! use imp_assm_global
	contains
	subroutine driving_presys()
		implicit none
        !scan & read Ny,Nflow
        call driving_input_read()		
		!alloc
		call driving_alloc_loop()
		!init
		call driving_init_loop()
		!read after alloc       
        call driving_input_read_after()
        !select the mas flowrate zone
        ! call select_maxflowzone()
	end subroutine driving_presys
	
	! subroutine driving_plain_scan()
		! implicit none
		
	! end subroutine driving_plain_scan
	
	subroutine driving_init_loop()
		implicit none
		call core%init()
		call IHX1%init()
		call PipeRI%init()
		call PipeIP%init()
		call PipePR%init()
		call pump1%init()
        decayheat=0.04
    end subroutine driving_init_loop
	
	subroutine driving_alloc_loop()
		implicit none
		call IHX1%free()
		call IHX1%alloc()
        call PipeRI%alloc()
		call PipeIP%alloc()
		call PipePR%alloc()      
        call core%alloc()
        if(pump1%is_table==.TRUE.) call pump1%alloc()
        if(is_THonly) then
            allocate(tpower1%pow(2,tpower1%Ntime))
        endif
        allocate(core%SAtable(core%Nzone))
        !用来临时代替core
        allocate(reInputdata%height(reInputdata%ny))
        allocate(reInputdata%sa(reInputdata%Ntype))
	end subroutine driving_alloc_loop
	
	subroutine driving_free_loop()
		implicit none
		call IHX1%free()
        if(pump1%is_table==.TRUE.) call pump1%free()
        if(is_THonly) then  
            if(allocated(tpower1%pow)) deallocate(tpower1%pow)
        endif
        if(allocated(reInputdata%sa)) deallocate(reInputdata%sa)
        call loopave%free()
	end subroutine driving_free_loop
    

end module Imp_driving_presys