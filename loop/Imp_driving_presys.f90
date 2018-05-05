module Imp_driving_presys
	use Imp_loop_global
	use Imp_inputcard
    use imp_re_input_global
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
    end subroutine driving_init_loop
	
	subroutine driving_alloc_loop()
		implicit none
		call IHX1%free()
		call IHX1%alloc()
        call PipeRI%alloc()
		call PipeIP%alloc()
		call PipePR%alloc()      
        call core%alloc()
        !用来临时代替core
        allocate(reInputdata%height(reInputdata%ny))
	end subroutine driving_alloc_loop
	
	subroutine driving_free_loop()
		implicit none
		call IHX1%free()
	end subroutine driving_free_loop
end module Imp_driving_presys