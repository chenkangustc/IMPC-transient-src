module Imp_inputcard
	use Imp_loop_global
	use Imp_timer_global
    use,intrinsic::ISO_FORTRAN_ENV
    use constants
	implicit none
	integer::file_i,file_o,file_t
	integer,parameter,private::N_keyword=8
    integer,parameter,private::MAX_REAL_PARAMETER=50
	integer,parameter,private::MAX_INT_PARAMETER=50
    character(len=MAX_WORD_LEN),parameter::FILE_IN='loopinput.case'
    character(len=MAX_WORD_LEN),parameter::FILE_OUT='loopoutput.txt'
	character(len=MAX_WORD_LEN),parameter::FILE_LTIME='looptimelist.txt'
    character(len=MAX_WORD_LEN)::INP_SECTION(N_keyword)
    
	
	contains
	subroutine driving_input_read()
		implicit none
		!local
		integer::io_error
		real::dummy_real(MAX_REAL_PARAMETER)
		integer::dummy_int(MAX_INT_PARAMETER)
		character(len=MAX_WORD_LEN)::aline
		character(len=MAX_WORD_LEN)::section_name,keyword
		! Variables
		call set_section_keyword()
        open(newunit=file_o,file=FILE_OUT,status='replace',action='write',iostat=io_error)
		open(newunit=file_t,file=FILE_LTIME,status='replace',action='write',iostat=io_error)   
		write(unit=file_t,fmt="('  time',' flowrate',' coreTin',' coreTout',' IHXTin',' IHXTout')")
		open(newunit=file_i,file=FILE_IN,status='old',action='read',iostat=io_error)   		
        !read(unit=file_i,fmt='(A)',iostat=io_error) aline
		do
			read(unit=file_i,fmt='(A)',iostat=io_error) aline
            if(io_error==IOSTAT_END) exit
			read(unit=aline,fmt=*,iostat=io_error) section_name       
			! if(is_keyword(INP_SECTION,section_name)) then
			!     backspace(file_i,iostat=io_error)
			!     exit
			! end if        
			if(is_keyword(INP_SECTION,section_name)) then
				select case(trim(adjustl(section_name)))
					case('pump')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:5)
					pump1%I=dummy_real(1)
					pump1%He=dummy_real(2)
					pump1%Qe=dummy_real(3)
					pump1%omegae=dummy_real(4)
					pump1%yita=dummy_real(5)
					
					case('pipePR')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:6),dummy_int(1)
					pipePR%ltotal=dummy_real(1)
					pipePR%Rtube=dummy_real(2)
					pipePR%thicks=dummy_real(3)
					pipePR%theta=dummy_real(4)
					pipePR%Q=dummy_real(5)
					pipePR%Ti=dummy_real(6)
					pipePR%Ny=dummy_int(1)
					
					case('reactor')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:6),dummy_int(1)
					core%ltotal=dummy_real(1)
					core%Rtube=dummy_real(2)
					core%thicks=dummy_real(3)
					core%theta=dummy_real(4)
					core%Q=dummy_real(5)
					core%Ti=dummy_real(6)
					core%Ny=dummy_int(1)
					
					case('pipeRI')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:6),dummy_int(1)
					pipeRI%ltotal=dummy_real(1)
					pipeRI%Rtube=dummy_real(2)
					pipeRI%thicks=dummy_real(3)
					pipeRI%theta=dummy_real(4)
					pipeRI%Q=dummy_real(5)
					pipeRI%Ti=dummy_real(6)
					pipeRI%Ny=dummy_int(1)
					
					case('IHX')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:9),dummy_int(1:2)
					IHX1%Lsingle=dummy_real(1)
					IHX1%Rtube=dummy_real(2)
					IHX1%thickt=dummy_real(3)
					IHX1%Plength=dummy_real(4)
					IHX1%AreaP=dummy_real(5)
					IHX1%thickv=dummy_real(6)
					
					IHX1%Qp=dummy_real(7)
					!IHX1%AreaS=dummy_real(7)
					IHX1%Qs=dummy_real(8)
					IHX1%Ti=dummy_real(9)
					IHX1%Ntube=dummy_int(1)
					IHX1%N=dummy_int(2)
					
					case('pipeIP')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:6),dummy_int(1)
					pipeIP%ltotal=dummy_real(1)
					pipeIP%Rtube=dummy_real(2)
					pipeIP%thicks=dummy_real(3)
					pipeIP%theta=dummy_real(4)
					pipeIP%Q=dummy_real(5)
					pipeIP%Ti=dummy_real(6)
					pipeIP%Ny=dummy_int(1)
					
					case('solid')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:4)
					IHX1%rhot=dummy_real(1)
					IHX1%shct=dummy_real(2)
					IHX1%rhov=dummy_real(3)
					IHX1%shcv=dummy_real(4)
					pipePR%rhos=dummy_real(3)
					pipePR%shcs=dummy_real(4)
					pipeRI%rhos=dummy_real(3)
					pipeRI%shcs=dummy_real(4)
					pipeIP%rhos=dummy_real(3)
					pipeIP%shcs=dummy_real(4)
					case('time')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1),dummy_int(1)
					timer1%ttotal=dummy_real(1)
					timer1%Nt=dummy_int(1)
				end select
			end if
			!print*,trim(aline)
			print*,trim(adjustl(section_name))
		end do
        close(file_i)
	end subroutine driving_input_read
	
    subroutine Set_section_keyword()
        implicit none
        INP_SECTION(1:N_keyword)=['pump   ',   &
                                & 'pipePR ',   &
                                & 'reactor',   &        
                                & 'pipeRI ',   &
                                & 'IHX    ',   &
                                & 'pipeIP ',   &
								& 'solid  ',   &
								& 'time   '     ]
    end subroutine Set_section_keyword
	
	function is_keyword(input,key) result(is_true)
        implicit none
        character(len=MAX_WORD_LEN),intent(in)::input(:)
        character(len=MAX_WORD_LEN),intent(in)::key
        !local
        logical::is_true
        integer::i,list
        list=size(input)
        is_true=.FALSE.
        do i=1,list,1
            if(trim(adjustl(input(i)))==trim(adjustl(key))) then
                is_true=.TRUE.
                exit
            endif
        end do
    end function is_keyword
end module Imp_inputcard