module Imp_inputcard
    use imp_re_input_global
	use Imp_loop_global
	use Imp_timer_global
    use,intrinsic::ISO_FORTRAN_ENV
    use constants
	implicit none
	integer::file_i,file_o,file_t,file_maxT,file_aveT,file_disT
	integer,parameter,private::N_keyword=13
    integer,parameter,private::MAX_REAL_PARAMETER=200
	integer,parameter,private::MAX_INT_PARAMETER=50
	integer,parameter,private::MAX_LOGICAL_PARAMETER=10
    character(len=MAX_WORD_LEN),parameter::FILE_IN='./src/loopinput.case'
    character(len=MAX_WORD_LEN),parameter::FILE_OUT='loopoutput.txt'
	character(len=MAX_WORD_LEN),parameter::FILE_LTIME='looptimelist.txt'
	character(len=MAX_WORD_LEN),parameter::FILE_MAXTIME='maxT.timelist'
	character(len=MAX_WORD_LEN),parameter::FILE_AVETIME='aveT.timelist'
	character(len=MAX_WORD_LEN),parameter::FILE_DIS='Tdis.txt'
    character(len=MAX_WORD_LEN)::INP_SECTION(N_keyword)
    
	
	contains
    subroutine Set_section_keyword()
        implicit none
        INP_SECTION(1:N_keyword)=['pump   ',   &
                                & 'control ',   &
                                & 'pipePR ',   &
                                & 'pinflow',   &        
                                & 'pipeRI ',   &
                                & 'IHX    ',   &
                                & 'pipeIP ',   &
                                & 'pingeom',  &
                                & 'pinmesh ',  &
                                & 'axil    ',  &
                                & 'height  ',  &
                                & 'rotate  ',  &
								& 'time   '     ]
    end subroutine Set_section_keyword
    
	subroutine driving_input_read()
		implicit none
		!local
		integer::io_error
		real::dummy_real(MAX_REAL_PARAMETER)
		integer::dummy_int(MAX_INT_PARAMETER)
        logical::dummy_logical(MAX_LOGICAL_PARAMETER)
		character(len=MAX_WORD_LEN)::aline
		character(len=MAX_WORD_LEN)::section_name,keyword
		! Variables
		call set_section_keyword()
        open(newunit=file_o,file=FILE_OUT,status='replace',action='write',iostat=io_error)
		open(newunit=file_t,file=FILE_LTIME,status='replace',action='write',iostat=io_error) 
		!write(unit=file_t,fmt="(F6.1,' ',F10.1,8F8.2)") current,powinput,Qloop,coreTin,coreTout,IHX1%Tpin,IHX1%Tpout,IHX1%Qs,IHX1%Tsin,IHX1%Tsout		
		write(unit=file_t,fmt="('   time','    pow','    flowrate',' coreTin',' coreTout',' IHXTin',' IHXTout',' IHXQs','  IHXTsin',' IHXTsout')")
		open(newunit=file_maxT,file=FILE_MAXTIME,status='replace',action='write',iostat=io_error)
        write(unit=file_maxT,fmt="('  time','  maxFuel','  maxCool',' maxTinner',' maxTouter')")
		open(newunit=file_aveT,file=FILE_AVETIME,status='replace',action='write',iostat=io_error)
		open(newunit=file_disT,file=FILE_DIS,status='replace',action='write',iostat=io_error)	
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
					case('control')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_logical(1)
                    pump1%is_table=dummy_logical(1)
					
                    case('pump')
                    if(pump1%is_table==.FALSE.)then
                        read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:5),dummy_int(1)
                        pump1%I=dummy_real(1)
                        pump1%He=dummy_real(2)
                        pump1%Qe=dummy_real(3)
                        pump1%omegae=dummy_real(4)
                        pump1%yita=dummy_real(5)
                        pump1%Nbranch=dummy_int(1)
                    else
                        read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1),dummy_int(1:2)
                        pump1%Qe=dummy_real(1)
                        pump1%Ntime=dummy_int(1)
                        pump1%Nbranch=dummy_int(2)
                    endif
					
					case('pipePR')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:6),dummy_int(1)
					pipePR%ltotal=dummy_real(1)
					pipePR%Rtube=dummy_real(2)
					pipePR%thicks=dummy_real(3)
					pipePR%theta=dummy_real(4)
					pipePR%Q=dummy_real(5)
					pipePR%Ti=dummy_real(6)
					pipePR%Ny=dummy_int(1)
					
					case('pinflow')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1:3),dummy_real(1)
                    core%Nflow=dummy_int(1)
                    core%Nflowsemi=dummy_int(2)
                    core%Nsplit=dummy_int(3)
                    core%sigmaPass=dummy_real(1)
					! core%ltotal=dummy_real(2)
					! core%Rtube=dummy_real(3)
					! core%thicks=dummy_real(4)
					! core%theta=dummy_real(5)
					! core%Q=dummy_real(6)
					! core%Ti=dummy_real(7)
					! core%Ny=dummy_int(4)
					
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
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:10),dummy_int(1:2)
					IHX1%Lsingle=dummy_real(1)
					IHX1%Rtube=dummy_real(2)
					IHX1%thickt=dummy_real(3)
					IHX1%Plength=dummy_real(4)
					IHX1%AreaP=dummy_real(5)
					IHX1%thickv=dummy_real(6)
					IHX1%Qp=dummy_real(7)
					IHX1%Qs=dummy_real(8)
                    IHX1%Tsin=dummy_real(9)
					IHX1%Ti=dummy_real(10)
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

                    case('pingeom')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:7),dummy_int(1:2)
                    reInputdata%npin = dummy_int(1)
                    reInputdata%nFuelPin = dummy_int(2)
                    reInputdata%xf = dummy_real(1)*0.001D0
                    reInputdata%xg = dummy_real(2)*0.001D0
                    reInputdata%xs = dummy_real(3)*0.001D0
                    reInputdata%pd = dummy_real(4)
                    reInputdata%f=dummy_real(5)
                    reInputdata%Tin=dummy_real(6)
                    reInputdata%Ti=dummy_real(7)
                    
                    case('pinmesh')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1:6)
                    ! this%acf,this%height,this%nf,this%ng,this%ns,this%f,this%Tin,this%pout,
                    ! this%uin,this%pin,this%Ti,this%ui,this%pi,this%alpha,this%sigma
                    reInputdata%nf=dummy_int(1)
                    reInputdata%ng=dummy_int(2)
                    reInputdata%ns=dummy_int(3)
                    reInputdata%ny=dummy_int(4)
                    reInputdata%ny_bottom=dummy_int(5)
                    reInputdata%ny_top=dummy_int(6)
					case('time')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1),dummy_int(1)
					timer1%ttotal=dummy_real(1)
					timer1%Nt=dummy_int(1)
				end select
			end if
			!print*,trim(aline)
			!print*,trim(adjustl(section_name))
		end do
        close(file_i)
	end subroutine driving_input_read
    
	subroutine driving_input_read_after()
        !local
        integer::io_error
		real::dummy_real(MAX_REAL_PARAMETER)
		integer::dummy_int(MAX_INT_PARAMETER)
		character(len=MAX_WORD_LEN)::aline
		character(len=MAX_WORD_LEN)::section_name,keyword
        
        open(newunit=file_i,file=FILE_IN,status='old',action='read',iostat=io_error)   		
        !read(unit=file_i,fmt='(A)',iostat=io_error) aline
		do
			read(unit=file_i,fmt='(A)',iostat=io_error) aline
            if(io_error==IOSTAT_END) exit
			read(unit=aline,fmt=*,iostat=io_error) section_name       
    
			if(is_keyword(INP_SECTION,section_name)) then
				select case(trim(adjustl(section_name)))
					case('pinflow')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1:3),dummy_real(1),dummy_int(4:core%Nflow+core%Nflowsemi+3)
                    core%fzone=dummy_int(4:core%Nflow+core%Nflowsemi+3)
                    case('height')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:reInputdata%ny)
                    reInputdata%height=dummy_real*0.01D0
                    case('rotate')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:2*pump1%Ntime)
                    pump1%rotate(1,:)=dummy_real(1:pump1%Ntime)!time
                    pump1%rotate(2,:)=dummy_real(pump1%Ntime+1:2*pump1%Ntime)!rotate
               end select
            endif
        end do
    end subroutine

	
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