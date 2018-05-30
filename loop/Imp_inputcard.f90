module Imp_inputcard
    use imp_re_input_global
	use Imp_loop_global
	use Imp_timer_global
    use,intrinsic::ISO_FORTRAN_ENV
    use constants
	implicit none
	integer::file_i,file_o,file_t,file_maxT,file_aveT,file_disT
	integer,parameter,private::N_keyword=25
    integer,parameter,private::MAX_REAL_PARAMETER=700
	integer,parameter,private::MAX_INT_PARAMETER=700
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
                                & 'assembly',   &        
                                & 'pipeRI ',   &
                                & 'IHX    ',   &
                                & 'pipeIP ',   &
                                & 'pingeom',  &
                                & 'pinmesh ',  &
                                & 'axil    ',  &
                                & 'height  ',  &
                                & 'rotate  ',  &
                                & 'power   ',  &
                                & 'IHXsTin ',  &
                                & 'IHXsflow',  &
                                & 'SAradiau',  &
                                & 'SAtype  ',  &
                                & 'Tb_RI  ',  &
                                & 'Tb_IP  ',  &
                                & 'Tb_PR  ',  &
                                & 'Bq_RI  ',  &
                                & 'Bq_IP  ',  &
                                & 'Bq_PR  ',  &
                                & 'IMPCpost',  &
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
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_logical(1:5)
                    is_THonly=dummy_logical(1)
                    is_tNK2TH=dummy_logical(2)
                    pump1%is_table=dummy_logical(3)
                    IHX1%is_flowtable=dummy_logical(4)
                    IHX1%is_Tintable=dummy_logical(5)
                    
					
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
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:8),dummy_logical(1),dummy_int(1)
					pipePR%ltotal=dummy_real(1)
					pipePR%Rtube=dummy_real(2)
					pipePR%thicks=dummy_real(3)
					pipePR%theta=dummy_real(4)
					pipePR%Q=dummy_real(5)
					pipePR%Ti=dummy_real(6)
					pipePR%fric=dummy_real(7)
					pipePR%K=dummy_real(8)
					pipePR%is_Tb=dummy_logical(1)
					pipePR%Ny=dummy_int(1)
					
					case('assembly')!径向描述
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1:4),dummy_real(1)
                    core%Nzone=dummy_int(1)
                    core%Nflow=dummy_int(2)
                    core%Nflowsemi=dummy_int(3)
                    core%Nsplit=dummy_int(4)
                    core%sigmaPass=dummy_real(1)
					! core%ltotal=dummy_real(2)
					! core%Rtube=dummy_real(3)
					! core%thicks=dummy_real(4)
					! core%theta=dummy_real(5)
					! core%Q=dummy_real(6)
					! core%Ti=dummy_real(7)
					! core%Ny=dummy_int(4)
					
					case('pipeRI')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:8),dummy_logical(1),dummy_int(1)
					pipeRI%ltotal=dummy_real(1)
					pipeRI%Rtube=dummy_real(2)
					pipeRI%thicks=dummy_real(3)
					pipeRI%theta=dummy_real(4)
					pipeRI%Q=dummy_real(5)
					pipeRI%Ti=dummy_real(6)
                    pipeRI%fric=dummy_real(7)
					pipeRI%K=dummy_real(8)
					pipeRI%is_Tb=dummy_real(1)
					pipeRI%Ny=dummy_int(1)
					
					case('IHX')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:12),dummy_int(1:2)
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
					IHX1%fricp=dummy_real(11)
					IHX1%kricp=dummy_real(12)
					IHX1%Ntube=dummy_int(1)
					IHX1%N=dummy_int(2)
					
					case('pipeIP')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:8),dummy_logical(1),dummy_int(1)
					pipeIP%ltotal=dummy_real(1)
					pipeIP%Rtube=dummy_real(2)
					pipeIP%thicks=dummy_real(3)
					pipeIP%theta=dummy_real(4)
					pipeIP%Q=dummy_real(5)
					pipeIP%Ti=dummy_real(6)
                    pipeIP%fric=dummy_real(7)
					pipeIP%K=dummy_real(8)
					pipeIP%is_Tb=dummy_logical(1)
					pipeIP%Ny=dummy_int(1)

                    case('pingeom')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:8),dummy_int(1:2)
                    reInputdata%npin = dummy_int(1)
                    reInputdata%nFuelPin = dummy_int(2)
                    reInputdata%xf = dummy_real(1)*0.001D0
                    reInputdata%xg = dummy_real(2)*0.001D0
                    reInputdata%xs = dummy_real(3)*0.001D0
                    reInputdata%pd = dummy_real(4)
                    reInputdata%f=dummy_real(5)
                    reInputdata%K=dummy_real(6)
                    reInputdata%Tin=dummy_real(7)
                    reInputdata%Ti=dummy_real(8)
                    
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
                    
					case('power')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1)
                    tpower1%Ntime=dummy_int(1)
                    
                    case('IHXsflow')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1)
                    IHX1%Nftime=dummy_int(1)
                    
                    case('IHXsTin')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1)
                    IHX1%NTtime=dummy_int(1)
                    
					case('SAtype')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1)
                    reInputdata%Ntype=dummy_int(1)
                    
					case('time')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1),dummy_int(1)
					timer1%ttotal=dummy_real(1)
					timer1%Nt=dummy_int(1)
                    
                    case('IMPCpost')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1)
                    izoneTdis=dummy_int(1)
				end select
			end if
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
					case('assembly')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1:4),dummy_real(1),dummy_int(5:core%Nflow+core%Nflowsemi+4)
                    core%fzone=dummy_int(5:core%Nflow+core%Nflowsemi+4)
                    case('height')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:reInputdata%ny)
                    reInputdata%height=dummy_real*0.01D0
                    case('rotate')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:2*pump1%Ntime)
                    pump1%rotate(1,:)=dummy_real(1:pump1%Ntime)!time
                    pump1%rotate(2,:)=dummy_real(pump1%Ntime+1:2*pump1%Ntime)!rotate
                    case('power')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1),dummy_real(1:2*tpower1%Ntime)
                    tpower1%pow(1,:)=dummy_real(1:tpower1%Ntime)
                    tpower1%pow(2,:)=dummy_real(tpower1%Ntime+1:2*tpower1%Ntime)
                    case('IHXsflow')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1),dummy_real(1:2*IHX1%Nftime)
                    IHX1%flowtable(1,:)=dummy_real(1:IHX1%Nftime)
                    IHX1%flowtable(2,:)=dummy_real(IHX1%Nftime+1:2*IHX1%Nftime)
                    case('IHXsTin')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1),dummy_real(1:2*IHX1%NTtime)
                    IHX1%Tintable(1,:)=dummy_real(1:IHX1%NTtime)
                    IHX1%Tintable(2,:)=dummy_real(IHX1%NTtime+1:2*IHX1%NTtime)
                    case('SAradiau')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1:core%Nzone)
                    core%SAtable(:)=dummy_int(1:core%Nzone)
                    case('SAtype')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_int(1),dummy_real(1:2*reInputdata%Ntype)
                    reInputdata%sa(:)%flowdis=dummy_real(1:reInputdata%Ntype)
                    reInputdata%sa(:)%powdis=dummy_real(reInputdata%Ntype+1:2*reInputdata%Ntype)
                    case('Tb_RI')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:PipeRI%Ny)
                    PipeRI%Tb(:)=dummy_real(1:PipeRI%Ny)
                    case('Tb_IP')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:PipeIP%Ny)                
                    PipeIP%Tb(:)=dummy_real(1:PipeIP%Ny)
                    case('Tb_PR')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:PipePR%Ny)
                    PipePR%Tb(:)=dummy_real(1:PipePR%Ny)
                    case('Bq_RI')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:PipeRI%Ny)
                    PipeRI%Bq(:)=dummy_real(1:PipeRI%Ny)
                    case('Bq_IP')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:PipeIP%Ny)                
                    PipeIP%Bq(:)=dummy_real(1:PipeIP%Ny)
                    case('Bq_PR')
                    read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:PipePR%Ny)
                    PipePR%Bq(:)=dummy_real(1:PipePR%Ny)
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