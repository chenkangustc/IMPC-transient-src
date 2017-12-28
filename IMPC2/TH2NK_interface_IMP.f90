!***************************************************************************************
! test if
!  module: debugassembly
!
!  PURPOSE:  Entry point for the console application.
!
!  pow(na,nr),fq_core(na,nr)     平均功率密度，功率峰因子
!  nr = SIZE(assembly, dim=1)    径向的组件数目
!  na = SIZE(assembly, dim=2)    轴向的节块数目，原输入变量不需要操作赋值的另外用局部变量表达
! 
!***************************************************************************************
    module TH2NK_interface_IMP
	 use constants
     use imp_re_input_global
     use imp_assm_global
     use imp_driving_pre_process
     use imp_driving_post_process
     use imp_power_header
     use imp_single_channel
    implicit none

     !real(KREAL),allocatable::power(:,:),fq_core(:,:)
     !integer M,N,i,j
    contains
    subroutine Perform_TH_imp(transient_flag, assembly, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)
		logical, intent(in)      :: transient_flag                              ! .TRUE. --transient
        real(KREAL), intent(in)  :: assembly(:, :)                              ! (zone, layer), in W, 各组件功率;
        real(KREAL), intent(in out)  :: Tfuel(:, :)                             ! (nr, na), in K, 各组件平均燃料温度;
        real(KREAL), intent(in out)  :: Tcoolant(:, :)                          ! (nr, na), in K, 各组件平均冷却剂温度;
        real(KREAL), intent(in out)  :: Rhocoolant(:, :)                        ! (nr, na), in Kg/m^3, 各组件平均冷却剂密度;
        real(KREAL), intent(in out)  :: max_Tfuel                               ! in K, 最热组件最大燃料温度;
        real(KREAL), intent(in out)  :: max_Tcoolant                            ! in K, 最热组件最大冷却剂温度;
        real(KREAL), intent(in out)  :: min_Rhocoolant                          ! in Kg/m^3, 最热组件最大冷却剂密度;
        real(KREAL), intent(in)  :: last                                        ! in s, 上一时间点
        real(KREAL), intent(in)  :: current                                     ! in s, 当前时间点
        real(KREAL), intent(in out)  :: toutlet                                 ! in K, 冷却剂出口平均温度
    
        !local
        REAL(KREAL)  :: last_, current_
        real(KREAL), allocatable  :: power(:, :)
        real(KREAL), allocatable  :: fq_core(:, :)
        integer  :: nr, na, npin,M,N
        !integer  :: ir, ia, ipin, itype
        !integer  :: i_allocate
		integer  :: i,j,k
        
        last_ = last
        current_ = current
        nr = SIZE(assembly, dim=1)                                              ! 径向的组件数目zone
        na = SIZE(assembly, dim=2)                                              ! 轴向的节块数目layer
        
    	M=size(assm1(1)%thermal%temperature,dim=1)
        N=size(assm1(1)%thermal%temperature,dim=2)
		!allocate(assembly(nr,na),Tfuel(nr,na),Tcoolant(nr,na),Rhocoolant(nr,na))
        allocate(power(M,N),fq_core(M,N))
		fq_core=1.0D0
        power=0.0 

		do i=1,nr,1
		  do j=1,na,1
		    imp_pow(i,j)=assembly(i,j)!W
		  enddo
		enddo
 
	 !do k=1,assm1%mesh%n_zone,1
	 !assm zone=1 组件1
	 !
	 do i=1,nr,1!zone
	   do j=1,assm1(i)%mesh%ny,1!dy
          do k=1,N,1
           !print*,'assembly=',assembly(i,j+assm1(i)%mesh%layer_bottom),'height=',assm1(i)%geom%height(j)
           if(k<=assm1(i)%mesh%Nf) power(j,k)=assembly(i,j+assm1(i)%mesh%layer_bottom)/(assm1(i)%geom%N_fuelpin*assm1(i)%geom%height(j)*3.14159*assm1(i)%geom%pellet**2)
          enddo
       enddo
	   if (transient_flag)  then
            call driving_imp_transient(assm1(i),power, fq_core,last_, current_)
        else
            call driving_imp_steady(assm1(i),power,fq_core)
        end if
 
	   do j=1,assm1(i)%mesh%Ny,1
	   Tfuel(i,j+assm1(i)%mesh%layer_bottom)=assm1(i)%thermal%temperature(j,1)
	   Tcoolant(i,j+assm1(i)%mesh%layer_bottom)=assm1(i)%thermal%temperature(j,N)
       !print*,'Tcoolant=',Tcoolant(i,j+assm1(i)%mesh%layer_bottom)
	   Rhocoolant(i,j+assm1(i)%mesh%layer_bottom)=assm1(i)%property%rho(j,N)
	   enddo
	 enddo
	  open(6,file='.\output\Tfuel.txt')
      write(6,*) Tfuel
      !write(2,*) assm1%pow%power
      close(6)  
	  open(7,file='.\output\Tcoolant.txt')
      write(7,*) Tcoolant
      !write(2,*) assm1%pow%power
      close(7) 
	  
	  if (allocated(power))       deallocate(power)
      if (allocated(fq_core))     deallocate(fq_core)
      
     !*********************************************
     !power come from other data,so it should be an interface in place with the data

     !do while(timer1%ctime<timer1%ttotal) 
     ! timer1%ctime=timer1%ctime+timer1%dt
     ! call update_power(power,fq_core,timer1%ltime,timer1%ctime)
     ! call driving_imp_transient(assm1,power, fq_core,timer1%ltime,timer1%ctime)
     ! call timer1%record(assm1%th_boundary%T%outlet,assm1%th_boundary%u%inlet,power(1,1))
     !print*,'ctime=',timer1%ctime
     ! timer1%ltime=timer1%ctime
     !enddo
   
     !call Run_output() 
     
	 !read(*,*)
     end subroutine Perform_TH_imp
    end module TH2NK_interface_IMP

