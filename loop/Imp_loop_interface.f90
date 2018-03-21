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
    module TH2NK_interface_loop
	 use constants 
     use imp_assm_global
     use Imp_cal_loop
    implicit none
    !real(KREAL),allocatable::power(:,:),fq_core(:,:)
    !integer M,N,i,j
    contains
    subroutine Perform_TH_loop(transient_flag, assembly, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)
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
		real(KREAL)  :: volumn,TVtotal,dr!used to calculate the average temperature of fuel
		real(KREAL)  :: TAoutlet_total,A_total!used to calculate the core average outlet temperature
		real(KREAL)	 :: xs,xg,xf !used to calculate the Tsurface
		real(KREAL)  :: max_T_inner,max_T_outer
		
		integer  :: nf,ng,ns,nRadial,ny
        integer  :: nr, na, npin,M,N
		integer  :: i,j,k
        
        last_ = last
        current_ = current
        nr = SIZE(assembly, dim=1)                                              ! 径向的组件数目zone
        na = SIZE(assembly, dim=2)                                              ! 轴向的节块数目layer     
		!=========================================================
		!if (transient_flag)  then
  !          call driving_loop_transient(assembly,last_, current_)
		!else
            call driving_loop_steady(assembly)
		!end if
	    !==========================================================
		!热工反馈
		dr=assm1(i)%geom%pellet/assm1(i)%mesh%Nf
		do j=1,assm1(i)%mesh%Ny,1
			volumn=0.0
			TVtotal=0.0
			do k=1,assm1(i)%mesh%Nf,1 !rod average		
				if (k==1) then
					TVtotal=TVtotal+assm1(i)%thermal%temperature(j,k)*3.14*(k*dr)**2*assm1(i)%geom%height(j)
					volumn=volumn+3.14*(k*dr)**2*assm1(i)%geom%height(j)
				else
					TVtotal=TVtotal+assm1(i)%thermal%temperature(j,k)*3.14*((k*dr)**2-((k-1)*dr)**2)*assm1(i)%geom%height(j)
					volumn=volumn+3.14*((k*dr)**2-((k-1)*dr)**2)*assm1(i)%geom%height(j)
				endif
			enddo
			Tfuel(i,j+assm1(i)%mesh%layer_bottom)=TVtotal/volumn
			Tcoolant(i,j+assm1(i)%mesh%layer_bottom)=assm1(i)%thermal%temperature(j,N)
			Rhocoolant(i,j+assm1(i)%mesh%layer_bottom)=assm1(i)%property%rho(j,N)
		enddo
			
		TAoutlet_total=0.0
		A_total=0.0
		do i=1,nr,1
			!volumn average
			TAoutlet_total=TAoutlet_total+assm1(i)%th_boundary%T%outlet*assm1(i)%hydrau%aflow*assm1(i)%th_boundary%u%outlet
			A_total=A_total+assm1(i)%hydrau%aflow*assm1(i)%th_boundary%u%outlet
		enddo
		toutlet=TAoutlet_total/A_total
		max_Tcoolant=0.0
		max_Tfuel=0.0
		!calculate max_Tcoolant max_Tfuel
		do i=1,nr,1
			do j=1,assm1(i)%mesh%Ny,1
				if(Tfuel(i,j)>max_Tfuel) 		max_Tfuel=Tfuel(i,j)
				if(Tcoolant(i,j)>max_Tcoolant)	max_Tcoolant=Tcoolant(i,j)
			enddo
		enddo
		!calculate the surface temperature
		do i=1,nr,1
			Nf=assm1(i)%mesh%Nf
			Ng=assm1(i)%mesh%Ng
			Ns=assm1(i)%mesh%Ns
			Ny=assm1(i)%mesh%Ny
			Nradial=Nf+Ng+Ns+1
			xf=assm1(i)%geom%pellet
			xg=assm1(i)%geom%bond
			xs=assm1(i)%geom%cladth		
			do j=1,assm1(i)%mesh%Ny,1
				assm1(i)%thermal%Tcoolant(j)=Tcoolant(i,j)
				assm1(i)%thermal%Tfuel(j)=Tfuel(i,j)
				assm1(i)%thermal%Tfuel_center(j)=assm1(i)%thermal%temperature(j,1)
				assm1(i)%thermal%Tfg(j)=(assm1(i)%property%ctc(j,Nf)*(Xg/Ng)*assm1(i)%thermal%temperature(j,Nf)+assm1(i)%property%ctc(j,Nf+1)*(Xf/Nf)*assm1(i)%thermal%temperature(j,Nf+1))/(assm1(i)%property%ctc(j,Nf)*(Xg/Ng)+assm1(i)%property%ctc(j,Nf+1)*(Xf/Nf))!芯块外边界
				assm1(i)%thermal%Tgs(j)=(assm1(i)%property%ctc(j,Nf+Ng)*(Xs/Ns)*assm1(i)%thermal%temperature(j,Nf+Ng)+assm1(i)%property%ctc(j,Nf+Ng+1)*(Xg/Ng)*assm1(i)%thermal%temperature(j,Nf+Ng+1))/(assm1(i)%property%ctc(j,Nf+Ng)*(Xs/Ns)+assm1(i)%property%ctc(j,Nf+Ng+1)*(Xg/Ng))!包壳内边界
				assm1(i)%thermal%Tsc(j)=(assm1(i)%property%htc(j)*assm1(i)%thermal%temperature(j,Nradial)+2*assm1(i)%property%ctc(j,Nradial-1)/(Xs/Ns)*assm1(i)%thermal%temperature(j,Nradial-1))/(assm1(i)%property%htc(j)+2*assm1(i)%property%ctc(j,Nradial-1)/(Xs/Ns))!包壳外边界
			enddo
		enddo!surface zone end
		!write current_,max_Tcoolant,toutlet,max_Tfuel,max_T_inner,max_T_outer
		max_T_inner=0.0
		max_T_outer=0.0
		do i=1,nr,1
			do j=1,assm1(i)%mesh%Ny,1
				if(max_T_inner<assm1(i)%thermal%Tgs(j)) max_T_inner=assm1(i)%thermal%Tgs(j)
				if(max_T_outer<assm1(i)%thermal%Tsc(j)) max_T_outer=assm1(i)%thermal%Tsc(j)
			enddo
        enddo
	end subroutine Perform_TH_loop
end module TH2NK_interface_loop

