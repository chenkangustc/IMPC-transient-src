module loop_reactor
    type reactor
        ! core
        integer Nzone
        real(KREAL)::Tcout
        real(KREAL)::Qcore
        !real(KREAL)::Qsa(:)         !Q(zone)组件流量
        !Geom
        real(KREAL):: pellet         !芯块半径  
        real(KREAL):: Bond           !元件气隙厚度  
        real(KREAL):: Cladth         !元件外壳厚度
        real(KREAL):: pitch          !组件外对边距（包含包壳厚度）
        real(KREAL):: pd             !燃料元件PD比
        real(KREAL):: rod            !元件半径
        real(KREAL):: area           !芯块横截面积
        integer Npin                 !燃料pin的个数
        real(KREAL),allocatable:: Height(:) !组件高度（活性区）
        !mesh
        integer nf
        integer ng
        integer ns
        integer ny
        integer layer_bottom
        integer layer_top
        real(KREAL),allocatable::r(:,:)
        real(KREAL),allocatable::z(:,:)
        !hydraulic
        real(KREAL):: fric
        real(KREAL):: areaf
        real(KREAL):: wet
        real(KREAL):: de
        ! property
        real(KREAL),allocatable::rho(:,:)!元件热物性
        real(KREAL),allocatable::shc(:,:)
        real(KREAL),allocatable::ctc(:,:)
        real(KREAL),allocatable::dvs(:,:)
        real(KREAL),allocatable::htc(:)
        !assembly thermal
        real(KREAL),allocatable::power(:,:,:)!W/m3 power(layer,nf+ns+nc+1)
        real(KREAL),allocatable::Temperature(:,:,:) !pvt
        real(KREAL),allocatable::Tfuel(:,:)
        real(KREAL),allocatable::Tcoolant(:,:)
        real(KREAL),allocatable::Tfuel_center(:,:)
        real(KREAL),allocatable::Tfg(:,:)
        real(KREAL),allocatable::Tgs(:,:)
        real(KREAL),allocatable::Tsc(:,:)
        real(KREAL),allocatable::Tfave(:)!Tfave(zone)
        real(KREAL),allocatable::Tcave(:)
        ! init
        real(KREAL):: Ti!初始温度
        !boundary
        real(KREAL),allocatable::Tin(:)
        real(KREAL),allocatable::Tout(:)
        real(KREAL):: Qpin(:)
    end type reactor
    contains
    subroutine driving_TH_core(this,transient_flag,Qin,Tin,assembly,last,current)
        class(reactor),intent(in out)::this
        logical,intent(in)::transient_flag
        real(KREAL),intent(in)::assembly(:,:)!power(zone,layer) W
        real(KREAL),intent(in)::Qin,Tin
        real(KREAL),intent(in),optional::last,current
        !local
        integer  :: i,j,k 
        integer  :: nr,na,M,N
        real(KREAL),allocatable::pow(:,:),fq_core(:,:)
        real(KREAL)::density,flowrate

        fq_core=1.0D0		
        nr = this%Nzone!zone                           
        na = this%Ny!layer                          
        M=size(this%temperature,dim=2)
        N=size(this%temperature,dim=3)		
        allocate(pow(na,N),fq_core(na,N))
        pow=0.0
        fq_core=1.0
        !Tin
        do i=1,nr,1
            this%Tin(i)=Tin
        enddo

        call driving_loop_flowAlloc(Qin)
        
        call powerDis_SAtoPin(assembly)
        
        do i=1,nr,1
            if (transient_flag)then
                call driving_imp_THtransient(i,last,current)
            else
                call driving_imp_THsteady(i)					
            endif
        enddo !nr
        !feedback values like Tout volum ave etc.
        call cal_core_Tave()

    end subroutine driving_TH_core
    subroutine powerDis_SAtoPin(this,SApower)
        class(reactor),intent(in out)::this
        real(KREAL),intent(in)::SApower(:,:)
        !local
        integer::i,j,k
        integer::Nzone,M,N
        Nzone=this%Nzone
        M=this%Ny
        N=this%Nf+this%Ng+this%Ns+1
        do i=1,Nzone,1
            do j=1,M,1!dy
                do k=1,N,1
                    if(k<=this%Nf) this%pow(i,j,k)=assembly(i,j+this%layer_bottom)/(this%Npin*this%height(j)*3.14159*this%pellet**2)
                enddo
            enddo
        enddo
    end subroutine powerDis_SAtoPin
	subroutine driving_loop_flowAlloc(this,Qin)
		class(reactor),intent(in out)::this
        real(KREAL),intent(in)::Qin
		!local
        this%Qcore=Qin
        this%Qpin=this%Qcore/(this%Nzone*this%Npin)
	end subroutine driving_loop_flowAlloc
    subroutine driving_imp_THsteady(this,izone)
        class(reactor),intent(in out)::this
        integer,intent(in)::izone
        ! local
        real(KREAL):: flag,dt
        flag=0.0
        dt=1.0!为保证方程求解dt不为0，无具体意义
        
        call cal_th_temperature_rhoi(izone,flag,dt)
        !cal Tfave & Tcave
        call cal_Tave(izone)
    end subroutine driving_imp_THsteady
    
    subroutine driving_imp_THtransient(this,izone,ltime,ctime)
        type(sys_assembly),intent(in out)::this
        integer,intent(in)::izone
        real(KREAL),intent(in)::ltime
        real(KREAL),intent(in)::ctime
        !local
        real(KREAL):: dt
        real(KREAL):: flag
        flag=1.0
        dt=ctime-ltime
        call cal_th_temperature_rhoi(izone,flag,dt)
        call cal_Tave(izone)
    end subroutine driving_imp_THtransient

    subroutine cal_th_temperature_rhoi(this,izone,flag,dt)
        implicit none
        class(reactor),intent(in out)::this
        integer,intent(in)::izone
        real(KREAL),intent(in)::flag
        real(KREAL),intent(in)::dt
        !local
        real(KREAL):: Xt,Df,Dg,Ds,Dy,erro
        integer  M,N,i,j,k,Num
        real(KREAL)::Qpin,Tin
        real(KREAL),allocatable::Ti(:,:),Tj(:,:),Tk(:,:),Td(:,:)
        real(KREAL),allocatable::rho(:,:),shc(:,:),ctc(:,:),dvs(:,:)
        real(KREAL),allocatable::aw(:,:),ae(:,:),ap(:,:),as(:,:),an(:,:),api(:,:),bs(:,:),q(:,:)
        real(KREAL)::kw,kn,ks 
        associate(Area=this%aflow,&
                  xf=this%pellet ,&
                  xg=this%Bond   ,&
                  xs=this%Cladth ,&
                  Nf=this%Nf     ,&
                  Ng=this%Ng     ,&
                  Ns=this%Ns     ,&
                  Ny=this%Ny                     
                  )
        Tin=this%Tin(izone)
        Qpin=this%Qpin(izone)
        rho=this%rho(izone,:,:)   
        shc=this%shc(izone,:,:)   
        ctc=this%ctc(izone,:,:)   
        dvs=this%dvs(izone,:,:)  
        Xt=Xf+Xg+Xs
        M=Ny+1          
        N=Nf+Ng+Ns+1
        Df=Xf/Nf        
        Dg=Xg/Ng        
        Ds=Xs/Ns      
        allocate(Ti(1:Ny,1:N),Tj(1:Ny,1:N),Tk(1:Ny,1:N),Td(1:Ny,1:N))
        allocate(aw(1:Ny,1:N),ae(1:Ny,1:N),ap(1:Ny,1:N),as(1:Ny,1:N),an(1:Ny,1:N),api(1:Ny,1:N),bs(1:Ny,1:N),q(1:Ny,1:N))
        !cal convection
        call cal_th_convection(izone)
        !set the pow
        Ti=this%Temperature(izone,:,:)
        q=this%power(izone,:,:)
        api=0.0
        Do i=1,Ny,1
            Dy=this%height(i+this%layer_bottom)
            Do j=1,N,1
                if (j==1)then!轴对称边界的控制体
                    aw(i,j)=0.0
                    ae(i,j)=(this%r(i,j)+Df/2.0)*CTC(i,j)/Df
                    as(i,j)=0.0
                    an(i,j)=0.0
                    if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*this%r(i,j)*Df/dt
                    ap(i,j)=ae(i,j)+api(i,j)
                    bs(i,j)=this%r(i,j)*Df*q(i,j)
                elseif (j>1.and.j<Nf)then!fuel内部控制体
                    aw(i,j)=(this%r(i,j)-Df/2.0)*CTC(i,j)/Df
                    ae(i,j)=(this%r(i,j)+Df/2.0)*CTC(i,j)/Df
                    as(i,j)=0.0
                    an(i,j)=0.0
                    if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*this%r(i,j)*Df/dt
                    ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
                    bs(i,j)=this%r(i,j)*Df*q(i,j)
                elseif(j==Nf)then!f-g边界左侧控制体
                    aw(i,j)=(this%r(i,j)-Df/2.0)*CTC(i,j)/Df
                    ae(i,j)=2*(this%r(i,j)+Df/2.0)*(Df+Dg)/(Df/CTC(i,j)+Dg/CTC(i,j+1))/(Df+Dg)
                    as(i,j)=0.0
                    an(i,j)=0.0
                    if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*this%r(i,j)*Df/dt
                    ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
                    bs(i,j)=this%r(i,j)*Df*q(i,j)  
                elseif (j==Nf+1)then!f-g边界右侧控制体
                    aw(i,j)=2*(this%r(i,j)-Dg/2.0)*(df+dg)/(df/CTC(i,j-1)+dg/CTC(i,j))/(df+dg)
                    ae(i,j)=(this%r(i,j)+Dg/2.0)*CTC(i,j)/Dg
                    as(i,j)=0.0
                    an(i,j)=0.0
                    if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*this%r(i,j)*Dg/dt
                    ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
                    bs(i,j)=this%r(i,j)*Dg*q(i,j)
                elseif(j>Nf+1.and.j<Nf+Ng)then!g气隙内部控制体
                    aw(i,j)=(this%r(i,j)-Dg/2.0)*CTC(i,j)/Dg
                    ae(i,j)=(this%r(i,j)+Dg/2.0)*CTC(i,j)/Dg
                    as(i,j)=0.0
                    an(i,j)=0.0
                    if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*this%r(i,j)*Dg/dt
                    ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
                    bs(i,j)=this%r(i,j)*Dg*q(i,j)
                elseif(j==Nf+Ng)then!g-c边界左侧控制体
                    aw(i,j)=(this%r(i,j)-Dg/2.0)*CTC(i,j)/Dg
                    ae(i,j)=2*(this%r(i,j)+Dg/2.0)*(Dg+Ds)/(Dg/CTC(i,j)+Ds/CTC(i,j+1))/(Dg+Ds)
                    as(i,j)=0.0
                    an(i,j)=0.0
                    if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*this%r(i,j)*Dg/dt
                    ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
                    bs(i,j)=this%r(i,j)*Dg*q(i,j)
                elseif(j==Nf+Ng+1)then!g-c边界右侧控制体
                    aw(i,j)=2*(this%r(i,j)-Ds/2.0)*(dg+ds)/(dg/CTC(i,j-1)+ds/CTC(i,j))/(dg+ds)
                    ae(i,j)=(this%r(i,j)+Ds/2.0)*CTC(i,j)/Ds
                    as(i,j)=0.0
                    an(i,j)=0.0
                    if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*this%r(i,j)*Ds/dt
                    ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
                    bs(i,j)=this%r(i,j)*Ds*q(i,j)
                elseif(j>Nf+Ng+1.and.j<Nf+Ng+Ns)then!c包壳内部控制体
                    aw(i,j)=(this%r(i,j)-Ds/2.0)*CTC(i,j)/Ds
                    ae(i,j)=(this%r(i,j)+Ds/2.0)*CTC(i,j)/Ds
                    as(i,j)=0.0
                    an(i,j)=0.0
                    if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*this%r(i,j)*Ds/dt
                    ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
                    bs(i,j)=this%r(i,j)*Ds*q(i,j)
                elseif(j==Nf+Ng+Ns)then!s-fluid边界左侧控制体
                    kw=2*CTC(i,j)*CTC(i,j-1)/(CTC(i,j)+CTC(i,j-1))
                    aw(i,j)=(this%r(i,j)-Ds/2.0)*kw/(this%r(i,j)*Ds**2)
                    ae(i,j)=this%htc(i)*(this%r(i,j)+Ds/2.0)/(this%r(i,j)*Ds)
                    as(i,j)=0.0
                    an(i,j)=0.0
                    if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)/dt
                    ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
                    bs(i,j)=q(i,j)
                elseif(j==Nf+Ng+Ns+1)then!fluid控制体
                    if(i==1)then!流体入口的控制体
                        aw(i,j)=this%htc(i)*2.0*PI*Xt*Dy
                        ae(i,j)=0.0
                        as(i,j)=0.0
                        an(i,j)=0.0
                        if(flag==1.0) api(i,j)=RHO(i,j)*Area*Dy*SHC(i,j)/dt
                        ap(i,j)=aw(i,j)+api(i,j)+Qpin*SHC(i,j)
                        bs(i,j)=Qpin*SHC(i,j)*this%Tin
                    elseif(i==Ny)then!流体出口
                        aw(i,j)=this%htc(i)*2.0*PI*Xt*Dy
                        ae(i,j)=0.0
                        as(i,j)=Qpin*SHC(i,j)
                        an(i,j)=0.0
                        if(flag==1.0) api(i,j)=RHO(i,j)*Area*Dy*SHC(i,j)/dt
                        ap(i,j)=as(i,j)+aw(i,j)+api(i,j)
                        bs(i,j)=0.0
                    else!流体内部，考虑导热
                        kn=2*CTC(i+1,j)*CTC(i,j)/(CTC(i+1,j)+CTC(i,j))
                        ks=2*CTC(i-1,j)*CTC(i,j)/(CTC(i-1,j)+CTC(i,j))
                        aw(i,j)=this%htc(i)*2.0*PI*Xt*Dy
                        ae(i,j)=0.0
                        as(i,j)=Qpin*SHC(i,j)+area*ks/Dy
                        an(i,j)=Area*kn/Dy
                        if(flag==1.0) api(i,j)=RHO(i,j)*Area*Dy*SHC(i,j)/dt
                        ap(i,j)=an(i,j)+as(i,j)+aw(i,j)+api(i,j)
                        bs(i,j)=0.0
                    endif
                endif              
            enddo
        enddo
        !do k=1,5000,1
        Tj=Ti
        erro=1.0
        num=1
        do while(erro.gt.1.0D-6)
        do i=1,Ny,1
            do j=1,N,1
                if(j==1)then
                 Tj(i,j)=(ae(i,j)*Tj(i,j+1)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)
                elseif(j>1.and.j<N)then
                 Tj(i,j)=(aw(i,j)*Tj(i,j-1)+ae(i,j)*Tj(i,j+1)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)
                elseif(j==N) then!fluid
                  if(i==1)then!fluid inlet
                    Tj(i,j)=(aw(i,j)*Tj(i,j-1)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)
                  elseif(i==Ny)then!fluid outlet
                    Tj(i,j)=(aw(i,j)*Tj(i,j-1)+as(i,j)*Tj(i-1,j)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)  
                  else!fluid inner
					Tj(i,j)=(aw(i,j)*Tj(i,j-1)+as(i,j)*Tj(i-1,j)+an(i,j)*Tj(i+1,j)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)
				 endif
                endif
				if(num.gt.1) Td(i,j)=abs((Tj(i,j)-Tk(i,j))/Tk(i,j))
		    enddo
        enddo
	    if(num.gt.1)  erro=maxval(Td)
	    Tk=Tj
	    num=num+1
        enddo
        do i=1,Ny,1!as for the solid,no need to know the inlet and outlet temperature
            do j=1,N,1
                    this%Temperature(izone,i,j)=Tj(i,j)
            enddo
        enddo
        this%Tout(izone)=this%Temperature(Ny,N)
        end associate
    end subroutine cal_th_temperature_RHOi
    subroutine cal_th_convection(this,izone)
        implicit none
        class(reactor),intent(in out)::this
        real(KREAL),intent(in)::izone
        !lcoal
        integer i,Ny,nr
        real(KREAL):: Qpin,velocity!单纯输入的变量可以用局部变量来替换
        real(KREAL),allocatable::rho(:),dvs(:),shc(:),ctc(:)
        Ny=this%Ny
        nr=this%nf+this%ng+this%ns
        allocate(rho(Ny),dvs(Ny),shc(Ny),ctc(Ny))
        associate(de=this%de        ,&
                  Area=this%aflow   ,&
                  wet=this%wet      
                  )
        Qpin=this%Qpin(izone)
        rho=this%rho(izone,:,nr+1)
        dvs=this%dvs(izone,:,nr+1)
        shc=this%shc(izone,:,nr+1)
        ctc=this%ctc(izone,:,nr+1)
        do i=1,Ny,1
            velocity=Qpin/Area
            call get_convection(De,Area,wet,RHO(i),velocity,DVS(i),SHC(i),CTC(i),this%htc(i))!DVS(i,N)动力粘度 Pa*s
        enddo
        end associate
    end subroutine cal_th_convection
    subroutine cal_Tave(this,izone)
        implicit none
        class(reactor),intent(in out)::this
        integer,intent(in)::izone
        integer::j,k,N
        real(KREAL)::volumn,TVtotal
        real(KREAL)::TLtotal,Ltotal
        real(KREAL)::dr
        !计算元件的芯块平均温度和冷却剂平均温度
        N=this%Nf+this%Ng+this%Ns+1
        dr=this%pellet/this%Nf
        !fuel
        volumn=0.0
        TVtotal=0.0
        do j=1,this%Ny,1
            do k=1,this%Nf,1 !rod average					
                if (k==1) then
                    TVtotal=TVtotal+this%temperature(izone,j,k)*PI*dr**2*this%height(j)
                    volumn=volumn+PI*dr**2*this%height(j)
                else
                    TVtotal=TVtotal+this%temperature(izone,j,k)*PI*((k*dr)**2-((k-1)*dr)**2)*this%height(j)
                    volumn=volumn+PI*((k*dr)**2-((k-1)*dr)**2)*this%height(j)
                endif
            enddo
        enddo
        this%Tfave(izone)=TVtotal/volumn
        !Tcave
        TLtotal=0.0
        Ltotal=0.0
        do j=1,this%Ny,1
            TLtotal=TLtotal+this%temperature(izone,j,N)*this%height(j)
            Ltotal=Ltotal+this%height(j)
        enddo
        this%Tcave(izone)=TLtotal/Ltotal
        !Tcoolant
        do j=1,this%Ny,1
            this%Tcoolant(izone,j)=this%temperature(izone,j,N)
        enddo
        !Tfuel
        do j=1,this%Ny,1
            TVtotal=0.0
            volumn=0.0
            do k=1,this%Nf,1 !rod average					
                if (k==1) then
                    TVtotal=TVtotal+this%temperature(izone,j,k)*PI*dr**2*this%height(j)
                    volumn=volumn+PI*dr**2*this%height(j)
                else
                    TVtotal=TVtotal+this%temperature(izone,j,k)*PI*((k*dr)**2-((k-1)*dr)**2)*this%height(j)
                    volumn=volumn+PI*((k*dr)**2-((k-1)*dr)**2)*this%height(j)
                endif
            enddo
            this%Tfuel(izone,j)=TVtotal/volumn
        enddo
    end subroutine cal_Tave
    
    subroutine cal_core_Tave(this)
        class(reactor),intent(in out)::this
        ! local
        integer::i
        integer::nr
        real(KREAL)::TQtotal,Qtotal
        nr=this%Nzone
        TQtotal=0.0
        Qtotal=0.0
        do i=1,nr,1
            TQtotal=TQtotal+this%Tout(i)*this%Qpin(i)
            Qtotal=Qtotal+this%Qpin(i)
        enddo
        this%Tcout=TQtotal/Qtotal
    end subroutine cal_core_Tave
end module loop_reactor