module  Imp_pump_header
	use constants
	use imp_property
	implicit none
    
    type pump
		real(KREAL)::I
		!OperatingCondition
		real(KREAL)::yita
		real(KREAL)::omega
		real(KREAL)::Q
		real(KREAL)::omegae
		real(KREAL)::Qe
		real(KREAL)::He!额定扬程
		!material
		real(KREAL)::rho
		!thermal
		real(KREAL)::T
        !branch
        integer::Nbranch
        !rotate
        integer::Ntime
        real,allocatable::rotate(:,:)
        !control
        logical::is_table
	  contains
		procedure,public::init=>init_pump
		procedure,public::alloc=>alloc_pump
		procedure,public::free=>free_pump
	end type pump
	private::init_pump
  contains
    subroutine init_pump(this)
		implicit none
		class(pump),intent(in out)::this
		!local
        real(KREAL)::temperature
		this%Q=this%Qe
		this%omega=this%omegae
		temperature=this%T
		this%rho=get_density_Na(temperature)
	end subroutine init_pump   
    
    subroutine free_pump(this)
        class(pump),intent(in out)::this
        if(allocated(this%rotate)) deallocate(this%rotate)
    end subroutine free_pump
    
    subroutine alloc_pump(this)
        class(pump),intent(in out)::this
        call this%free()
        allocate(this%rotate(2,this%Ntime))
    end subroutine alloc_pump
end module  Imp_pump_header