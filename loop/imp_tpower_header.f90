module imp_tpower_header
    use constants
    implicit none
    type tpower
        real(KREAL),allocatable::pow(:,:)
        integer::Ntime
    end type tpower
end module imp_tpower_header