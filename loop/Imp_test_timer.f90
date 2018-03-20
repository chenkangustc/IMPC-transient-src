module Imp_test_timer
    use constants
    type ImpTimer
        real(KREAL)::ttotal
        integer::Nt
        real(KREAL)::dt
        real(KREAL)::current
        real(KREAL)::last
    end type ImpTimer
  
end module Imp_test_timer