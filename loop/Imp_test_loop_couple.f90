module Imp_test_loop_couple
	use Imp_IHX_header
	use Imp_pipe_header
	use Imp_pump_header
    use Imp_coremodle_header
	use Imp_timer_global
	implicit none
		type(coremodle)::coret
		type(IHX)::IHXt
		type(pipe)::PipeRIt
		type(pipe)::PipeIPt
		type(pipe)::PipePRt
		type(pump)::pumpt
	contains
	subroutine test_loop_couple()
		implicit none
		integer::i
		real(KREAL)::last,current,dt
		coret=core
		IHXt=IHX1
		PipeRIt=PipeRI
		PipeIPt=PipeIP
		PipePRt=PipePR
		pumpt=pump1
		call
		
	end subroutine test_loop_couple
end module Imp_test_loop_couple