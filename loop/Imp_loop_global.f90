module Imp_loop_global
    use constants
	use Imp_IHX_header
	use Imp_pipe_header
	use Imp_pump_header
    use Imp_coremodle_header
    use imp_tpower_header
	implicit none
        logical::is_tNK2TH
        logical::is_THonly
        logical::is_natural
        real(KREAL)::Qloop
		type(coremodle)::core
		type(IHX)::IHX1
		type(pipe)::PipeRI
		type(pipe)::PipeIP
		type(pipe)::PipePR
		type(pump)::pump1
        !is_THonly
        type(tpower)::tpower1
        !post
        integer::izoneTdis
end module Imp_loop_global