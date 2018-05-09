module Imp_loop_global
	use Imp_IHX_header
	use Imp_pipe_header
	use Imp_pump_header
    use Imp_coremodle_header
	implicit none
        logical::is_tNK2TH
        logical::is_THonly
		type(coremodle)::core
		type(IHX)::IHX1
		type(pipe)::PipeRI
		type(pipe)::PipeIP
		type(pipe)::PipePR
		type(pump)::pump1
end module Imp_loop_global