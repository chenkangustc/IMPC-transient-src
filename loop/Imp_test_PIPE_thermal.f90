module Imp_test_PIPE_thermal
	use Imp_PIPE_header
	use Imp_loop_global
    use Imp_timer_global
	type(pipe)::pipet
	contains
	subroutine test_PIPE_thermal()
		implicit none
		integer::i
		real(KREAL)::Q,Tfin
		real(KREAL)::last,current,dt
		real(KREAL)::LL
        !=======================================================
        open(unit=11,file='./temperature.txt',status='replace')
        write(unit=11,fmt="('  height','  Tf  ',' Tshell')")
        open(unit=12,file='./timelist.txt',status='replace')
        write(unit=12,fmt="('  time','    Qp  ','  Tpin   ','  Tpout ')")
        !=======================================================
		pipet=pipeIP
		dt=timer1%ttotal/timer1%Nt
		last=0.0
		current=0.0
		
		call pipet%thcals()
        write(unit=12,fmt="(F6.1,' ',3F8.2)") current,pipet%Q,pipet%Tfin,pipet%Tfout
        
		do i=1,timer1%Nt,1
		   current=last+dt
		   !Q=pipet%Q*(1-real(i)/timer1%Nt)
           pipet%Q=pipet%Q
		   pipet%Tfin=pipet%Tfin*abs(sin(8*PI/timer1%ttotal*current))
		   call pipet%thcalt(last,current)
            write(unit=12,fmt="(F6.1,' ',3F8.2)") current,Q,pipet%Tfin,pipet%Tfout
		   last=current
        enddo
        
		LL=0.0
        do i=0,pipet%Ny+1,1
          if(i==0)then
              write(unit=11,fmt="('  0.0000',F6.1)") pipet%Tfin
          elseif(i==pipet%Ny+1)then
              write(unit=11,fmt="(F8.4,F6.1)") pipet%Ltotal,pipet%Tfout         
          else
			 if(i==1)then 
				LL=pipet%length(i)/2.0
			 else
				LL=LL+(pipet%length(i)+pipet%length(i-1))/2.0
             end if
              write(unit=11,fmt="(F8.4,2F6.1)") LL,pipet%Tf(i),pipet%Ts(i)
          end if
        enddo          
     close(11)
     close(12)
	end subroutine test_PIPE_thermal
end module Imp_test_PIPE_thermal