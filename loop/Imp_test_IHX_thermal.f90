module Imp_test_IHX_thermal
	use Imp_loop_global
	use Imp_timer_global
	use constants
	type(IHX)::IHXt
	contains
	subroutine test_IHX_thermal()
		implicit none
		integer::i
		real(KREAL)::Qp,Qs,Tpin,Tsin
		real(KREAL)::last,current,dt
        !+=======================================================
        open(unit=11,file='./temperature.txt',status='replace')
        write(unit=11,fmt="('  height','  Ts  ',' Ttube','  Tp ','  Tshell')")
        open(unit=12,file='./timelist.txt',status='replace')
        write(unit=12,fmt="('  time','   Qp  ','  Qs  ',' Tpout ','Tsout ')")
        !+=======================================================
		IHXt=IHX1
		dt=timer1%ttotal/timer1%Nt
		last=0.0
		current=0.0
		
		call IHXt%thcals()
		do i=1,timer1%Nt,1
		   current=last+dt
		   IHXt%Qp=IHXt%Qp*(1-real(i)/timer1%Nt)
		   IHXt%Qs=IHXt%Qs
		   IHXt%Tpin=IHXt%Tpin
		   IHXt%Tsin=IHXt%Tsin
		   call IHXt%thcalt(last,current)
            write(unit=12,fmt="(F6.1,' ',2F6.2,2F6.1)") current,Qp,Qs,IHXt%Tpout,IHXt%Tsout
		   last=current
        enddo
        
 
        do i=0,IHX1%N+1,1
          if(i==0)then
              write(unit=11,fmt="('  0.0000',F6.1,'      ',F6.1)") IHXt%Tsin,IHXt%Tpout
          elseif(i==IHX1%N+1)then
              write(unit=11,fmt="(F8.4,F6.1,'      ',F6.1)") IHXt%Lsingle,IHXt%Tsout,IHXt%Tpin          
          else
              write(unit=11,fmt="(F8.4,4F6.1)") IHXt%zz(i),IHXt%Ts(i),IHXt%Tt(i),IHXt%Tp(i),IHXt%Tv(i)
          end if
        enddo
           
     close(11)
     close(12)
	end subroutine test_IHX_thermal
end module Imp_test_IHX_thermal