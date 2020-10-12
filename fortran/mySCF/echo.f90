subroutine echo(nat, atom, chrg, xyz)

 implicit none

 character (len=2) ::atom
 integer ::nat
 real*8  ::chrg, xyz

 integer :: i

 dimension atom(nat)
 dimension chrg(nat)
 dimension xyz(3,nat)

 write(*,*)
 write(*,*) "****************************************INPUT****************************************"
 write(*,*)
 write(*,*) "atom      ", "charge  ", "coordinates"
 write(*,*) 

 do i=1,nat ! print nuclear charges and coordinates of atoms
  write(*,'(X,A,I3,X,F8.2,F10.4,F8.4,F8.4)') atom(i), i, chrg(i), xyz(:,i)
 end do
 
 write(*,*)
 write(*,*) "*************************************END*OF*INPUT************************************"


end subroutine echo
