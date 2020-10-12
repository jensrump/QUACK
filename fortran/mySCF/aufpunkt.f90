subroutine aufpunkt(nat, nao, slaterhere, xyz, xyzauf)

 implicit none

 real*8,  intent(in)  ::xyz
 integer, intent(in)  ::slaterhere
 real*8,  intent(out) ::xyzauf
 
 integer, intent(in)  ::nat, nao

 integer ::i, j, k

 dimension xyz(3,nat)
 dimension xyzauf(3,nao)
 dimension slaterhere(nat)

 k=0

 do i=1,nat
  do j=1,slaterhere(i)
   k=k+1
   xyzauf(:,k)=xyz(:,i)
  end do
 end do

end subroutine aufpunkt
