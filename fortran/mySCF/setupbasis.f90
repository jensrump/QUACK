subroutine setupbasis(nat, nao, gaussianzeta, gaussiancoeff, slaterhere, slaterzeta)

 implicit none

 real*8  ::gaussianzeta, gaussiancoeff, slaterzeta
 integer ::nat, nao, slaterhere

 real*8, allocatable ::z(:), C(:)
 integer ::i, j, k, l

 dimension gaussianzeta(6,nao)
 dimension gaussiancoeff(6,nao)
 dimension slaterhere(nat)
 dimension slaterzeta(nao)

 write(*,'(1X,A,13X,A,2X)',advance='NO') "setting up STO-6G basis", "..."

 k=0
 allocate( z(6), C(6) )

 do i=1,nat
  do j=1,slaterhere(i)
   k=k+1
   call slater(slaterzeta(k), z(:), C(:))
   do l=1,6
    gaussiancoeff(l,k)=C(l)
    gaussianzeta(l,k)=z(l)
   end do
  end do
 end do

 deallocate(z, C)

 write(*,*) "done"


end subroutine setupbasis
!*************************************************************************************************!
subroutine silentsetupbasis(nat, nao, gaussianzeta, gaussiancoeff, slaterhere, slaterzeta)

 implicit none

 real*8  ::gaussianzeta, gaussiancoeff, slaterzeta
 integer ::nat, nao, slaterhere

 real*8, allocatable ::z(:), C(:)
 integer ::i, j, k, l

 dimension gaussianzeta(6,nao)
 dimension gaussiancoeff(6,nao)
 dimension slaterhere(nat)
 dimension slaterzeta(nao)

 k=0
 allocate( z(6), C(6) )

 do i=1,nat
  do j=1,slaterhere(i)
   k=k+1
   call slater(slaterzeta(k), z(:), C(:))
   do l=1,6
    gaussiancoeff(l,k)=C(l)
    gaussianzeta(l,k)=z(l)
   end do
  end do
 end do

 deallocate(z, C)

end subroutine silentsetupbasis
