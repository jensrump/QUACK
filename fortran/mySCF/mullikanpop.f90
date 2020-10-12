subroutine mullikanpop(slaterhere, P, S, noe, nao, chrg, nat, atom)

 implicit none

!.:INPUT:.
 real*8,  intent(in)::P, S, chrg
 integer, intent(in)::noe, nao, nat, slaterhere
 character (len=2), intent(in)::atom

!.:INTERNAL:.
 real*8, allocatable::newnoe(:), PS(:,:), partialchrg(:)
 real*8 ::totalnoe
 integer::i, j, k, l

 dimension slaterhere(nat)
 dimension P(nao,nao)
 dimension S(nao,nao)
 dimension chrg(nat)
 dimension atom(nat)

 allocate( newnoe(nat), PS(nao,nao), partialchrg(nat) )
 newnoe=0
 totalnoe=0
 k=1

 PS=matmul(P,S)

 do i=1,nat
  if(i.gt.1) k=slaterhere(i-1)+1
  if(i.gt.1) l=slaterhere(i)+k-1
  if(i.eq.1) l=slaterhere(i)
  do j=k,l
   newnoe(i)=newnoe(i)+PS(j,j)
  end do
 end do

 partialchrg(:)=chrg(:)-newnoe(:)

 totalnoe=0

 do i=1,nao
  do j=1,nao
   totalnoe=totalnoe+P(i,j)*S(j,i)
  end do
 end do

 write(*,*)
 write(*,*) "Mullikan atomic population analysis:"
 write(*,*) "************************************"
 write(*,*)
 write(*,'(1X,A,2X,I5)') "number of electrons entered:", noe
 write(*,'(1X,A,5X,F10.5)') "number of electrons found:", totalnoe

 write(*,*)
 write(*,'(1X,A,5X,A)') "Atom", "partial charge"
 do i=1,nat
  write(*,'(1X,A2,I3,4X,F8.5)') atom(i), i, partialchrg(i)
 end do


 deallocate(newnoe, PS, partialchrg)

end subroutine mullikanpop
