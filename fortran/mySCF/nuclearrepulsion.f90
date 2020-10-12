subroutine nuclearrepulsion(nat, xyz, chrg, Enn)

 implicit none

 real*8  ::xyz, chrg
 integer ::atompairs, nat
 integer ::i, j, k

 real*8, intent(out)::Enn

 real*8, allocatable ::R(:)

 dimension xyz(3,nat)
 dimension chrg(nat)

! calculate number of pairs of atoms

 atompairs=0

 do i=1,nat-1
  atompairs=atompairs+i
 end do

 allocate( R(atompairs) )

! calculate interatomic distances

 k=0

 do i=1,nat-1
  do j=i+1,nat
   k=k+1
   R(k)=sqrt( (xyz(1,i)-xyz(1,j))**2 + (xyz(2,i)-xyz(2,j))**2 + (xyz(3,i)-xyz(3,j))**2 )
  end do
 end do

! calculate nuclear repulsion energy

 write(*,*)
 write(*,'(1X,A,A,2X)',advance='NO') "calculating nuclear repulsion energy", "..."

 k=0
 Enn=0

 do i=1,nat-1
  do j=i+1,nat
   k=k+1
   Enn=Enn+( (chrg(i)*chrg(j)) / R(k) )
  end do
 end do

 deallocate(R)

 write(*,*) "done"

!write(*,*)
!write(*,'(A,2X,F8.4,1X,A)') " nuclear repulsion energy:", Enn, "Hartree"

end subroutine nuclearrepulsion
!*************************************************************************************************!
subroutine silentnuclearrepulsion(nat, xyz, chrg, Enn)

 implicit none

 real*8  ::xyz, chrg
 integer ::atompairs, nat
 integer ::i, j, k

 real*8, intent(out)::Enn

 real*8, allocatable ::R(:)

 dimension xyz(3,nat)
 dimension chrg(nat)

! calculate number of pairs of atoms

 atompairs=0

 do i=1,nat-1
  atompairs=atompairs+i
 end do

 allocate( R(atompairs) )

! calculate interatomic distances

 k=0

 do i=1,nat-1
  do j=i+1,nat
   k=k+1
   R(k)=sqrt( (xyz(1,i)-xyz(1,j))**2 + (xyz(2,i)-xyz(2,j))**2 + (xyz(3,i)-xyz(3,j))**2 )
  end do
 end do

! calculate nuclear repulsion energy

 k=0
 Enn=0

 do i=1,nat-1
  do j=i+1,nat
   k=k+1
   Enn=Enn+( (chrg(i)*chrg(j)) / R(k) )
  end do
 end do

 deallocate(R)

end subroutine silentnuclearrepulsion
