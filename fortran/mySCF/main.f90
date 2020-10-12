program scf
! FORTRAN 90 !
implicit none                 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! system data
integer nat, noe, nao
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! nuclear data
! atom coordinates
real*8, allocatable ::xyz(:,:)
! atoms
character (len=2), allocatable ::atom(:)
! charges
real*8, allocatable ::chrg(:)            
! zetas of AOs 
real*8, allocatable ::slaterzeta(:)
! number of zetas on atom
integer, allocatable ::slaterhere(:)
! exponents of gaussians stored in matrix(number of exponent of slater function, exponents of the 6 primitives), same for coefficients)
real*8, allocatable ::gaussianzeta(:,:), gaussiancoeff(:,:)
! overlap, kinetic energy, potential energy matrices, T+V, two electron integrals
real*8, allocatable :: S(:,:), h(:,:), twoelint(:,:,:,:)
! transformation matrix
real*8, allocatable :: X(:,:)
! matrices
real*8, allocatable :: P(:,:)
! electron density
real*8, allocatable :: Rstart(:), Rfinish(:)
integer:: step
!potential surface scan
integer::atomnumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! other
character (len=50):: filename
character (len=3) :: method
character (len=4) :: input
integer i, j, k

!***RESULTS***!
real*8 ::Enn, Eel

!********************************!
!********** read input **********!
!********************************!

 call getarg(1,filename)

 open(unit=11,file=filename)
 read(11,*) nat, noe, nao, input, method
 if( input == 'rho' ) then
  allocate( Rstart(3), Rfinish(3) )
  read(11,*) Rstart(:), Rfinish(:), step
 end if
 if( input == 'scan' ) then
  allocate( Rfinish(3) )
  read(11,*) atomnumber, Rfinish(:), step
 end if
 allocate( xyz(3,nat), chrg(nat), slaterzeta(nao), slaterhere(nat) )
 k=0
 do i=1,nat !read coordinates, charges, slaterfunctions on corrsponding atom and exponents of slaterfunctions
  read(11,*) xyz(:,i), chrg(i), slaterhere(i)
  do j=1,slaterhere(i)
   k=k+1
   read(11,*) slaterzeta(k)
  end do
 end do
 close(11)

!**************************************************!
!********** match nuclear charge to atom **********!
!**************************************************!

 allocate( atom(nat) )
 call atomtype( chrg(:), atom(:), nat )

!***********************************!
!********** reprint input **********!
!***********************************!

 call echo( nat, atom(:), chrg(:), xyz(:,:) )

!********************************************************!
!********** calculate nuclear repulsion energy **********!
!********************************************************!

 call nuclearrepulsion( nat, xyz(:,:), chrg(:), Enn )

!*****************************************!
!********** set up STO-6G Basis **********!
!*****************************************!

 allocate( gaussianzeta(6,nao), gaussiancoeff(6,nao) )

 call setupbasis( nat, nao, gaussianzeta(:,:) , gaussiancoeff(:,:), slaterhere(:), slaterzeta(:) )

!*************************************!
!********** set up matrices **********!
!*************************************!

 allocate( S(nao,nao), h(nao,nao), twoelint(nao,nao,nao,nao), X(nao,nao) )

 call matrices( nat, nao, xyz(:,:), slaterhere(:), gaussianzeta(:,:), gaussiancoeff(:,:), chrg(:),&
               S(:,:), h(:,:), twoelint(:,:,:,:), X(:,:) )

!*************************************!
!********** enter SCF cycle **********!
!*************************************!

 allocate( P(nao,nao) )

 call scfmodule( nao, noe, Enn, h(:,:), twoelint(:,:,:,:), X(:,:), Eel, P(:,:), method )

!*****************************************************************!
!********** perform Mullikan atomic population analysis **********!
!*****************************************************************!

 call mullikanpop( slaterhere(:), P(:,:), S(:,:), noe, nao, chrg(:), nat, atom(:) )

!****************************************************!
!********** perform potential surface scan **********!
!****************************************************!

 if( input == 'scan' ) then
  call surfacescan( nat, xyz(:,:), chrg(:), nao, slaterhere(:), gaussianzeta(:,:), gaussiancoeff(:,:),&
                   noe, atom(:), atomnumber, Rfinish(:), step, Enn, Eel, method )
  deallocate(Rfinish)
 end if

!************************************************!
!********** calculate electron density **********!
!************************************************!

 if( input == 'rho' ) then
   call electrondensity( nat, nao, slaterhere(:), xyz(:,:), gaussianzeta(:,:), gaussiancoeff(:,:),&
                      P(:,:), step, Rstart(:), Rfinish(:), step )
   deallocate(Rstart, Rfinish)
 end if

 deallocate(S, h, twoelint, P, X)

!******************************************!
!********** perform optimization **********!
!******************************************!

 if( input == 'opt' .OR. input == 'sopt' .OR. input == 'xopt' ) then
  call optimization( nat, nao, noe, xyz(:,:), chrg(:), slaterzeta(:), slaterhere(:), gaussianzeta(:,:),&
                    gaussiancoeff(:,:), input, Enn, Eel, atom(:), method )
 end if

 deallocate(xyz, chrg, slaterzeta, slaterhere, gaussianzeta, gaussiancoeff, atom)

end program scf
