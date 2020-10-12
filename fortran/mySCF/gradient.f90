subroutine optimization( nat, nao, noe, xyz, chrg, slaterzeta, slaterhere, gaussianzeta,&
                        gaussiancoeff, input, Enn, Eel, atom, method )

 implicit none

!.:INPUT:.
 real*8,            intent(in)::chrg, gaussianzeta, gaussiancoeff
 integer,           intent(in)::nat, nao, noe, slaterhere
 character (len=2), intent(in)::atom
 character (len=3), intent(in)::method

!.:INOUT:.
 real*8,            intent(inout)::xyz, slaterzeta, Enn, Eel
 character (len=4), intent(inout)::input

!.:INTERNAL:.
 real*8 ::oldE, deltaE, totalE
 integer::i, j, k, counter

 dimension xyz(3,nat)
 dimension chrg(nat)
 dimension slaterzeta(nao)
 dimension slaterhere(nat)
 dimension gaussianzeta(6,nao)
 dimension gaussiancoeff(6,nao)
 dimension atom(nat)

 counter=0

 if( input == 'opt' .OR. input == 'sopt') then
  call optimizationroutine( nat, nao, noe, xyz(:,:), chrg(:), slaterzeta(:), slaterhere(:), gaussianzeta(:,:),&
                    gaussiancoeff(:,:), input, Enn, Eel, atom(:), counter, method )
  write(*,*)
  write(*,*) "************************************************"
  write(*,'(1X,A,1X,I3,1X,A)') "********** CONVERGED AFTER", counter, "cycles **********"
  write(*,*) "************************************************"
  if( input == 'opt' ) then
   write(*,*)
   write(*,*) "Optimized structure:"
   write(*,*) "********************"
   write(*,*)
   write(*,'(1X,I3)') nat
   write(*,'(1X,A,1X,F15.10)') "FINAL SCF ENERGY =", Enn+Eel
   do i=1,nat
    write(*,'(1X,A,2X,F11.7,1X,F11.7,1X,F11.7)') atom(i), xyz(:,i)
   end do
  end if
 k=0
  if( input == 'sopt' ) then
   write(*,*)
   write(*,*) "Optimized slater exponents:"
   write(*,*) "***************************"
   write(*,*)
   write(*,'(1X,A,1X,F15.10)') "FINAL SCF ENERGY =", Enn+Eel
   write(*,*)
   write(*,'(1X,A,5X,A)') "Atom", "new exponent"
   do i=1,nat
    write(*,'(1X,A2,I3)',advance='YES') atom(i), i
    do j=1,slaterhere(i)
     k=k+1
     write(*,'(8X,F14.8,5X,F13.8)') slaterzeta(k)
    end do
   end do
  end if
 end if

 oldE=0
 deltaE=-1

 if( input == 'xopt' ) then
  do while( deltaE <= -5.0d-6 )
   input='opt'
   call optimizationroutine( nat, nao, noe, xyz(:,:), chrg(:), slaterzeta(:), slaterhere(:), gaussianzeta(:,:),&
                     gaussiancoeff(:,:), input, Enn, Eel, atom(:), counter, method )
   input='sopt'
   call optimizationroutine( nat, nao, noe, xyz(:,:), chrg(:), slaterzeta(:), slaterhere(:), gaussianzeta(:,:),&
                     gaussiancoeff(:,:), input, Enn, Eel, atom(:), counter, method )
   totalE=Enn+Eel
   deltaE=totalE-oldE 
   oldE=totalE
  end do
  write(*,*)
  write(*,*) "************************************************"
  write(*,'(1X,A,1X,I3,1X,A)') "********** CONVERGED AFTER", counter, "cycles **********"
  write(*,*) "************************************************"
  write(*,*)
  write(*,*) "Optimized structure:"
  write(*,*) "********************"
  write(*,*)
  write(*,'(1X,I3)') nat
  write(*,'(1X,A,1X,F15.10)') "FINAL SCF ENERGY =", totalE
  do i=1,nat
   write(*,'(1X,A,2X,F11.7,1X,F11.7,1X,F11.7)') atom(i), xyz(:,i)
  end do
 k=0
  write(*,*)
  write(*,*) "Optimized slater exponents:"
  write(*,*) "***************************"
  write(*,*)
  write(*,'(1X,A,5X,A)') "Atom", "new exponent"
  do i=1,nat
   write(*,'(1X,A2,I3)',advance='YES') atom(i), i
   do j=1,slaterhere(i)
    k=k+1
    write(*,'(8X,F14.8,5X,F13.8)') slaterzeta(k)
   end do
  end do
 end if

end subroutine optimization
!*************************************************************************************************!
subroutine optimizationroutine( nat, nao, noe, xyz, chrg, slaterzeta, slaterhere, gaussianzeta,&
                        gaussiancoeff, input, inputEnn, inputEel, atom, counter, method )

 implicit none

!.:INPUT:.
 real*8,            intent(in)::chrg, gaussianzeta, gaussiancoeff
 integer,           intent(in)::nat, nao, noe, slaterhere
 character (len=4), intent(in)::input
 character (len=2), intent(in)::atom
 character (len=3), intent(in)::method

!.:INUOT:.
 real*8,  intent(inout)::xyz, slaterzeta, inputEnn, inputEel
 integer, intent(inout)::counter

!.:INTERNAL:.
 real*8, allocatable::thisxyz(:,:), thiszeta(:), S(:,:), h(:,:), twoelint(:,:,:,:), X(:,:)
 real*8, allocatable::P(:,:), oldxyz(:,:), oldzeta(:), xyzgrad(:,:), zetagrad(:)
 real*8, allocatable::previousxyz(:,:), previouszeta(:)
 real*8 ::oldE, deltaE, Enn, Eel, newE
 integer::i, j, k, converged

 dimension xyz(3,nat)
 dimension chrg(nat)
 dimension slaterzeta(nao)
 dimension slaterhere(nat)
 dimension gaussianzeta(6,nao)
 dimension gaussiancoeff(6,nao)
 dimension atom(nat)

 allocate( S(nao,nao), P(nao,nao), previousxyz(3,nat), previouszeta(nao) )
 allocate( oldxyz(3,nat), oldzeta(nao), thisxyz(3,nat), thiszeta(nao), xyzgrad(3,nat), zetagrad(nao) )
 newE=inputEnn+inputEel
 deltaE=-1
 converged=0
 oldxyz(:,:)=xyz(:,:)
 oldzeta(:)=slaterzeta(:)
 xyzgrad(:,:)=0
 zetagrad(:)=0
 previousxyz(:,:)=xyz(:,:)
 previouszeta(:)=slaterzeta(:)

 do while(converged.eq.0 .OR. deltaE <= -5.0d-6)
  counter=counter+1
  oldE=newE
  write(*,*)
  write(*,*) "*******************************************"
  write(*,'(1X,A,1X,I3,1X,A)') "********** Optimization step", counter, "**********"
  write(*,*) "*******************************************"
  write(*,*)
  call gradient( nat, nao, noe, oldxyz(:,:), chrg(:), oldzeta(:), slaterhere(:), gaussianzeta(:,:),&
                gaussiancoeff(:,:), input, thisxyz(:,:), thiszeta(:), converged, xyzgrad(:,:),&
                zetagrad(:), previousxyz(:,:), previouszeta(:), method )
  oldxyz(:,:)=thisxyz(:,:)
  oldzeta(:)=thiszeta(:)
  if(input == 'opt' ) then
   call nuclearrepulsion( nat, thisxyz(:,:), chrg(:), Enn )
   allocate( h(nao,nao), twoelint(nao,nao,nao,nao), X(nao,nao) )
   call matrices( nat, nao, thisxyz(:,:), slaterhere(:), gaussianzeta(:,:), gaussiancoeff(:,:), chrg(:),&
                 S(:,:), h(:,:), twoelint(:,:,:,:), X(:,:) )
   call scfmodule( nao, noe, Enn, h(:,:), twoelint(:,:,:,:), X(:,:), Eel, P(:,:), method )
   deallocate(h, twoelint, X)
   newE=Enn+Eel
   deltaE=newE-oldE
   inputEnn=Enn
   inputEel=Eel
   write(*,*)
   write(*,*) "New structure:"
   write(*,*) "**************"
   write(*,*)
   write(*,'(1X,I3)') nat
   write(*,'(1X,A,1X,F15.10)') "FINAL SCF ENERGY =", newE
   do i=1,nat
    write(*,'(1X,A,2X,F11.7,1X,F11.7,1X,F11.7)') atom(i), thisxyz(:,i)
   end do
  end if
  if(input == 'sopt') then
   call setupbasis( nat, nao, gaussianzeta(:,:), gaussiancoeff(:,:), slaterhere(:), thiszeta(:) )
   allocate( h(nao,nao), twoelint(nao,nao,nao,nao), X(nao,nao) )
   call matrices( nat, nao, xyz(:,:), slaterhere(:), gaussianzeta(:,:), gaussiancoeff(:,:), chrg(:),&
                 S(:,:), h(:,:), twoelint(:,:,:,:), X(:,:) )
   call scfmodule( nao, noe, inputEnn, h(:,:), twoelint(:,:,:,:), X(:,:), Eel, P(:,:), method )
   deallocate(h, twoelint, X)
   newE=inputEnn+Eel
   deltaE=newE-oldE
   inputEel=Eel
   write(*,*)
   write(*,*) "New slater exponents:"
   write(*,*) "*********************"
   write(*,*)
   write(*,'(1X,A,5X,A,5X,A)') "Atom", "old exponent", "new exponent"
   k=0
   do i=1,nat
    write(*,'(1X,A2,I3)',advance='YES') atom(i), i
    do j=1,slaterhere(i)
     k=k+1
     write(*,'(8X,F14.9,5X,F12.9)') slaterzeta(k), thiszeta(k)
    end do
   end do
  end if
  write(*,*)
  write(*,'(1X,A,1X,F12.8,1X,A)') "change in energy after optimization step:", deltaE, "Hartree"
 end do

 deallocate(xyzgrad, zetagrad, previousxyz, previouszeta)

 call mullikanpop( slaterhere(:), P(:,:), S(:,:), noe, nao, chrg(:), nat, atom(:) )
 deallocate(P, S)

 slaterzeta(:)=thiszeta(:)
 xyz(:,:)=thisxyz(:,:)

 deallocate(thisxyz, thiszeta)

end subroutine optimizationroutine
!*************************************************************************************************!
subroutine gradient( nat, nao, noe, oldxyz, chrg, oldzeta, slaterhere, oldgaussianzeta, oldgaussiancoeff,&
                    input, nextxyz, nextzeta, converged, oldxyzgrad, oldzetagrad,&
                    previousxyz, previouszeta, method )

 implicit none

!.:INPUT:.
 real*8,            intent(in)::oldxyz, chrg, oldzeta, oldgaussianzeta, oldgaussiancoeff
 integer,           intent(in)::nat, nao, noe, slaterhere
 character (len=4), intent(in)::input
 character (len=3), intent(in)::method

!.:OUTPUT:.
 real*8,  intent(out)::nextxyz, nextzeta
 integer, intent(out)::converged

!.:INOUT:.
 real*8, intent(inout)::oldxyzgrad, oldzetagrad, previousxyz, previouszeta

!.:INTERNAL:.
 real*8,  allocatable::newxyz(:,:), gaussianzeta(:,:), gaussiancoeff(:,:), S(:,:), h(:,:), twoelint(:,:,:,:)
 real*8,  allocatable::X(:,:), P(:,:), xyzgrad(:,:), newzeta(:), zetagrad(:)
 real*8 ::delta, Enn, Eel, Energy1, Energy2, two, limit
 integer::i, j

 dimension oldxyz(3,nat)
 dimension chrg(nat)
 dimension oldzeta(nao)
 dimension slaterhere(nat)
 dimension oldgaussianzeta(6,nao)
 dimension oldgaussiancoeff(6,nao)
 dimension nextxyz(3,nat)
 dimension nextzeta(nao)
 dimension oldxyzgrad(3,nat)
 dimension oldzetagrad(nao)
 dimension previousxyz(3,nat)
 dimension previouszeta(nao)

 two=2.0d0
 converged=0
 delta=5.0d-6
 limit=5.0d-6

!********************************************!
!********** calculate xyz-gradient **********!
!********************************************!

 if(input == 'opt') then
  write(*,*) "calculating gradient of atomic coordinates"
  nextxyz(:,:)=oldxyz(:,:)
  allocate( newxyz(3,nat), xyzgrad(3,nat) )
  do i=2,nat
   write(*,'(1X)',advance='NO')
   do j=1,3
    newxyz(:,:)=nextxyz(:,:) 
    newxyz(j,i)=nextxyz(j,i)+delta
    call silentnuclearrepulsion( nat, newxyz(:,:), chrg(:), Enn )
    allocate( S(nao,nao), h(nao,nao), twoelint(nao,nao,nao,nao), X(nao,nao) )
    call silentmatrices( nat, nao, newxyz(:,:), slaterhere(:), oldgaussianzeta(:,:), oldgaussiancoeff(:,:), chrg(:),&
                  S(:,:), h(:,:), twoelint(:,:,:,:), X(:,:) )
    deallocate(S)
    allocate( P(nao,nao) )
    call silentscfmodule( nao, noe, Enn, h(:,:), twoelint(:,:,:,:), X(:,:), Eel, P(:,:), method )
    deallocate(h, twoelint, X, P)
    Energy1=Eel+Enn
    newxyz(j,i)=oldxyz(j,i)-delta
    call silentnuclearrepulsion( nat, newxyz(:,:), chrg(:), Enn )
    allocate( S(nao,nao), h(nao,nao), twoelint(nao,nao,nao,nao), X(nao,nao) )
    call silentmatrices( nat, nao, newxyz(:,:), slaterhere(:), oldgaussianzeta(:,:), oldgaussiancoeff(:,:), chrg(:),&
                  S(:,:), h(:,:), twoelint(:,:,:,:), X(:,:) )
    deallocate(S)
    allocate( P(nao,nao) )
    call silentscfmodule( nao, noe, Enn, h(:,:), twoelint(:,:,:,:), X(:,:), Eel, P(:,:), method )
    deallocate(h, twoelint, X, P)
    Energy2=Eel+Enn
    xyzgrad(j,i)=Energy1-Energy2
    xyzgrad(j,i)=xyzgrad(j,i)/(two*delta)
    if( previousxyz(j,i) .ne. oldxyz(j,i) ) nextxyz(j,i)=oldxyz(j,i)-xyzgrad(j,i)*abs(&
                                            ((oldxyz(j,i)-previousxyz(j,i))/(xyzgrad(j,i)-oldxyzgrad(j,i))))
    if( previousxyz(j,i) .eq. oldxyz(j,i) ) nextxyz(j,i)=oldxyz(j,i)-xyzgrad(j,i)
    oldxyzgrad(j,i)=xyzgrad(j,i)
    if( abs(xyzgrad(j,i)) .le. limit ) converged=converged+1
    write(*,'(A)',advance='NO') "."
   end do
   write(*,'(1X,I3,1X,A,1X,I3,1X,A)') 3*(i-1), "of", 3*nat-3, "gradients done"
  end do
  write(*,'(1X,I3,1X,A,1X,I3,1X,A)') converged, "out of", 3*nat-3, "coordinates converged"
  if( converged .eq. 3*nat-3 ) then
   converged=1
  else
   converged=0
  end if
  nextzeta(:)=oldzeta(:)
  deallocate(newxyz, xyzgrad)
  previousxyz(:,:)=oldxyz(:,:)
 end if

!********************************************************!
!********** calculate slater-exponent-gradient **********!
!********************************************************!

 if(input =='sopt') then
  write(*,*) "calculating gradient of slater exponents"
  Enn=0
  nextzeta(:)=oldzeta(:)
  allocate( newzeta(nao), zetagrad(nao) )
  do i=1,nao
   newzeta(:)=nextzeta(:)
   newzeta(i)=nextzeta(i)+delta
   allocate( gaussianzeta(6,nao), gaussiancoeff(6,nao) )
   call silentsetupbasis( nat, nao, gaussianzeta(:,:), gaussiancoeff(:,:), slaterhere(:), newzeta(:) )
   allocate( S(nao,nao), h(nao,nao), twoelint(nao,nao,nao,nao), X(nao,nao) )
   call silentmatrices( nat, nao, oldxyz(:,:), slaterhere(:), gaussianzeta(:,:), gaussiancoeff(:,:), chrg(:),&
                 S(:,:), h(:,:), twoelint(:,:,:,:), X(:,:) )
   deallocate(gaussianzeta, gaussiancoeff, S)
   allocate( P(nao,nao) )
   call silentscfmodule( nao, noe, Enn, h(:,:), twoelint(:,:,:,:), X(:,:), Eel, P(:,:), method )
   deallocate(h, twoelint, X, P)
   Energy1=Eel
   newzeta(i)=oldzeta(i)-delta
   allocate( gaussianzeta(6,nao), gaussiancoeff(6,nao) )
   call silentsetupbasis( nat, nao, gaussianzeta(:,:), gaussiancoeff(:,:), slaterhere(:), newzeta(:) )
   allocate( S(nao,nao), h(nao,nao), twoelint(nao,nao,nao,nao), X(nao,nao) )
   call silentmatrices( nat, nao, oldxyz(:,:), slaterhere(:), gaussianzeta(:,:), gaussiancoeff(:,:), chrg(:),&
                 S(:,:), h(:,:), twoelint(:,:,:,:), X(:,:) )
   deallocate(gaussianzeta, gaussiancoeff, S)
   allocate( P(nao,nao) )
   call silentscfmodule( nao, noe, Enn, h(:,:), twoelint(:,:,:,:), X(:,:), Eel, P(:,:), method )
   deallocate(h, twoelint, X, P)
   Energy2=Eel
   zetagrad(i)=Energy1-Energy2
   zetagrad(i)=zetagrad(i)/(two*delta)
   if( previouszeta(i) .ne. oldzeta(i) ) nextzeta(i)=oldzeta(i)-zetagrad(i)*abs(&
                                           ((oldzeta(i)-previouszeta(i))/(zetagrad(i)-oldzetagrad(i))))
   if( previouszeta(i) .eq. oldzeta(i) ) nextzeta(i)=oldzeta(i)-zetagrad(i)
   oldzetagrad(i)=zetagrad(i)
   if( abs(zetagrad(i)) .le. limit ) converged=converged+1
   write(*,'(1X,I3,1X,A,1X,I3,1X,A)') i, "of", nao, "gradients done"
  end do
  nextxyz(:,:)=oldxyz(:,:)
  deallocate(newzeta, zetagrad)
  write(*,'(1X,I3,1X,A,1X,I3,1X,A)') converged, "out of", nao, "exponents converged"
  write(*,*)
  if( converged .eq. nao ) then
   converged=1
  else
   converged=0
  end if
  previouszeta(:)=oldzeta(:)
 end if

end subroutine gradient
