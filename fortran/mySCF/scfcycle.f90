subroutine scfmodule(nao, noe, Enn, h, twoelint, X, Eel, oldP, method)

 implicit none
!.:INPUT:.
 real*8,            intent(in)::h, twoelint, X, Enn
 integer,           intent(in)::nao, noe
 character (len=3), intent(in)::method
!.:OUTPUT:.
 real*8, intent(out)::oldP, Eel
!.:INTERNAL:.
 real*8, allocatable::newP(:,:), C(:,:), orbitalenergies(:)
 real*8 ::deltaE, oldEel, totalE, mptwoE
 integer::iteration

 dimension h(nao,nao)
 dimension twoelint(nao,nao,nao,nao)
 dimension X(nao,nao)
 dimension oldP(nao,nao)

 allocate( newP(nao,nao), C(nao,nao), orbitalenergies(nao) )
 oldEel=0
 oldP=0
 iteration=0
 deltaE=-1
 Eel=0

 call scfcycle( nao, noe, h(:,:), twoelint(:,:,:,:), newP(:,:), X(:,:), Eel, oldP(:,:), C(:,:),&
               orbitalenergies(:) )

 write(*,*)
 write(*,*) "Beginning SCF-procedure:"
 write(*,*) "************************"
 write(*,*)
 write(*,'(1X,A9,2X,A15,2X,A20,2X,A19)') "iteration", "total energy", "electronic energy", "change in energy"
 do while(deltaE <= -5.0d-6)
  iteration=iteration+1
  if(iteration.ge.128) then
   write(*,*) "FAILURE TO CONVERGE IN 128 CYCLES, EXITING SCF-MODULE"
   exit
  end if
  oldP=newP
  call scfcycle( nao, noe, h(:,:), twoelint(:,:,:,:), newP(:,:), X(:,:), Eel, oldP(:,:), C(:,:),&
                orbitalenergies(:) )
  deltaE=Eel-oldEel
  totalE=Eel+Enn
  write(*,'(1X,I3,8X,F15.8,2X,F20.8,2X,F19.8)') iteration, totalE, Eel, deltaE
  oldEel=Eel
 end do
 
 oldP=newP

 if(method == 'MP2') then
  write(*,*)
  write(*,'(1X,A)',advance='NO') "Calculating MP2 Energy ..."
  call mptwo(nao, noe, C(:,:), twoelint(:,:,:,:), orbitalenergies(:), mptwoE)
  write(*,*) "done"
 else
  mptwoE=0
 end if

 totalE=totalE+mptwoE

 write(*,*)
 write(*,'(1X,A,2X,F20.8,2X,A)') "Final single point energy:", totalE, "Hartree"
 write(*,*)
 write(*,*) "Components:"
 write(*,'(1X,A18,F20.8,2X,A)') "Nuclear Repulsion:", Enn, "Hartree"
 write(*,'(1X,A18,F20.8,2X,A)') "Electronic Energy:", Eel, "Hartree"
 write(*,'(1X,A18,F20.8,2X,A)') "HF Energy:", Enn+Eel, "Hartree"
 if(method == 'MP2') write(*,'(1X,A18,F20.8,2X,A)') "MP2 Correction:", mptwoE, "Hartree"

 Eel=Eel+mptwoE

 deallocate(newP, C, orbitalenergies)

end subroutine scfmodule
!*************************************************************************************************!
subroutine scfcycle(nao, noe, h, twoelint, P, X, Eel, oldP, C, orbitalenergies)

 implicit none
!.:INPUT:.
 real*8,  intent(in)::h, twoelint, X, oldP
 integer, intent(in)::nao, noe
!.:OUTPUT:.
 real*8, intent(out)::P, Eel, C, orbitalenergies
!.:INTERNAL:.
 real*8, allocatable ::vecX(:)
 real*8, allocatable ::G(:,:), F(:,:)
 real*8, allocatable ::work(:), onehalf, two
 integer ::i, j, k, l

 dimension h(nao,nao)
 dimension twoelint(nao,nao,nao,nao)
 dimension P(nao,nao)
 dimension X(nao,nao)
 dimension oldP(nao,nao)
 dimension C(nao,nao)
 dimension orbitalenergies(nao)

 allocate( G(nao,nao), F(nao,nao) )
 G=0
 F=0
 C=0
 orbitalenergies=0
 onehalf=0.5d0
 two=2.0d0

!*************************************!
!********** begin scf cycle **********!
!*************************************!

  do k=1,nao
   do l=1,nao
    G(:,:)=G(:,:)+oldP(k,l)*(twoelint(:,:,l,k)-onehalf*twoelint(:,k,l,:))
   end do
  end do

  F(:,:)=h(:,:)+G(:,:)


!*************************************************!
!********** calculate electronic energy **********!
!*************************************************!

 Eel=0

 do i=1,nao
  do j=1,nao
   Eel=Eel+onehalf*oldP(i,j)*(h(i,j)+F(i,j))
  end do
 end do

!****************************************************!
!********** calculate new density matrix ************!
!****************************************************!

! transform F to F' = X**T*F*X
 
 F(:,:)=matmul(X,F)
 F(:,:)=matmul(F,X)

! diagonalize F' to obtain C' and orbital energies e from F'*C'=C'*e

 allocate( vecX(nao*(nao+1)/2) )
 call matrixpacker( nao, F(:,:), vecX(:) )
 allocate( work(3*nao) )
 call dspev( 'V', 'U', nao, vecX(:), orbitalenergies(:), C(:,:), nao, work(:), i )
 if(i.lt.0.OR.i.gt.0) then
  write(*,*) "ERROR WHILE DIAGONALIZING F"
 end if
 deallocate(work, vecX)

! calculate C=X*C'

 C(:,:)=matmul(X,C)

! calculating a new density matrix P

 P=0

 do i=1,nao
  do j=1,nao
   do k=1,noe/2
    P(i,j)=P(i,j)+two*C(i,k)*C(j,k)
   end do
  end do
 end do

 deallocate(F, G)

end subroutine scfcycle
!*************************************************************************************************!
subroutine mptwo(nao, noe, C, twoelint, orbitalenergies, mptwoE)

 implicit none

!.:INPUT:.
 real*8,  intent(in)::C, twoelint, orbitalenergies
 integer, intent(in)::nao, noe

!.:OUTPUT:.
 real*8, intent(out)::mptwoE

!.:INTERNAL:.
 real*8, allocatable::motwoelint(:,:,:,:), motwoelinta(:,:,:,:), motwoelintb(:,:,:,:), motwoelintc(:,:,:,:)
 real*8             ::two
 integer            ::i, j, k, l, p, q, r, s

 dimension C(nao,nao)
 dimension twoelint(nao,nao,nao,nao)
 dimension orbitalenergies(nao)

 two=2.0d0
 mptwoE=0

!********************************************************!
!********** Transformation from AO to MO basis **********!
!********************************************************!

 allocate( motwoelint(nao,nao,nao,nao) )
 motwoelint=0

 !do p=1,nao
 ! do q=1,nao
 !  do r=1,nao
 !   do s=1,nao
 !    do i=1,nao
 !     do j=1,nao
 !      do k=1,nao
 !       do l=1,nao
 !        motwoelint(p,q,r,s)=motwoelint(p,q,r,s)+C(i,p)*C(j,q)*twoelint(i,j,k,l)*C(k,r)*C(l,s)
 !       end do
 !      end do
 !     end do
 !    end do
 !   end do
 !  end do
 ! end do
 !end do

 allocate( motwoelinta(nao,nao,nao,nao), motwoelintb(nao,nao,nao,nao), motwoelintc(nao,nao,nao,nao) )
 motwoelinta=0
 motwoelintb=0
 motwoelintc=0

 do p=1,nao
  do i=1,nao
   motwoelinta(p,:,:,:)=motwoelinta(p,:,:,:)+C(i,p)*twoelint(i,:,:,:)
  end do
  do q=1,nao
   do j=1,nao
    motwoelintb(p,q,:,:)=motwoelintb(p,q,:,:)+C(j,q)*motwoelinta(p,j,:,:)
   end do
   do r=1,nao
    do k=1,nao
     motwoelintc(p,q,r,:)=motwoelintc(p,q,r,:)+C(k,r)*motwoelintb(p,q,k,:)
    end do
    do s=1,nao 
     do l=1,nao
      motwoelint(p,q,r,s)=motwoelint(p,q,r,s)+C(l,s)*motwoelintc(p,q,r,l)
     end do
    end do
   end do
  end do
 end do

 deallocate(motwoelinta, motwoelintb, motwoelintc)

!**********************************************!
!********** Caculation of MP2 Energy **********!
!**********************************************!

 do p=1,noe/2
  do q=1,noe/2
   do r=(noe/2)+1,nao
    do s=(noe/2)+1,nao
     mptwoE=mptwoE+( motwoelint(p,r,q,s)*(two*motwoelint(p,r,q,s)-motwoelint(p,s,q,r)) ) /&
                   ( orbitalenergies(p)+orbitalenergies(q)-orbitalenergies(r)-orbitalenergies(s) )
    end do
   end do
  end do
 end do

 deallocate(motwoelint)

end subroutine mptwo
!*************************************************************************************************!
subroutine silentscfmodule(nao, noe, Enn, h, twoelint, X, Eel, oldP, method)

 implicit none
!.:INPUT:.
 real*8,    intent(in)::h, twoelint, X, Enn
 integer,   intent(in)::nao, noe
 character, intent(in)::method
!.:OUTPUT:.
 real*8, intent(out)::oldP, Eel
!.:INTERNAL:.
 real*8, allocatable::newP(:,:), C(:,:), orbitalenergies(:)
 real*8 ::deltaE, oldEel, totalE, mptwoE
 integer::iteration

 dimension h(nao,nao)
 dimension twoelint(nao,nao,nao,nao)
 dimension X(nao,nao)
 dimension oldP(nao,nao)

 allocate( newP(nao,nao), C(nao,nao), orbitalenergies(nao) )
 oldEel=0
 oldP=0
 iteration=0
 deltaE=-1
 Eel=0

 call scfcycle( nao, noe, h(:,:), twoelint(:,:,:,:), newP(:,:), X(:,:), Eel, oldP(:,:), C(:,:),&
               orbitalenergies(:) )

 do while(deltaE <= -5.0d-6)
  iteration=iteration+1
  if(iteration.ge.128) then
   write(*,*) "FAILURE TO CONVERGE IN 128 CYCLES, EXITING SCF-MODULE"
   exit
  end if
  oldP=newP
  call scfcycle( nao, noe, h(:,:), twoelint(:,:,:,:), newP(:,:), X(:,:), Eel, oldP(:,:), C(:,:),&
                orbitalenergies(:) )
  deltaE=Eel-oldEel
  totalE=Eel+Enn
  oldEel=Eel
 end do
 
 oldP=newP

 if(method == "MP2") then
  call mptwo(nao, noe, C(:,:), twoelint(:,:,:,:), orbitalenergies(:), mptwoE)
 else
  mptwoE=0
 end if

 totalE=totalE+mptwoE

 Eel=Eel+mptwoE

 deallocate(newP, C, orbitalenergies)

end subroutine silentscfmodule
