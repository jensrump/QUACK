subroutine matrices(nat, nao, xyz, slaterhere, gaussianzeta, gaussiancoeff, chrg, S, h, twoelint, X)

 implicit none

!.:INPUT:.
 real*8,  intent(in) ::xyz, gaussianzeta, gaussiancoeff, chrg
 integer, intent(in) ::nat, nao, slaterhere
!.:OUTPUT:.
 real*8,  intent(out)::S, h, twoelint, X
!.:INTERNAL:.
 integer::npa, npb, npc, npd
 integer::i, j, k, l
 real*8 ::Gijkl, Sij, Tij, Vij
 real*8, allocatable ::vecS(:), eigenvalues(:), eigenmat(:,:), eigenvectors(:,:)
 real*8, allocatable ::xyzauf(:,:), XT(:,:)
 real*8, allocatable ::work(:)
 real*8, allocatable ::xyza(:), xyzb(:), xyzc(:), xyzd(:)
 real*8, allocatable ::za(:), zb(:), zc(:), zd(:)
 real*8, allocatable ::ca(:), cb(:), cc(:), cd(:)

 dimension xyz(3,nat)
 dimension slaterhere(nat)
 dimension gaussianzeta(6,nao)
 dimension gaussiancoeff(6,nao)
 dimension chrg(nat)
 dimension S(nao,nao)
 dimension h(nao,nao)
 dimension twoelint(nao,nao,nao,nao)
 dimension X(nao,nao)

!******************************************************!
!********** calculate one electron integrals **********!
!******************************************************!

 write(*,'(1X,A,2X,A,2X)',advance='NO') "calculating one electron integrals", "..."

 npa=6
 npb=6
 allocate( xyzauf(3,nao) )

 call aufpunkt( nat, nao, slaterhere(:), xyz(:,:), xyzauf(:,:) )

 allocate( xyza(3), xyzb(3), za(6), zb(6), ca(6), cb(6) )
 do i=1,nao
  xyza(:)=xyzauf(:,i)
  za(:)=gaussianzeta(:,i)
  ca(:)=gaussiancoeff(:,i)
  do j=1,nao
   xyzb(:)=xyzauf(:,j)
   zb(:)=gaussianzeta(:,j)
   cb(:)=gaussiancoeff(:,j)
   call stvint( npa, npb, nat, xyz(:,:), chrg(:), xyza(:), xyzb(:), za(:), zb(:), ca(:), cb(:),&
               Sij, Tij, Vij )
   S(i,j)=Sij
   h(i,j)=Tij+Vij
  end do
 end do

 write(*,*) "done"

!******************************************************!
!********** calculate two electron integrals **********!
!******************************************************!

 write(*,'(1X,A,2X,A,2X)',advance='NO') "calculating two electron integrals", "..."

 npc=6
 npd=6

 allocate( xyzc(3), xyzd(3), zc(6), zd(6), cc(6), cd(6) )
 do i=1,nao
  xyza(:)=xyzauf(:,i)
  za(:)=gaussianzeta(:,i)
  ca(:)=gaussiancoeff(:,i)
  do j=i,nao
   xyzb(:)=xyzauf(:,j)
   zb(:)=gaussianzeta(:,j)
   cb(:)=gaussiancoeff(:,j)
   do k=i,nao
    xyzc(:)=xyzauf(:,k)
    zc(:)=gaussianzeta(:,k)
    cc(:)=gaussiancoeff(:,k)
    do l=k,nao
     if(j.eq.i.AND.k.eq.j.AND.l.eq.k+1) cycle
     xyzd(:)=xyzauf(:,l)
     zd(:)=gaussianzeta(:,l)
     cd(:)=gaussiancoeff(:,l)
     call twoint( npa, npb, npc, npd, xyza(:), xyzb(:), xyzc(:), xyzd(:), za(:), zb(:), zc(:), zd(:),&
                 ca(:), cb(:), cc(:), cd(:), Gijkl )
     twoelint(i,j,k,l)=Gijkl
     twoelint(i,j,l,k)=Gijkl
     twoelint(j,i,k,l)=Gijkl
     twoelint(j,i,l,k)=Gijkl
     twoelint(k,l,i,j)=Gijkl
     twoelint(k,l,j,i)=Gijkl
     twoelint(l,k,i,j)=Gijkl
     twoelint(l,k,j,i)=Gijkl
    end do
   end do
  end do
 end do

 write(*,*) "done"

 deallocate(xyzauf, xyza, xyzb, xyzc, xyzd, za, zb, zc, zd, ca, cb, cc, cd)

!*******************************************************!
!********** calculate transformation matrix X **********!
!*******************************************************!

 write(*,'(1X,A,3X,A,2X)',advance='NO') "calculating transformation matrix", "..."

! pack density matrix S and calculate its eigenvalues and eigenvectors

 allocate( vecS(nao*(nao+1)/2) )
 call matrixpacker( nao, S(:,:), vecS(:) )
 allocate( eigenvalues(nao), eigenvectors(nao,nao), work(3*nao) )
 call dspev( 'V', 'U', nao, vecS(:), eigenvalues(:), eigenvectors(:,:), nao, work(:), i )
 if(i.gt.0.OR.i.lt.0) then
  write(*,*) "ERROR WHILE DIAGONALIZING S"
 end if
 deallocate(work)

! calculate square root of eigenvalues

 eigenvalues(:)=1.0d0/sqrt(eigenvalues(:))

! put eigenvalues in matrix SQRT(E)

 allocate( eigenmat(nao,nao) )
 eigenmat(:,:)=0.0d0
 do i=1,nao
  eigenmat(i,i)=eigenvalues(i)
 end do
 deallocate(eigenvalues)

! transpose eigenvector matrix V

 allocate( XT(nao,nao) )
 XT=transpose(eigenvectors)

! calculate S**1/2 = V*SQRT(E)*V**-1

 X=matmul(eigenvectors,eigenmat)
 X=matmul(X,XT)

 write(*,*) "done"

end subroutine matrices
!*************************************************************************************************!
subroutine silentmatrices(nat, nao, xyz, slaterhere, gaussianzeta, gaussiancoeff, chrg, S, h, twoelint, X)

 implicit none

!.:INPUT:.
 real*8,  intent(in) ::xyz, gaussianzeta, gaussiancoeff, chrg
 integer, intent(in) ::nat, nao, slaterhere
!.:OUTPUT:.
 real*8,  intent(out)::S, h, twoelint, X
!.:INTERNAL:.
 integer::npa, npb, npc, npd
 integer::i, j, k, l
 real*8 ::Gijkl, Sij, Tij, Vij
 real*8, allocatable ::vecS(:), eigenvalues(:), eigenmat(:,:), eigenvectors(:,:)
 real*8, allocatable ::xyzauf(:,:), XT(:,:)
 real*8, allocatable ::work(:)
 real*8, allocatable ::xyza(:), xyzb(:), xyzc(:), xyzd(:)
 real*8, allocatable ::za(:), zb(:), zc(:), zd(:)
 real*8, allocatable ::ca(:), cb(:), cc(:), cd(:)

 dimension xyz(3,nat)
 dimension slaterhere(nat)
 dimension gaussianzeta(6,nao)
 dimension gaussiancoeff(6,nao)
 dimension chrg(nat)
 dimension S(nao,nao)
 dimension h(nao,nao)
 dimension twoelint(nao,nao,nao,nao)
 dimension X(nao,nao)

!******************************************************!
!********** calculate one electron integrals **********!
!******************************************************!

 npa=6
 npb=6
 allocate( xyzauf(3,nao) )

 call aufpunkt( nat, nao, slaterhere(:), xyz(:,:), xyzauf(:,:) )

 allocate( xyza(3), xyzb(3), za(6), zb(6), ca(6), cb(6) )
 do i=1,nao
  xyza(:)=xyzauf(:,i)
  za(:)=gaussianzeta(:,i)
  ca(:)=gaussiancoeff(:,i)
  do j=1,nao
   xyzb(:)=xyzauf(:,j)
   zb(:)=gaussianzeta(:,j)
   cb(:)=gaussiancoeff(:,j)
   call stvint( npa, npb, nat, xyz(:,:), chrg(:), xyza(:), xyzb(:), za(:), zb(:), ca(:), cb(:),&
               Sij, Tij, Vij )
   S(i,j)=Sij
   h(i,j)=Tij+Vij
  end do
 end do

!******************************************************!
!********** calculate two electron integrals **********!
!******************************************************!

 npc=6
 npd=6

 allocate( xyzc(3), xyzd(3), zc(6), zd(6), cc(6), cd(6) )
 do i=1,nao
  xyza(:)=xyzauf(:,i)
  za(:)=gaussianzeta(:,i)
  ca(:)=gaussiancoeff(:,i)
  do j=i,nao
   xyzb(:)=xyzauf(:,j)
   zb(:)=gaussianzeta(:,j)
   cb(:)=gaussiancoeff(:,j)
   do k=i,nao
    xyzc(:)=xyzauf(:,k)
    zc(:)=gaussianzeta(:,k)
    cc(:)=gaussiancoeff(:,k)
    do l=k,nao
     if(j.eq.i.AND.k.eq.j.AND.l.eq.k+1) cycle
     xyzd(:)=xyzauf(:,l)
     zd(:)=gaussianzeta(:,l)
     cd(:)=gaussiancoeff(:,l)
     call twoint( npa, npb, npc, npd, xyza(:), xyzb(:), xyzc(:), xyzd(:), za(:), zb(:), zc(:), zd(:),&
                 ca(:), cb(:), cc(:), cd(:), Gijkl )
     twoelint(i,j,k,l)=Gijkl
     twoelint(i,j,l,k)=Gijkl
     twoelint(j,i,k,l)=Gijkl
     twoelint(j,i,l,k)=Gijkl
     twoelint(k,l,i,j)=Gijkl
     twoelint(k,l,j,i)=Gijkl
     twoelint(l,k,i,j)=Gijkl
     twoelint(l,k,j,i)=Gijkl
    end do
   end do
  end do
 end do

 deallocate(xyzauf, xyza, xyzb, xyzc, xyzd, za, zb, zc, zd, ca, cb, cc, cd)

!*******************************************************!
!********** calculate transformation matrix X **********!
!*******************************************************!

! pack density matrix S and calculate its eigenvalues and eigenvectors

 allocate( vecS(nao*(nao+1)/2) )
 call matrixpacker( nao, S(:,:), vecS(:) )
 allocate( eigenvalues(nao), eigenvectors(nao,nao), work(3*nao) )
 call dspev( 'V', 'U', nao, vecS(:), eigenvalues(:), eigenvectors(:,:), nao, work(:), i )
 if(i.gt.0.OR.i.lt.0) then
  write(*,*) "ERROR WHILE DIAGONALIZING S"
 end if
 deallocate(work)

! calculate square root of eigenvalues

 eigenvalues(:)=1.0d0/sqrt(eigenvalues(:))

! put eigenvalues in matrix SQRT(E)

 allocate( eigenmat(nao,nao) )
 eigenmat(:,:)=0.0d0
 do i=1,nao
  eigenmat(i,i)=eigenvalues(i)
 end do
 deallocate(eigenvalues)

! transpose eigenvector matrix V

 allocate( XT(nao,nao) )
 XT=transpose(eigenvectors)

! calculate S**1/2 = V*SQRT(E)*V**-1

 X=matmul(eigenvectors,eigenmat)
 X=matmul(X,XT)

end subroutine silentmatrices
