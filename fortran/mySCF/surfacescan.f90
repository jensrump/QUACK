subroutine surfacescan( nat, xyz, chrg, nao, slaterhere, gaussianzeta, gaussiancoeff, noe, atom,&
                       atomnumber, Rfinish, step, oldEnn, oldEel, method )

 implicit none

!.:INPUT:.
 real*8,            intent(in)::xyz, chrg, gaussianzeta, gaussiancoeff, Rfinish, oldEnn, oldEel
 integer,           intent(in)::nat, nao, slaterhere, noe, atomnumber, step
 character (len=2), intent(in)::atom
 character (len=3), intent(in)::method

!.:INTERNAL:.
 real*8, allocatable::S(:,:), h(:,:), twoelint(:,:,:,:), X(:,:), P(:,:)
 real*8, allocatable::xyznew(:,:), d(:), Energies(:)
 real*8 ::Enn, Eel
 integer::i

 dimension xyz(3,nat)
 dimension chrg(nat)
 dimension slaterhere(nat)
 dimension gaussianzeta(6,nao)
 dimension gaussiancoeff(6,nao)
 dimension atom(nat)
 dimension Rfinish(3)

 write(*,*)
 write(*,*) "Beginning potential surface scan:"
 write(*,*) "*********************************"

 allocate( xyznew(3,nat), S(nao,nao), h(nao,nao), twoelint(nao,nao,nao,nao), X(nao,nao), P(nao,nao),&
          d(step), Energies(step) )

 xyznew(:,:)=xyz(:,:)
 
 do i=1, step
  write(*,*)
  write(*,'(1X,A,1X,I3,A)') "**********STEP", i, "**********"
  xyznew(:,atomnumber)=xyz(:,atomnumber)+dble(i)*(Rfinish(:)-xyz(:,atomnumber))/step
  d(i)=sqrt(dot_product(xyz(:,atomnumber)-xyznew(:,atomnumber),xyz(:,atomnumber)-xyznew(:,atomnumber)))
  call nuclearrepulsion( nat, xyznew(:,:), chrg(:), Enn )
  call matrices( nat, nao, xyznew(:,:), slaterhere(:), gaussianzeta(:,:), gaussiancoeff(:,:), chrg(:),&
               S(:,:), h(:,:), twoelint(:,:,:,:), X(:,:) )
  call scfmodule( nao, noe, Enn, h(:,:), twoelint(:,:,:,:), X(:,:), Eel, P(:,:), method )
  call mullikanpop( slaterhere(:), P(:,:), S(:,:), noe, nao, chrg(:), nat, atom(:) )
  Energies(i)=Enn+Eel
 end do

 write(*,*)
 write(*,*) "Scan summary:"
 write(*,*) "*************"
 write(*,*)
 write(*,'(1X,A,2X,A,2X,A)') "step", "distance from origin", "single point energy"
 write(*,'(1X,I4,2X,F17.8,2X,F13.8)') 0, 0.0d0, oldEnn+oldEel
 do i=1,step
  write(*,'(1X,I4,2X,F17.8,2X,F13.8)') i, d(i), Energies(i)
 end do

 deallocate(xyznew, S, h, twoelint, X, P, d, Energies)

end subroutine surfacescan
