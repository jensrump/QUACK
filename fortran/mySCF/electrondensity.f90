subroutine electrondensity( nat, nao, slaterhere, xyz, gaussianzeta, gaussiancoeff, P, step,&
                           Rstart, Rfinish)

 implicit none

!.:INPUT:.
 real*8,  intent(in)::xyz, gaussianzeta, gaussiancoeff, P, Rstart, Rfinish
 integer, intent(in)::nat, nao, slaterhere, step

!.:INTERNAL:.
 real*8, allocatable::xyzauf(:,:), xyza(:), xyzb(:), xyzc(:), R(:), rP(:)
 real*8, allocatable::za(:), zb(:), Ca(:), Cb(:)
 real*8 ::d, rho
 integer::i, j, k, l, m

 dimension xyz(3,nat)
 dimension gaussianzeta(6,nao)
 dimension gaussiancoeff(6,nao)
 dimension slaterhere(nao)
 dimension P(nao,nao)
 dimension Rstart(3)
 dimension Rfinish(3)

 allocate( xyzauf(3,nao), xyza(3), xyzb(3), xyzc(3), R(3), rP(3), za(6), zb(6), Ca(6), Cb(6) )

 call aufpunkt( nat, nao, slaterhere(:), xyz(:,:), xyzauf(:,:) )

 write(*,*)
 write(*,*) "Electron density along specified path:"
 write(*,*) "**************************************"
 write(*,*)
 write(*,'(1X,A,2X,A)') "distance from origin", "electron density"

 do m=0,step
  rho=0
  R(:)=Rstart(:)+dble(m)*(Rfinish(:)-Rstart(:))/step
  d=sqrt(dot_product(Rstart-R,Rstart-R))
  do i=1,nao
   xyza(:)=xyzauf(:,i)
   za(:)=gaussianzeta(:,i)
   Ca(:)=gaussiancoeff(:,i)
   do j=1,nao
    xyzb(:)=xyzauf(:,j)
    xyzc(:)=xyza(:)-xyzb(:)
    zb(:)=gaussianzeta(:,j)
    Cb(:)=gaussiancoeff(:,j)
    do k=1,6
     do l=1,6
      rP(:)=( (za(k)*xyza(:) + zb(l)*xyzb(:)) / (za(k)+zb(l)) )
      rP(:)=R(:)-rP(:)
      rho=rho + P(i,j) * Ca(k)*Cb(l) * EXP( -( za(k)+zb(l) ) * dot_product(rP,rP))&
      * EXP( -( za(k)*zb(l)/(za(k)+zb(l)))*(dot_product(xyzc,xyzc)) )
     end do
    end do
   end do
  end do
  write(*,'(1X,F20.5,2X,F16.5)') d, rho
 end do

 deallocate(xyzauf, xyza, xyzb, xyzc, R, rP, za, zb, Ca, Cb)

end subroutine electrondensity
