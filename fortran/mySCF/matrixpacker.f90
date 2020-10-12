subroutine matrixpacker(n, matrix, vector)

 implicit none
 
 real*8  ::matrix, vector
 integer ::i, j, k, n

 dimension matrix(n,n)
 dimension vector(n*(n+1)/2)

 k=0

 do i=1,n
  do j=1,i
   k=k+1
   vector(k)=matrix(j,i)
  end do
 end do

end subroutine matrixpacker
