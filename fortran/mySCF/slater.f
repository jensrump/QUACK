
***********************************************************************
* set up a sto-6g basis for a 1s function with slater exponent zeff
* and normalize the primitives
* return two vectors of order 6
***********************************************************************

      subroutine slater(zeff,z,c)
      implicit real*8 (a-h,o-z)
      dimension allz(6),allc(6)
      dimension z(6),c(6)
      DATA PI2/6.36619772367582D-1/
      
      ALLZ(1) =2.310303149D01
      ALLZ(2) =4.235915534D00
      ALLZ(3) =1.185056519D00
      ALLZ(4) =4.070988982D-01
      ALLZ(5) =1.580884151D-01
      ALLZ(6) =6.510953954D-02

      ALLC(1) =9.163596280D-03
      ALLC(2) =4.936149294D-02
      ALLC(3) =1.685383049D-01
      ALLC(4) =3.705627997D-01
      ALLC(5) =4.164915298D-01
      ALLC(6) =1.303340841D-01

      do i=1,6
         z(i)=allz(i)*zeff**2
         C(I)=allC(I)*(PI2*Z(I))**0.750D0
      enddo

      return
      end

