***********************************************************************
* overlap  integral [s:s]
* kinetic  energy integral [s:t:s]
* potential energy integral [s:sum z/r:s]
*
* s.grimme, aug 1993 (quantenchemie vol 5, page 592...)
*
* npa     : # of primitives in <a|
* npb     : # "      "      "  |b>
* nat     : # of atoms in entire system
* xyz(3,*): atomic coordinates 
* chrg(*) : atomic charges
* xyza(3) : aufpunkt of <a|     
* xyzb(3) : aufpunkt of |b>     
* za()  : exponents of primitives in <a| 
* zb()  : exponents of primitives in |b> 
* ca()  : contraction coeffs. of primitives in <a| 
* cb()  : contraction coeffs. of primitives in |b> 
*
* s00 : overlap
* t00 : kinetic energy
* v00 : potential energy
***********************************************************************

      subroutine stvint(npa,npb,nat,xyz,chrg,xyza,xyzb,za,zb,ca,cb,
     .                  s00,t00,v00)
      implicit real*8 (a-h,o-z)   

      dimension xyz(3,nat)
      dimension xyza(3),xyzb(3),za(npa),zb(npb),ca(npa),cb(npb)
      dimension chrg(nat)

      data pi /3.14159265358979d0/
      data tpi/6.283185307179580/
      data one/1.0d0/,two/2.0d0/,three/3.0d0/,p15/1.50d0/

      absq=(xyza(1)-xyzb(1))**2
     .    +(xyza(2)-xyzb(2))**2
     .    +(xyza(3)-xyzb(3))**2 

      s00=0.0d0
      t00=0.0d0
      v00=0.0d0

      do i=1,npa     
         alp=za(i)
         pqx=alp*xyza(1)
         pqy=alp*xyza(2)
         pqz=alp*xyza(3)
         do 10 j=1,npb     
            bet=zb(j)
            ab=alp*bet
            eab=alp+bet
            eabo=one/eab
            fn=ca(i)*cb(j)
            eabm=ab*eabo
            sqm=absq*eabm
            abab=dexp(-sqm)	
            dab=fn*(pi*eabo)**p15*abab            

ccccc
c s c
ccccc
            s00=s00+dab

ccccc
c t c
ccccc
            t00=t00+eabm*(three-two*sqm)*dab             
ccccc
c v c
ccccc
            f=fn*tpi*eabo*abab            
            pq1=(pqx+bet*xyzb(1))*eabo
            pq2=(pqy+bet*xyzb(2))*eabo
            pq3=(pqz+bet*xyzb(3))*eabo

            do 20 ic1=1,nat
               dx=pq1-xyz(1,ic1)
               dy=pq2-xyz(2,ic1)
               dz=pq3-xyz(3,ic1)
               cpsq=dx*dx+dy*dy+dz*dz 
               v00=v00-chrg(ic1)*f02(eab*cpsq)*f
  20        continue

  10     continue
      enddo

      return
      end


***********************************************************************
* two-electron repulsion integral [ab|cd] over s-functions
* quantity is given in chemist's notation
*
* s.grimme, aug 1993 (quantenchemie vol 5, page 592...)
* npa     : # of primitives in a (same for b,c,d)
* xyza(3) : aufpunkt of a (same for b,c,d)     
* za()  : exponents of primitives in a (same for b,c,d)
* ca()  : contraction coeffs. of primitives in a (same for b,c,d)
*
* g : two-electron integral
***********************************************************************

      subroutine twoint(npa,npb,npc,npd,xyza,xyzb,xyzc,xyzd,
     .                  za,zb,zc,zd,ca,cb,cc,cd,g)

      implicit real*8 (a-h,o-z) 
      dimension xyza(3),xyzb(3),xyzc(3),xyzd(3)
      dimension za(npa),zb(npb),zc(npc),zd(npd)
      dimension ca(npa),cb(npb),cc(npc),cd(npd)

      data twopi25/34.98683665524963d0/
      data pi2/6.36619772367582d-1/

      dimension p(3), q(3)

c R(a-b)
      rab=(xyza(1)-xyzb(1))**2
     .   +(xyza(2)-xyzb(2))**2
     .   +(xyza(3)-xyzb(3))**2
c R(c-d)
      rcd=(xyzc(1)-xyzd(1))**2
     .   +(xyzc(2)-xyzd(2))**2
     .   +(xyzc(3)-xyzd(3))**2 

      g=0.0d0

      do 10 i=1,npa
         alp=za(i)
         c1 =ca(i)
         do 10 j=1,npb 
            bet=zb(j)
            c12=cb(j)*c1
            eab=alp+bet
            eabo=1.0d0/eab
c new gaussian at p
            do m=1,3
               p(m)=(alp*xyza(m)+bet*xyzb(m))*eabo
            enddo
            ab2=dexp(-alp*bet*rab*eabo)

            do 20 k=1,npc
               gam=zc(k)
               c3 =cc(k)
               do 20 l=1,npd 
                  del=zd(l)
                  c34=cd(l)*c3                  
                  ecd=gam+del
                  ecdo=1.0d0/ecd
c new gaussian at q
                  do m=1,3
                     q(m)=(gam*xyzc(m)+del*xyzd(m))*ecdo
                  enddo
                  cd2=dexp(-gam*del*rcd*ecdo)
                  abab=ab2*cd2

                  rpq=(p(1)-q(1))**2+
     .                (p(2)-q(2))**2+
     .                (p(3)-q(3))**2

                  eabcd=eab*ecd
                  abcd=eab+ecd
 
                  t=rpq*eabcd/abcd
                  g=g+twopi25/(eabcd*dsqrt(abcd))*f02(t)*
     .                abab*c12*c34


 20         continue
 10         continue

      return
      end



************************************************************************
* F0 FUNCTION 
************************************************************************

      REAL*8 FUNCTION F02(ARG)
      IMPLICIT REAL*8 (A-H,O-Z)        
      DATA PI/3.14159265358979D0/

      IF(ARG.LT.1.0D-10) GOTO 1

      F02=0.50D0*DSQRT(PI/ARG)*DERF(DSQRT(ARG))
      RETURN

    1 F02=1.0D0-0.33333333333333D0*ARG
      RETURN

      END

