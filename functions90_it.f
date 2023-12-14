!      MODULE FUNCTIONS
!      CONTAINS
      function isamax (n,sx,incx)
c  ******************************************************************
c
c 1. function
c  isamax returns the index of element absolute mumimum value 
c 2.  description
c  n     ; in   : Number of elements to process in the vector to be searched
c                  (n=vector length if incx=1;n=vector length/2 if incx=2;
c                    and so on).
c                 if n<=0,isamax returns 0 .
c  sx    ; in   : Real vector to be searched
c  incx  ; in   : Increment between elements of sx
c  ******************************************************************
c
      real :: sx(*), vmax
      integer :: n, incx, ind, j,  i, isamax
      intrinsic abs
c
      ind=0
!      if(n.le.0) goto 20
      if(n.gt.0) then
        if(incx.ge.0) then
          j=1
        else
          j=1-incx*(n-1)
        endif
        vmax=abs(sx(j))
        ind=1
        do 10 i=2,n
          j=j+incx
!        if(vmax.ge.abs(sx(j))) goto 10
        if(vmax.lt.abs(sx(j))) then
          vmax=abs(sx(j))
          ind=i
        end if
   10   end do
        end if
   20   isamax=ind
      return
      end function isamax
!
!
      real function ssum (dim,ary,dummy)                                    1
C...Translated by  fopp     4.02E36 17:59:18  12/14/93
C...Switches: -eadejlpuvx18 -dchimr07 -e1 -gi
      integer   k,dummy,dim                                                 2
      real, dimension(dim) ::  ary                                          3
c     ssum = sum(ary,dim)                                                   4
      ssum=0.                                                               5
*VDIR NODEP
      do k=1,dim,dummy                                                      6
         ssum=ssum+ary(k)                                                   7
      enddo                                                                 8
      return                                                               10
      end function ssum                                                    10
                                                                           11
      real function sdot (dim,v1,d1,v2,d2)
C...Translated by  fopp     4.02E36 17:59:18  12/14/93
C...Switches: -eadejlpuvx18 -dchimr07 -e1 -gi
      integer j1
      integer   dim,d1,d2 
      real :: v1(dim),v2(dim)
      real*8 d3
      d3 = 0
*VDIR NODEP
      do j1 = 1, dim 
         d3 = d3 + v1(1+(j1-1)*d1)*v2(1+(j1-1)*d2)
      end do
      sdot = d3
      return
      end function sdot
 
      subroutine scopy(nn,x,incx,y,incy)
      real :: x(1),y(1)
      integer nn, incx, incy, i
      do i = 1, nn
      y(1+(i-1)*incy) = x(1+(i-1)*incx)
      end do
      return
      end subroutine scopy

      SUBROUTINE SAXPY(n,sa,sx,incx,sy,incy)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      REAL, DIMENSION (1) ::  sx, sy
      REAL :: sa
      INTEGER :: I,incx,incy,IX,IY,M,MP1,n
C
      IF(n.LE.0)RETURN
      IF (sa .EQ. 0.0) RETURN
      IF(incx.NE.1.OR.incy.NE.1)THEN
C   
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(incx.LT.0)IX = (-N+1)*incx + 1
      IF(incy.LT.0)IY = (-N+1)*incy + 1
      DO 10 I = 1,n
        sy(IY) = sy(IY) + sa*sx(IX)
        IX = IX + incx
        IY = IY + incy
   10 END DO
      RETURN
      END if
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(n,4)
      IF( M .NE. 0 )THEN
      DO 30 I = 1,M
        sy(I) = sy(I) + sa*sx(I)
   30 END DO
      IF( n .LT. 4 ) RETURN
      END if
   40 MP1 = M + 1
      DO 50 I = MP1,n,4
        sy(I) = sy(I) + sa*sx(I)
        sy(I + 1) = sy(I + 1) + sa*sx(I + 1)
        sy(I + 2) = sy(I + 2) + sa*sx(I + 2)
        sy(I + 3) = sy(I + 3) + sa*sx(I + 3)
   50 END DO
      RETURN
      END SUBROUTINE SAXPY

      SUBROUTINE SSCAL(n,sa,sx,incx)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      real :: sa
      real, dimension (1) :: sx
      integer :: I,incx,M,MP1,Nincx,n
C
      IF(N.LE.0)return
!      IF(incx.EQ.1)GO TO 20
      IF(incx.ne.1)then
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      Nincx = n*incx
      DO 10 I = 1,Nincx,incx
        sx(I) = sa*sx(I)
   10 END DO
      RETURN
      end if
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(n,5)
!      IF( M .EQ. 0 ) GO TO 40
      IF( M .ne. 0 ) then
      DO 30 I = 1,M
        sx(I) = sa*sx(I)
   30 end do
      IF( n .LT. 5 ) return
      end if
   40 MP1 = M + 1
      DO 50 I = MP1,n,5
        sx(I) = sa*sx(I)
        sx(I + 1) = sa*sx(I + 1)
        sx(I + 2) = sa*sx(I + 2)
        sx(I + 3) = sa*sx(I + 3)
        sx(I + 4) = sa*sx(I + 4)
   50 END DO
      RETURN
      END SUBROUTINE SSCAL
      SUBROUTINE SSWAP (n,sx,incx,sy,incy)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      real :: sx(1),sy(1),STEMP
      integer :: I,incx,incy,IX,IY,M,MP1,n
C
      IF(n.LE.0)return
!      IF(incx.EQ.1.AND.incy.EQ.1)GO TO 20
      IF(incx.ne.1.or.incy.ne.1)then
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(incx.LT.0)IX = (-N+1)*incx + 1
      IF(incy.LT.0)IY = (-N+1)*incy + 1
      DO 10 I = 1,n
        STEMP = sx(IX)
        sx(IX) = sy(IY)
        sy(IY) = STEMP
        IX = IX + incx
        IY = IY + incy
   10 end do
      RETURN
      end if
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(n,3)
!      IF( M .EQ. 0 ) GO TO 40
      IF( M .ne. 0 ) then
      DO 30 I = 1,M
        STEMP = sx(I)
        sx(I) = sy(I)
        sy(I) = STEMP
   30 end do
      IF( n .LT. 3 ) return
      end if
   40 MP1 = M + 1
      DO 50 I = MP1,n,3
        STEMP = sx(I)
        sx(I) = sy(I)
        sy(I) = STEMP
        STEMP = sx(I + 1)
        sx(I + 1) = sy(I + 1)
        sy(I + 1) = STEMP
        STEMP = sx(I + 2)
        sx(I + 2) = sy(I + 2)
        sy(I + 2) = STEMP
   50 end do
      RETURN
      END SUBROUTINE SSWAP
      real function cvmgp (TR,FR,MR)
      real :: TR, FR, MR
      cvmgp = merge(TR,FR,MR>=0.0)
      return
      end function cvmgp
      real function cvmgz (TR,FR,MR)
      real :: TR, FR, MR
      cvmgz = merge(TR,FR,MR==0.0)
      return
      end function cvmgz
C****************************************************************************A
      subroutine second(tt)
C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
C.. In/Out Status: Not Read, Overwritten ..
      real tt
C
C.. Local Scalars ..
!      integer i
!      real tt1
C
C.. Local Arrays ..
      real ra(2)
C
C.. External Functions ..
!      real etime
!      external etime
C
C ... Executable Statements ...
C
!      do i=1,2
!         ra(i) = 0.0
!      end do
C      call etime(ra)
      call cpu_time(tt)
!      tt1 = etime(ra)
!      tt = ra(1)
      end subroutine second
C****************************************************************************A
!      function random()
!      real :: random
!      external rand
!      call random_number(random)
!      return
!      end function random
C****************************************************************************A
      subroutine fdate(t)
C
C.. Implicits ..
      implicit none
C
C.. External Functions..
      external date
      character t*24
c$$$      call date(t)
      end subroutine fdate
C****************************************************************************A
ccc---------------------------------------------------------------------
      subroutine spline (n, x, y, b, c, d)                              00002230
      integer n                                                         00002240
      real x(n), y(n), b(n), c(n), d(n)                                 00002250
c                                                                       00002260
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed      00002270
c  for a cubic interpolating spline                                     00002280
c                                                                       00002290
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3  00002300
c                                                                       00002310
c    for  x(i) .le. x .le. x(i+1)                                       00002320
c                                                                       00002330
c  input..                                                              00002340
c                                                                       00002350
c    n = the number of data points or knots (n.ge.2)                    00002360
c    x = the abscissas of the knots in strictly increasing order        00002370
c    y = the ordinates of the knots                                     00002380
c                                                                       00002390
c  output..                                                             00002400
c                                                                       00002410
c    b, c, d  = arrays of spline coefficients as defined above.         00002420
c                                                                       00002430
c  using  p  to denote differentiation,                                 00002440
c                                                                       00002450
c    y(i) = s(x(i))                                                     00002460
c    b(i) = sp(x(i))                                                    00002470
c    c(i) = spp(x(i))/2                                                 00002480
c    d(i) = sppp(x(i))/6  (derivative from the right)                   00002490
c                                                                       00002500
c  the accompanying function subprogram  seval  can be used             00002510
c  to evaluate the spline.                                              00002520
c                                                                       00002530
c                                                                       00002540
      integer nm1, ib, i                                                00002550
      real t                                                            00002560
c                                                                       00002570
      nm1 = n-1                                                         00002580
      if ( n .lt. 2 ) return                                            00002590
      if ( n .lt. 3 ) go to 50                                          00002600
c                                                                       00002610
c  set up tridiagonal system                                            00002620
c                                                                       00002630
c  b = diagonal, d = offdiagonal, c = right hand side.                  00002640
c                                                                       00002650
      d(1) = x(2) - x(1)                                                00002660
      c(2) = (y(2) - y(1))/d(1)                                         00002670
      do 10 i = 2, nm1                                                  00002680
         d(i) = x(i+1) - x(i)                                           00002690
         b(i) = 2.*(d(i-1) + d(i))                                      00002700
         c(i+1) = (y(i+1) - y(i))/d(i)                                  00002710
         c(i) = c(i+1) - c(i)                                           00002720
   10 continue                                                          00002730
c                                                                       00002740
c  end conditions.  third derivatives at  x(1)  and  x(n)               00002750
c  obtained from divided differences                                    00002760
c                                                                       00002770
      b(1) = -d(1)                                                      00002780
      b(n) = -d(n-1)                                                    00002790
      c(1) = 0.                                                         00002800
      c(n) = 0.                                                         00002810
      if ( n .eq. 3 ) go to 15                                          00002820
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))                        00002830
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))              00002840
      c(1) = c(1)*d(1)**2/(x(4)-x(1))                                   00002850
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))                              00002860
c                                                                       00002870
c  forward elimination                                                  00002880
c                                                                       00002890
   15 do 20 i = 2, n                                                    00002900
         t = d(i-1)/b(i-1)                                              00002910
         b(i) = b(i) - t*d(i-1)                                         00002920
         c(i) = c(i) - t*c(i-1)                                         00002930
   20 continue                                                          00002940
c                                                                       00002950
c  back substitution                                                    00002960
c                                                                       00002970
      c(n) = c(n)/b(n)                                                  00002980
      do 30 ib = 1, nm1                                                 00002990
         i = n-ib                                                       00003000
         c(i) = (c(i) - d(i)*c(i+1))/b(i)                               00003010
   30 continue                                                          00003020
c                                                                       00003030
c  c(i) is now the sigma(i) of the text                                 00003040
c                                                                       00003050
c  compute polynomial coefficients                                      00003060
c                                                                       00003070
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))         00003080
      do 40 i = 1, nm1                                                  00003090
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))          00003100
         d(i) = (c(i+1) - c(i))/d(i)                                    00003110
         c(i) = 3.*c(i)                                                 00003120
   40 continue                                                          00003130
      c(n) = 3.*c(n)                                                    00003140
      d(n) = d(n-1)                                                     00003150
      return                                                            00003160
c                                                                       00003170
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))                                    00003180
      c(1) = 0.                                                         00003190
      d(1) = 0.                                                         00003200
      b(2) = b(1)                                                       00003210
      c(2) = 0.                                                         00003220
      d(2) = 0.                                                         00003230
      return                                                            00003240
      end                                                               00003250
      real function seval(n, u, x, y, b, c, d)                          00003260
      integer n                                                         00003270
      real  u, x(n), y(n), b(n), c(n), d(n)                             00003280
c                                                                       00003290
c  this subroutine evaluates the cubic spline function                  00003300
c                                                                       00003310
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3 00003320
c                                                                       00003330
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule                00003340
c                                                                       00003350
c  if  u .lt. x(1) then  i = 1  is used.                                00003360
c  if  u .ge. x(n) then  i = n  is used.                                00003370
c                                                                       00003380
c  input..                                                              00003390
c                                                                       00003400
c    n = the number of data points                                      00003410
c    u = the abscissa at which the spline is to be evaluated            00003420
c    x,y = the arrays of data abscissas and ordinates                   00003430
c    b,c,d = arrays of spline coefficients computed by spline           00003440
c                                                                       00003450
c  if  u  is not in the same interval as the previous call, then a      00003460
c  binary search is performed to determine the proper interval.         00003470
c                                                                       00003480
      integer i, j, k                                                   00003490
      real dx                                                           00003500
      data i/1/                                                         00003510
      if ( i .ge. n ) i = 1                                             00003520
      if ( u .lt. x(i) ) go to 10                                       00003530
      if ( u .le. x(i+1) ) go to 30                                     00003540
c                                                                       00003550
c  binary search                                                        00003560
c                                                                       00003570
   10 i = 1                                                             00003580
      j = n+1                                                           00003590
   20 k = (i+j)/2                                                       00003600
      if ( u .lt. x(k) ) j = k                                          00003610
      if ( u .ge. x(k) ) i = k                                          00003620
      if ( j .gt. i+1 ) go to 20                                        00003630
c                                                                       00003640
c  evaluate spline                                                      00003650
c                                                                       00003660
   30 dx = u - x(i)                                                     00003670
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))                    00003680
      return                                                            00003690
      end                                                               00003700
ccc---------------------------------------------------------------------
!
!      END MODULE FUNCTIONS
