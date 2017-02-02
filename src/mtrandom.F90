! A Fortran-program for MT19937: Real number version

! Code converted using TO_F90 by Alan Miller
! Date: 1999-11-26  Time: 17:09:23
! This version is compatible with Lahey's ELF90 compiler,
! and should be compatible with any full Fortran 90 or 95 compiler.
! Notice the strange way in which umask is specified for ELF90.

!   genrand() generates one pseudorandom real number (double) which is
! uniformly distributed on [0,1]-interval, for each call.
! sgenrand(seed) set initial values to the working area of 624 words.
! Before genrand(), sgenrand(seed) must be called once.  (seed is any 32-bit
! integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.

! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Library General Public License as published by
! the Free Software Foundation; either version 2 of the License, or (at your
! option) any later version.   This library is distributed in the hope that
! it will be useful, but WITHOUT ANY WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General Public License
! along with this library; if not, write to the Free Foundation, Inc.,
! 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.

! N.B. Subroutine sgrnd has been modified in this version to prevent
!      integer overflows which occur with some compilers, particularly
!      when compiler check options are switched on.
!      It is slower that MT19937 because of this, but should yield
!      exactly the same sequence of random numbers.
! Modification by H.D. (Skip) Knoble

!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.

!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed

! This program uses the following standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.

!***********************************************************************

module mtrandom
IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

! Period parameters
INTEGER, PARAMETER :: n = 624, n1 = n+1, m = 397, mata = -1727483681, nrnd=n
!                                    constant vector a
INTEGER, PARAMETER :: umask = -2147483647 - 1
!                                    most significant w-r bits
INTEGER, PARAMETER :: lmask =  2147483647
!                                    least significant r bits
! Tempering parameters
INTEGER, PARAMETER :: tmaskb= -1658038656, tmaskc= -272236544

!                     the array for the state vector
INTEGER, SAVE, target  :: mt(0:n-1), mti = n1
!                     mti==N+1 means mt[N] is not initialized

PRIVATE
PUBLIC :: dp, sgrnd, grnd, nrnd

CONTAINS


SUBROUTINE sgrnd(seed,mt_tgt,mti_tgt)

INTEGER, INTENT(IN)   :: seed
INTEGER, DIMENSION(0:n-1), target,optional:: mt_tgt
INTEGER, target, optional:: mti_tgt
INTEGER, DIMENSION(:), pointer:: mt_ptr
INTEGER, pointer:: mti_ptr

REAL (dp)  :: two31, temp
INTEGER    :: itemp, itemp2
!    setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
!    [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]

if (present(mt_tgt)) then 
   mt_ptr=>mt_tgt
else 
   mt_ptr=>mt
endif
if (present(mti_tgt)) then 
   mti_ptr=>mti_tgt
else 
   mti_ptr=>mti
endif
   

mt_ptr(0)= IAND(seed, -1)
two31 = 2.D0**31
DO  mti_ptr=1,n-1
! The following code in this loop is equivalent to:
! the single line of code:
! mt(mti) = IAND(69069 * mt(mti-1), -1)
!
! The code here is used instead to prevent integer overflow.
   temp=69069.D0*DBLE(mt_ptr(mti_ptr-1))
   itemp=mod(temp,two31)
   itemp2=temp/two31
   if (mod(itemp2,2).ne.0) then
     if (itemp.gt.0) then
       itemp=itemp-two31
     else
       itemp=itemp+two31
     endif
   endif
   mt_ptr(mti_ptr)=itemp
END DO


RETURN
END SUBROUTINE sgrnd
!***********************************************************************

FUNCTION grnd(mt_tgt,mti_tgt) RESULT(fn_val)

REAL (dp) :: fn_val

INTEGER, DIMENSION(0:n-1), target,optional:: mt_tgt
INTEGER, target, optional:: mti_tgt
INTEGER, DIMENSION(:), pointer:: mt_ptr
INTEGER, pointer:: mti_ptr

INTEGER, SAVE :: mag01(0:1) = (/ 0, mata /)
!                        mag01(x) = x * MATA for x=0,1
INTEGER       :: kk, y

! These statement functions have been replaced with separate functions
! tshftu(y) = ISHFT(y,-11)
! tshfts(y) = ISHFT(y,7)
! tshftt(y) = ISHFT(y,15)
! tshftl(y) = ISHFT(y,-18)
if (present(mt_tgt)) then 
   mt_ptr=>mt_tgt
else 
   mt_ptr=>mt
endif
if (present(mti_tgt)) then 
   mti_ptr=>mti_tgt
else 
   mti_ptr=>mti
endif

IF(mti_ptr >= n) THEN
!                       generate N words at one time
  IF(mti_ptr == n+1) THEN
!                            if sgrnd() has not been called,
    CALL sgrnd(4357,mt_ptr,mti_ptr)
!                              a default initial seed is used
  END IF

  DO  kk = 0, n-m-1
    y = IOR(IAND(mt_ptr(kk),umask), IAND(mt_ptr(kk+1),lmask))
    mt_ptr(kk) = IEOR(IEOR(mt_ptr(kk+m), ISHFT(y,-1)),mag01(IAND(y,1)))
  END DO
  DO  kk = n-m, n-2
    y = IOR(IAND(mt_ptr(kk),umask), IAND(mt_ptr(kk+1),lmask))
    mt_ptr(kk) = IEOR(IEOR(mt_ptr(kk+(m-n)), ISHFT(y,-1)),mag01(IAND(y,1)))
  END DO
  y = IOR(IAND(mt_ptr(n-1),umask), IAND(mt_ptr(0),lmask))
  mt_ptr(n-1) = IEOR(IEOR(mt_ptr(m-1), ISHFT(y,-1)),mag01(IAND(y,1)))
  mti_ptr = 0
END IF

y = mt_ptr(mti_ptr)
mti_ptr = mti_ptr + 1
y = IEOR(y, tshftu(y))
y = IEOR(y, IAND(tshfts(y),tmaskb))
y = IEOR(y, IAND(tshftt(y),tmaskc))
y = IEOR(y, tshftl(y))

IF(y < 0) THEN
  fn_val = (DBLE(y) + 2.0D0**32) / (2.0D0**32 - 1.0D0)
ELSE
  fn_val = DBLE(y) / (2.0D0**32 - 1.0D0)
END IF

RETURN
END FUNCTION grnd


FUNCTION tshftu(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,-11)
RETURN
END FUNCTION tshftu


FUNCTION tshfts(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,7)
RETURN
END FUNCTION tshfts


FUNCTION tshftt(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,15)
RETURN
END FUNCTION tshftt


FUNCTION tshftl(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,-18)
RETURN
END FUNCTION tshftl


end module mtrandom

