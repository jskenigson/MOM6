!> A module with intrinsic functions written for reproducible answers.
!! Originally this module existed because some intrinsic functions were
!! not supported by some compilers. We now also supply intrinsic functions
!! calculated from series or other means. These functions are not as efficient
!! as those in the math libraries but can be used during non-time critical
!! initialization of data.
module MOM_intrinsic_functions

use iso_fortran_env, only : stdout=>output_unit
use iso_fortran_env, only : stderr=>error_unit

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public :: invcosh
public :: sin_m6
public :: cos_m6
public :: sind_m6
public :: cosd_m6
public :: pi
public :: pi_180
public :: intrinsics_unit_tests

!> The ratio of a circle's circumference to its diameter, approximately
!! 22/7, or 355/113, ...
!! Some people can recite hundreds of digits (base 10). Here are just
!! 41 digits: 3.141592653589793238462643383279502884197169399...
!! In IEEE 754 single precision floating point is 3.141593
!! (24 mantissa bits or 7 digits).
!! In IEEE 754 double precision floating point representation pi is
!! 3.14159265358979323846 (52 mantissa bits or 21 digits) which is the
!! value found in the C library math.h. We provide more digits (100)
!! here for no better reason than the compilers handle it.
real, parameter :: pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
real, parameter :: pi_180 = 0.0174532925199432957692369076848861271344287188854172545609719144017100911460344944368224156963450948216

logical :: use_fortran_intrinsics = .false.

contains

!> Evaluate the inverse cosh, either using a math library or an
!! equivalent expression
function invcosh(x)
  real, intent(in) :: x !< The argument of the inverse of cosh.  NaNs will
                        !! occur if x<1, but there is no error checking
  real :: invcosh

#ifdef __INTEL_COMPILER
  invcosh = acosh(x)
#else
  invcosh = log(x+sqrt(x*x-1))
#endif

end function invcosh

!> Returns sin(x) where x is in radians
real elemental function sin_m6(x)
  real, intent(in) :: x !< Argument of sin [radians]
  real :: a ! abs(x) or multiples thereof
  real :: s ! 1 or -1 depending on quadrant

  s = 1.0
  a = abs(x)
  if (a>2.0*pi) then
    ! Reduce range to 0..2pi
    a = mod(a, 2.0*pi)
  endif
  if (a>pi) then
    ! Reduce range to 0..pi
    a = a - pi
    s = -1.0
  endif
  if (a<=0.5*pi) then
    sin_m6 = sign(sin_Taylor(a), s*x)
  elseif (a<=pi) then
    sin_m6 = sign(sin_Taylor(pi-a), s*x)
  endif

  if (use_fortran_intrinsics) sin_m6 = sin(x)

end function sin_m6

!> Returns sin(x) where x is in degrees
real elemental function sind_m6(x)
  real, intent(in) :: x !< Argument of sin [degrees]
  real :: a ! abs(x) or multiples thereof
  real :: s ! 1 or -1 depending on quadrant
  real :: d2r ! degrees to radians

  d2r = pi / 180.
  s = 1.0
  a = abs(x)
  if (a>360.) then
    ! Reduce range to 0..360
    a = mod(a, 360.)
  endif
  if (a>180.) then
    ! Reduce range to 0..180
    a = a - 180.
    s = -1.0
  endif
  if (a<=90.) then
    sind_m6 = sign(sin_Taylor(d2r*a), s*x)
  elseif (a<=180.) then
    sind_m6 = sign(sin_Taylor(d2r*(180.-a)), s*x)
  endif

  if (use_fortran_intrinsics) sind_m6 = sin(d2r*x)

end function sind_m6

!> Returns sin(x) if x is in range -pi/2..pi/2 calculated using Taylor
!! series. This approach adds Taylor series terms from smallest to largest
!! and is thus as accurate as the underlying f.p. representation can be.
real elemental function sin_Taylor(x)
  real, intent(in) :: x !< Argument of sin in range -pi/2..pi/2 [radians]
  ! Local variables
  integer, parameter :: n = 16 ! N-1 number of terms in series
  ! Coefficients in Taylor series
  ! https://en.wikipedia.org/wiki/Sine#Series_definition
  ! 15 terms of the Taylor series are needed for 64 bit precision when x=pi
  ! 12 terms of the Taylor series are needed for 64 bit precision when x=pi/2
  !  9 terms of the Taylor series are needed for 64 bit precision when x=pi/4
  real, parameter :: C(16) = (/0.16666666666666666, 8.3333333333333333E-003, &
                           1.9841269841269839E-004, 2.7557319223985884E-006, &
                           2.5052108385441710E-008, 1.6059043836821608E-010, &
                           7.6471637318198144E-013, 2.8114572543455198E-015, &
                           8.2206352466243264E-018, 1.9572941063391257E-020, &
                           3.8681701706306830E-023, 6.4469502843844724E-026, &
                           9.1836898637955449E-029, 1.1309962886447716E-031, &
                           1.2161250415535179E-034, 1.1516335620771951E-037/)
  real :: x2 ! x**2
  real :: xxx(n) ! computed powers of x**2
  real :: r ! accumulated terms
  integer :: j ! term number

  x2 = x*x
  xxx(1) = -x2
  do j = 2, n
    xxx(j) = -xxx(j-1) * x2
  enddo
  r = 0.0
  do j = n, 1, -1
    r = r + C(j) * xxx(j)
  enddo

  sin_Taylor = ( 1.0 + r ) * x

end function sin_Taylor

! !> Returns sin(x) if x is in range -pi/2..pi/2 calculated using Horner's
! !! method applied to the Taylor series polynomial.
! real elemental function sin_Horner(x)
!   real, intent(in) :: x !< Argument of sin in range -pi/2..pi/2 [radians]
!   ! Local variables
!   integer, parameter :: n = 16 ! N-1 number of terms in series
!   ! Coefficients in Taylor series, divided by coefficient of previous term
!   ! 15 terms of the Taylor series are needed for 64 bit precision when x=pi
!   ! 12 terms of the Taylor series are needed for 64 bit precision when x=pi/2
!   !  9 terms of the Taylor series are needed for 64 bit precision when x=pi/4
!   real, parameter :: C(19) = (/0.16666666666666667,0.05,0.0238095238095238081, &
!                    0.013888888888888889,0.00909090909090909,0.00641025641025641, &
!                    0.004761904761904762,0.003676470588235294,0.0029239766081871343, &
!                    0.002380952380952381,0.001976284584980237,0.0016666666666666667, &
!                    0.0014245014245014246,0.0012315270935960591,0.001075268817204301, &
!                    0.000946969696969697,0.0008403361344537816,0.0007507507507507508, &
!                    0.0006747638326585695/) ! https://en.wikipedia.org/wiki/Sine#Series_definition
!   real :: x2 ! x**2
!   real :: r ! accumulated terms
!   integer :: j ! term number

!   x2 = x*x
!   r = 1.0
!   do j = n, 1, -1
!     r = 1.0 - ( C(j) * x2 ) * r
!   enddo

!   sin_Horner = r * x

! end function sin_Horner

!> Returns cos(x) where x is in radians
real elemental function cos_m6(x)
  real, intent(in) :: x !< Argument of cos [radians]
  real :: a ! abs(x) or multiples thereof

  a = abs(x)
  if (a>=2.0*pi) then
    ! Reduce range to 0..2pi
    a = mod(a, 2.0*pi)
  endif
  if (a>pi) then
    ! Reduce range to 0..pi
    a = abs(pi - a)
  endif
  if (a<=0.5*pi) then
    cos_m6 = cos_Taylor(a)
  elseif (a<=pi) then
    cos_m6 = -cos_Taylor(pi-a)
  endif

  if (use_fortran_intrinsics) cos_m6 = cos(x)

end function cos_m6

!> Returns cos(x) where x is in degrees
real elemental function cosd_m6(x)
  real, intent(in) :: x !< Argument of cos [degrees]
  real :: a ! abs(x) or multiples thereof
  real :: d2r ! degrees to radians

  d2r = pi / 180.
  a = abs(x)
  if (a>=360.) then
    ! Reduce range to 0..2pi
    a = mod(a, 360.)
  endif
  if (a>180.) then
    ! Reduce range to 0..pi
    a = abs(180. - a)
  endif
  if (a<=90.) then
    cosd_m6 = cos_Taylor(d2r*a)
  elseif (a<=180.) then
    cosd_m6 = -cos_Taylor(d2r*(180.-a))
  endif

  if (use_fortran_intrinsics) cosd_m6 = cos(d2r*x)

end function cosd_m6

!> Returns cos(x) if x is in range -pi/2..pi/2 calculated using Taylor
!! series. This approach adds Taylor series terms from smallest to largest
!! and is thus as accurate as the underlying f.p. representation can be.
real elemental function cos_Taylor(x)
  real, intent(in) :: x !< Argument of sin in range -pi/2..pi/2 [radians]
  ! Local variables
  integer, parameter :: n = 16 ! N-1 number of terms in series
  ! Coefficients in Taylor series
  ! https://en.wikipedia.org/wiki/Sine#Series_definition
  real, parameter :: C(20) = (/ 0.50000000000000000, 4.1666666666666664E-002, &
                            1.3888888888888887E-003, 2.4801587301587298E-005, &
                            2.7557319223985888E-007, 2.0876756987868096E-009, &
                            1.1470745597729723E-011, 4.7794773323873846E-014, &
                            1.5619206968586225E-016, 4.1103176233121644E-019, &
                            8.8967913924505722E-022, 1.6117375710961182E-024, &
                            2.4795962632247972E-027, 3.2798892370698378E-030, &
                            3.7699876288159054E-033, 3.8003907548547434E-036, &
                            3.3871575355211618E-039, 2.6882202662866363E-042, &
                            1.9119632050402820E-045, 1.2256174391283858E-048/)
  real :: x2 ! x**2
  real :: xxx(n) ! computed powers of x**2
  real :: r ! accumulated terms
  integer :: j ! term number

  x2 = x*x
  xxx(1) = -x2
  do j = 2, n
    xxx(j) = -xxx(j-1) * x2
  enddo
  r = 0.0
  do j = n, 1, -1
    r = r + C(j) * xxx(j)
  enddo

  cos_Taylor = 1.0 + r

end function cos_Taylor

!> Runs unit tests for MOM_intrinsic_functions.
!! Returns .true. if a test fails, otherwise returns .false.
logical function intrinsics_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables
  logical :: fail ! True if a test failed
  real :: x ! Temporary argument

  if (verbose) write(stdout,*) '==== MOM_intrinsic_functions: intrinsics_unit_tests ==='
  fail = .false.  ! Start with no fails

  if (verbose) write(stdout,'(a21,1pe24.16)') 'module pi:',pi
  if (verbose) write(stdout,'(a21,2a24,x,a)') '','result','correct result','err/epsilon'

  fail = test(fail, pi, 4.0 * atan( 1.0 ), 'module pi (v. library)')

  ! Sine tests
  if (verbose) write(stdout,*) 'Tests of sin()'
  fail = test(fail, sin_m6(0.0), 0., 'sin(0)')
  fail = test(fail, sin_m6(pi/12.), 0.25*(sqrt(6.)-sqrt(2.)), 'sin(pi/12)=0.2588...', inexact=1.)
  fail = test(fail, sin_m6(pi/6.), .5, 'sin(pi/6)=0.5', inexact=1.)
  fail = test(fail, sin_m6(0.25*pi), 0.5*sqrt(2.), 'sin(pi/4)=sqrt(0.5)', inexact=1.)
  fail = test(fail, sin_m6(pi/3.), 0.5*sqrt(3.), 'sin(pi/3)=sqrt(3/4)', inexact=1.)
  fail = test(fail, sin_m6(0.5*pi), 1.0, 'sin(pi/2)=1')
  fail = test(fail, sin_m6(pi), 0., 'sin(pi)')
  fail = test(fail, sin_m6(1.5*pi), -1.0, 'sin(3/2 pi)')
  fail = test(fail, sin_m6(2.5*pi), 1.0, 'sin(5/2 pi)')
  fail = test(fail, sin_m6(-2.5*pi), -1.0, 'sin(-5/2 pi)')

  ! Cosine tests
  if (verbose) write(stdout,*) 'Tests of cos()'
  fail = test(fail, cos_m6(0.), 1., 'cos(0)=1')
  fail = test(fail, cos_m6(0.25*pi), sqrt(0.5), 'cos(pi/4)=sqrt(0.5)',inexact=1.)
  fail = test(fail, cos_m6(0.5*pi), 0., 'cos(pi/2)=0')
  fail = test(fail, cos_m6(pi), -1., 'cos(pi)=-1')
  fail = test(fail, cos_m6(1.5*pi), 0., 'cos(3/2 pi)=0')
  fail = test(fail, cos_m6(2.0*pi), 1., 'cos(2pi)=-1')

  ! Tests that sin(x)**2 + cos(x)**2 = 1 (or less within a bit)
  if (verbose) write(stdout,*) 'Tests of sin(x)**2 + cos(x)**2'
  x = 0.5*pi
  fail = test(fail, cos_m6(x)**2+sin_m6(x)**2, 1., 'cos^2+sin^2, x=pi/2')
  x = pi/3.
  fail = test(fail, cos_m6(x)**2+sin_m6(x)**2, 1., 'cos^2+sin^2, x=pi/3')
  x = 0.25
  fail = test(fail, cos_m6(x)**2+sin_m6(x)**2, 1., 'cos^2+sin^2, x=1/4', inexact=1.)
  x = 0.5
  fail = test(fail, cos_m6(x)**2+sin_m6(x)**2, 1., 'cos^2+sin^2, x=1/2')

  ! Sine tests in degrees
  if (verbose) write(stdout,*) 'Tests of sind() in degrees'
  fail = test(fail, sind_m6(0.0), 0., 'sin(0)')
  fail = test(fail, sind_m6(30.), .5, 'sin(30)=0.5', inexact=1.)
  fail = test(fail, sind_m6(45.), 0.5*sqrt(2.), 'sin(45)=sqrt(0.5)', inexact=1.)
  fail = test(fail, sind_m6(60.), 0.5*sqrt(3.), 'sin(60)=sqrt(3/4)', inexact=1.)
  fail = test(fail, sind_m6(90.), 1.0, 'sin(90)=1')
  fail = test(fail, sind_m6(180.), 0., 'sin(180)')
  fail = test(fail, sind_m6(270.), -1.0, 'sin(270)')

  ! Cosine tests in degrees
  if (verbose) write(stdout,*) 'Tests of cosd() in degrees'
  fail = test(fail, cosd_m6(0.), 1., 'cos(0)=1')
  fail = test(fail, cosd_m6(45.), sqrt(0.5), 'cos(45)=sqrt(0.5)',inexact=1.)
  fail = test(fail, cosd_m6(90.), 0., 'cos(90)=0')
  fail = test(fail, cosd_m6(180.), -1., 'cos(180)=-1')
  fail = test(fail, cosd_m6(270.), 0., 'cos(270)=0')
  fail = test(fail, cosd_m6(360.), 1., 'cos(360)=-1')

  if (verbose .and. .not. fail) write(stdout,*) 'Pass'
  intrinsics_unit_tests = fail

  contains

  !> Returns true if a==b or previous failed
  logical function test(previous_test, val, correctval, msg, inexact)
    logical, intent(in) :: previous_test !< True if an earlier test failed
    real, intent(in) :: val !< Value to test against correct value
    real, intent(in) :: correctval !< Correct value
    character(len=*), intent(in) :: msg !< Label
    real, optional :: inexact !< If present allow a relative error of inexact*epsilon
    ! Local variables
    integer :: chan ! Stream to print to
    real :: err ! error
    err = abs( val - correctval )
    test = err > 0. ! Comparison must be exact
    if (present(inexact)) then
      test = err > inexact*val*epsilon(val) ! Allow a round off difference
    endif
    chan = stdout
    if (test.and..not.verbose) chan = stderr
    if (verbose.or.test) then
      if (test) then
        err = err / ( abs(val) * epsilon(val) )
        write(chan,'(a20,":",2(1pe24.16),g9.2,a)') msg,val,correctval,err,' <--- FAIL!'
      elseif (err<=0.) then
        write(chan,'(a20,":",2(1pe24.16))') msg,val,correctval
      else
        err = err / ( abs(val) * epsilon(val) )
        if (err>1.) then
          write(chan,'(a20,":",2(1pe24.16),g9.2,a)') msg,val,correctval,err,' <- differ in last bits'
        else
          write(chan,'(a20,":",2(1pe24.16),g9.2,a)') msg,val,correctval,err,' <- differs in last bit'
        endif
      endif
    endif
    test = test .or. previous_test
  end function test

end function intrinsics_unit_tests

end module MOM_intrinsic_functions
