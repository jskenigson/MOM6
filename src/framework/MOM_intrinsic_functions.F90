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
public :: tan_m6
public :: atan_m6
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

!> Returns sin(x).
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

  contains

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

  !> Returns sin(x) if x is in range -pi/2..pi/2 calculated using Horner's
  !! method applied to the Taylor series polynomial.
  real elemental function sin_Horner(x)
    real, intent(in) :: x !< Argument of sin in range -pi/2..pi/2 [radians]
    ! Local variables
    integer, parameter :: n = 16 ! N-1 number of terms in series
    ! Coefficients in Taylor series, divided by coefficient of previous term
    ! 15 terms of the Taylor series are needed for 64 bit precision when x=pi
    ! 12 terms of the Taylor series are needed for 64 bit precision when x=pi/2
    !  9 terms of the Taylor series are needed for 64 bit precision when x=pi/4
    real, parameter :: C(19) = (/0.16666666666666667,0.05,0.0238095238095238081, &
                     0.013888888888888889,0.00909090909090909,0.00641025641025641, &
                     0.004761904761904762,0.003676470588235294,0.0029239766081871343, &
                     0.002380952380952381,0.001976284584980237,0.0016666666666666667, &
                     0.0014245014245014246,0.0012315270935960591,0.001075268817204301, &
                     0.000946969696969697,0.0008403361344537816,0.0007507507507507508, &
                     0.0006747638326585695/) ! https://en.wikipedia.org/wiki/Sine#Series_definition
    real :: x2 ! x**2
    real :: r ! accumulated terms
    integer :: j ! term number

    x2 = x*x
    r = 1.0
    do j = n, 1, -1
      r = 1.0 - ( C(j) * x2 ) * r
    enddo

    sin_Horner = r * x

  end function sin_Horner

end function sin_m6

!> Returns cos(x).
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

  contains

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

end function cos_m6

!> Returns tan(x). Input assumed to be in range -pi/2..pi/2.
real elemental function tan_m6(x)
  real, intent(in) :: x !< Argument of tan(x) [radians]
  real :: a ! abs(x) or multiples thereof
  real :: d ! denominator in division expressions

  a = abs(x)
  if (a<=0.125*pi) then
    ! a is in range 0..pi/8 for which tan_series is accurate
    ! (or 0..22.5 in degrees)
    tan_m6 = tan_series( a )
    tan_m6 = sign( tan_m6, x)
  elseif (a<=0.25*pi) then
    ! Reduce range from pi/8..pi/4 to pi/16..pi/8
    ! (or 22.5...45 to 11.25...22.5 in degrees)
    a = 0.5 * a
    tan_m6 = tan_series( a )
    !d = 1.0 / ( 1.0 - tan_m6**2 )
    d = 1.0 / ( ( 1.0 - tan_m6 ) * ( 1.0 + tan_m6 ) )
    tan_m6 = sign( 2.0 * tan_m6 * d, x )
  else ! a>pi/4
    ! Reduce range from
    ! (or 45...90 to 11.25...22.5 in degrees)
    ! or pi/4..pi/2 to pi/8..pi/8
    a = 0.25 * a
    tan_m6 = tan_series( a )
    !d = 1.0 / ( 1.0 - tan_m6**2 )
    d = 1.0 / ( ( 1.0 - tan_m6 ) * ( 1.0 + tan_m6 ) )
    tan_m6 = 2.0 * tan_m6 * d
    !d = 1.0 / ( 1.0 - tan_m6**2 )
    d = 1.0 / ( ( 1.0 - tan_m6 ) * ( 1.0 + tan_m6 ) )
    tan_m6 = sign( 2.0 * tan_m6 * d, x )
  endif

  if (use_fortran_intrinsics) tan_m6 = tan(x)

  contains

  !> Returns tangent of x if x is in range -pi/6..pi/6.
  real elemental function tan_series(x)
    real, intent(in) :: x !< Argument of tan(x) in range -pi/6..pi/6 [radians]
    ! Local variables
    !integer, parameter :: N(17) =(/ 1, 1, 2, 17, 62, 1382, 21844, 929569, 6404582, 443861162, 18888466084, &
    !  113927491862, 58870668456604, 8374643517010684, 689005380505609448, 129848163681107301953, &
    !  1736640792209901647222/) ! Numerators from http://oeis.org/A002430
    !integer, parameter :: D(17) = (/1, 3, 15, 315, 2835, 155925, 6081075, 638512875, 10854718875, 1856156927625, &
    !  194896477400625, 2900518163668125, 3698160658676859375, 1298054391195577640625, 263505041412702261046875, &
    !  122529844256906551386796875, 4043484860477916195764296875/) ! Denominators from http://oeis.org/A036279
    real, parameter :: C(17) = (/ 1.0, 0.33333333333333333, 0.13333333333333333, 5.3968253968253971E-002, &
      2.1869488536155203E-002, 8.8632355299021973E-003, 3.5921280365724811E-003, 1.4558343870513183E-003, &
      5.9002744094558595E-004, 2.3912911424355248E-004, 9.6915379569294509E-005, 3.9278323883316833E-005, &
      1.5918905069328964E-005, 6.4516892156554306E-006, 2.6147711512907546E-006, 1.0597268320104654E-006, &
      4.2949110782738063E-007 /) ! Coefficients in Taylor series (N/D)
    real :: x2 ! x**2
    real :: r ! accumulated terms
    integer :: j ! term number
    real :: term(17) ! Actual coefficient in Taylor series

    x2 = x*x
    r = 1.0
    do j = 1, 17
      term(j) = C(j) * r
      r = r * x2
    enddo
    r = 0.0
    do j = 17, 1, -1
      r = r + term(j)
    enddo

    tan_series = r * x

  end function tan_series

end function tan_m6

!> Returns arctan(x). Results in in range -pi/2..pi/2. [radians]
real elemental function atan_m6(x)
  real, intent(in) :: x !< Argument of atan(x)
  real :: a ! abs(x)

  a = abs(x)
  if (a<=0.5) then
    atan_m6 = sign( atan_series( a ), x)
  elseif (a<2.0) then
    atan_m6 = sign( p4atx( a - 1.0 ), x)
  else ! a>=2
    atan_m6 = atan_series( 1.0 / a )
    atan_m6 = sign( 0.5*pi - atan_m6, x)
  endif

  if (use_fortran_intrinsics) atan_m6 = atan(x)

  contains

  !> Returns the arctangent of x if x is in range -3/4..3/4
  !! This is calculated using Horner's method for simplicity.
  real elemental function atan_series(x)
    real, intent(in) :: x !< Argument of arctangent in range -3/4..3/4
    ! Local variables
    integer, parameter :: n = 56 ! Number of terms (could be dynamically determined)
    real :: x2 ! x**2
    real :: r ! accumulated terms
    real :: d ! polynomial coefficient
    integer :: j ! term number

    x2 = min(1.0, x*x)
    r = 1.0 / float(2*n-1)
    do j = n-1, 1, -1
      d = 1.0 / float(2*j-1)
      r = d - x2 * r
    enddo

    atan_series = r * x

  end function atan_series

  !> Returns pi/4 + arctan ( x / ( 2 + x) )
  real elemental function p4atx(x)
    real, intent(in) :: x !< Argument of function
    real :: d ! reciprocal of denominator

    d = 1.0 / ( 2.0 + x )
    p4atx = 0.25*pi + atan_series( x*d )

  end function p4atx

end function atan_m6

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

  ! Tangent tests
  if (verbose) write(stdout,*) 'Tests of tan()'
  fail = test(fail, tan_m6(0.0), 0., 'tan(0)')
  fail = test(fail, tan_m6(0.0625*pi), tan(0.0625*pi), 'tan(pi/16) v. lib')
  fail = test(fail, tan_m6(0.125*pi), tan(0.125*pi), 'tan(pi/8) v. lib')
  fail = test(fail, tan_m6(0.1875*pi), tan(0.1875*pi), 'tan(3/16*pi) v. lib', inexact=1.)
  fail = test(fail, tan_m6(pi/3.0), sqrt(3.0), 'tan(pi/6)=1/sqrt(3)', inexact=1.)
  fail = test(fail, tan_m6(0.25*pi), 1.0, 'tan(pi/4)=1', inexact=1.)
  fail = test(fail, tan_m6(pi/3.0), sqrt(3.0), 'tan(pi/3)=sqrt(3)', inexact=1.)
  fail = test(fail, tan_m6(0.375*pi), tan(0.375*pi), 'tan(3/8 pi) v. lib', inexact=2.)
  fail = test(fail, tan_m6(0.49*pi), tan(0.49*pi), 'tan(49/50 pi) v. lib', inexact=36.)

  ! Arc-tangent tests
  if (verbose) write(stdout,*) 'Tests of atan()'
  fail = test(fail, atan_m6(1.0 / sqrt(3.0)), pi/6.0, 'atan(3**-0.5)')
  fail = test(fail, atan_m6(0.0), 0., 'atan(0)')
  fail = test(fail, atan_m6(0.125), atan(0.125), 'atan(1/8)')
  fail = test(fail, atan_m6(0.75), atan(0.75), 'atan(3/4)')
  fail = test(fail, atan_m6(1.0), 0.25*pi, 'atan(1)')
  fail = test(fail, atan_m6(1.5), atan(1.5), 'atan(3/2)')
  fail = test(fail, atan_m6(4.0), atan(4.0), 'atan(4)', inexact=1.)
  fail = test(fail, atan_m6(-1.0), -0.25*pi, 'atan(-1)')

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
