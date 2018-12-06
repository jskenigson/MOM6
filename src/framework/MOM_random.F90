!> Provides gridded random number capability
module MOM_random

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_domains,       only : MOM_domain_type
use MOM_time_manager,  only : time_type, set_date, get_date

use MersenneTwister_mod, only : randomNumberSequence ! Random number class from FMS
use MersenneTwister_mod, only : new_RandomNumberSequence ! Constructor/initializer
use MersenneTwister_mod, only : getRandomReal ! Generates a random number

implicit none ; private

public :: random_construct
public :: random_construct_2d
public :: random_01
public :: random_norm
public :: random_01_2d
public :: random_norm_2d
public :: random_unit_tests

#include <MOM_memory.h>

!> Container for pseudo-random number generators
type, public :: PRNG ; private

  type(randomNumberSequence) :: stream0d !< Scalar random number generator for whole model
  type(randomNumberSequence), dimension(:,:), allocatable :: stream2d !< Random number generator for each cell on horizontal grid

end type PRNG

contains

!> Returns a real random number between 0 and 1
real function random_01(CS)
  type(PRNG), intent(inout) :: CS !< Container for pseudo-random number generators

  random_01 = getRandomReal(CS%stream0d)

end function random_01

!> Returns a normally distributed (approx) real random number between 0 and 1
real function random_norm(CS)
  type(PRNG), intent(inout) :: CS !< Container for pseudo-random number generators
  ! Local variables
  integer :: i

  random_norm = getRandomReal(CS%stream0d)
  do i = 1,11
    random_norm = random_norm + getRandomReal(CS%stream0d)
  enddo
  random_norm = random_norm - 6.

end function random_norm

!> Generates random numbers for each cell of the model grid
subroutine random_01_2d(CS, G, rand)
  type(PRNG),            intent(inout) :: CS   !< Container for pseudo-random number generators
  type(ocean_grid_type), intent(in)    :: G    !< Ocean grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: rand !< Random numbers between 0 and 1
  ! Local variables
  integer :: i,j

  do j = G%jsd,G%jed
    do i = G%isd,G%ied
      rand(i,j) = getRandomReal( CS%stream2d(i,j) )
    enddo
  enddo

end subroutine random_01_2d

!> Generates normally distributed random numbers for each cell of the model grid
subroutine random_norm_2d(CS, G, rand)
  type(PRNG),            intent(inout) :: CS   !< Container for pseudo-random number generators
  type(ocean_grid_type), intent(in)    :: G    !< Ocean grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: rand !< Random numbers between 0 and 1
  ! Local variables
  integer :: i,j,n

  do j = G%jsd,G%jed
    do i = G%isd,G%ied
      rand(i,j) = getRandomReal( CS%stream2d(i,j) ) - 0.5
    enddo
    do n = 1,11
      do i = G%isd,G%ied
        rand(i,j) = rand(i,j) + ( getRandomReal( CS%stream2d(i,j) ) - 0.5 )
      enddo
    enddo
  enddo

end subroutine random_norm_2d

!> Constructor for scalar PRNG. Also can be used to reset to seed.
subroutine random_construct(CS, Time, seed)
  type(PRNG),      intent(inout) :: CS !< Container for pseudo-random number generators
  type(time_type), intent(in)    :: Time !< Current model time
  integer,         intent(in)    :: seed !< Seed for PRNG
  ! Local variables
  integer :: tseed

  tseed = seed_from_time(Time)
  tseed = ieor(tseed, seed)
  CS%stream0d = new_RandomNumberSequence(tseed)

end subroutine random_construct

!> Constructor for gridded PRNG
subroutine random_construct_2d(CS, G, Time, seed)
  type(PRNG),            intent(inout) :: CS   !< Container for pseudo-random number generators
  type(ocean_grid_type), intent(in)    :: G    !< Ocean grid structure
  type(time_type),       intent(in)    :: Time !< Current model time
  integer,               intent(in)    :: seed !< Seed for PRNG
  ! Local variables
  integer :: i,j,sseed,tseed

  if (.not. allocated(CS%stream2d)) allocate( CS%stream2d(G%isd:G%ied,G%jsd:G%jed) )

  tseed = seed_from_time(Time)
  tseed = ieor(tseed*9007, seed)
  do j = G%jsd,G%jed
    do i = G%isd,G%ied
      sseed = seed_from_space(G, i, j)
      sseed = ieor(tseed, sseed*7993)
      CS%stream2d(i,j) = new_RandomNumberSequence(sseed)
    enddo
  enddo

end subroutine random_construct_2d

!> Create seed from time
integer function seed_from_time(Time)
  type(time_type), intent(in)    :: Time !< Current model time
  ! Local variables
  integer :: yr,mo,dy,hr,mn,sc

  call get_date(Time,yr,mo,dy,hr,mn,sc)
  seed_from_time = sc + 61*(mn + 61*hr) + 379
  seed_from_time = mod(dy + 32*(mo + 13*yr), 2147483648) + int(43200*(1.+sin(real(seed_from_time))))
  
end function seed_from_time

!> Create seed from position index
!!
!! \note It would be cleaner to use the hor_index_type instead of the ocean_grod_type
!! but unfortunately the former does not have all the information about periodicity and
!! global size.
integer function seed_from_space(G, ic, jc)
  type(ocean_grid_type), intent(in) :: G  !< Ocean grid structure
  integer,               intent(in) :: ic !< Computational i-index
  integer,               intent(in) :: jc !< Computational j-index
  ! Local variables
  integer :: ig, jg, ni, nj, ij

  ni = G%Domain%niglobal
  nj = G%Domain%njglobal
  ! Periodicity is assumed here but does not break non-periodic models
  ig = mod(G%HI%idg_offset + ic - 1 + ni, ni)+1
  jg = max(G%HI%jdg_offset + jc, 0)
  if (jg>nj) then ! Tri-polar hard-coded until we put needed info in HI
    jg = 2*nj+1-jg
    ig = ni+1-ig
  endif
  seed_from_space = ig + ni*(jg-1)

end function seed_from_space

!> Destructor for PRNG
subroutine random_end(CS)
  type(PRNG), pointer :: CS !< Container for pseudo-random number generators
  if (allocated(CS%stream2d)) deallocate(CS%stream2d)
  deallocate(CS)
end subroutine random_end

!> Runs some statistical tests on the PRNG
logical function random_unit_tests(verbose)
  logical :: verbose !< True if results should be written to stdout
  ! Local variables
  type(PRNG) :: test ! Generator
  type(time_type) :: Time ! Model time
  real :: r1, r2, r3 ! Some random numbers
  real :: mean, var, ar1, std ! Some statistics
  integer :: stdunit ! For messages
  integer, parameter :: n = 600
  integer :: i,j,ni,nj
  ! Fake being on a decomposed domain
  type(ocean_grid_type), pointer :: G => null() !< Not the real grid
  type(MOM_domain_type), pointer :: Domain => null() !< Not the real domain
  real, dimension(:,:), allocatable :: r2d ! Random numbers

  random_unit_tests = .false.
  stdunit = 6
  write(stdunit,'(1x,a)') '==== MOM_random: randum_unit_tests ======================='

  ! Check time-based seed
  Time = set_date(1903, 11, 21, 13, 47, 29)
  i = seed_from_time(Time)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, i>0, 'time seed', ivalue=i)
  Time = set_date(1903, 11, 22, 13, 47, 29)
  i = seed_from_time(Time)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, i>0, 'time seed', ivalue=i)
  Time = set_date(1903, 11, 21, 13, 47, 30)
  i = seed_from_time(Time)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, i>0, 'time seed', ivalue=i)

  ! Generate a random number, r1
  call random_construct(test, Time, 1)
  r1 = random_01(test) 
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r1)/=0., 'first call', r1)

  ! Check that we get a different number on a second call
  r2 = random_01(test) 
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r2-r1)>0., 'consecutive test', r2)

  ! Check that we can reproduce r1
  call random_construct(test, Time, 1)
  r3 = random_01(test) 
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r3-r1)==0., 'reproduce test', r2)

  ! Check that we get a different number for a different date
  Time = set_date(1903, 11, 21, 13, 0, 29)
  call random_construct(test, Time, 1)
  r2 = random_01(test) 
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r2-r1)>0., 'date test', r2)

  ! Check statistics of large samples for uniform generator
  mean = 0. ; std = 0. ; r2 = 0.
  do i = 1, n 
    r1 = random_01(test)
    mean = mean + (r1 - 0.5)
    var = var + (r1 - 0.5)**2
    ar1 = ar1 + (r1 - 0.5)*(r2 - 0.5)
    r2 = r1 ! Keep copy of last value
  enddo
  mean = mean / real(n)
  var = var / real(n)
  ar1 = ar1 / real(n)
  std = sqrt(var)
  r3 = 1./sqrt(real(12*n)) ! Standard error of mean
  r2 = mean*sqrt(real(12*n)) ! Non-dimensional error in mean
  r3 = std*sqrt(12.) ! Non-dimensional standard deviation
  r1 = ar1 / var * sqrt(real(n))
  if (verbose) then
    write(stdunit,'(1x,"random: ",a)') 'Uniform 0..1 generator'
    write(stdunit,'(1x,"random: ",a,f12.9)') 'mean =',mean,'std =',std,'AR1 =',ar1
    write(stdunit,'(1x,"random: ",a,f12.9)') 'non-dim. error in mean =',r2,'non-dim. standard deviation =',r3,'non-dim. AR1 =',r1
  endif
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r2)<2., 'n>>1, mean within 2 sigma [uniform]', r2)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r3-1.)<1./sqrt(real(n)), 'n>>1, std ~ 1/sqrt(12) [uniform]', r3-1.)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r1)<2., 'n>>1, AR1 < std/sqrt(n) [uniform]', r1)

  ! Check statistics of large samples for normal generator
  mean = 0. ; std = 0. ; r2 = 0.
  do i = 1, n 
    r1 = random_norm(test)
    mean = mean + r1
    var = var + r1**2
    ar1 = ar1 + r1*r2
    r2 = r1 ! Keep copy of last value
  enddo
  mean = mean / real(n)
  var = var / real(n)
  ar1 = ar1 / real(n)
  std = sqrt(var)
  r3 = 1./sqrt(real(n)) ! Standard error of mean
  r2 = mean*sqrt(real(n)) ! Non-dimensional error in mean
  r3 = std ! Non-dimensional standard deviation
  r1 = ar1 / var * sqrt(real(n))
  if (verbose) then
    write(stdunit,'(1x,"random: ",a)') 'Normal generator'
    write(stdunit,'(1x,"random: ",a,f12.9)') 'mean =',mean,'std =',std,'AR1 =',ar1
    write(stdunit,'(1x,"random: ",a,f12.9)') 'non-dim. error in mean =',r2,'non-dim. standard deviation =',r3,'non-dim. AR1 =',r1
  endif
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r2)<2., 'n>>1, mean within 2 sigma [norm]', r2)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r3-1.)<1./sqrt(real(n)), 'n>>1, std ~ 1/sqrt(12) [norm]', r3-1.)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r1)<2., 'n>>1, AR1 < std/sqrt(n) [norm]', r1)

  ! Now fake being on a decomposed domain
  ni = 6
  nj = 9
  allocate(Domain)
  Domain%niglobal = ni
  Domain%njglobal = nj
  allocate(G)
  G%Domain => Domain
  G%isd = 0
  G%ied = ni+1
  G%jsd = 0
  G%jed = nj+1
  G%HI%idg_offset = 0
  G%HI%jdg_offset = 0

  ! Check space-based seed
  i = seed_from_space(G,1,1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, i==1, 'index seed (1,1)', ivalue=i)
  j = seed_from_space(G,ni+1,1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, j==i, 'index seed (n+1,1)', ivalue=j)
  i = seed_from_space(G,ni,1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, i==ni, 'index seed (n,1)', ivalue=i)
  j = seed_from_space(G,0,1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, j==i, 'index seed (0,1)', ivalue=j)
  i = seed_from_space(G,1,nj)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, i==ni*nj-ni+1, 'index seed (1,n)', ivalue=i)
  j = seed_from_space(G,ni,nj+1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, j==i, 'index seed (n,n+1)', ivalue=j)
  i = seed_from_space(G,ni,nj)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, i==ni*nj, 'index seed (1,n)', ivalue=i)
  j = seed_from_space(G,1,nj+1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, j==i, 'index seed (n,n+1)', ivalue=j)

  ! Check 2d random number generator 0..1
  allocate( r2d(G%isd:G%ied,G%jsd:G%jed) )
  call random_construct_2d(test, G, Time, 111)
  r2d(:,:) = -9. ! Use -9. to detect unset values
  call random_01_2d(test, G, r2d)
  r1 = minval(r2d)
  r2 = maxval(r2d)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, r1>=0., '2d all set', r1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, r2<=1., '2d all valid', r2)
  mean = sum( r2d(1:ni,1:nj) - 0.5 )/real(ni*nj)
  var = sum( (r2d(1:ni,1:nj) - 0.5 - mean)**2 )/real(ni*nj)
  std = sqrt(var)
  r3 = 1./sqrt(real(12*ni*nj)) ! Standard error of mean
  r2 = mean*sqrt(real(12*ni*nj)) ! Non-dimensional error in mean
  r3 = std*sqrt(12.) ! Non-dimensional standard deviation
  if (verbose) then
    write(stdunit,'(1x,"random: ",a)') '2D uniform 0..1 generator'
    write(stdunit,'(1x,"random: ",a,f12.9)') 'mean =',mean,'std =',std
    write(stdunit,'(1x,"random: ",a,f12.9)') 'non-dim. error in mean =',r2,'non-dim. standard deviation =',r3
  endif
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r2)<2., '2d, mean within 2 sigma [uniform]', r2)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r3-1.)<1./sqrt(real(n)), '2d, std ~ 1/sqrt(12) [uniform]', r3-1.)
  if (verbose) then
    write(stdunit,'(1x,"random:")')
    write(stdunit,'(1x,"random:",8f8.5)') r2d
    write(stdunit,'(1x,"random:")')
  endif

  ! Check 2d normal random number generator
  call random_norm_2d(test, G, r2d)
  mean = sum( r2d(1:ni,1:nj) )/real(ni*nj)
  var = sum( r2d(1:ni,1:nj)**2 )/real(ni*nj)
  std = sqrt(var)
  r3 = 1./sqrt(real(ni*nj)) ! Standard error of mean
  r2 = mean*sqrt(real(ni*nj)) ! Non-dimensional error in mean
  r3 = std ! Non-dimensional standard deviation
  if (verbose) then
    write(stdunit,'(1x,"random: ",a)') '2D normal generator'
    write(stdunit,'(1x,"random: ",a,f12.9)') 'mean =',mean,'std =',std
    write(stdunit,'(1x,"random: ",a,f12.9)') 'non-dim. error in mean =',r2,'non-dim. standard deviation =',r3
  endif
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r2)<2., '2d, mean within 2 sigma [norm]', r2)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, abs(r3-1.)<1./sqrt(real(n)), '2d, std ~ 1/sqrt(12) [norm]', r3-1.)

  ! Clean up
  deallocate(r2d)
  deallocate(Domain)
  deallocate(G)

  if (.not.random_unit_tests) write(stdunit,'(1x,a)') 'Pass'

end function random_unit_tests

!> Convenience function for reporting result of test
logical function test_fn(verbose, good, label, rvalue, ivalue)
  logical, intent(in) :: verbose !< Verbosity
  logical, intent(in) :: good !< True if pass, false otherwise
  character(len=*), intent(in) :: label !< Label for messages
  real,    intent(in) :: rvalue !< Result of calculation
  integer, intent(in) :: ivalue !< Result of calculation
  optional :: rvalue, ivalue

  if (present(ivalue)) then
    if (.not. good) then
      write(6,'(1x,a,i9,1x,a,a)') 'random: result =',ivalue,label,' <------- FAIL!'
    elseif (verbose) then
      write(6,'(1x,a,i9,1x,a)') 'random: result =',ivalue,label
    endif
  else
    if (.not. good) then
      write(6,'(1x,a,1pe15.8,1x,a,a)') 'random: result =',rvalue,label,' <------- FAIL!'
    elseif (verbose) then
      write(6,'(1x,a,1pe15.8,1x,a)') 'random: result =',rvalue,label
    endif
  endif
  test_fn = .not. good
  
end function test_fn

end module MOM_random
