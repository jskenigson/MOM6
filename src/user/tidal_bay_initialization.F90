module tidal_bay_initialization
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_NONE, OBC_segment_type
use MOM_open_boundary,  only : OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary,  only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, set_time, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

public tidal_bay_set_OBC_data

contains

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine tidal_bay_set_OBC_data(OBC, G, GV, h, eta, Time)
  type(ocean_OBC_type),   pointer    :: OBC  !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used.
  type(ocean_grid_type),  intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV  !<  Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h !< layer thickness.
  real, dimension(SZI_(G),SZJ_(G)),         intent(in) :: eta !< free surface height (m)
  type(time_type),        intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the tidal_bay example.
  real :: time_sec, cff, cff2, tide_flow
  real :: my_area, my_flux
  real :: PI
  character(len=40)  :: mod = "tidal_bay_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, n, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment
  integer :: ni_seg, nj_seg   ! number of src gridpoints along the segments
  integer :: i2, j2           ! indices for referencing local domain array
  integer :: ishift, jshift   ! offsets for staggered locations

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  PI = 4.0*atan(1.0) ;

  if (.not.associated(OBC)) return

  ! Should be doing sum reduction if open boundary on more than one PE.
  time_sec = time_type_to_real(Time)
  cff = 0.1*sin(2.0*PI*time_sec/(12.0*3600.0))
  tide_flow = 3.0e6
  my_area=0.0
  my_flux=0.0
  ! need to rewrite along segment
  do j=jsd,jed ; do I=IsdB,IedB
    if (OBC%OBC_segment_u(I,j) /= OBC_NONE) then
      do k=1,nz
        cff2 = h(I,j,k)*G%dyCu(I,j)
        my_area = my_area + cff2
      enddo
    endif
  enddo ; enddo
  my_flux = -tide_flow*SIN(2.0*PI*time_sec/(12.0*3600.0))

  do n = 1, OBC%number_of_segments
    segment => OBC%OBC_segment_number(n)

    if (.not. segment%on_pe) cycle ! continue to next segment if not in computational domain

    ! Segment indices are on q points - segments are one wide and long enough
    ! for all the internal velocity points. This give something that's actually
    ! one too long in the long direction.
    ni_seg = segment%Ie_obc-segment%Is_obc + 1
    nj_seg = segment%Je_obc-segment%Js_obc + 1

    if (.not. ASSOCIATED(segment%Cg)) then ! finishing allocating storage for segments
      allocate(segment%Cg(ni_seg,nj_seg));       segment%Cg(:,:) = 0.0
      allocate(segment%Htot(ni_seg,nj_seg));     segment%Htot(:,:) = 0.0
      allocate(segment%h(ni_seg,nj_seg,G%ke));   segment%h(:,:,:) = 0.0
      allocate(segment%unh(ni_seg,nj_seg,G%ke)); segment%unh(:,:,:) = 0.0
      allocate(segment%unhbt(ni_seg,nj_seg));    segment%unhbt(:,:) = 0.0
      allocate(segment%unbt(ni_seg,nj_seg));     segment%unbt(:,:) = my_flux/my_area
      allocate(segment%eta(ni_seg,nj_seg));      segment%eta(:,:) = cff
    endif

    ishift = 0; jshift = 0
    if (segment%direction == OBC_DIRECTION_W .or. segment%direction == OBC_DIRECTION_E) &
        ishift = 1
    if (segment%direction == OBC_DIRECTION_S .or. segment%direction == OBC_DIRECTION_N) &
        jshift = 1

    do j=1,nj_seg
      do i=1,ni_seg
        i2 =  segment%Is_obc + i - 1
        j2 =  segment%Js_obc + j - 1
        if ((i2 .gt. ied .or. i2 .lt. isd) .or. (j2 .gt. jed .or. j2 .lt. jsd)) cycle
        segment%Cg(i,j) = sqrt(GV%g_prime(1)*(0.5*(G%bathyT(i2,j2) + G%bathyT(i2+ishift,j2+jshift))))
        if (GV%Boussinesq) then
           segment%Htot(i,j) = 0.5*((G%bathyT(i2,j2)*GV%m_to_H + eta(i2,j2)) + &
                              (G%bathyT(i2+ishift,j2+jshift)*GV%m_to_H + eta(i2+ishift,j2+jshift)))
        else
           segment%Htot(i,j) = 0.5*(eta(i2,j2)+eta(i2+ishift,j2+jshift))
        endif
        do k=1,G%ke
           segment%h(i,j,k) = 0.5*(h(i2,j2,k)+h(i2+ishift,j2+jshift,k))
        enddo
      enddo
    enddo
  enddo ! end segment loop

end subroutine tidal_bay_set_OBC_data

!> \class tidal_bay_Initialization
!!
!! The module configures the model for the "tidal_bay" experiment.
!! tidal_bay = Tidally resonant bay from Zygmunt Kowalik's class on tides.
end module tidal_bay_initialization
