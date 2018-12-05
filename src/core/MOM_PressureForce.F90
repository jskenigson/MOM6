!> A thin wrapper for Boussinesq/non-Boussinesq forms of the pressure force calculation.
module MOM_PressureForce

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : diag_ctrl, time_type, post_data, register_diag_field
use MOM_EOS, only : calculate_density
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_PressureForce_AFV, only : PressureForce_AFV_Bouss, PressureForce_AFV_nonBouss
use MOM_PressureForce_AFV, only : PressureForce_AFV_init, PressureForce_AFV_end
use MOM_PressureForce_AFV, only : PressureForce_AFV_CS
use MOM_PressureForce_blk_AFV, only : PressureForce_blk_AFV_Bouss, PressureForce_blk_AFV_nonBouss
use MOM_PressureForce_blk_AFV, only : PressureForce_blk_AFV_init, PressureForce_blk_AFV_end
use MOM_PressureForce_blk_AFV, only : PressureForce_blk_AFV_CS
use MOM_PressureForce_Mont, only : PressureForce_Mont_Bouss, PressureForce_Mont_nonBouss
use MOM_PressureForce_Mont, only : PressureForce_Mont_init, PressureForce_Mont_end
use MOM_PressureForce_Mont, only : PressureForce_Mont_CS
use MOM_tidal_forcing, only : tidal_forcing_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_ALE, only: ALE_CS
implicit none ; private

#include <MOM_memory.h>

public PressureForce, PressureForce_init, PressureForce_end

!> Pressure force control structure
type, public :: PressureForce_CS ; private
  logical :: Analytic_FV_PGF !< If true, use the analytic finite volume form
                             !! (Adcroft et al., Ocean Mod. 2008) of the PGF.
  logical :: blocked_AFV     !< If true, used the blocked version of the ANALYTIC_FV_PGF
                             !! code.  The value of this parameter should not change answers.
  real :: Brankart_factor    !< Scaling for mean contribution to Brankart effect.
  real :: Brankart_noise     !< Scaling for random contribution to Brankart effect.
  !> Control structure for the analytically integrated finite volume pressure force
  type(PressureForce_AFV_CS), pointer :: PressureForce_AFV_CSp => NULL()
  !> Control structure for the analytically integrated finite volume pressure force
  type(PressureForce_blk_AFV_CS), pointer :: PressureForce_blk_AFV_CSp => NULL()
  !> Control structure for the Montgomery potential form of pressure force
  type(PressureForce_Mont_CS), pointer :: PressureForce_Mont_CSp => NULL()
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the timing of diagnostic output.

  !>@{ Diagnostic IDs
  integer :: id_brankart_anom = -1, id_brankart_pfu = -1, id_brankart_pfv = -1
  !!@}
end type PressureForce_CS

contains

!> A thin layer between the model and the Boussinesq and non-Boussinesq pressure force routines.
subroutine PressureForce(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, p_atm, pbce, eta)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)  :: h    !< Layer thicknesses, in H (usually m or kg m-2)
  type(thermo_var_ptrs),   intent(in)  :: tv   !< A structure pointing to various thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(out) :: PFu  !< Zonal pressure force acceleration (m/s2)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(out) :: PFv  !< Meridional pressure force acceleration (m/s2)
  type(PressureForce_CS),  pointer     :: CS   !< Pressure force control structure
  type(ALE_CS),            pointer     :: ALE_CSp !< ALE control structure
  real, dimension(:,:), &
                 optional, pointer     :: p_atm !< The pressure at the ice-ocean or
                                               !! atmosphere-ocean interface in Pa.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                 optional, intent(out) :: pbce !< The baroclinic pressure anomaly in each layer
                                               !! due to eta anomalies, in m2 s-2 H-1.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(out) :: eta  !< The bottom mass used to calculate PFu and PFv,
                                               !! in H, with any tidal contributions.
  ! Local variables
  type(thermo_var_ptrs) :: tv_tmp ! A temporary work space for spoofing T and S
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), target :: Ttmp, Stmp ! Alternatives for T and S
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: delta_T, delta_S ! Differences for T and S
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: PFu2 ! Alternative PFu
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: PFv2 ! Alternative PFv
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: pbce2 ! Alternative pbce
  real, dimension(SZI_(G),SZJ_(G)) :: eta2 ! Alternative eta
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: rho1,rho2 ! For diagnostics
  logical :: do_Brankart

  tv_tmp%eqn_of_state => tv%eqn_of_state
  tv_tmp%P_ref = tv%P_ref
  do_Brankart = .false.
  if (CS%Brankart_factor>0.) do_Brankart = .true.
  if (do_Brankart) then
    Ttmp(:,:,:) = 0.
    Stmp(:,:,:) = 0.
    PFu2(:,:,:) = 0.
    PFv2(:,:,:) = 0.
    pbce2(:,:,:) = 0.
    eta2(:,:) = 0.
    rho1(:,:,:) = 0.
    rho2(:,:,:) = 0.
    tv_tmp%T => Ttmp
    tv_tmp%S => Stmp
  else
    tv_tmp%T => tv%T
    tv_tmp%S => tv%S
  endif

  if (do_Brankart) then
    call lsfit_TS(CS, G, tv, h, tv%T, tv%S, delta_T, delta_S)
    Stmp(:,:,:) = max(0., tv%S(:,:,:) + delta_S(:,:,:) )
    Ttmp(:,:,:) = tv%T(:,:,:) + delta_T(:,:,:)
    if (CS%id_brankart_anom>0) then
      call density_from_TS(G, GV, tv, h, Ttmp, Stmp, rho1)
    endif
    if (CS%Analytic_FV_PGF .and. CS%blocked_AFV) then
      if (GV%Boussinesq) then
        call PressureForce_blk_AFV_Bouss(h, tv_tmp, PFu2, PFv2, G, GV, US, &
                 CS%PressureForce_blk_AFV_CSp, ALE_CSp, p_atm, pbce2, eta2)
      else
        call PressureForce_blk_AFV_nonBouss(h, tv_tmp, PFu2, PFv2, G, GV, US, &
                 CS%PressureForce_blk_AFV_CSp, p_atm, pbce2, eta2)
      endif
    elseif (CS%Analytic_FV_PGF) then
      if (GV%Boussinesq) then
        call PressureForce_AFV_Bouss(h, tv_tmp, PFu2, PFv2, G, GV, US, CS%PressureForce_AFV_CSp, &
                                     ALE_CSp, p_atm, pbce2, eta2)
      else
        call PressureForce_AFV_nonBouss(h, tv_tmp, PFu2, PFv2, G, GV, US, CS%PressureForce_AFV_CSp, &
                                        ALE_CSp, p_atm, pbce2, eta2)
      endif
    else
      if (GV%Boussinesq) then
        call PressureForce_Mont_Bouss(h, tv_tmp, PFu2, PFv2, G, GV, US, CS%PressureForce_Mont_CSp, &
                                      p_atm, pbce2, eta2)
      else
        call PressureForce_Mont_nonBouss(h, tv_tmp, PFu2, PFv2, G, GV, US, CS%PressureForce_Mont_CSp, &
                                         p_atm, pbce2, eta2)
      endif
    endif
    Stmp(:,:,:) = max(0., tv%S(:,:,:) - delta_S(:,:,:) )
    Ttmp(:,:,:) = tv%T(:,:,:) - delta_T(:,:,:)
    if (CS%id_brankart_anom>0) then
      call density_from_TS(G, GV, tv, h, Ttmp, Stmp, rho2)
      rho1 = 0.5 * ( rho1(:,:,:) + rho2(:,:,:) ) ! < Brankart density, BAR[ rho(T,S) ]
      call density_from_TS(G, GV, tv, h, tv%T, tv%S, rho2) ! Normal density rho(BAR[T],BAR[S])
      rho2 = rho1(:,:,:) - rho2(:,:,:)
      call post_data(CS%id_brankart_anom, rho2, CS%diag)
    endif
  endif

  if (CS%Analytic_FV_PGF .and. CS%blocked_AFV) then
    if (GV%Boussinesq) then
      call PressureForce_blk_AFV_Bouss(h, tv_tmp, PFu, PFv, G, GV, US, &
               CS%PressureForce_blk_AFV_CSp, ALE_CSp, p_atm, pbce, eta)
    else
      call PressureForce_blk_AFV_nonBouss(h, tv_tmp, PFu, PFv, G, GV, US, &
               CS%PressureForce_blk_AFV_CSp, p_atm, pbce, eta)
    endif
  elseif (CS%Analytic_FV_PGF) then
    if (GV%Boussinesq) then
      call PressureForce_AFV_Bouss(h, tv_tmp, PFu, PFv, G, GV, US, CS%PressureForce_AFV_CSp, &
                                   ALE_CSp, p_atm, pbce, eta)
    else
      call PressureForce_AFV_nonBouss(h, tv_tmp, PFu, PFv, G, GV, US, CS%PressureForce_AFV_CSp, &
                                      ALE_CSp, p_atm, pbce, eta)
    endif
  else
    if (GV%Boussinesq) then
      call PressureForce_Mont_Bouss(h, tv_tmp, PFu, PFv, G, GV, US, CS%PressureForce_Mont_CSp, &
                                    p_atm, pbce, eta)
    else
      call PressureForce_Mont_nonBouss(h, tv_tmp, PFu, PFv, G, GV, US, CS%PressureForce_Mont_CSp, &
                                       p_atm, pbce, eta)
    endif
  endif

  if (do_Brankart) then
    PFu(:,:,:)= 0.5 * ( PFu(:,:,:) + PFu2(:,:,:) )
    PFv(:,:,:)= 0.5 * ( PFv(:,:,:) + PFv2(:,:,:) )
    if (present(pbce)) pbce(:,:,:)= 0.5 * ( pbce(:,:,:) + pbce2(:,:,:) )
    if (present(eta)) eta(:,:)= 0.5 * ( eta(:,:) + eta2(:,:) )
    if (CS%id_brankart_pfu>0 .or. CS%id_brankart_pfv>0) then
      if (CS%Analytic_FV_PGF .and. CS%blocked_AFV) then
        if (GV%Boussinesq) then
          call PressureForce_blk_AFV_Bouss(h, tv, PFu2, PFv2, G, GV, US, &
                   CS%PressureForce_blk_AFV_CSp, ALE_CSp, p_atm, pbce2, eta2)
        else
          call PressureForce_blk_AFV_nonBouss(h, tv, PFu2, PFv2, G, GV, US, &
                   CS%PressureForce_blk_AFV_CSp, p_atm, pbce2, eta2)
        endif
      elseif (CS%Analytic_FV_PGF) then
        if (GV%Boussinesq) then
          call PressureForce_AFV_Bouss(h, tv, PFu2, PFv2, G, GV, US, CS%PressureForce_AFV_CSp, &
                                       ALE_CSp, p_atm, pbce2, eta2)
        else
          call PressureForce_AFV_nonBouss(h, tv, PFu2, PFv2, G, GV, US, CS%PressureForce_AFV_CSp, &
                                          ALE_CSp, p_atm, pbce2, eta2)
        endif
      else
        if (GV%Boussinesq) then
          call PressureForce_Mont_Bouss(h, tv, PFu2, PFv2, G, GV, US, CS%PressureForce_Mont_CSp, &
                                        p_atm, pbce2, eta2)
        else
          call PressureForce_Mont_nonBouss(h, tv, PFu2, PFv2, G, GV, US, CS%PressureForce_Mont_CSp, &
                                           p_atm, pbce2, eta2)
        endif
      endif
      if (CS%id_brankart_pfu>0) then
        PFu2(:,:,:) = PFu(:,:,:) - PFu2(:,:,:)
        call post_data(CS%id_brankart_pfu, PFu2, CS%diag)
      endif
      if (CS%id_brankart_pfv>0) then
        PFv2(:,:,:) = PFv(:,:,:) - PFv2(:,:,:)
        call post_data(CS%id_brankart_pfv, PFv2, CS%diag)
      endif
    endif
  endif

end subroutine Pressureforce

!> Calculate a T/S difference for each cell to use in Brankart's density correction
subroutine lsfit_TS(CS, G, tv, h, T, S, delta_T, delta_S)
  type(PressureForce_CS), pointer     :: CS !< Pressure force control structure
  type(ocean_grid_type),  intent(in)  :: G  !< The ocean's grid structure
  type(thermo_var_ptrs),  intent(in)  :: tv !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h       !< Layer thicknesses, in H (usually m or kg m-2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: T       !< Temperature (degC)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: S       !< Salinity (1e-3)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: delta_T !< Temperature difference (degC)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: delta_S !< Salinity difference (1e-3)
  ! Local variables
  integer :: i,j,k
  real :: Tl(5), Sl(5), hl(5) ! Copies of local stencil
  real :: mn_S, mn_T, mn_S2, mn_ST, mn_T2, vr_S, vr_T, mn_H, cv_ST, sd_S, sd_T
  real :: dSdT, dTdS, del_T, del_S, scl, Tf

  Tf = -1.9

  do k = 1, G%ke
    do j = G%jsc-2,G%jec+2 ; do i = G%isc-2,G%iec+2
      ! Least squares fit to data in five-point stencil
      Tl(1) = T(i,j,k) ; Tl(2) = T(i-1,j,k) ; Tl(3) = T(i+1,j,k) ; Tl(4) = T(i,j-1,k) ; Tl(5) = T(i,j+1,k)
      Sl(1) = S(i,j,k) ; Sl(2) = S(i-1,j,k) ; Sl(3) = S(i+1,j,k) ; Sl(4) = S(i,j-1,k) ; Sl(5) = S(i,j+1,k)
      hl(1) = h(i,j,k) ; hl(2) = h(i-1,j,k) ; hl(3) = h(i+1,j,k) ; hl(4) = h(i,j-1,k) ; hl(5) = h(i,j+1,k)
      mn_H = hl(1) + ( ( hl(2) + hl(3) ) + ( hl(4) + hl(5) ) )
      if (mn_H>0.) mn_H = 1. / mn_H ! Hereafter, mn_H is the reciprocal of mean h for the stencil
      ! Mean of S,T
      mn_S = ( hl(1)*Sl(1) + ( ( hl(2)*Sl(2) + hl(3)*Sl(3) ) + ( hl(4)*Sl(4) + hl(5)*Sl(5) ) ) ) * mn_H
      mn_T = ( hl(1)*Tl(1) + ( ( hl(2)*Tl(2) + hl(3)*Tl(3) ) + ( hl(4)*Tl(4) + hl(5)*Tl(5) ) ) ) * mn_H
      ! Adjust S,T vectors to have zero mean
      Sl(:) = Sl(:) - mn_S ; mn_S = 0.
      Tl(:) = Tl(:) - mn_T ; mn_T = 0.
      ! Variance of S,T (or mean square error if mean is removed above)
      mn_S2 = ( hl(1)*Sl(1)*Sl(1) + ( ( hl(2)*Sl(2)*Sl(2) + hl(3)*Sl(3)*Sl(3) ) &
                                    + ( hl(4)*Sl(4)*Sl(4) + hl(5)*Sl(5)*Sl(5) ) ) ) * mn_H
      mn_T2 = ( hl(1)*Tl(1)*Tl(1) + ( ( hl(2)*Tl(2)*Tl(2) + hl(3)*Tl(3)*Tl(3) ) &
                                    + ( hl(4)*Tl(4)*Tl(4) + hl(5)*Tl(5)*Tl(5) ) ) ) * mn_H
      vr_S = max(0., mn_S2 - mn_S*mn_S) ! Variance should be positive but round-off can violate this. Calculating
      vr_T = max(0., mn_T2 - mn_T*mn_T) ! variance directly would fix this but require more operations.
      sd_S = sqrt( vr_S )
      sd_T = sqrt( vr_T )
      ! Covariance of S,T
      mn_ST = ( hl(1)*Sl(1)*Tl(1) + ( ( hl(2)*Sl(2)*Tl(2) + hl(3)*Sl(3)*Tl(3) ) &
                                    + ( hl(4)*Sl(4)*Tl(4) + hl(5)*Sl(5)*Tl(5) ) ) ) * mn_H
      cv_ST = mn_ST - mn_S*mn_T

      del_T = sd_T
      del_S = 0.
      if (sd_T > 0.) del_S = cv_ST/sd_T
      ! Scale down results to avoid spuriously large deltas
      scl = 1.
      if (abs(del_S)>sd_S) scl = sd_S/abs(del_S)
      ! Bound for fresh water
      if (CS%Brankart_factor*abs(del_S)>S(i,j,k)) scl = min( scl, S(i,j,k)/(CS%Brankart_factor*abs(del_S)) )
      ! Bound for freezing water (approximately because using Tf=0). Should really offset T by actual Tfreeze...)
      if (CS%Brankart_factor*del_T>T(i,j,k)-Tf.and.del_T>0.) scl = min( scl, max(T(i,j,k)-Tf,0.)/(CS%Brankart_factor*del_T) )
      del_S = scl * del_S
      del_T = scl * del_T
      delta_S(i,j,k) = del_S * CS%Brankart_factor
      delta_T(i,j,k) = del_T * CS%Brankart_factor
    enddo ; enddo
  enddo
end subroutine lsfit_TS

!> Calculate a density from an anomalous T/S
subroutine density_from_TS(G, GV, tv, h, T, S, rho)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV !< Vertical grid structure
  type(thermo_var_ptrs),   intent(in)  :: tv !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h   !< Layer thicknesses, in H (usually m or kg m-2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: T   !< Temperature (degC)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: S   !< Salinity (1e-3)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: rho !< Density (kg m-3)
  ! Local variables
  integer :: i,j,k,i0,ni
  real, dimension(SZI_(G)) :: p ! Pressure (Pa)

  i0 = G%isc - G%isd + 1
  ni = G%iec - G%isc + 1
  do j = G%jsc,G%jec
    p(:) = 0.
    do k = 1, G%ke
      ! Calculate pressure mid-layer
      do i = G%isc,G%iec
        p(i) = p(i) + h(i,j,k) * GV%H_to_Pa * 0.5
      enddo
      call calculate_density(T(:,j,k), S(:,j,k), p, rho(:,j,k), i0, ni, tv%eqn_of_state)
      ! Calculate pressure at bottom of layer to use as the top of the next layer
      do i = G%isc,G%iec
        p(i) = p(i) + h(i,j,k) * GV%H_to_Pa * 0.5
      enddo
    enddo
  enddo
end subroutine density_from_TS

!> Initialize the pressure force control structure
subroutine PressureForce_init(Time, G, GV, US, param_file, diag, CS, tides_CSp)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(PressureForce_CS),  pointer       :: CS   !< Pressure force control structure
  type(tidal_forcing_CS), optional, pointer :: tides_CSp !< Tide control structure
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_PressureForce" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "PressureForce_init called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ANALYTIC_FV_PGF", CS%Analytic_FV_PGF, &
                 "If true the pressure gradient forces are calculated \n"//&
                 "with a finite volume form that analytically integrates \n"//&
                 "the equations of state in pressure to avoid any \n"//&
                 "possibility of numerical thermobaric instability, as \n"//&
                 "described in Adcroft et al., O. Mod. (2008).", default=.true.)
  call get_param(param_file, mdl, "BLOCKED_ANALYTIC_FV_PGF", CS%blocked_AFV, &
                 "If true, used the blocked version of the ANALYTIC_FV_PGF \n"//&
                 "code.  The value of this parameter should not change answers.", &
                 default=.false., do_not_log=.true., debuggingParam=.true.)
  call get_param(param_file, mdl, "BRANKART_FACTOR", CS%Brankart_factor, &
                 "A non-dimensional coefficient to multiply local\n"//&
                 "T/S differences when incorporating unresolved\n"//&
                 "perturbations in the mean density.", default=0., units='nondom', do_not_log=.true.)
  call get_param(param_file, mdl, "BRANKART_NOISE", CS%Brankart_noise, &
                 "Amplitude of stochastic contribution to\n"//&
                 "T/S differences representing unresolved\n"//&
                 "perturbations in the mean density.", default=0., units='nondom', do_not_log=.true.)
  if (CS%Brankart_noise /= 0.) call MOM_error(FATAL, "interpret_eos_selection: "//&
                 "Brankart noise not implemented yet!")
  if (CS%Brankart_factor /= 0.) then
    CS%id_brankart_anom = register_diag_field('ocean_model', 'Brankart_anom', diag%axesTL, Time, &
                 'Brankart density anomoly', 'kg m-3')
    CS%id_brankart_pfu = register_diag_field('ocean_model', 'Brankart_PFu', diag%axesCuL, Time, &
                 'Brankart zonal pressure force anomoly', 'kg m-3')
    CS%id_brankart_pfv = register_diag_field('ocean_model', 'Brankart_PFv', diag%axesCvL, Time, &
                 'Brankart zonal pressure force anomoly', 'kg m-3')
  endif

  if (CS%Analytic_FV_PGF .and. CS%blocked_AFV) then
    call PressureForce_blk_AFV_init(Time, G, GV, US, param_file, diag, &
             CS%PressureForce_blk_AFV_CSp, tides_CSp)
  elseif (CS%Analytic_FV_PGF) then
    call PressureForce_AFV_init(Time, G, GV, US, param_file, diag, &
             CS%PressureForce_AFV_CSp, tides_CSp)
  else
    call PressureForce_Mont_init(Time, G, GV, US, param_file, diag, &
             CS%PressureForce_Mont_CSp, tides_CSp)
  endif

end subroutine PressureForce_init

!> Deallocate the pressure force control structure
subroutine PressureForce_end(CS)
  type(PressureForce_CS), pointer :: CS !< Pressure force control structure

  if (CS%Analytic_FV_PGF .and. CS%blocked_AFV) then
    call PressureForce_blk_AFV_end(CS%PressureForce_blk_AFV_CSp)
  elseif (CS%Analytic_FV_PGF) then
    call PressureForce_AFV_end(CS%PressureForce_AFV_CSp)
  else
    call PressureForce_Mont_end(CS%PressureForce_Mont_CSp)
  endif

  if (associated(CS)) deallocate(CS)
end subroutine PressureForce_end

!> \namespace mom_pressureforce
!!
!! This thin module provides a branch to two forms of the horizontal accelerations
!! due to pressure gradients. The two options currently available are a
!! Montgomery potential form (used in traditional isopycnal layer models), and the
!! analytic finite volume form.

end module MOM_PressureForce
