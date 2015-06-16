!> Applies counter-gradient re-stratification fluxes in the vertical
module MOM_counter_grad_restrat

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_EOS, only : calculate_density

implicit none ; private

#include <MOM_memory.h>

public counter_grad_restrat, counter_grad_restrat_init

type, public :: counter_grad_restrat_CS ; private
  real :: ml_restrat_coef  !<   A nondimensional factor by which the 
                           !! instability is enhanced over what would be
                           !! predicted based on the resolved  gradients.  This
                           !! increases with grid spacing^2, up to something
                           !! of order 500.
  real :: MLE_density_diff !< Density difference used in detecting mixed-layer
                           !! depth (kg/m3)
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                           !! timing of diagnostic output.
  integer :: id_MLD = -1, id_deltaT = -1, id_heat = -1
end type counter_grad_restrat_CS

contains

subroutine counter_grad_restrat(tv, h, dt, G, CS)
  type(thermo_var_ptrs),                  intent(inout) :: tv  !< Thermodynamics structure.
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h   !< Layer thicknesses, in m or kg/m2.
  real,                                   intent(in)    :: dt  !< Time-step, in s.
  type(ocean_grid_type),                  intent(in)    :: G   !< Grid structure.
  type(counter_grad_restrat_CS),          pointer       :: CS  !< Control structure for this module.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    MLD,       & ! The diagnosed MLD
    Rml_av,    & !
    htot,      & !
    deltaT,    & ! 
    heat,      & ! 
    T_av,      & ! 
    T_min,     & ! 
    T_max,     & ! 
    S_av,      & ! 
    S_min,     & ! 
    S_max        ! 
  real, dimension(SZI_(G)) :: dK, dKm1, deltaRhoAtK, deltaRhoAtKm1, p0, Rho0, rhoSurf, pref_MLD
  real :: Tlowest, Thighest, depthT, depthB, df, aFac, ddRho
  integer, dimension(SZI_(G),SZJ_(G)) :: kMLD ! Last level within the mixed-layer (0 for land)
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.associated(tv%eqn_of_state)) call MOM_error(FATAL, "MOM_counter_grad_restrat: "// &
         "An equation of state must be used with this module.")

  ! Calculate mixed-layer depth
  !! TODO: use derivatives and mid-MLD pressure. Currently this is sigma-0. -AJA
  pRef_MLD(:) = 0.
  do j = js-1, je+1
    dK(:) = 0.5 * h(:,j,1) ! Depth of center of surface layer
    call calculate_density(tv%T(:,j,1), tv%S(:,j,1), pRef_MLD, rhoSurf, is-1, ie-is+3, tv%eqn_of_state)
    deltaRhoAtK(:) = 0.
    MLD(:,j) = 0.
    kMLD(:,j) = 0
    do k = 2, nz
      dKm1(:) = dK(:) ! Depth of center of layer K-1
      dK(:) = dK(:) + 0.5 * ( h(:,j,k) + h(:,j,k-1) ) ! Depth of center of layer K

      ! Mixed-layer depth, using sigma-0 (surface reference pressure)
      deltaRhoAtKm1(:) = deltaRhoAtK(:) ! Store value from previous iteration of K
      call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pRef_MLD, deltaRhoAtK, is-1, ie-is+3, tv%eqn_of_state)
      deltaRhoAtK(:) = deltaRhoAtK(:) - rhoSurf(:) ! Density difference between layer K and surface
      do i = is-1, ie+1
        ddRho = deltaRhoAtK(i) - deltaRhoAtKm1(i)
        if ((MLD(i,j)==0.) .and. (ddRho>0.) .and. &
            (deltaRhoAtKm1(i)<CS%MLE_density_diff) .and. (deltaRhoAtK(i)>=CS%MLE_density_diff)) then
          aFac = ( CS%MLE_density_diff - deltaRhoAtKm1(i) ) / ddRho
          MLD(i,j) = dK(i) * aFac + dKm1(i) * (1. - aFac)
          kMLD(i,j) = k
        endif
      enddo ! i-loop

    enddo ! k-loop
    do i = is-1, ie+1
      if ((MLD(i,j)==0.) .and. (deltaRhoAtK(i)<CS%MLE_density_diff)) then
        MLD(i,j) = dK(i) ! Assume mixing to the bottom
        kMLD(i,j) = nz
      endif
    enddo
  enddo ! j-loop

  ! Column-wise bounds and averages
  do j=js-1,je+1
    Rml_av(:,j) = 0.
    T_av(:,j) = 0.
    T_min(:,j) = tv%T(:,j,1)
    T_max(:,j) = tv%T(:,j,1)
    S_av(:,j) = 0.
    S_min(:,j) = tv%S(:,j,1)
    S_max(:,j) = tv%S(:,j,1)
    htot(i,j) = 0.
    do k=1,nz
      call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p0,Rho0(:),is-1,ie-is+3,tv%eqn_of_state)
      do i=is-1,ie+1
        if (htot(i,j) < MLD(i,j)) then
          Rml_av(i,j) = Rml_av(i,j) + h(i,j,k)*Rho0(i)
          T_av(i,j) = T_av(i,j) + h(i,j,k)*tv%T(i,j,k)
          S_av(i,j) = S_av(i,j) + h(i,j,k)*tv%S(i,j,k)
          htot(i,j) = htot(i,j) + h(i,j,k)
          T_min(i,j) = min(T_min(i,j), tv%T(i,j,k))
          T_max(i,j) = max(T_max(i,j), tv%T(i,j,k))
          S_min(i,j) = min(S_min(i,j), tv%S(i,j,k))
          S_max(i,j) = max(S_max(i,j), tv%S(i,j,k))
        endif
      enddo
    enddo
    do i=is-1,ie+1
      Rml_av(i,j) = Rml_av(i,j) / (htot(i,j) + G%H_subroundoff)
      T_av(i,j) = T_av(i,j) / (htot(i,j) + G%H_subroundoff)
      S_av(i,j) = S_av(i,j) / (htot(i,j) + G%H_subroundoff)
    enddo
  enddo

  do j=js,je 
    do i=is,ie
      Tlowest = T_min(i,j)
      Thighest = T_max(i,j)
      if (kMLD(i-1,j)>0) then
        Tlowest = min(Tlowest, T_min(i-1,j))
        Thighest = max(Thighest, T_max(i-1,j))
      endif
      if (kMLD(i+1,j)>0) then
        Tlowest = min(Tlowest, T_min(i+1,j))
        Thighest = max(Thighest, T_max(i+1,j))
      endif
      if (kMLD(i,j-1)>0) then
        Tlowest = min(Tlowest, T_min(i,j-1))
        Thighest = max(Thighest, T_max(i,j-1))
      endif
      if (kMLD(i,j+1)>0) then
        Tlowest = min(Tlowest, T_min(i,j+1))
        Thighest = max(Thighest, T_max(i,j+1))
      endif
      if (kMLD(i,j)<nz) then
        Tlowest = min(Tlowest, tv%T(i,j,kMLD(i,j)+1))
        Thighest = max(Thighest, tv%T(i,j,kMLD(i,j)+1))
      endif
      ! Determine largest deltaT that will not create new extrema
      deltaT(i,j) = max(0., T_min(i,j) - Tlowest)
      deltaT(i,j) = min(deltaT(i,j), max(0., Thighest - T_max(i,j)))
      depthT = 0.
      heat(i,j) = 0.
      do k=1,kMLD(i,j)
        depthB = min(depthT + h(i,j,k), MLD(i,j))
        df = 1. - (depthT + depthB) / (MLD(i,j) + G%H_subroundoff)
        df = deltaT(i,j) * df * (depthB - depthT)/(h(i,j,k) + G%H_subroundoff)
        tv%T(i,j,k) = tv%T(i,j,k) + df
        depthT = depthB ! Cycle for next loop
        heat(i,j) = heat(i,j) + df * (depthB - depthT)
      enddo
      heat(i,j) = heat(i,j) / (MLD(i,j) + G%H_subroundoff)
    enddo
  enddo

! Offer fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_MLD > 0) call post_data(CS%id_MLD, MLD*G%H_to_m, CS%diag)
    if (CS%id_deltaT > 0) call post_data(CS%id_deltaT, deltaT, CS%diag)
    if (CS%id_heat > 0) call post_data(CS%id_heat, heat, CS%diag)
  endif

end subroutine counter_grad_restrat

logical function counter_grad_restrat_init(Time, G, param_file, diag, CS)
  type(time_type),             intent(in)    :: Time       !< Time structure
  type(ocean_grid_type),       intent(in)    :: G          !< Grid structure
  type(param_file_type),       intent(in)    :: param_file !< Parameter file handle
  type(diag_ctrl), target,     intent(inout) :: diag       !< Diagnostics structure
  type(counter_grad_restrat_CS), pointer     :: CS         !< Control structure for this module
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_counter_grad_restrat"  ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "counter_grad_restrat_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "COUNTER_GRAD_RESTRAT", counter_grad_restrat_init, &
             "If true, a counter-gradient re-stratifying flux is imposed\n"//&
             "in the mixed layer.", default=.false.)
  if (.not. counter_grad_restrat_init) return

  call get_param(param_file, mod, "FOX_KEMPER_ML_RESTRAT_COEF", CS%ml_restrat_coef, &
             "A nondimensional coefficient that is proportional to \n"//&
             "the ratio of the deformation radius to the dominant \n"//&
             "lengthscale of the submesoscale mixed layer \n"//&
             "instabilities, times the minimum of the ratio of the \n"//&
             "mesoscale eddy kinetic energy to the large-scale \n"//&
             "geostrophic kinetic energy or 1 plus the square of the \n"//&
             "grid spacing over the deformation radius, as detailed \n"//&
             "by Fox-Kemper et al. (2010)", units="nondim", default=0.0)
  if (G%nkml==0) then
    call get_param(param_file, mod, "MLE_DENSITY_DIFF", CS%MLE_density_diff, &
             "Density difference used to detect the mixed-layer\n"//&
             "depth used for the mixed-layer eddy parameterization\n"//&
             "by Fox-Kemper et al. (2010)", units="kg/m3", default=0.03)
  else
    ! Nonsense values to cause problems when these parameters are not used
    CS%MLE_density_diff = -9.e9
  endif
  CS%id_deltaT = register_diag_field('ocean_model', 'dt_counter_grad', diag%axesT1, Time, &
      'Mixed Layer Depth as used in the mixed-layer counter-gradient parameterization', 'meter')
  CS%id_heat = register_diag_field('ocean_model', 'heat_counter_grad', diag%axesT1, Time, &
      'Mixed Layer Depth as used in the mixed-layer counter-gradient parameterization', 'meter')
  CS%id_MLD = register_diag_field('ocean_model', 'MLD_counter_grad', diag%axesT1, Time, &
      'Mixed Layer Depth as used in the mixed-layer counter-gradient parameterization', 'meter')

end function counter_grad_restrat_init

end module MOM_counter_grad_restrat
