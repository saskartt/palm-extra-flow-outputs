!> @file user_module.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2021 Leibniz Universitaet Hannover
! Copyright 2024 University of Helsinki
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> User module for various additional outputs from the model. See README.md.
!> Author: Sasu Karttunen <sasu.karttunen@helsinki.fi>
!--------------------------------------------------------------------------------------------------!
 MODULE user

    USE arrays_3d

    USE control_parameters

    USE cpulog

    USE indices

    USE kinds

#if defined( __parallel )
    USE MPI
#endif

    USE pegrid

    USE statistics

    USE surface_mod

    IMPLICIT NONE

!-- Please change the revision number in case that the current revision does not match with previous
!-- revisions (e.g. if routines have been added/deleted or if parameter lists in subroutines have
!-- been changed).
    CHARACTER (LEN=12) ::  user_interface_current_revision = '1.1'

    INTEGER(iwp) ::  dots_num_palm      !<
    INTEGER(iwp) ::  dots_num_user = 0  !<
    INTEGER(iwp) ::  user_idummy        !<

    INTEGER(iwp) ::  total_atm_gp

    LOGICAL ::  user_module_enabled = .FALSE.  !<

    REAL(wp) ::  user_rdummy  !<

!
!-- Sample for user-defined output
!    REAL(wp) :: global_parameter  !< user defined global parameter
!

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pres_drag_norm_x        !< user defined array
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pres_drag_norm_x_av     !< user defined array
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pres_drag_norm_y        !< user defined array
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pres_drag_norm_y_av     !< user defined array

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  theta_product      !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  theta_product_av   !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  utheta_product     !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  utheta_product_av  !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vtheta_product     !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vtheta_product_av  !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  uq_product         !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  uq_product_av      !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vq_product         !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vq_product_av      !< user defined array

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  wu_sgs             !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  wu_sgs_av          !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  wv_sgs             !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  wv_sgs_av          !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  wtheta_sgs         !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  wtheta_sgs_av      !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  wthetav_sgs        !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  wthetav_sgs_av     !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  wq_sgs             !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  wq_sgs_av          !< user defined array

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  uv_sgs             !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  uv_sgs_av          !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  uw_sgs             !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  uw_sgs_av          !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  utheta_sgs         !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  utheta_sgs_av      !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  uq_sgs             !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  uq_sgs_av          !< user defined array

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vu_sgs             !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vu_sgs_av          !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vw_sgs             !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vw_sgs_av          !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vtheta_sgs         !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vtheta_sgs_av      !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vq_sgs             !< user defined array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vq_sgs_av          !< user defined array
!    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  ustvst  !< user defined array

    SAVE

    PRIVATE

!
!- Public functions
    PUBLIC                                                                                         &
       user_actions,                                                                               &
       user_boundary_conditions,                                                                   &
       user_check_data_output,                                                                     &
       user_check_data_output_pr,                                                                  &
       user_check_data_output_ts,                                                                  &
       user_check_parameters,                                                                      &
       user_data_output_2d,                                                                        &
       user_data_output_3d,                                                                        &
       user_define_netcdf_grid,                                                                    &
       user_exchange_horiz,                                                                        &
       user_header,                                                                                &
       user_init,                                                                                  &
       user_init_arrays,                                                                           &
       user_last_actions,                                                                          &
       user_parin,                                                                                 &
       user_prognostic_equations,                                                                  &
       user_rrd_global,                                                                            &
       user_rrd_local,                                                                             &
       user_statistics,                                                                            &
       user_swap_timelevel,                                                                        &
       user_3d_data_averaging,                                                                     &
       user_wrd_global,                                                                            &
       user_wrd_local

!
!- Public parameters, constants and initial values
   PUBLIC                                                                                          &
      user_module_enabled


    INTERFACE user_parin
       MODULE PROCEDURE user_parin
    END INTERFACE user_parin

    INTERFACE user_boundary_conditions
       MODULE PROCEDURE user_boundary_conditions
    END INTERFACE user_boundary_conditions

    INTERFACE user_check_parameters
       MODULE PROCEDURE user_check_parameters
    END INTERFACE user_check_parameters

    INTERFACE user_check_data_output_ts
       MODULE PROCEDURE user_check_data_output_ts
    END INTERFACE user_check_data_output_ts

    INTERFACE user_check_data_output_pr
       MODULE PROCEDURE user_check_data_output_pr
    END INTERFACE user_check_data_output_pr

    INTERFACE user_check_data_output
       MODULE PROCEDURE user_check_data_output
    END INTERFACE user_check_data_output

    INTERFACE user_define_netcdf_grid
       MODULE PROCEDURE user_define_netcdf_grid
    END INTERFACE user_define_netcdf_grid

    INTERFACE user_exchange_horiz
       MODULE PROCEDURE user_exchange_horiz
    END INTERFACE user_exchange_horiz

    INTERFACE user_init
       MODULE PROCEDURE user_init
    END INTERFACE user_init

    INTERFACE user_init_arrays
       MODULE PROCEDURE user_init_arrays
    END INTERFACE user_init_arrays

    INTERFACE user_header
       MODULE PROCEDURE user_header
    END INTERFACE user_header

    INTERFACE user_actions
       MODULE PROCEDURE user_actions
       MODULE PROCEDURE user_actions_ij
    END INTERFACE user_actions

    INTERFACE user_prognostic_equations
       MODULE PROCEDURE user_prognostic_equations
       MODULE PROCEDURE user_prognostic_equations_ij
    END INTERFACE user_prognostic_equations

    INTERFACE user_3d_data_averaging
       MODULE PROCEDURE user_3d_data_averaging
    END INTERFACE user_3d_data_averaging

    INTERFACE user_data_output_2d
       MODULE PROCEDURE user_data_output_2d
    END INTERFACE user_data_output_2d

    INTERFACE user_data_output_3d
       MODULE PROCEDURE user_data_output_3d
    END INTERFACE user_data_output_3d

    INTERFACE user_statistics
       MODULE PROCEDURE user_statistics
    END INTERFACE user_statistics

    INTERFACE user_swap_timelevel
       MODULE PROCEDURE user_swap_timelevel
    END INTERFACE user_swap_timelevel

    INTERFACE user_rrd_global
       MODULE PROCEDURE user_rrd_global_ftn
       MODULE PROCEDURE user_rrd_global_mpi
    END INTERFACE user_rrd_global

    INTERFACE user_rrd_local
       MODULE PROCEDURE user_rrd_local_ftn
       MODULE PROCEDURE user_rrd_local_mpi
    END INTERFACE user_rrd_local

    INTERFACE user_wrd_global
       MODULE PROCEDURE user_wrd_global
    END INTERFACE user_wrd_global

    INTERFACE user_wrd_local
       MODULE PROCEDURE user_wrd_local
    END INTERFACE user_wrd_local

    INTERFACE user_last_actions
       MODULE PROCEDURE user_last_actions
    END INTERFACE user_last_actions


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &user_parameters for user module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_parin

    CHARACTER (LEN=80) ::  line  !< string containing the last line read from namelist file

    INTEGER(iwp) ::  cmi        !< end position of major index in user_interface_current_revision
    INTEGER(iwp) ::  io_status  !< status after reading the namelist file
    INTEGER(iwp) ::  rmi        !< end position of major index in user_interface_required_revision

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /user_parameters/  data_output_masks_user,                                            &
                                data_output_pr_user,                                               &
                                data_output_user,                                                  &
                                region,                                                            &
                                switch_off_module

!
!-- Next statement is to avoid compiler warnings about unused variables. Please remove in case
!-- that you are using them.
    IF ( dots_num_palm == 0  .OR.  dots_num_user == 0  .OR.  user_idummy == 0  .OR.                &
         user_rdummy == 0.0_wp )  CONTINUE

!
!-- Check for the user's interface version. Only major numbers must match.
    cmi = INDEX( user_interface_current_revision, '.' ) - 1
    rmi = INDEX( user_interface_required_revision, '.' ) - 1
    IF ( user_interface_current_revision(1:cmi) /= user_interface_required_revision(1:rmi) )  THEN
       message_string = 'current user-interface revision "' //                                     &
                        TRIM( user_interface_current_revision ) // '" does ' //                    &
                        'not match the required revision "' //                                     &
                        TRIM( user_interface_required_revision ) // '"'
        CALL message( 'user_parin', 'USR0001', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Position the namelist-file at the beginning (it has already been opened in parin), and try to
!-- read (find) a namelist named "user_parameters".
    REWIND ( 11 )
    READ( 11, user_parameters, IOSTAT=io_status )

!
!-- Actions depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    User namelist found and correctly read. Set default module switch to true. This activates
!--    calls of the user-interface subroutines.
       IF ( .NOT. switch_off_module )  user_module_enabled = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    User namelist was found, but contained errors. Print an error message containing the line
!--    that caused the problem
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'user_parameters', line )

    ENDIF

 END SUBROUTINE user_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set boundary conditions for a user-defined prognostic quantity.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_boundary_conditions

!
!-- Example code for a user-defined prognostic quantity s1.
!    IMPLICIT NONE

!    INTEGER(iwp) ::  i                            !< grid index x direction.
!    INTEGER(iwp) ::  j                            !< grid index y direction.
!    INTEGER(iwp) ::  k                            !< grid index z direction.
!    INTEGER(iwp) ::  m                            !< running index surface elements.

!!
!!-- Boundary conditions for s1.
!!-- Bottom boundary: Neumann condition because s1 flux is always given.
!    !$OMP PARALLEL DO PRIVATE( i, j, k )
!    DO  m = 1, bc_hv%ns
!       i = bc_hv%i(m)
!       j = bc_hv%j(m)
!       k = bc_hv%k(m)
!       s1_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = s1_p(k,j,i)
!    ENDDO
!!
!!-- Top boundary: Dirichlet or Neumann
!    IF ( ibc_s1_t == 0 )  THEN
!        s1_p(nzt+1,:,:) = s1(nzt+1,:,:)
!    ELSEIF ( ibc_sa_t == 1 )  THEN
!        s1_p(nzt+1,:,:) = s1_p(nzt,:,:)
!    ENDIF

 END SUBROUTINE user_boundary_conditions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check &user_parameters control parameters and deduce further quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_check_parameters

!
!-- Here the user may add code to check the validity of further &userpar control parameters or
!-- deduce further quantities.


 END SUBROUTINE user_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set module-specific timeseries units and labels
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )

    INTEGER(iwp),      INTENT(IN)     ::  dots_max  !<
    INTEGER(iwp),      INTENT(INOUT)  ::  dots_num  !<

    CHARACTER(LEN=*), DIMENSION(dots_max), INTENT(INOUT)  ::  dots_label  !<
    CHARACTER(LEN=*), DIMENSION(dots_max), INTENT(INOUT)  ::  dots_unit   !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( dots_num == 0  .OR.  dots_label(1)(1:1) == ' '  .OR.  dots_unit(1)(1:1) == ' ' )  CONTINUE

!
!-- Sample for user-defined time series:
!-- For each time series quantity you have to give a label and a unit, which will be used for the
!-- NetCDF file. They must not contain more than seven characters. The value of dots_num has to be
!-- increased by the number of new time series quantities. Its old value has to be stored in
!-- dots_num_palm. See routine user_statistics on how to calculate and output these quantities.

!    dots_num_palm = dots_num

!    dots_num = dots_num + 1
!    dots_num_user = dots_num_user + 1
!    dots_label(dots_num) = 'abs_umx'
!    dots_unit(dots_num)  = 'm/s'

!    dots_num = dots_num + 1
!    dots_num_user = dots_num_user + 1
!    dots_label(dots_num) = 'abs_vmx'
!    dots_unit(dots_num)  = 'm/s'


 END SUBROUTINE user_check_data_output_ts


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of user defined profile output quantities. For those variables not recognized by the
!> user, the parameter unit is set to "illegal", which tells the calling routine that the
!> output variable is not defined and leads to a program abort.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_check_data_output_pr( variable, var_count, unit, dopr_unit )


    USE profil_parameter


    CHARACTER (LEN=*) ::  unit      !<
    CHARACTER (LEN=*) ::  variable  !<
    CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit

!    INTEGER(iwp) ::  user_pr_index  !<
    INTEGER(iwp) ::  var_count      !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( unit(1:1) == ' '  .OR.  dopr_unit(1:1) == ' '  .OR.  var_count == 0 )  CONTINUE

    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary.
!--    Add additional CASE statements depending on the number of quantities for which profiles are
!--    to be calculated. The respective calculations to be performed have to be added in routine
!--    user_statistics. The quantities are (internally) identified by a user-profile-number
!--    (see variable "user_pr_index" below). The first user-profile must be assigned the number
!--    "pr_palm+max_pr_cs+max_pr_det+max_pr_salsa+1", the second one "pr_palm+max_pr_cs+max_pr_det+
!--    max_pr_salsa+2", etc. The respective user-profile-numbers have also to be used in routine 
!--    user_statistics!
!       CASE ( 'u*v*' )                      ! quantity string as given in data_output_pr_user
!          user_pr_index = pr_palm + max_pr_cs + max_pr_det + max_pr_salsa + 1
!          dopr_index(var_count)  = user_pr_index    ! quantities' user-profile-number
!          dopr_unit = 'm2/s2'  ! quantity unit
!          unit = dopr_unit
!          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 )
!                                            ! grid on which the quantity is defined (use zu or zw)
!

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE user_check_data_output_pr


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of user defined output quantities. For those variables not recognized by the user,
!> the parameter unit is set to "illegal", which tells the calling routine that the output variable
!> is not defined and leads to a program abort.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_check_data_output( variable, unit )


    CHARACTER (LEN=*) ::  unit      !<
    CHARACTER (LEN=*) ::  variable  !<


    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary
       CASE ( 'pres_drag_norm_x*', 'pres_drag_norm_y*' )
          unit = 'kgm/s2'

       CASE ( 'theta_product' )
          unit = 'K2'

       CASE ( 'wu_sgs', 'wv_sgs', 'uv_sgs', 'uw_sgs', 'vu_sgs', 'vw_sgs' )
          unit = 'm2/s2'

       CASE ( 'utheta_product', 'vtheta_product' )
         unit = 'Km/s'

       CASE ( 'uq_product', 'vq_product' )
         unit = 'm/s'

       CASE ( 'wtheta_sgs', 'utheta_sgs', 'vtheta_sgs' )
          unit = 'W/m2'

       CASE ( 'wthetav_sgs' )
          unit = 'W/m2'

       CASE ( 'wq_sgs', 'uq_sgs', 'vq_sgs' )
          unit = 'W/m2'

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE user_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange ghost points for a user-defined prognostic quantity.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_exchange_horiz( location )

    CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string

    SELECT CASE ( location )

       CASE ( 'before_prognostic_equation' )

       CASE ( 'after_prognostic_equation' )

!
!--        Example for an additional passive scalar quantity s1.
!          CALL exchange_horiz( s1_p, nbgp )

    END SELECT

 END SUBROUTINE user_exchange_horiz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize user-defined arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_init_arrays


!    INTEGER(iwp) :: i       !< loop index
!    INTEGER(iwp) :: j       !< loop index
!    INTEGER(iwp) :: region  !< index for loop over statistic regions

!
!-- Allocate user-defined arrays and set flags for statistic regions.
!-- Sample for user-defined output
    ALLOCATE( pres_drag_norm_x(nysg:nyng,nxlg:nxrg) )
    ALLOCATE( pres_drag_norm_y(nysg:nyng,nxlg:nxrg) )

    ALLOCATE( theta_product(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( utheta_product(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( vtheta_product(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( uq_product(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( vq_product(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    ALLOCATE( wu_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( wv_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( wtheta_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( wthetav_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( wq_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    ALLOCATE( uv_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( uw_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( utheta_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( uq_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    ALLOCATE( vu_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( vw_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( vtheta_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( vq_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!   ALLOCATE( ustvst(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

!
!-- Example for defining a statistic region:
!-- ATTENTION: rmask = 0 is required at the ghost boundaries to guarantee correct statistic
!--            evaluations (otherwise ghost points would be counted twice). This setting has
!--            already been cared for in routine init_3d_model. Please don't set the ghost points
!--            /= 0. i.e. run the following loop only over nxl,nxr and nys,nyn.
!     IF ( statistic_regions >= 1 )  THEN
!        region = 1
!
!        rmask(:,:,region) = 0.0_wp
!        DO  i = nxl, nxr
!           IF ( i >= INT( 0.25 * nx ) .AND. i <= INT( 0.75 * nx ) )  THEN
!              DO  j = nys, nyn
!                 IF ( i >= INT( 0.25 * ny ) .AND. i <= INT( 0.75 * ny ) )  THEN
!                    rmask(j,i,region) = 1.0_wp
!                 ENDIF
!              ENDDO
!           ENDIF
!        ENDDO
!
!     ENDIF

 END SUBROUTINE user_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execution of user-defined initializing actions
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_init

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  total_atm_gp_l

!    CHARACTER(LEN=20) :: field_char  !<
!
!-- Here the user-defined initializing actions follow:
!-- Sample for user-defined output
!    ustvst = 0.0_wp

!
!-- Count total atmospheric grid points.
    total_atm_gp = 0
    total_atm_gp_l = 0
    DO  i = nxl, nxr
       DO  j =  nys, nyn
          DO  k = nzb, nzt+1
             total_atm_gp_l = total_atm_gp_l + MERGE( 1, 0, BTEST( topo_flags(k,j,i), 22 ) )
          ENDDO
       ENDDO
    ENDDO

    CALL MPI_ALLREDUCE( total_atm_gp_l, total_atm_gp, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )


 END SUBROUTINE user_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the grids on which user-defined output quantities are defined. Allowed values for grid_x are
!> "x" and "xu", for grid_y "y" and "yv", and for grid_z "zu" and "zw".
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )


    CHARACTER (LEN=*) ::  grid_x     !<
    CHARACTER (LEN=*) ::  grid_y     !<
    CHARACTER (LEN=*) ::  grid_z     !<
    CHARACTER (LEN=*) ::  variable   !<

    LOGICAL ::  found   !<


    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary
       CASE ( 'theta_product', 'theta_product_xy', 'theta_product_xz', 'theta_product_yz' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'utheta_product', 'utheta_product_xy', 'utheta_product_xz', 'utheta_product_yz' )
          found = .TRUE.
          grid_x = 'xu'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'vtheta_product', 'vtheta_product_xy', 'vtheta_product_xz', 'vtheta_product_yz' )
          found = .TRUE.
          grid_x = 'x'
          grid_y = 'yv'
          grid_z = 'zu'

       CASE ( 'uq_product', 'uq_product_xy', 'uq_product_xz', 'uq_product_yz' )
          found = .TRUE.
          grid_x = 'xu'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'vq_product', 'vq_product_xy', 'vq_product_xz', 'vq_product_yz' )
          found = .TRUE.
          grid_x = 'x'
          grid_y = 'yv'
          grid_z = 'zu'

       CASE ( 'wu_sgs', 'wu_sgs_xy', 'wu_sgs_xz', 'wu_sgs_yz' )
          found  = .TRUE.
          grid_x = 'xu'
          grid_y = 'y'
          grid_z = 'zw'

       CASE ( 'wv_sgs', 'wv_sgs_xy', 'wv_sgs_xz', 'wv_sgs_yz' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'yv'
          grid_z = 'zw'

       CASE ( 'wtheta_sgs', 'wtheta_sgs_xy', 'wtheta_sgs_xz', 'wtheta_sgs_yz' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'

       CASE ( 'wthetav_sgs', 'wthetav_sgs_xy', 'wthetav_sgs_xz', 'wthetav_sgs_yz' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'

       CASE ( 'wq_sgs', 'wq_sgs_xy', 'wq_sgs_xz', 'wq_sgs_yz' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'

       CASE ( 'uv_sgs', 'uv_sgs_xy', 'uv_sgs_xz', 'uv_sgs_yz' )
          found  = .TRUE.
          grid_x = 'xu'
          grid_y = 'yv'
          grid_z = 'zu'

       CASE ( 'uw_sgs', 'uw_sgs_xy', 'uw_sgs_xz', 'uw_sgs_yz' )
          found  = .TRUE.
          grid_x = 'xu'
          grid_y = 'y'
          grid_z = 'zw'

       CASE ( 'utheta_sgs', 'utheta_sgs_xy', 'utheta_sgs_xz', 'utheta_sgs_yz' )
          found  = .TRUE.
          grid_x = 'xu'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'uq_sgs', 'uq_sgs_xy', 'uq_sgs_xz', 'uq_sgs_yz' )
          found  = .TRUE.
          grid_x = 'xu'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'vu_sgs', 'vu_sgs_xy', 'vu_sgs_xz', 'vu_sgs_yz' )
          found  = .TRUE.
          grid_x = 'xu'
          grid_y = 'yv'
          grid_z = 'zu'

       CASE ( 'vw_sgs', 'vw_sgs_xy', 'vw_sgs_xz', 'vw_sgs_yz' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'yv'
          grid_z = 'zw'

       CASE ( 'vtheta_sgs', 'vtheta_sgs_xy', 'vtheta_sgs_xz', 'vtheta_sgs_yz' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'yv'
          grid_z = 'zu'

       CASE ( 'vq_sgs', 'vq_sgs_xy', 'vq_sgs_xz', 'vq_sgs_yz' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'yv'
          grid_z = 'zu'

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT


 END SUBROUTINE user_define_netcdf_grid




!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Print a header with user-defined information.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_header( io )


    INTEGER(iwp) ::  i   !<
    INTEGER(iwp) ::  io  !<
!
!-- If no user-defined variables are read from the namelist-file, no information will be printed.
    IF ( .NOT. user_module_enabled )  THEN
       WRITE ( io, 100 )
       RETURN
    ENDIF

!
!-- Printing the information.
    WRITE ( io, 110 )

    IF ( statistic_regions /= 0 )  THEN
       WRITE ( io, 200 )
       DO  i = 0, statistic_regions
          WRITE ( io, 201 )  i, region(i)
       ENDDO
    ENDIF

!
!-- Format-descriptors
100 FORMAT (//' *** no user-defined variables found'/)
110 FORMAT (//1X,78('#') // ' User-defined variables and actions:' /                               &
            ' -----------------------------------'//)
200 FORMAT (' Output of profiles and time series for following regions:' /)
201 FORMAT (4X,'Region ',I1,':   ',A)


 END SUBROUTINE user_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_actions( location )

    USE grid_variables,                                                                            &
        ONLY:   dx,                                                                                &
                dy,                                                                                &
                ddx,                                                                               &
                ddy

    CHARACTER(LEN=*) ::  location  !<

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  m  !<

    REAL(wp) ::  flag

    REAL(wp) ::  total_pres
    REAL(wp) ::  total_pres_l

    CALL cpu_log( log_point(24), 'user_actions', 'start' )

!
!-- Here the user-defined actions follow. No calls for single grid points are allowed at locations
!-- before and after the timestep, since these calls are not within an i,j-loop
    SELECT CASE ( location )

       CASE ( 'before_timestep' )
!
!--       Enter actions to be done before every timestep here

       CASE ( 'before_prognostic_equations' )
!
!--       Enter actions to be done before all prognostic equations here

       CASE ( 'after_integration' )
!
!--       Enter actions to be done after every time integration (before data output)

!
!--       Compute average perturbation pressure, which is not necessarily zero for nest domains due
!--       due to all-neumann boundary conditions.
          total_pres_l = 0.0_wp
          DO  i = nxl, nxr
             DO  j =  nys, nyn
                DO  k = nzb, nzt+1
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 22 ) )
                   total_pres_l = total_pres_l + p(k,j,i) * flag
                ENDDO
             ENDDO
          ENDDO

          CALL MPI_ALLREDUCE( total_pres_l, total_pres, 1, MPI_REAL, MPI_SUM, comm2d, ierr )

!
!--       Calculate the mean atmospheric perturbation pressure.
          total_pres = total_pres / total_atm_gp

!
!--       Compute diagnostic user outputs.

!
!--       Pressure drag on right (east) and left (west) faces from normalized perturbation pressure.
          pres_drag_norm_x = 0.0_wp
          DO  m = 1, surf_def%ns
             i = surf_def%i(m)
             j = surf_def%j(m)
             k = surf_def%k(m)
             pres_drag_norm_x(j,i) = pres_drag_norm_x(j,i)                                         &
                                + MERGE( - ( p(k,j,i) - total_pres ) * dy * dz(1),                 &
                                         0.0_wp, surf_def%eastward(m) )
             pres_drag_norm_x(j,i) = pres_drag_norm_x(j,i)                                         &
                                + MERGE( ( p(k,j,i) - total_pres ) * dy * dz(1),                   &
                                         0.0_wp, surf_def%westward(m) )
          ENDDO
          DO  m = 1, surf_lsm%ns
             i = surf_lsm%i(m)
             j = surf_lsm%j(m)
             k = surf_lsm%k(m)
             pres_drag_norm_x(j,i) = pres_drag_norm_x(j,i)                                         &
                                      + MERGE( - ( p(k,j,i) - total_pres ) * dy * dz(1),           &
                                               0.0_wp, surf_lsm%eastward(m) )
             pres_drag_norm_x(j,i) = pres_drag_norm_x(j,i)                                         &
                                + MERGE( ( p(k,j,i) - total_pres ) * dy * dz(1),                   &
                                         0.0_wp, surf_lsm%westward(m) )
          ENDDO
          DO  m = 1, surf_usm%ns
             i = surf_usm%i(m)
             j = surf_usm%j(m)
             k = surf_usm%k(m)
             pres_drag_norm_x(j,i) = pres_drag_norm_x(j,i)                                         &
                                + MERGE( - ( p(k,j,i) - total_pres ) * dy * dz(1),                 &
                                         0.0_wp, surf_usm%eastward(m) )
             pres_drag_norm_x(j,i) = pres_drag_norm_x(j,i)                                         &
                                + MERGE( ( p(k,j,i) - total_pres ) * dy * dz(1),                   &
                                         0.0_wp, surf_usm%westward(m) )
          ENDDO

!
!--       Pressure drag on north and south faces from normalized perturbation pressure.
          pres_drag_norm_y = 0.0_wp
          DO  m = 1, surf_def%ns
             i = surf_def%i(m)
             j = surf_def%j(m)
             k = surf_def%k(m)
             pres_drag_norm_y(j,i) = pres_drag_norm_y(j,i)                                         &
                                + MERGE( - ( p(k,j,i) - total_pres ) * dx * dz(1),                 &
                                         0.0_wp, surf_def%northward(m) )
             pres_drag_norm_y(j,i) = pres_drag_norm_y(j,i)                                         &
                                + MERGE( ( p(k,j,i) - total_pres ) * dx * dz(1),                   &
                                         0.0_wp, surf_def%southward(m) )
          ENDDO
          DO  m = 1, surf_lsm%ns
             i = surf_lsm%i(m)
             j = surf_lsm%j(m)
             k = surf_lsm%k(m)
             pres_drag_norm_y(j,i) = pres_drag_norm_y(j,i)                                         &
                                      + MERGE( - ( p(k,j,i) - total_pres ) * dx * dz(1),           &
                                               0.0_wp, surf_lsm%northward(m) )
             pres_drag_norm_y(j,i) = pres_drag_norm_y(j,i)                                         &
                                + MERGE( ( p(k,j,i) - total_pres ) * dx * dz(1),                   &
                                         0.0_wp, surf_lsm%southward(m) )
          ENDDO
          DO  m = 1, surf_usm%ns
             i = surf_usm%i(m)
             j = surf_usm%j(m)
             k = surf_usm%k(m)
             pres_drag_norm_y(j,i) = pres_drag_norm_y(j,i)                                         &
                                + MERGE( - ( p(k,j,i) - total_pres ) * dx * dz(1),                 &
                                         0.0_wp, surf_usm%northward(m) )
             pres_drag_norm_y(j,i) = pres_drag_norm_y(j,i)                                         &
                                + MERGE( ( p(k,j,i) - total_pres ) * dx * dz(1),                   &
                                         0.0_wp, surf_usm%southward(m) )
          ENDDO



          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb, nzt
                   theta_product(k,j,i) = pt(k,j,i)**2                                             &
                                       * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                   utheta_product(k,j,i) = u(k,j,i) * 0.5_wp * ( pt(k,j,i) + pt(k,j,i+1) )         &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                   vtheta_product(k,j,i) = v(k,j,i) * 0.5_wp * ( pt(k,j,i) + pt(k,j+1,i) )         &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )

                   uq_product(k,j,i) = u(k,j,i) * 0.5_wp * ( q(k,j,i) + q(k,j,i+1) )               &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                   vq_product(k,j,i) = v(k,j,i) * 0.5_wp * ( q(k,j,i) + q(k,j+1,i) )               &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )

                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 23 ) )                  &
                          * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 9  ) )
                                     
                   wu_sgs(k,j,i) = - 0.25_wp * (                                                   &
                                        km(k,j,i) + km(k+1,j,i) + km(k,j,i-1) + km(k+1,j,i-1) )    &
                                            *( ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)             &
                                              + ( w(k,j,i)   - w(k,j,i-1) ) * ddx                  &
                                         ) * rho_air_zw(k)                                         &
                                           * momentumflux_output_conversion(k)                     &
                                           * flag
                                          

                   wv_sgs(k,j,i) = - 0.25_wp * (                                                   &
                                        km(k,j,i) + km(k+1,j,i) + km(k,j-1,i) + km(k+1,j-1,i) )    &
                                            * ( ( v(k+1,j,i) - v(k,j,i)   ) * ddzu(k+1)            &
                                              + ( w(k,j,i)   - w(k,j-1,i) ) * ddy                  &
                                         ) * rho_air_zw(k)                                         &
                                           * momentumflux_output_conversion(k)                     &
                                           * flag

                   wtheta_sgs(k,j,i) = - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )                      &
                                                * ( pt(k+1,j,i) - pt(k,j,i) )                      &
                                                * rho_air_zw(k)                                    &
                                                * heatflux_output_conversion(k)                    &
                                                * ddzu(k+1)                                        &
                                                * flag

                   wthetav_sgs(k,j,i) = - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )                     &
                                                 * ( vpt(k+1,j,i) - vpt(k,j,i) )                   &
                                                 * rho_air_zw(k)                                   &
                                                 * heatflux_output_conversion(k)                   &
                                                 * ddzu(k+1)                                       &
                                                 * flag
                   
                   wq_sgs(k,j,i) = - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )                          &
                                            * ( q(k+1,j,i) - q(k,j,i) )                            &
                                            * rho_air_zw(k)                                        &
                                            * heatflux_output_conversion(k)                        &
                                            * ddzu(k+1)                                            &
                                            * flag


                   uv_sgs(k,j,i) = - 0.25_wp * (                                                   &
                                        km(k,j,i) + km(k,j+1,i) + km(k,j,i-1) + km(k,j+1,i-1) )    &
                                            *( ( u(k,j+1,i) - u(k,j,i)   ) * ddx                   &
                                              + ( v(k,j,i)   - v(k,j,i-1) ) * ddy                  &
                                         ) * rho_air_zw(k)                                         &
                                           * momentumflux_output_conversion(k)                     &
                                           * flag
                                          

                   uw_sgs(k,j,i) = - 0.25_wp * (                                                   &
                                        km(k,j,i) + km(k+1,j,i) + km(k,j,i-1) + km(k+1,j,i-1) )    &
                                            * ( ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)            &
                                              + ( w(k,j,i)   - w(k,j,i-1) ) * ddx )                &
                                           * rho_air_zw(k)                                         &
                                           * momentumflux_output_conversion(k)                     &
                                           * flag

                   utheta_sgs(k,j,i) = - 0.5_wp * ( kh(k,j,i) + kh(k,j,i+1) )                      &
                                                * ( pt(k,j,i+1) - pt(k,j,i) )                      &
                                                * rho_air_zw(k)                                    &
                                                * heatflux_output_conversion(k)                    &
                                                * ddx                                              &
                                                * flag
                   
                   wq_sgs(k,j,i) = - 0.5_wp * ( kh(k,j,i) + kh(k,j,i+1) )                          &
                                            * ( q(k,j,i+1) - q(k,j,i) )                            &
                                            * rho_air_zw(k)                                        &
                                            * heatflux_output_conversion(k)                        &
                                            * ddx                                                  &
                                            * flag

                   vu_sgs(k,j,i) = - 0.25_wp * (                                                   &
                                        km(k,j,i) + km(k,j+1,i) + km(k,j,i-1) + km(k,j+1,i-1) )    &
                                            *( ( v(k,j,i+1) - v(k,j,i)   ) * ddy                   &
                                              + ( u(k,j,i)   - v(k,j-1,i) ) * ddy )                &
                                           * rho_air_zw(k)                                         &
                                           * momentumflux_output_conversion(k)                     &
                                           * flag
                                          

                   vw_sgs(k,j,i) = - 0.25_wp * (                                                   &
                                        km(k,j,i) + km(k+1,j,i) + km(k,j-1,i) + km(k+1,j-1,i) )    &
                                            * ( ( v(k+1,j,i) - v(k,j,i)   ) * ddzu(k+1)            &
                                              + ( w(k,j,i)   - w(k,j-1,i) ) * ddy )                &
                                           * rho_air_zw(k)                                         &
                                           * momentumflux_output_conversion(k)                     &
                                           * flag

                   vtheta_sgs(k,j,i) = - 0.5_wp * ( kh(k,j,i) + kh(k,j+1,i) )                      &
                                                * ( pt(k,j+1,i) - pt(k,j,i) )                      &
                                                * rho_air_zw(k)                                    &
                                                * heatflux_output_conversion(k)                    &
                                                * ddy                                              &
                                                * flag
                   
                   vq_sgs(k,j,i) = - 0.5_wp * ( kh(k,j,i) + kh(k,j+1,i) )                          &
                                            * ( q(k,j+1,i) - q(k,j,i) )                            &
                                            * rho_air_zw(k)                                        &
                                            * heatflux_output_conversion(k)                        &
                                            * ddy                                                  &
                                            * flag
                ENDDO ! k
!
!--              Add contribution from surfaces
                DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                   IF ( surf_lsm%upward(m)  .OR.  surf_lsm%downward(m) )  THEN
                      k = surf_lsm%k(m) + MERGE( surf_lsm%koff(m), 0, surf_lsm%upward(m) )
      
                      wu_sgs(k,j,i) = wu_sgs(k,j,i) +                                              &
                                     momentumflux_output_conversion(k) *                           &
                                     surf_lsm%usws(m)
                      wv_sgs(k,j,i) = wv_sgs(k,j,i) +                                              &
                                     momentumflux_output_conversion(k) *                           &
                                     surf_lsm%vsws(m)
                      wtheta_sgs(k,j,i) = wtheta_sgs(k,j,i) +                                      &
                                        heatflux_output_conversion(k) *                            &
                                        surf_lsm%shf_agg(m)
      
                      wthetav_sgs(k,j,i) = wthetav_sgs(k,j,i) +  (                                 &
                                           ( 1.0_wp + 0.61_wp * q(k,j,i) ) *                       &
                                           surf_lsm%shf_agg(m) + 0.61_wp * pt(k,j,i) *             &
                                           surf_lsm%qsws_agg(m) ) * heatflux_output_conversion(k)
      
                      wq_sgs(k,j,i) = wq_sgs(k,j,i) +                                              &
                                      waterflux_output_conversion(k) *                             &
                                      surf_lsm%qsws_agg(m)
                   ENDIF
                ENDDO ! m
             ENDDO ! j
          ENDDO ! i


!          DO  i = nxlg, nxr
!             DO  j = nysg, nyn
!                DO  k = nzb, nzt+1
!                   ustvst(k,j,i) =  &
!                      ( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) - hom(k,1,1,0) ) *                      &
!                      ( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) - hom(k,1,2,0) )
!                ENDDO
!             ENDDO
!          ENDDO


       CASE ( 'after_timestep' )
!
!--       Enter actions to be done after every timestep here


       CASE ( 'u-tendency' )
!
!--       Enter actions to be done in the u-tendency term here


       CASE ( 'v-tendency' )


       CASE ( 'w-tendency' )


       CASE ( 'pt-tendency' )


       CASE ( 'sa-tendency' )


       CASE ( 'e-tendency' )


       CASE ( 'q-tendency' )


       CASE ( 's-tendency' )


       CASE DEFAULT
          CONTINUE

    END SELECT

    CALL cpu_log( log_point(24), 'user_actions', 'stop' )

 END SUBROUTINE user_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_actions_ij( i, j, location )


    CHARACTER(LEN=*) ::  location  !<

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<

!
!-- Here the user-defined actions follow
    SELECT CASE ( location )

       CASE ( 'u-tendency' )

!
!--       Next line is to avoid compiler warning about unused variables. Please remove.
          IF ( i == 0  .OR.  j == 0 )  CONTINUE

!
!--       Enter actions to be done in the u-tendency term here


       CASE ( 'v-tendency' )


       CASE ( 'w-tendency' )


       CASE ( 'pt-tendency' )


       CASE ( 'sa-tendency' )


       CASE ( 'e-tendency' )


       CASE ( 'q-tendency' )


       CASE ( 's-tendency' )


       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE user_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execution of user-defined prognostic equations.
!> Vector optimized version (call for all grid points).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_prognostic_equations
!
!-- What follows is an example about how to add an equation for an additional passive scalar s1.
!-- Keep in mind, that this scalar has to be introduced in several other user-interface routines,
!-- too, e.g. in user_init or the routines reading/writing restart data.
!    USE advec_s_bc_mod,                                                                            &
!        ONLY:  advec_s_bc
!
!    USE advec_s_pw_mod,                                                                            &
!        ONLY:  advec_s_pw
!
!    USE advec_s_up_mod,                                                                            &
!        ONLY:  advec_s_up
!
!    USE advec_ws,                                                                                  &
!        ONLY:  advec_s_ws
!
!    USE arrays_3d,                                                                                 &
!        ONLY:  rdf_sc,                                                                             &
!               tend,                                                                               &
!
!    USE control_parameters,                                                                        &
!        ONLY:  bc_dirichlet_l,                                                                     &
!               bc_dirichlet_n,                                                                     &
!               bc_dirichlet_r,                                                                     &
!               bc_dirichlet_s,                                                                     &
!               bc_radiation_l,                                                                     &
!               bc_radiation_n,                                                                     &
!               bc_radiation_r,                                                                     &
!               bc_radiation_s,                                                                     &
!               dt_3d,                                                                              &
!               intermediate_timestep_count,                                                        &
!               intermediate_timestep_count_max,                                                    &
!               scalar_advec,                                                                       &
!               simulated_time,                                                                     &
!               timestep_scheme,                                                                    &
!               tsc,                                                                                &
!               ws_scheme_sca
!
!    USE diffusion_s_mod,                                                                           &
!        ONLY:  diffusion_s
!
!    IMPLICIT NONE
!
!    INTEGER(iwp) ::  i       !< loop index
!    INTEGER(iwp) ::  j       !< loop index
!    INTEGER(iwp) ::  k       !< loop index
!
!    REAL(wp)     ::  sbt     !< weighting factor for sub-time step
!
!
!!
!!-- Compute a prognostic equations for an additional scalar s1.
!
!!
!!-- tendency terms with communication
!    sbt = tsc(2)
!    IF ( scalar_advec == 'bc-scheme' )  THEN
!
!       IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!!
!!--       Bott-Chlond scheme always uses Euler time step. Thus:
!          sbt = 1.0_wp
!       ENDIF
!       tend = 0.0_wp
!       CALL advec_s_bc( s1, 's1' )
!
!    ENDIF
!
!!
!!-- tendency terms with no communication
!    IF ( scalar_advec /= 'bc-scheme' )  THEN
!       tend = 0.0_wp
!       IF ( timestep_scheme(1:5) == 'runge' )  THEN
!          IF ( ws_scheme_sca )  THEN
!             CALL advec_s_ws( advc_flags_s, s1, 's1',                                           &
!                              bc_dirichlet_l  .OR.  bc_radiation_l,                             &
!                              bc_dirichlet_n  .OR.  bc_radiation_n,                             &
!                              bc_dirichlet_r  .OR.  bc_radiation_r,                             &
!                              bc_dirichlet_s  .OR.  bc_radiation_s )
!          ELSE
!             CALL advec_s_pw( s1 )
!          ENDIF
!       ELSE
!          CALL advec_s_up( s1 )
!       ENDIF
!    ENDIF
!
!    CALL diffusion_s( s1, surf_top%s1sws, surf_def%s1sws, surf_lsm%s1sws, surf_usm%s1sws )
!
!!
!!-- Prognostic equation for s1
!    DO  i = nxl, nxr
!       DO  j = nys, nyn
!          !following directive is required to vectorize on Intel19
!          !DIR$ IVDEP
!          DO  k = nzb+1, nzt
!             s1_p(k,j,i) = s1(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) + tsc(3) * ts1_m(k,j,i) )  &
!                                         - tsc(5) * rdf_sc(k) * ( s1(k,j,i) - s1_init(k) )      &
!                                       )                                                        &
!                                         * MERGE( 1.0_wp, 0.0_wp,                               &
!                                                  BTEST( topo_flags(k,j,i), 0 )                 &
!                                                )
!             IF ( s1_p(k,j,i) < 0.0_wp )  s1_p(k,j,i) = 0.1_wp * s1(k,j,i)
!          ENDDO
!       ENDDO
!    ENDDO
!
!!
!!-- Calculate tendencies for the next Runge-Kutta step
!    IF ( timestep_scheme(1:5) == 'runge' )  THEN
!       IF ( intermediate_timestep_count == 1 )  THEN
!          DO  i = nxl, nxr
!             DO  j = nys, nyn
!                DO  k = nzb+1, nzt
!                  ts1_m(k,j,i) = tend(k,j,i)
!                ENDDO
!             ENDDO
!          ENDDO
!       ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
!          DO  i = nxl, nxr
!             DO  j = nys, nyn
!                DO  k = nzb+1, nzt
!                   ts1_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * ts1_m(k,j,i)
!                ENDDO
!             ENDDO
!          ENDDO
!       ENDIF
!    ENDIF

 END SUBROUTINE user_prognostic_equations


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execution of user-defined prognostic equations.
!> Cache optimized version (call for all k along grid point (j,i)).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_prognostic_equations_ij( i, j, i_omp_start, tn )
!
!-- Use code based on the vector optimized routine above.
!-- See e.g. routine ocean_prognostic_equations_ij on how to setup for a specific i,j index.
    INTEGER(iwp) ::  i             !< loop index x direction
    INTEGER(iwp) ::  i_omp_start   !< first loop index of i-loop in calling routine prognostic_equations
    INTEGER(iwp) ::  j             !< loop index y direction
!    INTEGER(iwp) ::  k             !< loop index z direction
    INTEGER(iwp) ::  tn            !< task number of openmp task


!
!-- Dummy code to avoid compiler warnings. Please remove
    IF ( i == -100 )  tn = i_omp_start + j + tn

 END SUBROUTINE user_prognostic_equations_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels for a user-defined prognostic quantity.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_swap_timelevel( mod_count )

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  mod_count  !< flag defining where pointers point to


!
!-- Example code for a user-defined prognostic quantity s1.
    SELECT CASE ( mod_count )

       CASE ( 0 )
!          s1 => s1_1;    s1_p => s1_2

       CASE ( 1 )
!          s1 => s1_2;    s1_p => s1_1

    END SELECT


 END SUBROUTINE user_swap_timelevel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average user-defined output quantities as well as allocate the array necessary
!> for storing the average.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_3d_data_averaging( mode, variable )


    CHARACTER(LEN=*) ::  mode      !<
    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

!
!--       Uncomment and extend the following lines, if necessary.
!--       The arrays for storing the user defined quantities (here u2_av) have to be declared and
!--       defined by the user!
!--       Sample for user-defined output:
          CASE ( 'pres_drag_norm_x*' )
             IF ( .NOT. ALLOCATED( pres_drag_norm_x_av ) )  THEN
                ALLOCATE( pres_drag_norm_x_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             pres_drag_norm_x_av = 0.0_wp

          CASE ( 'pres_drag_norm_y*' )
             IF ( .NOT. ALLOCATED( pres_drag_norm_y_av ) )  THEN
                ALLOCATE( pres_drag_norm_y_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             pres_drag_norm_y_av = 0.0_wp

          CASE ( 'theta_product' )
             IF ( .NOT. ALLOCATED( theta_product_av ) )  THEN
                ALLOCATE( theta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             theta_product_av = 0.0_wp

          CASE ( 'utheta_product' )
             IF ( .NOT. ALLOCATED( utheta_product_av ) )  THEN
                ALLOCATE( utheta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             utheta_product_av = 0.0_wp

          CASE ( 'vtheta_product' )
             IF ( .NOT. ALLOCATED( vtheta_product_av ) )  THEN
                ALLOCATE( vtheta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             vtheta_product_av = 0.0_wp

          CASE ( 'uq_product' )
             IF ( .NOT. ALLOCATED( uq_product_av ) )  THEN
                ALLOCATE( uq_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             uq_product_av = 0.0_wp

          CASE ( 'vq_product' )
             IF ( .NOT. ALLOCATED( vq_product_av ) )  THEN
                ALLOCATE( vq_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             vq_product_av = 0.0_wp

          CASE ( 'wu_sgs' )
             IF ( .NOT. ALLOCATED( wu_sgs_av ) )  THEN
                ALLOCATE( wu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             wu_sgs_av = 0.0_wp

          CASE ( 'wv_sgs' )
             IF ( .NOT. ALLOCATED( wv_sgs_av ) )  THEN
                ALLOCATE( wv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             wv_sgs_av = 0.0_wp

          CASE ( 'wtheta_sgs' )
             IF ( .NOT. ALLOCATED( wtheta_sgs_av ) )  THEN
                ALLOCATE( wtheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             wtheta_sgs_av = 0.0_wp

          CASE ( 'wthetav_sgs' )
             IF ( .NOT. ALLOCATED( wthetav_sgs_av ) )  THEN
                ALLOCATE( wthetav_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             wthetav_sgs_av = 0.0_wp

          CASE ( 'wq_sgs' )
             IF ( .NOT. ALLOCATED( wq_sgs_av ) )  THEN
                ALLOCATE( wq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             wq_sgs_av = 0.0_wp


          CASE ( 'uv_sgs' )
             IF ( .NOT. ALLOCATED( uv_sgs_av ) )  THEN
                ALLOCATE( uv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             uv_sgs_av = 0.0_wp

          CASE ( 'uw_sgs' )
             IF ( .NOT. ALLOCATED( uw_sgs_av ) )  THEN
                ALLOCATE( uw_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             uw_sgs_av = 0.0_wp

          CASE ( 'utheta_sgs' )
             IF ( .NOT. ALLOCATED( utheta_sgs_av ) )  THEN
                ALLOCATE( utheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             utheta_sgs_av = 0.0_wp

          CASE ( 'uq_sgs' )
             IF ( .NOT. ALLOCATED( uq_sgs_av ) )  THEN
                ALLOCATE( uq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             uq_sgs_av = 0.0_wp

          CASE ( 'vu_sgs' )
             IF ( .NOT. ALLOCATED( vu_sgs_av ) )  THEN
                ALLOCATE( vu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             vu_sgs_av = 0.0_wp

          CASE ( 'vw_sgs' )
             IF ( .NOT. ALLOCATED( vw_sgs_av ) )  THEN
                ALLOCATE( vw_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             vw_sgs_av = 0.0_wp

          CASE ( 'vtheta_sgs' )
             IF ( .NOT. ALLOCATED( vtheta_sgs_av ) )  THEN
                ALLOCATE( vtheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             vtheta_sgs_av = 0.0_wp

          CASE ( 'vq_sgs' )
             IF ( .NOT. ALLOCATED( vq_sgs_av ) )  THEN
                ALLOCATE( vq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             vq_sgs_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

!
!--       Uncomment and extend the following lines, if necessary.
!--       The arrays for storing the user defined quantities (here u2 and u2_av) have to be declared
!--       and defined by the user!
!--       Sample for user-defined output:
          CASE ( 'pres_drag_norm_x*' )
             IF ( ALLOCATED( pres_drag_norm_x ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      pres_drag_norm_x_av(j,i) = pres_drag_norm_x_av(j,i) + pres_drag_norm_x(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pres_drag_norm_y*' )
             IF ( ALLOCATED( pres_drag_norm_y ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      pres_drag_norm_y_av(j,i) = pres_drag_norm_y_av(j,i) + pres_drag_norm_y(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'theta_product' )
             IF ( ALLOCATED( theta_product ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         theta_product_av(k,j,i) = theta_product_av(k,j,i) + theta_product(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'utheta_product' )
             IF ( ALLOCATED( utheta_product ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         utheta_product_av(k,j,i) = utheta_product_av(k,j,i) + utheta_product(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vtheta_product' )
             IF ( ALLOCATED( vtheta_product ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vtheta_product_av(k,j,i) = vtheta_product_av(k,j,i) + vtheta_product(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uq_product' )
             IF ( ALLOCATED( uq_product ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         uq_product_av(k,j,i) = uq_product_av(k,j,i) + uq_product(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vq_product' )
             IF ( ALLOCATED( vq_product ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vq_product_av(k,j,i) = vq_product_av(k,j,i) + vq_product(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wu_sgs' )
             IF ( ALLOCATED( wu_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wu_sgs_av(k,j,i) = wu_sgs_av(k,j,i) + wu_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wv_sgs' )
             IF ( ALLOCATED( wv_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wv_sgs_av(k,j,i) = wv_sgs_av(k,j,i) + wv_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wtheta_sgs' )
             IF ( ALLOCATED( wtheta_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wtheta_sgs_av(k,j,i) = wtheta_sgs_av(k,j,i) + wtheta_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wthetav_sgs' )
             IF ( ALLOCATED( wthetav_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wthetav_sgs_av(k,j,i) = wthetav_sgs_av(k,j,i) + wthetav_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wq_sgs' )
             IF ( ALLOCATED( wq_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wq_sgs_av(k,j,i) = wq_sgs_av(k,j,i) + wq_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF


          CASE ( 'uv_sgs' )
             IF ( ALLOCATED( uv_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         uv_sgs_av(k,j,i) = uv_sgs_av(k,j,i) + uv_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uw_sgs' )
             IF ( ALLOCATED( uw_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         uw_sgs_av(k,j,i) = uw_sgs_av(k,j,i) + uw_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'utheta_sgs' )
             IF ( ALLOCATED( utheta_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         utheta_sgs_av(k,j,i) = utheta_sgs_av(k,j,i) + utheta_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uq_sgs' )
             IF ( ALLOCATED( uq_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         uq_sgs_av(k,j,i) = uq_sgs_av(k,j,i) + uq_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vu_sgs' )
             IF ( ALLOCATED( vu_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vu_sgs_av(k,j,i) = vu_sgs_av(k,j,i) + vu_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vw_sgs' )
             IF ( ALLOCATED( vw_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vw_sgs_av(k,j,i) = vw_sgs_av(k,j,i) + vw_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vtheta_sgs' )
             IF ( ALLOCATED( vtheta_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vtheta_sgs_av(k,j,i) = vtheta_sgs_av(k,j,i) + vtheta_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vq_sgs' )
             IF ( ALLOCATED( vq_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vq_sgs_av(k,j,i) = vq_sgs_av(k,j,i) + vq_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

!
!--       Uncomment and extend the following lines, if necessary.
!--       The arrays for storing the user defined quantities (here u2_av) have to be declared and
!--       defined by the user!
!--       Sample for user-defined output:
          CASE ( 'pres_drag_norm_x*' )
             IF ( ALLOCATED( pres_drag_norm_x ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                     pres_drag_norm_x_av(j,i) = pres_drag_norm_x_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pres_drag_norm_y*' )
             IF ( ALLOCATED( pres_drag_norm_y ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                     pres_drag_norm_y_av(j,i) = pres_drag_norm_y_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'theta_product' )
             IF ( ALLOCATED( theta_product ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         theta_product_av(k,j,i) = theta_product_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF


          CASE ( 'utheta_product' )
             IF ( ALLOCATED( utheta_product ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         utheta_product_av(k,j,i) = utheta_product_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
         
          CASE ( 'vtheta_product' )
             IF ( ALLOCATED( vtheta_product ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vtheta_product_av(k,j,i) = vtheta_product_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uq_product' )
             IF ( ALLOCATED( uq_product ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         uq_product_av(k,j,i) = uq_product_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
         
          CASE ( 'vq_product' )
             IF ( ALLOCATED( vq_product ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vq_product_av(k,j,i) = vq_product_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF


          CASE ( 'wu_sgs' )
             IF ( ALLOCATED( wu_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wu_sgs_av(k,j,i) = wu_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wv_sgs' )
             IF ( ALLOCATED( wv_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wv_sgs_av(k,j,i) = wv_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
 
          CASE ( 'wtheta_sgs' )
             IF ( ALLOCATED( wtheta_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wtheta_sgs_av(k,j,i) = wtheta_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wthetav_sgs' )
             IF ( ALLOCATED( wthetav_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wthetav_sgs_av(k,j,i) = wthetav_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wq_sgs' )
             IF ( ALLOCATED( wq_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wq_sgs_av(k,j,i) = wq_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uv_sgs' )
             IF ( ALLOCATED( uv_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         uv_sgs_av(k,j,i) = uv_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uw_sgs' )
             IF ( ALLOCATED( uw_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         uw_sgs_av(k,j,i) = uw_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
 
          CASE ( 'utheta_sgs' )
             IF ( ALLOCATED( utheta_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         utheta_sgs_av(k,j,i) = utheta_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uq_sgs' )
             IF ( ALLOCATED( uq_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         uq_sgs_av(k,j,i) = uq_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vu_sgs' )
             IF ( ALLOCATED( vu_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vu_sgs_av(k,j,i) = vu_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vw_sgs' )
             IF ( ALLOCATED( vw_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vw_sgs_av(k,j,i) = vw_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
 
          CASE ( 'vtheta_sgs' )
             IF ( ALLOCATED( vtheta_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vtheta_sgs_av(k,j,i) = vtheta_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vq_sgs' )
             IF ( ALLOCATED( vq_sgs ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vq_sgs_av(k,j,i) = vq_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

       END SELECT

    ENDIF


 END SUBROUTINE user_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a temporary array with indices
!> (i,j,k) and sets the grid on which it is defined. Allowed values for grid are "zu" and "zw".
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_data_output_2d( av, variable, found, grid, local_pf, two_d, nzb_do, nzt_do )


    CHARACTER(LEN=*) ::  grid      !<
    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp) ::  av      !< flag to control data output of instantaneous or time-averaged data
    INTEGER(iwp) ::  i       !< grid index along x-direction
    INTEGER(iwp) ::  j       !< grid index along y-direction
    INTEGER(iwp) ::  k       !< grid index along z-direction
    !INTEGER(iwp) ::  m       !< running index surface elements
    INTEGER(iwp) ::  nzb_do  !< lower limit of the domain (usually nzb)
    INTEGER(iwp) ::  nzt_do  !< upper limit of the domain (usually nzt+1)

    LOGICAL      ::  found  !<
    LOGICAL      ::  two_d  !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !<



    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary.
!--    The arrays for storing the user defined quantities (here u2 and u2_av) have to be declared
!--    and defined by the user!
!--    Sample for user-defined output:
       CASE ( 'pres_drag_norm_x*' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = pres_drag_norm_x(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( pres_drag_norm_x_av ) )  THEN
                ALLOCATE( pres_drag_norm_x_av(nysg:nyng,nxlg:nxrg) )
                pres_drag_norm_x_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = pres_drag_norm_x_av(j,i)
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu1'
          two_d = .TRUE.

       CASE ( 'pres_drag_norm_y*' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = pres_drag_norm_y(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( pres_drag_norm_y_av ) )  THEN
                ALLOCATE( pres_drag_norm_y_av(nysg:nyng,nxlg:nxrg) )
                pres_drag_norm_y_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = pres_drag_norm_y_av(j,i)
                ENDDO
             ENDDO
          ENDIF
  
          grid = 'zu1'
          two_d = .TRUE.

       CASE ( 'theta_product_xy', 'theta_product_xz', 'theta_product_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = theta_product(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( theta_product_av ) )  THEN
                ALLOCATE( theta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                theta_product_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = theta_product_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu'

       CASE ( 'utheta_product_xy', 'utheta_product_xz', 'utheta_product_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = utheta_product(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( utheta_product_av ) )  THEN
                ALLOCATE( utheta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                utheta_product_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = utheta_product_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu'

       CASE ( 'vtheta_product_xy', 'vtheta_product_xz', 'vtheta_product_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vtheta_product(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vtheta_product_av ) )  THEN
                ALLOCATE( vtheta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vtheta_product_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vtheta_product_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu'

       CASE ( 'uq_product_xy', 'uq_product_xz', 'uq_product_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uq_product(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uq_product_av ) )  THEN
                ALLOCATE( uq_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uq_product_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uq_product_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu'

       CASE ( 'vq_product_xy', 'vq_product_xz', 'vq_product_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vq_product(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vq_product_av ) )  THEN
                ALLOCATE( vq_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vq_product_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vq_product_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu'

       CASE ( 'wu_sgs_xy', 'wu_sgs_xz', 'wu_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wu_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wu_sgs_av ) )  THEN
                ALLOCATE( wu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wu_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wu_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zw'

       CASE ( 'wv_sgs_xy', 'wv_sgs_xz', 'wv_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wv_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wv_sgs_av ) )  THEN
                ALLOCATE( wv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wv_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wv_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zw'

       CASE ( 'wtheta_sgs_xy', 'wtheta_sgs_xz', 'wtheta_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wtheta_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wtheta_sgs_av ) )  THEN
                ALLOCATE( wtheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wtheta_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wtheta_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zw'

       CASE ( 'wthetav_sgs_xy', 'wthetav_sgs_xz', 'wthetav_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wthetav_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wthetav_sgs_av ) )  THEN
                ALLOCATE( wthetav_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wthetav_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wthetav_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zw'

       CASE ( 'wq_sgs_xy', 'wq_sgs_xz', 'wq_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wq_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wq_sgs_av ) )  THEN
                ALLOCATE( wq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wq_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wq_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zw'

       CASE ( 'uv_sgs_xy', 'uv_sgs_xz', 'uv_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uv_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uv_sgs_av ) )  THEN
                ALLOCATE( uv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uv_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uv_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu'

       CASE ( 'uw_sgs_xy', 'uw_sgs_xz', 'uw_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uw_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uw_sgs_av ) )  THEN
                ALLOCATE( uw_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uw_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uw_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zw'

       CASE ( 'utheta_sgs_xy', 'utheta_sgs_xz', 'utheta_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = utheta_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( utheta_sgs_av ) )  THEN
                ALLOCATE( utheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                utheta_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = utheta_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu'

       CASE ( 'uq_sgs_xy', 'uq_sgs_xz', 'uq_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uq_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uq_sgs_av ) )  THEN
                ALLOCATE( uq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uq_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uq_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu'

       CASE ( 'vu_sgs_xy', 'vu_sgs_xz', 'vu_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vu_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vu_sgs_av ) )  THEN
                ALLOCATE( vu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vu_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vu_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu'

       CASE ( 'vw_sgs_xy', 'vw_sgs_xz', 'vw_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vw_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vw_sgs_av ) )  THEN
                ALLOCATE( vw_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vw_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vw_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zw'

       CASE ( 'vtheta_sgs_xy', 'vtheta_sgs_xz', 'vtheta_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vtheta_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vtheta_sgs_av ) )  THEN
                ALLOCATE( vtheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vtheta_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vtheta_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu'

       CASE ( 'vq_sgs_xy', 'vq_sgs_xz', 'vq_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vq_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vq_sgs_av ) )  THEN
                ALLOCATE( vq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vq_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vq_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zu'
!
!--    In case two-dimensional surface variables are output, the user has to access related
!--    surface-type. Uncomment and extend following lines appropriately (example output of vertical
!--    surface momentum flux of u-component). Please note, surface elements can be distributed over
!--    several data types, depending on their respective surface properties. The example shows how
!--    to sample data on horizontally-upward facing surfaces
!       CASE ( 'usws_xy' )
!          IF ( av == 0 )  THEN
!
!--           Horizontal default-type surfaces
!             DO  m = 1, surf_def%ns
!                i = surf_def%i(m)
!                j = surf_def%j(m)
!                local_pf(i,j,1) = MERGE( surf_def%usws(m),                                        &
!                                         local_pf(i,j,1),                                         &
!                                         surf_def%upward(m) )
!             ENDDO
!
!--           Horizontal natural-type surfaces
!             DO  m = 1, surf_lsm%ns
!                i = surf_lsm%i(m)
!                j = surf_lsmj(m)
!                local_pf(i,j,1) = MERGE( surf_lsm%usws(m),                                        &
!                                         local_pf(i,j,1),                                         &
!                                         surf_lsm%upward(m) )
!             ENDDO
!
!--           Horizontal urban-type surfaces
!             DO  m = 1, surf_usm%ns
!                i = surf_usm%i(m)
!                j = surf_usm%j(m)
!                local_pf(i,j,1) = MERGE( surf_usm%usws(m),                                        &
!                                         local_pf(i,j,1),                                         &
!                                         surf_usm%upward(m) )
!             ENDDO
!          ENDIF
!
!          grid = 'zu'
!--


       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT


 END SUBROUTINE user_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a temporary array with indices
!> (i,j,k).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )


    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp) ::  av     !<
    INTEGER(iwp) ::  i      !<
    INTEGER(iwp) ::  j      !<
    INTEGER(iwp) ::  k      !<
    INTEGER(iwp) ::  nzb_do !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL      ::  found  !<

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !<


!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( av == 0  .OR.  local_pf(nxl,nys,nzb_do) == 0.0_wp )  CONTINUE

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary.
!--    The arrays for storing the user defined quantities (here u2 and u2_av) have to be declared
!--    and defined by the user!
!--    Sample for user-defined output:
       CASE ( 'theta_product' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = theta_product(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( theta_product_av ) )  THEN
                ALLOCATE( theta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                theta_product_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = theta_product_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF


       CASE ( 'utheta_product' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = utheta_product(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( utheta_product_av ) )  THEN
                ALLOCATE( utheta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                utheta_product_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = utheta_product_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'vtheta_product' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vtheta_product(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vtheta_product_av ) )  THEN
                ALLOCATE( vtheta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vtheta_product_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vtheta_product_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'uq_product' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uq_product(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uq_product_av ) )  THEN
                ALLOCATE( uq_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uq_product_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uq_product_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'vq_product' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vq_product(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vq_product_av ) )  THEN
                ALLOCATE( vq_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vq_product_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vq_product_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'wu_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wu_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wu_sgs_av ) )  THEN
                ALLOCATE( wu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wu_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wu_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'wv_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wv_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wv_sgs_av ) )  THEN
                ALLOCATE( wv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wv_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wv_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'wtheta_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wtheta_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wtheta_sgs_av ) )  THEN
                ALLOCATE( wtheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wtheta_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wtheta_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'wthetav_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wthetav_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wthetav_sgs_av ) )  THEN
                ALLOCATE( wthetav_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wthetav_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wthetav_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'wq_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wq_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wq_sgs_av ) )  THEN
                ALLOCATE( wq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wq_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wq_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          
       CASE ( 'uv_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uv_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uv_sgs_av ) )  THEN
                ALLOCATE( uv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uv_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uv_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'uw_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uw_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uw_sgs_av ) )  THEN
                ALLOCATE( uw_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uw_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uw_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'utheta_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = utheta_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( utheta_sgs_av ) )  THEN
                ALLOCATE( utheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                utheta_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = utheta_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'uq_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uq_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uq_sgs_av ) )  THEN
                ALLOCATE( uq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uq_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = uq_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
!
       CASE ( 'vu_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vu_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vu_sgs_av ) )  THEN
                ALLOCATE( vu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vu_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vu_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'vw_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vw_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vw_sgs_av ) )  THEN
                ALLOCATE( vw_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vw_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vw_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'vtheta_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vtheta_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vtheta_sgs_av ) )  THEN
                ALLOCATE( vtheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vtheta_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vtheta_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'vq_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vq_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vq_sgs_av ) )  THEN
                ALLOCATE( vq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vq_sgs_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = vq_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
!

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE user_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of user-defined statistics, i.e. horizontally averaged profiles and time series.
!> This routine is called for every statistic region sr defined by the user, but at least for the
!> region "total domain" (sr=0). See section 3.5.4 on how to define, calculate, and output user
!> defined quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_statistics( mode, sr, tn )


    CHARACTER(LEN=*) ::  mode  !<
!    INTEGER(iwp) ::  i   !<
!    INTEGER(iwp) ::  j   !<
!    INTEGER(iwp) ::  k   !<
    INTEGER(iwp) ::  sr  !<
    INTEGER(iwp) ::  tn  !<

!    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ts_value_l  !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( sr == 0  .OR.  tn == 0 )  CONTINUE

    IF ( mode == 'profiles' )  THEN

!
!--    Sample on how to calculate horizontally averaged profiles of user-defined quantities. Each
!--    quantity is identified by the index "pr_palm+#" where "#" is an integer starting from 1.
!--    These user-profile-numbers must also be assigned to the respective strings given by
!--    data_output_pr_user in routine user_check_data_output_pr.
!       !$OMP DO
!       DO  i = nxl, nxr
!          DO  j = nys, nyn
!             DO  k = nzb+1, nzt
!!
!!--             Sample on how to calculate the profile of the resolved-scale horizontal momentum
!!--             flux u*v*
!                ind = pr_palm + max_pr_cs + max_pr_det + max_pr_salsa + 1
!                sums_l(k,ind,tn) = sums_l(k,pr_ind,tn) +                                           &
!                                   ( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) - hom(k,1,1,sr) ) *        &
!                                   ( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) - hom(k,1,2,sr) ) *        &
!                                   rmask(j,i,sr) * MERGE( 1.0_wp, 0.0_wp,                          &
!                                   BTEST( topo_flags(k,j,i), 0 ) )
!!
!!--             Further profiles can be defined and calculated by increasing the second index of
!!--             array sums_l (replace ... appropriately)
!                ind = pr_palm + max_pr_cs + max_pr_det + max_pr_salsa + 2
!                sums_l(k,ind,tn) = sums_l(k,ind,tn) + ...   * rmask(j,i,sr)
!             ENDDO
!          ENDDO
!       ENDDO

    ELSEIF ( mode == 'time_series' )  THEN


!       ALLOCATE ( ts_value_l(dots_num_user) )
!
!--    Sample on how to add values for the user-defined time series quantities.
!--    These have to be defined before in routine user_init. This sample creates two time series for
!--    the absolut values of the horizontal velocities u and v.
!       ts_value_l = 0.0_wp
!       ts_value_l(1) = ABS( u_max )
!       ts_value_l(2) = ABS( v_max )
!
!--     Collect / send values to PE0, because only PE0 outputs the time series.
!--     CAUTION: Collection is done by taking the sum over all processors. You may have to normalize
!--              this sum, depending on the quantity that you like to calculate. For serial runs,
!--              nothing has to be done.
!--     HINT: If the time series value that you are calculating has the same value on all PEs, you
!--           can omit the MPI_ALLREDUCE call and assign ts_value(dots_num_palm+1:,sr) = ts_value_l directly.
!#if defined( __parallel )
!       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
!       CALL MPI_ALLREDUCE( ts_value_l(1), ts_value(dots_num_palm+1,sr), dots_num_user, MPI_REAL,   &
!                           MPI_MAX, comm2d, ierr )
!#else
!       ts_value(dots_num_palm+1:dots_num_palm+dots_num_user,sr) = ts_value_l
!#endif

    ENDIF

 END SUBROUTINE user_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_rrd_global_ftn( found )


    LOGICAL, INTENT(OUT)  ::  found  !<


    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'global_paramter' )
!          READ ( 13 )  global_parameter

       CASE DEFAULT

          found = .FALSE.

    END SELECT


 END SUBROUTINE user_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_rrd_global_mpi

!    USE restart_data_mpi_io_mod,                                                                   &
!        ONLY:  rrd_mpi_io

!    CALL rrd_mpi_io( 'global_parameter', global_parameter )
    CONTINUE

 END SUBROUTINE user_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!> Subdomain
!> index limits on file are given by nxl_on_file, etc. Indices nxlc, etc. indicate the range of
!> gridpoints to be mapped from the subdomain on file (f) to the subdomain of the current PE (c).
!> They have been calculated in routine rrd_local.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf, nync,   &
                                nyn_on_file, nysf, nysc, nys_on_file, tmp_3d, found )


    INTEGER(iwp) ::  idum            !<
    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  nxlc            !<
    INTEGER(iwp) ::  nxlf            !<
    INTEGER(iwp) ::  nxl_on_file     !<
    INTEGER(iwp) ::  nxrc            !<
    INTEGER(iwp) ::  nxrf            !<
    INTEGER(iwp) ::  nxr_on_file     !<
    INTEGER(iwp) ::  nync            !<
    INTEGER(iwp) ::  nynf            !<
    INTEGER(iwp) ::  nyn_on_file     !<
    INTEGER(iwp) ::  nysc            !<
    INTEGER(iwp) ::  nysf            !<
    INTEGER(iwp) ::  nys_on_file     !<

    LOGICAL, INTENT(OUT)  ::  found  !<

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d  !<

    REAL(wp), DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_2d

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    idum = k + nxlc + nxlf + nxrc + nxrf + nync + nynf + nysc + nysf +                             &
           INT( tmp_3d(nzb,nys_on_file,nxl_on_file) )

!
!-- Here the reading of user-defined restart data follows:
!-- Sample for user-defined output

    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'pres_drag_norm_x_av' )
          IF ( .NOT. ALLOCATED( pres_drag_norm_x_av ) )  THEN
             ALLOCATE( pres_drag_norm_x_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          pres_drag_norm_x_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
             tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'pres_drag_norm_y_av' )
          IF ( .NOT. ALLOCATED( pres_drag_norm_y_av ) )  THEN
             ALLOCATE( pres_drag_norm_y_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          pres_drag_norm_y_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
             tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'theta_product_av' )
          IF ( .NOT. ALLOCATED( theta_product_av ) )  THEN
               ALLOCATE( theta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             theta_product_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                         &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)


       CASE ( 'utheta_product_av' )
          IF ( .NOT. ALLOCATED( utheta_product_av ) )  THEN
               ALLOCATE( utheta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             utheta_product_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                        &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'vtheta_product_av' )
          IF ( .NOT. ALLOCATED( vtheta_product_av ) )  THEN
               ALLOCATE( vtheta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             vtheta_product_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                        &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'uq_product_av' )
          IF ( .NOT. ALLOCATED( uq_product_av ) )  THEN
               ALLOCATE( uq_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             uq_product_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'vq_product_av' )
          IF ( .NOT. ALLOCATED( vq_product_av ) )  THEN
               ALLOCATE( vq_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             vq_product_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)  

       CASE ( 'wu_sgs_av' )
          IF ( .NOT. ALLOCATED( wu_sgs_av ) )  THEN
               ALLOCATE( wu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             wu_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wv_sgs_av' )
          IF ( .NOT. ALLOCATED( wv_sgs_av ) )  THEN
               ALLOCATE( wv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             wv_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wtheta_sgs_av' )
          IF ( .NOT. ALLOCATED( wtheta_sgs_av ) )  THEN
               ALLOCATE( wtheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             wtheta_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wthetav_sgs_av' )
          IF ( .NOT. ALLOCATED( wthetav_sgs_av ) )  THEN
               ALLOCATE( wthetav_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             wthetav_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wq_sgs_av' )
          IF ( .NOT. ALLOCATED( wq_sgs_av ) )  THEN
               ALLOCATE( wq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             wq_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'uv_sgs_av' )
          IF ( .NOT. ALLOCATED( uv_sgs_av ) )  THEN
               ALLOCATE( uv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             uv_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'uw_sgs_av' )
          IF ( .NOT. ALLOCATED( uw_sgs_av ) )  THEN
               ALLOCATE( uw_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             uw_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'utheta_sgs_av' )
          IF ( .NOT. ALLOCATED( utheta_sgs_av ) )  THEN
               ALLOCATE( utheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             utheta_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'uq_sgs_av' )
          IF ( .NOT. ALLOCATED( uq_sgs_av ) )  THEN
               ALLOCATE( uq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             uq_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'vu_sgs_av' )
          IF ( .NOT. ALLOCATED( vu_sgs_av ) )  THEN
               ALLOCATE( vu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             vu_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'vw_sgs_av' )
          IF ( .NOT. ALLOCATED( vw_sgs_av ) )  THEN
               ALLOCATE( vw_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             vw_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'vtheta_sgs_av' )
          IF ( .NOT. ALLOCATED( vtheta_sgs_av ) )  THEN
               ALLOCATE( vtheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             vtheta_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'vq_sgs_av' )
          IF ( .NOT. ALLOCATED( vq_sgs_av ) )  THEN
               ALLOCATE( vq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             vq_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)


       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE user_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_rrd_local_mpi

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array, rrd_mpi_io

    LOGICAL  ::  array_found     !<

!
!--  Restart input of time-averaged quantities is skipped in case of cyclic-fill initialization.
!--  This case, input of time-averaged data is useless and can lead to faulty averaging.
    IF ( .NOT. cyclic_fill_initialization )  THEN
       CALL rd_mpi_io_check_array( 'pres_drag_norm_x_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( pres_drag_norm_x_av ) )  THEN
             ALLOCATE( pres_drag_norm_x_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          CALL rrd_mpi_io( 'pres_drag_norm_x_av', pres_drag_norm_x_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'pres_drag_norm_y_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( pres_drag_norm_y_av ) )  THEN
             ALLOCATE( pres_drag_norm_y_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          CALL rrd_mpi_io( 'pres_drag_norm_y_av', pres_drag_norm_y_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'theta_product_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( theta_product_av ) )  THEN
             ALLOCATE( theta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          CALL rrd_mpi_io( 'theta_product_av', theta_product_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'utheta_product_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( utheta_product_av ) )  THEN
             ALLOCATE( utheta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          CALL rrd_mpi_io( 'utheta_product_av', utheta_product_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vtheta_product_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vtheta_product_av ) )  THEN
             ALLOCATE( vtheta_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          CALL rrd_mpi_io( 'vtheta_product_av', vtheta_product_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'uq_product_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( uq_product_av ) )  THEN
             ALLOCATE( uq_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          CALL rrd_mpi_io( 'uq_product_av', uq_product_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vq_product_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vq_product_av ) )  THEN
             ALLOCATE( vq_product_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          CALL rrd_mpi_io( 'vq_product_av', vq_product_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'wu_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( wu_sgs_av ) )  ALLOCATE( wu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'wu_sgs_av', wu_sgs_av )
       ENDIF
       
       CALL rd_mpi_io_check_array( 'wv_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( wv_sgs_av ) )  ALLOCATE( wv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'wv_sgs_av', wv_sgs_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'wtheta_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( wtheta_sgs_av ) )  ALLOCATE( wtheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'wtheta_sgs_av', wtheta_sgs_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'wthetav_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( wthetav_sgs_av ) )  ALLOCATE( wthetav_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'wthetav_sgs_av', wthetav_sgs_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'wq_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( wq_sgs_av ) )  ALLOCATE( wq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'wq_sgs_av', wq_sgs_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'uv_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( uv_sgs_av ) )  ALLOCATE( uv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'uv_sgs_av', uv_sgs_av )
       ENDIF
       
       CALL rd_mpi_io_check_array( 'uw_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( uw_sgs_av ) )  ALLOCATE( uw_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'uw_sgs_av', uw_sgs_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'utheta_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( utheta_sgs_av ) )  ALLOCATE( utheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'utheta_sgs_av', utheta_sgs_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'uq_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( uq_sgs_av ) )  ALLOCATE( uq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'uq_sgs_av', uq_sgs_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vu_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vu_sgs_av ) )  ALLOCATE( vu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vu_sgs_av', vu_sgs_av )
       ENDIF
       
       CALL rd_mpi_io_check_array( 'vw_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vw_sgs_av ) )  ALLOCATE( vw_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vw_sgs_av', vw_sgs_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vtheta_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vtheta_sgs_av ) )  ALLOCATE( vtheta_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vtheta_sgs_av', vtheta_sgs_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vq_sgs_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vq_sgs_av ) )  ALLOCATE( vq_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vq_sgs_av', vq_sgs_av )
       ENDIF
    ENDIF

    CONTINUE

 END SUBROUTINE user_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes global and user-defined restart data into binary file(s) for restart runs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_wrd_global

!    USE restart_data_mpi_io_mod,                                                                   &
!        ONLY:  wrd_mpi_io

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

!       CALL wrd_write_string( 'global_parameter' )
!       WRITE ( 14 )  global_parameter

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

!    CALL rrd_mpi_io( 'global_parameter', global_parameter )

    ENDIF

 END SUBROUTINE user_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes processor specific and user-defined restart data into binary file(s) for restart runs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_wrd_local

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  wrd_mpi_io

!
!-- Here the user-defined actions at the end of a job follow.
!-- Sample for user-defined output:
    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       IF ( ALLOCATED( pres_drag_norm_x_av ) )  THEN
          CALL wrd_write_string( 'pres_drag_norm_x_av' )
          WRITE ( 14 )  pres_drag_norm_x_av
       ENDIF

       IF ( ALLOCATED( pres_drag_norm_y_av ) )  THEN
         CALL wrd_write_string( 'pres_drag_norm_y_av' )
         WRITE ( 14 )  pres_drag_norm_y_av
       ENDIF

       IF ( ALLOCATED( theta_product_av ) )  THEN
          CALL wrd_write_string( 'theta_product_av' )
          WRITE ( 14 )  theta_product_av
       ENDIF

       IF ( ALLOCATED( utheta_product_av ) )  THEN
          CALL wrd_write_string( 'utheta_product_av' )
          WRITE ( 14 )  utheta_product_av
       ENDIF

       IF ( ALLOCATED( vtheta_product_av ) )  THEN
          CALL wrd_write_string( 'vtheta_product_av' )
          WRITE ( 14 )  vtheta_product_av
       ENDIF

       IF ( ALLOCATED( uq_product_av ) )  THEN
          CALL wrd_write_string( 'uq_product_av' )
          WRITE ( 14 )  uq_product_av
       ENDIF

       IF ( ALLOCATED( vq_product_av ) )  THEN
          CALL wrd_write_string( 'vq_product_av' )
          WRITE ( 14 )  vq_product_av
       ENDIF

       IF ( ALLOCATED( wu_sgs_av ) )  THEN
          CALL wrd_write_string( 'wu_sgs_av' )
          WRITE ( 14 )  wu_sgs_av
       ENDIF

       IF ( ALLOCATED( wv_sgs_av ) )  THEN
          CALL wrd_write_string( 'wv_sgs_av' )
          WRITE ( 14 )  wv_sgs_av
       ENDIF

       IF ( ALLOCATED( wtheta_sgs_av ) )  THEN
          CALL wrd_write_string( 'wtheta_sgs_av' )
          WRITE ( 14 )  wtheta_sgs_av
       ENDIF

       IF ( ALLOCATED( wthetav_sgs_av ) )  THEN
          CALL wrd_write_string( 'wthetav_sgs_av' )
          WRITE ( 14 )  wthetav_sgs_av
       ENDIF

       IF ( ALLOCATED( wq_sgs_av ) )  THEN
          CALL wrd_write_string( 'wq_sgs_av' )
          WRITE ( 14 )  wq_sgs_av
       ENDIF

       IF ( ALLOCATED( uv_sgs_av ) )  THEN
          CALL wrd_write_string( 'uv_sgs_av' )
          WRITE ( 14 )  uv_sgs_av
       ENDIF

       IF ( ALLOCATED( uw_sgs_av ) )  THEN
          CALL wrd_write_string( 'uw_sgs_av' )
          WRITE ( 14 )  uw_sgs_av
       ENDIF

       IF ( ALLOCATED( utheta_sgs_av ) )  THEN
          CALL wrd_write_string( 'utheta_sgs_av' )
          WRITE ( 14 )  utheta_sgs_av
       ENDIF

       IF ( ALLOCATED( uq_sgs_av ) )  THEN
          CALL wrd_write_string( 'uq_sgs_av' )
          WRITE ( 14 )  uq_sgs_av
       ENDIF

       IF ( ALLOCATED( vu_sgs_av ) )  THEN
          CALL wrd_write_string( 'vu_sgs_av' )
          WRITE ( 14 )  vu_sgs_av
       ENDIF

       IF ( ALLOCATED( vw_sgs_av ) )  THEN
          CALL wrd_write_string( 'vw_sgs_av' )
          WRITE ( 14 )  vw_sgs_av
       ENDIF

       IF ( ALLOCATED( vtheta_sgs_av ) )  THEN
          CALL wrd_write_string( 'vtheta_sgs_av' )
          WRITE ( 14 )  vtheta_sgs_av
       ENDIF
       
       IF ( ALLOCATED( vq_sgs_av ) )  THEN
          CALL wrd_write_string( 'vq_sgs_av' )
          WRITE ( 14 )  vq_sgs_av
       ENDIF

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       IF ( ALLOCATED( pres_drag_norm_x_av ) )  THEN
          CALL wrd_mpi_io( 'pres_drag_norm_x_av', pres_drag_norm_x_av )
       ENDIF

       IF ( ALLOCATED( pres_drag_norm_y_av ) )  THEN
          CALL wrd_mpi_io( 'pres_drag_norm_y_av', pres_drag_norm_y_av )
       ENDIF

       IF ( ALLOCATED( theta_product_av ) )  CALL wrd_mpi_io( 'theta_product_av', theta_product_av )
       IF ( ALLOCATED( utheta_product_av ) )  CALL wrd_mpi_io( 'utheta_product_av', utheta_product_av )
       IF ( ALLOCATED( vtheta_product_av ) )  CALL wrd_mpi_io( 'vtheta_product_av', vtheta_product_av )
       IF ( ALLOCATED( uq_product_av ) )  CALL wrd_mpi_io( 'uq_product_av', uq_product_av )
       IF ( ALLOCATED( vq_product_av ) )  CALL wrd_mpi_io( 'vq_product_av', vq_product_av )
       IF ( ALLOCATED( wu_sgs_av ) )  CALL wrd_mpi_io( 'wu_sgs_av', wu_sgs_av )
       IF ( ALLOCATED( wv_sgs_av ) )  CALL wrd_mpi_io( 'wv_sgs_av', wv_sgs_av )
       IF ( ALLOCATED( wtheta_sgs_av ) )  CALL wrd_mpi_io( 'wtheta_sgs_av', wtheta_sgs_av )
       IF ( ALLOCATED( wthetav_sgs_av ) )  CALL wrd_mpi_io( 'wthetav_sgs_av', wthetav_sgs_av )
       IF ( ALLOCATED( wq_sgs_av ) )  CALL wrd_mpi_io( 'wq_sgs_av', wq_sgs_av )


       IF ( ALLOCATED( uv_sgs_av ) )  CALL wrd_mpi_io( 'uv_sgs_av', uv_sgs_av )
       IF ( ALLOCATED( uw_sgs_av ) )  CALL wrd_mpi_io( 'uw_sgs_av', uw_sgs_av )
       IF ( ALLOCATED( utheta_sgs_av ) )  CALL wrd_mpi_io( 'utheta_sgs_av', utheta_sgs_av )
       IF ( ALLOCATED( uq_sgs_av ) )  CALL wrd_mpi_io( 'uq_sgs_av', uq_sgs_av )

       IF ( ALLOCATED( vu_sgs_av ) )  CALL wrd_mpi_io( 'vu_sgs_av', vu_sgs_av )
       IF ( ALLOCATED( vw_sgs_av ) )  CALL wrd_mpi_io( 'vw_sgs_av', vw_sgs_av )
       IF ( ALLOCATED( vtheta_sgs_av ) )  CALL wrd_mpi_io( 'vtheta_sgs_av', vtheta_sgs_av )
       IF ( ALLOCATED( vq_sgs_av ) )  CALL wrd_mpi_io( 'vq_sgs_av', vq_sgs_av )


    ENDIF

 END SUBROUTINE user_wrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execution of user-defined actions at the end of a job.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_last_actions

!
!-- Here the user-defined actions at the end of a job might follow.


 END SUBROUTINE user_last_actions

 END MODULE user
