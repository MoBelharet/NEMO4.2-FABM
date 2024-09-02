MODULE trcsms_fabm
   !!======================================================================
   !!                         ***  MODULE trcsms_fabm  ***
   !! TOP :   Main module of the FABM tracers
   !!======================================================================
   !! History :   1.0  !  2015-04  (PML) Original code
   !!----------------------------------------------------------------------
#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                               FABM tracers
   !!----------------------------------------------------------------------
   !! trc_sms_fabm       : FABM model main routine
   !! trc_sms_fabm_alloc : allocate arrays specific to FABM sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trcbc
   USE trd_oce
   USE trdtrc
#if defined key_trdtrc && defined key_iomput
   USE trdtrc_oce
#endif

   USE oce, only: ts,uu,vv  ! ocean state variables
   USE sbc_oce, only: lk_oasis,fr_i
   USE dom_oce !  depth and thicknesses
   USE zdf_oce !  diffusivity
   USE sbcapr  !  airpressure
   USE zdfdrg
   USE iom
   USE lib_mpp
   USE xios
   USE cpl_oasis3
   USE st2D_fabm
   USE inputs_fabm
   USE vertical_movement_fabm
   USE fabm
   USE fabm_types,only:type_interior_standard_variable

   IMPLICIT NONE

!#  include "vectopt_loop_substitute.h90"
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   PRIVATE

   PUBLIC   trc_sms_fabm       ! called by trcsms.F90 module
   PUBLIC   trc_sms_fabm_alloc ! called by trcini_fabm.F90 module
   PUBLIC   trc_sms_fabm_check_mass
   PUBLIC   st2d_fabm_nxt ! 2D state intergration
   PUBLIC   compute_fabm ! Compute FABM sources, sinks and diagnostics
   PUBLIC   check_state

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: flux    ! Cross-interface flux of pelagic variables (# m-2 s-1)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: current_total   ! Work array for mass aggregation

   ! Arrays for environmental variables
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:) :: prn,rho
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)   :: taubot
#if defined key_qco
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:) :: gdept_dummy, e3t_dummy
#endif
  !--------Mokrane ----------------------------------------------
  REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:) :: zcmask
  !-------------------------------------------------------------

   REAL(wp), PUBLIC, TARGET :: daynumber_in_year

   ! state check type
   TYPE, PUBLIC :: type_state
      LOGICAL             :: valid
      LOGICAL             :: repaired
   END TYPE

   ! State repair counters
   INTEGER, SAVE :: repair_interior_count = 0
   INTEGER, SAVE :: repair_surface_count  = 0
   INTEGER, SAVE :: repair_bottom_count   = 0

   ! Coupler parameters
   INTEGER, PUBLIC :: nn_adv  ! Vertical advection scheme for sinking/floating/movement
                              ! (1: 1st order upwind, 3: 3rd order TVD)

   ! Flag indicating whether model%start has been called (will be done on-demand)
!   LOGICAL, SAVE :: started = .false.

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL licence     (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_fabm( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_fabm  ***
      !!
      !! ** Purpose :   main routine of FABM model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kbb, Kmm, Krhs ! time level indices      
      INTEGER :: ji, jj, jn, jk
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztrfabm 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: pdat
      REAL(wp), DIMENSION(jpi,jpj)    :: vint

!!----------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('trc_sms_fabm')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,'(a,i0,a,i4.4,a,i2.2,a,i2.2,a,i5,a)') &
          ' trc_sms_fabm:  FABM model, iteration ',kt,' ', &
          nyear,'-',nmonth,'-',nday,' ',nsec_day," secs"
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      IF (kt == nittrc000) CALL nemo_fabm_start


      ! update links to FABM to ensure that FABM is pointing to the correct time index !
      ! Send pointers to state data to FABM
      do jn=1,jp_fabm
         CALL model%link_interior_state_data(jn,tr(:,:,:,jp_fabm_m1+jn,Kmm))
      end do
      DO jn=1,jp_fabm_surface
         CALL model%link_surface_state_data(jn,fabm_st2Dn(:,:,jn))
      END DO
      DO jn=1,jp_fabm_bottom
         CALL model%link_bottom_state_data(jn,fabm_st2Dn(:,:,jp_fabm_surface+jn))
      END DO

      ! Send pointers to environmental data to FABM
      CALL model%link_interior_data(fabm_standard_variables%temperature, ts(:,:,:,jp_tem,Kmm))
      CALL model%link_interior_data(fabm_standard_variables%practical_salinity, ts(:,:,:,jp_sal,Kmm))
#if defined key_qco
      DO_3D(0,0,0,0,1,jpkm1)
         gdept_dummy(ji,jj,jk) = (gdept_0(ji,jj,jk)*(1._wp+r3t(ji,jj,Kmm)))
         e3t_dummy(ji,jj,jk) = (e3t_0(ji,jj,jk)*(1._wp+r3t(ji,jj,Kmm)*tmask(ji,jj,jk)))
      END_3D

      CALL model%link_interior_data(fabm_standard_variables%depth, gdept_dummy(:,:,:))
      CALL model%link_interior_data(fabm_standard_variables%cell_thickness, e3t_dummy(:,:,:))
#else
      CALL model%link_interior_data(fabm_standard_variables%depth, gdept(:,:,:,Kmm))
      CALL model%link_interior_data(fabm_standard_variables%cell_thickness, e3t(:,:,:, Kmm))
#endif

      CALL update_inputs( kt )

      CALL compute_fabm( kt, Kbb, Kmm , Krhs )

      CALL compute_vertical_movement( kt, nn_adv, Kbb, Kmm, Krhs  )

      CALL st2d_fabm_nxt( kt, Kbb, Kmm, Krhs  )

      IF( l_trdtrc )  ztrfabm(:,:,:) = 0._wp

      ! CALL trc_bc       ( kt )        from NEMO4 boundary conditions now called with TRP
      ! CALL trc_rnf_fabm ( kt, Kmm, Krhs  ) ! River forcings 

      ! Send 3D diagnostics to output (these apply to time "n")
      DO jn = 1, size(model%interior_diagnostic_variables)
         IF (model%interior_diagnostic_variables(jn)%save) THEN
            ! Save 3D field
            pdat => model%get_interior_diagnostic_data(jn)
            CALL iom_put(model%interior_diagnostic_variables(jn)%name, pdat)

            ! Save depth integral if selected for output in XIOS
            IF (iom_use(TRIM(model%interior_diagnostic_variables(jn)%name)//'_VINT')) THEN
               vint = 0._wp
               DO jk = 1, jpkm1
                  vint = vint + pdat(:,:,jk) * e3t(:,:,jk,Kmm) * tmask(:,:,jk)
               END DO
               CALL iom_put(TRIM(model%interior_diagnostic_variables(jn)%name)//'_VINT', vint)
            END IF
         END IF
      END DO

      ! Send 2D diagnostics to output (these apply to time "n")
      DO jn = 1, size(model%horizontal_diagnostic_variables)
         IF (model%horizontal_diagnostic_variables(jn)%save) &
             CALL iom_put( model%horizontal_diagnostic_variables(jn)%name, model%get_horizontal_diagnostic_data(jn))
      END DO
      IF( l_trdtrc ) THEN      ! Save the trends in the mixed layer
          DO jn = jp_fabm0, jp_fabm1
            ztrfabm(:,:,:) = tr(:,:,:,jn, Krhs)
            CALL trd_trc( ztrfabm, jn, jptra_sms, kt, Kmm )   ! save trends
          END DO
      END IF
 
      IF( ln_timing )  CALL timing_stop('trc_sms_fabm')

   END SUBROUTINE trc_sms_fabm

   SUBROUTINE compute_fabm( kt , Kbb, Kmm , Krhs )
      INTEGER, INTENT(in) :: kt   ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kbb, Kmm, Krhs ! time level indices
      INTEGER :: ji,jj,jk,jn
!      TYPE(type_state) :: valid_state
      REAL(wp) :: zalfg, zztmp, zztmp2

!      ! Validate current model state (setting argument to .TRUE. enables repair=clipping)
!      valid_state%valid = .FALSE.
!      valid_state%repaired = .FALSE.
!      valid_state = check_state(.TRUE.)
!      IF (.NOT. valid_state%valid) THEN
!         WRITE(numout,*) "Invalid value in FABM encountered in area ",narea,"!!!"
!#if defined key_iomput
!         CALL xios_finalize                ! end mpp communications with xios
!         IF( lk_oasis ) CALL cpl_finalize    ! end coupling and mpp communications with OASIS
!#else
!         IF( lk_oasis ) THEN
!            CALL cpl_finalize              ! end coupling and mpp communications with OASIS
!         ELSE
!            IF( lk_mpp )   CALL mppstop    ! end mpp communications
!         ENDIF
!#endif
!      END IF

!      IF (valid_state%repaired) THEN
!         WRITE(numout,*) "Total interior repairs up to now on process",narea,":",repair_interior_count
!         WRITE(numout,*) "Total surface repairs up to now on process",narea,":",repair_surface_count
!         WRITE(numout,*) "Total bottom repairs up to now on process",narea,":",repair_bottom_count
!      ENDIF

      daynumber_in_year = fjulday - fjulstartyear + 1

      ! Update the dummy depths if key_qco is defined
      ! ------------------------------------
#if defined key_qco
      DO_3D(0,0,0,0,1,jpkm1)
         gdept_dummy(ji,jj,jk) = (gdept_0(ji,jj,jk)*(1._wp+r3t(ji,jj,Kmm)))
         e3t_dummy(ji,jj,jk) = (e3t_0(ji,jj,jk)*(1._wp+r3t(ji,jj,Kmm)*tmask(ji,jj,jk)))
      END_3D
#endif

      ! Compute the now hydrostatic pressure (copied from istate.F90 NEMO3.6)
      ! ------------------------------------

      IF (ALLOCATED(rho)) rho = rho0 * ( 1._wp + rhd )      
      IF (ALLOCATED(prn)) THEN
         zalfg = 0.5e-4_wp * grav ! FABM wants dbar, convert from Pa (and multiply with 0.5 to average 2 cell thicknesses below)
         prn(:,:,1) = 10.1325_wp + zalfg * e3t(:,:,1, Kmm) * rho(:,:,1)
         DO jk = 2, jpkm1                                              ! Vertical integration from the surface
            prn(:,:,jk) = prn(:,:,jk-1) + zalfg * ( &
                        e3t(:,:,jk-1, Kmm) * rho(:,:,jk-1)  &
                        + e3t(:,:,jk, Kmm) * rho(:,:,jk) )
         END DO
      END IF

      ! Compute the bottom stress (copied from diawri.F90)
      ! ------------------------------------
      IF (ALLOCATED(taubot)) THEN
         CALL zdf_drg( kt, Kmm, mbkt , r_Cdmin_bot, r_Cdmax_bot,   &   ! <<== in
            &              r_z0_bot,   r_ke0_bot,    rCd0_bot,   &
            &                                        rz0_bot,rCdU_bot )     ! ==>> out : bottom drag [m/s]

         zztmp = rho0 * 0.25_wp
         taubot(:,:) = 0._wp
         DO_2D( 0, 0, 0, 0 )
            zztmp2 = (  ( rCdU_bot(ji+1,jj)+rCdU_bot(ji  ,jj) ) * uu(ji  ,jj,mbku(ji  ,jj),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji  ,jj)+rCdU_bot(ji-1,jj) ) * uu(ji-1,jj,mbku(ji-1,jj),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj  ) ) * vv(ji,jj  ,mbkv(ji,jj  ),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji,jj  )+rCdU_bot(ji,jj-1) ) * vv(ji,jj-1,mbkv(ji,jj-1),Kmm)  )**2
            taubot(ji,jj) = zztmp * SQRT( zztmp2 ) * tmask(ji,jj,1)
            !
         END_2D

      END IF

      CALL model%prepare_inputs(real(kt, wp),nyear,nmonth,nday,REAL(nsec_day,wp))

      ! Zero rate array of interface-attached state variables
      fabm_st2Da(:,:,:) = 0._wp

      ! Compute interfacial source terms and fluxes
      DO jj=ntsj, ntej  !2,jpjm1  TODO cleanup old loop indexes
         ! Process bottom (get_bottom_sources increments rather than sets, so zero flux array first)
         flux(:,:) = 0._wp
         CALL model%get_bottom_sources(ntsi,ntei,jj,flux,fabm_st2Da(ntsi:ntei,jj,jp_fabm_surface+1:))
         DO jn=1,jp_fabm
            ! Divide bottom fluxes by height of bottom layer and add to source terms.
            DO ji=ntsi,ntei !fs_2,fs_jpim1
               tr(ji,jj,mbkt(ji,jj),jp_fabm_m1+jn, Krhs) = tr(ji,jj,mbkt(ji,jj),jp_fabm_m1+jn, Krhs) + flux(ji,jn)/e3t(ji,jj,mbkt(ji,jj), Kmm)
            END DO
         END DO
         ! Process surface (fabm_do_surface increments rather than sets, so zero flux array first)
         flux(:,:) = 0._wp
         CALL model%get_surface_sources(ntsi,ntei,jj,flux,fabm_st2Da(ntsi:ntei,jj,1:jp_fabm_surface))
         ! Divide surface fluxes by height of surface layer and add to source terms.
         DO jn=1,jp_fabm
            DO ji=ntsi,ntei
          tr(ji,jj,1,jp_fabm_m1+jn, Krhs) = tr(ji,jj,1,jp_fabm_m1+jn, Krhs) + flux(ji,jn)/e3t(ji,jj,1, Kmm)
            END DO ! ji
         END DO ! jn
      END DO ! jj

      ! Compute interior source terms (NB fabm_do increments rather than sets)
      DO jk=1,jpkm1
          DO jj=ntsj, ntej
            CALL model%get_interior_sources(ntsi,ntei,jj,jk,tr(ntsi:ntei,jj,jk,jp_fabm0:jp_fabm1, Krhs))
          END DO
      END DO

      CALL model%finalize_outputs()

   END SUBROUTINE compute_fabm

   FUNCTION check_state(repair) RESULT(exit_state)
      LOGICAL, INTENT(IN) :: repair
      TYPE(type_state) :: exit_state

      INTEGER             :: jj,jk
      LOGICAL             :: valid_int,valid_sf,valid_bt

      exit_state%valid = .TRUE.
      exit_state%repaired =.FALSE.
      valid_int = .FALSE.
      valid_sf = .FALSE.
      valid_bt = .FALSE.
      DO jk=1,jpkm1
         DO jj=ntsj, ntej
            CALL model%check_interior_state(ntsi,ntei,jj,jk,repair,valid_int)
            IF (repair.AND..NOT.valid_int) THEN
               repair_interior_count = repair_interior_count + 1
               exit_state%repaired = .TRUE.
            END IF
            IF (.NOT.(valid_int.OR.repair)) exit_state%valid = .FALSE.
         END DO
      END DO
      DO jj=ntsj, ntej
         CALL model%check_surface_state(ntsi,ntei,jj,repair,valid_sf)
         IF (repair.AND..NOT.valid_sf) THEN
            repair_surface_count = repair_surface_count + 1
            exit_state%repaired = .TRUE.
         END IF
         IF (.NOT.(valid_sf.AND.valid_bt).AND..NOT.repair) exit_state%valid = .FALSE.
         CALL model%check_bottom_state(ntsi,ntei,jj,repair,valid_bt)
         IF (repair.AND..NOT.valid_bt) THEN
            repair_bottom_count = repair_bottom_count + 1
            exit_state%repaired = .TRUE.
         END IF
         IF (.NOT.(valid_sf.AND.valid_bt).AND..NOT.repair) exit_state%valid = .FALSE.
      END DO
   END FUNCTION



   SUBROUTINE trc_sms_fabm_check_mass()
      REAL(wp) :: total(SIZE(model%conserved_quantities))
      REAL(wp) :: fullTotal(SIZE(model%conserved_quantities))
      INTEGER :: ji,jk,jj,jn

      total(:) = 0._wp
      !IF (.NOT. started) CALL nemo_fabm_start

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_fabm_check_mass:  Total conserved quantities'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      DO jk=1,jpkm1
         DO jj=ntsj, ntej
            CALL model%get_interior_conserved_quantities(ntsi,ntei,jj,jk,current_total)
            DO jn=1,SIZE(model%conserved_quantities)
               DO ji=ntsi,ntei
                  total(jn) = total(jn) + cvol(ji,jj,jk) * current_total(ji,jn) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
      END DO

      DO jj=ntsj, ntej
         CALL model%get_horizontal_conserved_quantities(ntsi,ntei,jj,current_total)
         DO jn=1,SIZE(model%conserved_quantities)
            DO ji=ntsi,ntei
               total(jn) = total(jn) + e1e2t(ji,jj) * current_total(ji,jn) * tmask_i(ji,jj)
            END DO
         END DO
      END DO

      IF( lk_mpp ) CALL mpp_sum('trcsms_fabm',total)
      fullTotal = REAL(total,wp)

      DO jn=1,SIZE(model%conserved_quantities)
         IF(lwp) WRITE(numout,*) 'FABM '//TRIM(model%conserved_quantities(jn)%name),fullTotal(jn),TRIM(model%conserved_quantities(jn)%units)//'*m3'
      END DO

   END SUBROUTINE trc_sms_fabm_check_mass

   SUBROUTINE st2d_fabm_nxt( kt, Kbb, Kmm, Krhs   )
      !!----------------------------------------------------------------------
      !!                     ***  st2d_fabm_nxt  ***
      !!
      !! ** Purpose :   routine to integrate 2d states in time
      !!
      !! ** Method  :   based on integration of 3D passive tracer fields
      !!                implemented in TOP/TRP/trcnxt.F90, plus
      !!                tra_nxt_fix in OCE/TRA/tranxt.F90. Similar to
      !!                time integration of sea surface height in
      !!                OCE/DYN/sshwzv.F90.
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kbb, Kmm, Krhs ! time level indices        
      REAL(wp) :: z2dt
      INTEGER :: jn

!!----------------------------------------------------------------------
      !
! TODO check EULER stuff
!      IF( neuler == 0 .AND. kt == nittrc000 ) THEN
      IF ( l_1st_euler .OR. ln_top_euler ) THEN
          z2dt = rdt                  ! set time step size (Euler)
      ELSE
          z2dt = 2._wp * rdt          ! set time step size (Leapfrog)
      ENDIF
!          z2dt = rdt  
      ! Forward Euler time step to compute "now"
      DO jn=1,jp_fabm_surface+jp_fabm_bottom
         fabm_st2Da(:,:,jn) = (fabm_st2Db(:,:,jn) + z2dt * fabm_st2Da(:,:,jn)) * tmask(:,:,1)
      ENDDO
      
!      IF( neuler == 0 .AND. kt == nittrc000 ) THEN        ! Euler time-stepping at first time-step
      IF ( l_1st_euler .OR. ln_top_euler ) THEN     
         !                                                ! (only swap)
         fabm_st2Dn(:,:,:) = fabm_st2Da(:,:,:)
! TODO   rewrite  fabm_st2Db  using  time indexes
!          as they get swaped in OCE/step.F90    sets "before" fields = "now" fields.
         fabm_st2Db(:,:,:)  = fabm_st2Dn(:,:,:) ; !dows this fix this issue?

         !
      ELSE
         ! Update now state + Asselin filter time stepping
         fabm_st2Db(:,:,:) = (1._wp - 2._wp*rn_atfp) * fabm_st2Dn(:,:,:) + &
             rn_atfp * ( fabm_st2Db(:,:,:) + fabm_st2Da(:,:,:) )
         fabm_st2Dn(:,:,:) = fabm_st2Da(:,:,:)
      ENDIF

   END SUBROUTINE st2d_fabm_nxt

   INTEGER FUNCTION trc_sms_fabm_alloc(Kmm)
      INTEGER, INTENT( in ) ::    Kmm  ! time level indices    
      INTEGER :: jn
      INTEGER :: ji,jj,jk
      !------- Mokrane ---------
      INTEGER  :: ik50                !  last level where depth less than 50 m
      REAL(wp) :: zsurfc, zsurfp,ze3t, ze3t2, zcslp
      REAL(wp), PARAMETER :: distcoast = 5.e3
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_fabm_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to FABM
      ALLOCATE( lk_rad_fabm(jp_fabm))
      IF (model%variable_needs_values(fabm_standard_variables%pressure)) ALLOCATE(prn(jpi, jpj, jpk))
      IF (ALLOCATED(prn) .or. model%variable_needs_values(fabm_standard_variables%density)) ALLOCATE(rho(jpi, jpj, jpk))
      IF (model%variable_needs_values(fabm_standard_variables%bottom_stress)) ALLOCATE(taubot(jpi, jpj))

      ! Allocate arrays to hold state for surface-attached and bottom-attached state variables
      ALLOCATE(fabm_st2Dn(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))
      ALLOCATE(fabm_st2Da(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))
      ALLOCATE(fabm_st2Db(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))

      ! Work array to hold surface and bottom fluxes
      ALLOCATE(flux(ntsi:ntei,jp_fabm))

      ! Allocate work arrays for vertical movement
      ALLOCATE(w_ct(ntsi:ntei,1:jpkm1,jp_fabm))
      ALLOCATE(current_total(ntsi:ntei,SIZE(model%conserved_quantities)))
#if defined key_trdtrc && defined key_iomput
      IF( lk_trdtrc ) ALLOCATE(tr_vmv(jpi,jpj,jpk,jp_fabm))
      IF( lk_trdtrc ) ALLOCATE(tr_inp(jpi,jpj,jpk))
#endif
#if defined key_qco
      ALLOCATE(gdept_dummy(jpi,jpj,jpk), e3t_dummy(jpi,jpj,jpk))
#endif 
      
      !------------ Mokrane --------------
      ALLOCATE(zcmask(jpi,jpj,jpk) )
      !-----------------------------------

      trc_sms_fabm_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_fabm_alloc /= 0 ) CALL ctl_warn('trc_sms_fabm_alloc : failed to allocate arrays')
      !

      ! Provide FABM with domain extents
      CALL model%set_domain(jpi, jpj, jpk, rdt)
      CALL model%set_domain_start(ntsi, ntsj, 1)
      CALL model%set_domain_stop(ntei, ntej, jpkm1)

      ! Provide FABM with the vertical indices of the surface and bottom, and the land-sea mask.
      call model%set_bottom_index(mbkt)  ! NB mbkt extents should match dimension lengths provided to model%set_domain
      call model%set_mask(tmask,tmask(:,:,1)) ! NB tmask extents should match dimension lengths provided to model%set_domain

      ! Send pointers to state data to FABM
      ! TODO Isn't it initialization ???       
      do jn=1,jp_fabm
        tr(:,:,:,jp_fabm_m1+jn,Kmm) = model%interior_state_variables(jn)%initial_value * tmask
        call model%link_interior_state_data(jn,tr(:,:,:,jp_fabm_m1+jn,Kmm))
      end do
      DO jn=1,jp_fabm_surface
        fabm_st2Dn(:,:,jn) = model%surface_state_variables(jn)%initial_value * tmask(:,:,1)
        CALL model%link_surface_state_data(jn,fabm_st2Dn(:,:,jn))
      END DO
      DO jn=1,jp_fabm_bottom
         fabm_st2Dn(:,:,jp_fabm_surface+jn) = model%bottom_state_variables(jn)%initial_value * tmask(:,:,1)
         CALL model%link_bottom_state_data(jn,fabm_st2Dn(:,:,jp_fabm_surface+jn))
      END DO

      ! Send pointers to environmental data to FABM
      
      CALL model%link_interior_data(fabm_standard_variables%temperature, ts(:,:,:,jp_tem,Kmm))
      CALL model%link_horizontal_data(fabm_standard_variables%surface_temperature, ts(:,:,1,jp_tem,Kmm)) ! Mokrane
      CALL model%link_interior_data(fabm_standard_variables%practical_salinity, ts(:,:,:,jp_sal,Kmm))
      IF (ALLOCATED(rho)) CALL model%link_interior_data(fabm_standard_variables%density, rho(:,:,:))
      IF (ALLOCATED(prn)) CALL model%link_interior_data(fabm_standard_variables%pressure, prn)
      IF (ALLOCATED(taubot)) CALL model%link_horizontal_data(fabm_standard_variables%bottom_stress, taubot(:,:))
      CALL model%link_horizontal_data(fabm_standard_variables%latitude, gphit)
      CALL model%link_horizontal_data(fabm_standard_variables%longitude, glamt)
      CALL model%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year, daynumber_in_year)
      CALL model%link_horizontal_data(fabm_standard_variables%wind_speed, wndm(:,:))
      CALL model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux, qsr(:,:))
      CALL model%link_horizontal_data(fabm_standard_variables%surface_mean_downwelling_shortwave_flux,qsr_mean(:,:))
      CALL model%link_horizontal_data(fabm_standard_variables%bottom_depth_below_geoid, ht_0(:,:)) ! Mokrane
      CALL model%link_horizontal_data(fabm_standard_variables%ice_area_fraction, fr_i(:,:))
! TODO  rewrite FABM models such that they wouldn't include temporal integration!
      !CALL model%link_scalar(fabm_standard_variables%timestep, rdt)
! TODO replace dummy values if you are using IOW ERGOM v2 and CDOM addon
      CALL model%link_horizontal_data(fabm_standard_variables%bottom_depth, ht_0(:,:))
      ! use avs because temperature diffusivity (avt) is/might be affected by doublediffusivity
      !CALL model%link_interior_data(fabm_standard_variables%diffusivity, avs(:,:,:))
      CALL model%link_interior_data(type_interior_standard_variable(name='vertical_tracer_diffusivity', units='m2 s-1'), avt(:,:,:)) !Mokrane    
      ! in some FABM versions, vertical_tracer_diffusivity is a standard variable
      ! CALL model%link_interior_data(fabm_standard_variables%vertical_tracer_diffusivity, avs(:,:,:))
      CALL model%link_horizontal_data(fabm_standard_variables%surface_air_pressure, apr(:,:))
      CALL model%link_horizontal_data(fabm_standard_variables%mixed_layer_thickness_defined_by_vertical_tracer_diffusivity, hmld(:,:)) ! Mokrane
      !CALL model%link_scalar(type_global_standard_variable(name='physical_time_step', units='s'), rDt_trc)
      !CALL model%link_scalar(fabm_standard_variables%physical_time_step, rDt_trc)
      
      CALL model%link_horizontal_data(fabm_standard_variables%cell_area,  e1e2t(:,:)) ! Mokrane: grid cell area
      CALL model%link_horizontal_data(model%get_horizontal_variable_id('fmmflx'), fmmflx(:,:))
      CALL model%link_interior_data(type_interior_standard_variable(name='ref_cell_thickness' , units='m'), e3t_0(:,:,:))
      !---------- Mokrane --------------------

      ik50 = 5        !  last level where depth less than 50 m
      DO jk = jpkm1, 1, -1
         IF( gdept_1d(jk) > 50. )   ik50 = jk - 1
      END DO

      zcmask(:,:,:) = 0.

      DO_3D( 0, 0, 0, 0, 1, ik50 )
         ze3t   = e3t_0(ji,jj,jk)
         zsurfc =  e1u(ji,jj) * ( 1. - umask(ji  ,jj  ,jk) )   &
                 + e1u(ji,jj) * ( 1. - umask(ji-1,jj  ,jk) )   &
                 + e2v(ji,jj) * ( 1. - vmask(ji  ,jj  ,jk) )   &
                 + e2v(ji,jj) * ( 1. - vmask(ji  ,jj-1,jk) )
         zsurfp = zsurfc * ze3t / e1e2t(ji,jj)
         ! estimation of the coastal slope : 5 km off the coast
         ze3t2 = ze3t * ze3t
         zcslp = SQRT( ( distcoast*distcoast + ze3t2 ) / ze3t2 )
         zcmask(ji,jj,jk) = zcslp * zsurfp
       
      END_3D


      CALL model%link_interior_data(type_interior_standard_variable(name='coastal_island_mask', units='1'), zcmask(:,:,:) )


      !----------------------------------------
      

#if defined key_qco
      DO_3D(0,0,0,0,1,jpkm1)
         gdept_dummy(ji,jj,jk) = (gdept_0(ji,jj,jk)*(1._wp+r3t(ji,jj,Kmm)))
         e3t_dummy(ji,jj,jk) = (e3t_0(ji,jj,jk)*(1._wp+r3t(ji,jj,Kmm)*tmask(ji,jj,jk)))
      END_3D

     CALL model%link_interior_data(fabm_standard_variables%depth, gdept_dummy(:,:,:))
     CALL model%link_interior_data(fabm_standard_variables%cell_thickness, e3t_dummy(:,:,:)) 
#else
     CALL model%link_interior_data(fabm_standard_variables%depth, gdept(:,:,:,Kmm))
     CALL model%link_interior_data(fabm_standard_variables%cell_thickness, e3t(:,:,:, Kmm))
#endif

      ! Obtain user-specified input variables (read from NetCDF file)
      call link_inputs
      call update_inputs( nit000, .false. )

      ! Set mask for negativity corrections to the relevant states
      lk_rad_fabm(:) = .FALSE.
      DO jn=1,jp_fabm
        IF (model%interior_state_variables(jn)%minimum >= 0._wp) THEN
          lk_rad_fabm(jn) = .TRUE.
          IF(lwp) WRITE(numout,*) 'FABM clipping for '//TRIM(model%interior_state_variables(jn)%name)//' activated.'
        END IF
      END DO

      ! Copy initial condition for interface-attached state variables to "previous" state field
      ! NB NEMO does this itself for pelagic state variables (trb) in TOP/trcini.F90.
      fabm_st2Db = fabm_st2Dn

   END FUNCTION trc_sms_fabm_alloc

   SUBROUTINE nemo_fabm_start()
      INTEGER :: jn

      ! Make FABM aware of diagnostics that are not needed [not included in output]
      ! This works only after iom has completely initialised, because it depends on iom_use
      DO jn=1,size(model%interior_diagnostic_variables)
         model%interior_diagnostic_variables(jn)%save = iom_use(model%interior_diagnostic_variables(jn)%name) &
            .or. iom_use(TRIM(model%interior_diagnostic_variables(jn)%name)//'_VINT')
      END DO
      DO jn=1,size(model%horizontal_diagnostic_variables)
         model%horizontal_diagnostic_variables(jn)%save = iom_use(model%horizontal_diagnostic_variables(jn)%name)
      END DO

      ! Check whether FABM has all required data
      ! [after this, the save attribute of diagnostic variables can no longer change!]
      CALL model%start()

      !started = .TRUE.
   END SUBROUTINE

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No FABM model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_fabm( kt, Kbb, Kmm, Krhs )              ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      INTEGER, INTENT( in ) ::   Kbb, Kmm, Krhs ! time level indices
      WRITE(*,*) 'trc_sms_fabm: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_fabm
#endif

   !!======================================================================
END MODULE trcsms_fabm
