MODULE trcwri_fabm
   !!======================================================================
   !!                       *** MODULE trcwri_fabm ***
   !!    fabm :   Output of FABM tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && key_fabm && defined key_xios
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                           FABM model
   !!----------------------------------------------------------------------
   !! trc_wri_fabm   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   USE sbc_oce, only: lk_oasis
   USE trcsms_fabm!, only: trc_sms_fabm_check_mass, check_state
   USE par_fabm
   USE lib_mpp
   USE xios
   USE cpl_oasis3
   USE st2d_fabm
   USE,INTRINSIC :: iso_fortran_env, only: output_unit

   IMPLICIT NONE
   PRIVATE

#if defined key_tracer_budget
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: tr_temp
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: fabm_st2d_temp
#endif

   ! State repair counters
   INTEGER, SAVE :: repair_interior_count = 0
   INTEGER, SAVE :: repair_surface_count  = 0
   INTEGER, SAVE :: repair_bottom_count   = 0

   INTERFACE trc_wri_fabm
       MODULE PROCEDURE wri_fabm,wri_fabm_fl
   END INTERFACE trc_wri_fabm

   PUBLIC trc_wri_fabm 

CONTAINS

   SUBROUTINE wri_fabm_fl (kt, Kmm, fl)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )               :: fl
      INTEGER, INTENT( in )               :: kt
      INTEGER, INTENT( in )               :: Kmm  ! time level indices 
#if defined key_tracer_budget
      INTEGER              :: jn
      CHARACTER (len=20)   :: cltra
      REAL(wp), DIMENSION(jpi,jpj,jpk)    :: trpool !temporary storage pool 3D
      REAL(wp), DIMENSION(jpi,jpj)        :: st2dpool !temporary storage pool 2D
      !!---------------------------------------------------------------------

      ! write the tracer concentrations in the file
      ! ---------------------------------------
! depth integrated
! for strict budgetting write this out at end of timestep as an average between 'now' and 'after' at kt
      DO jn = 1, jp_fabm1
        IF(ln_trdtrc (jp_fabm_m1+jn))THEN
         trpool(:,:,:) = 0.5 * ( tr(:,:,:,jp_fabm_m1+jn, Kmm)*e3t(:,:,:, Kmm) + &
                                 tr_temp(:,:,:,jn)*e3t(:,:,:, Kmm) )
         cltra = TRIM( model%interior_state_variables(jn)%name )//"_e3t"     ! depth integrated output
         IF( kt == nittrc000 ) write(output_unit,*)'output pool ',cltra
         CALL iom_put( cltra, trpool)
        ENDIF
      END DO
#else
      CONTINUE
#endif

   END SUBROUTINE wri_fabm_fl


   SUBROUTINE wri_fabm (kt, Kmm)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )               :: kt
      INTEGER, INTENT( in )   ::    Kmm  ! time level indices       
      INTEGER                 ::    jn, jk
      REAL(wp), DIMENSION(jpi,jpj)    :: vint
      TYPE(type_state) :: valid_state
      !!---------------------------------------------------------------------

      ! ----- Mokrane --------
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
      !--------------------------------------
      
      IF ( kt > nittrc000 ) THEN !can only check state if FABM model has been started!
         ! Validate current model state (setting argument to .TRUE. enables repair=clipping)
         valid_state%valid = .FALSE.
         valid_state%repaired = .FALSE.
         valid_state = check_state(.TRUE.)
         IF (.NOT. valid_state%valid) THEN
            WRITE(numout,*) "Invalid value in FABM encountered in area ",narea,"!!!"
#if defined key_iomput
            CALL xios_finalize                ! end mpp communications with xios
            IF( lk_oasis ) CALL cpl_finalize    ! end coupling and mpp communications with OASIS
#else
            IF( lk_oasis ) THEN
               CALL cpl_finalize              ! end coupling and mpp communications with OASIS
            ELSE
               IF( lk_mpp )   CALL mppstop    ! end mpp communications
            ENDIF
#endif
         END IF
      END IF

!      IF (valid_state%repaired) THEN
!         WRITE(numout,*) "Total interior repairs up to now on process",narea,":",repair_interior_count
!         WRITE(numout,*) "Total surface repairs up to now on process",narea,":",repair_surface_count
!         WRITE(numout,*) "Total bottom repairs up to now on process",narea,":",repair_bottom_count
!      ENDIF

#if defined key_tracer_budget
      IF( kt == nittrc000 ) THEN
         ALLOCATE(tr_temp(jpi,jpj,jpk,jp_fabm),fabm_st2d_temp(jpi,jpj,jp_fabm_surface+jp_fabm_bottom))
      ENDIF
      tr_temp(:,:,:,:)      = tr(:,:,:,jp_fabm0:jp_fabm1, Kmm) ! slwa save for tracer budget (unfiltered trn)
      fabm_st2d_temp(:,:,:) = fabm_st2dn(:,:,:)
#endif
      DO jn = 1, jp_fabm
         ! Save 3D field
         CALL iom_put(model%interior_state_variables(jn)%name, tr(:,:,:,jp_fabm_m1+jn,Kmm))

         ! Save depth integral if selected for output in XIOS
         IF (iom_use(TRIM(model%interior_state_variables(jn)%name)//'_VINT')) THEN
            vint = 0._wp
            DO jk = 1, jpkm1
            ! TODO  fix the vertical integration here; something with the e3t(:,:,jk,Kmm)
               vint = vint !+ tr(:,:,jk,jp_fabm_m1+jn,Kmm) * e3t(:,:,jk,Kmm) * tmask(:,:,jk)
            END DO
            CALL iom_put(TRIM(model%interior_state_variables(jn)%name)//'_VINT', vint)
         END IF
      END DO
      DO jn = 1, jp_fabm_surface
         CALL iom_put( model%surface_state_variables(jn)%name, fabm_st2dn(:,:,jn) )
      END DO
      DO jn = 1, jp_fabm_bottom
         CALL iom_put( model%bottom_state_variables(jn)%name, fabm_st2dn(:,:,jp_fabm_surface+jn) )
      END DO

      ! write 3D diagnostics in the file
      ! ---------------------------------------
!      DO jn = 1, size(model%interior_diagnostic_variables)
!         IF (model%interior_diagnostic_variables(jn)%save) &
!             CALL iom_put( model%interior_diagnostic_variables(jn)%name, model%get_interior_diagnostic_data(jn))
!      END DO

      ! write 2D diagnostics in the file
      ! ---------------------------------------
!      DO jn = 1, size(model%horizontal_diagnostic_variables)
!         IF (model%horizontal_diagnostic_variables(jn)%save) &
!             CALL iom_put( model%horizontal_diagnostic_variables(jn)%name, model%get_horizontal_diagnostic_data(jn))
!      END DO
      !

      IF ( kt > nittrc000 ) CALL trc_sms_fabm_check_mass
   END SUBROUTINE wri_fabm

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   INTERFACE trc_wri_fabm
       MODULE PROCEDURE wri_fabm,wri_fabm_fl
   END INTERFACE trc_wri_fabm

   PUBLIC trc_wri_fabm

   CONTAINS

   SUBROUTINE wri_fabm_fl (kt, Kmm, fl)
      INTEGER, INTENT( in )               :: fl
      INTEGER, INTENT( in )               :: kt
      INTEGER, INTENT( in )               :: Kmm       
   END SUBROUTINE wri_fabm_fl

   SUBROUTINE wri_fabm (kt, Kmm)          ! Empty routine  
      INTEGER, INTENT( in )               :: kt
      INTEGER, INTENT( in )               :: Kmm       
   END SUBROUTINE wri_fabm
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcwri_fabm.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (see ./LICENSE)
   !!======================================================================
END MODULE trcwri_fabm
