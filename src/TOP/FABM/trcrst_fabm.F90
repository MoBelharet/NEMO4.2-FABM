MODULE trcrst_fabm
   !!======================================================================
   !!                      ***  MODULE trcrst_fabm  ***
   !! Read and write additional restart fields used by FABM
   !!======================================================================
   !! History :
   !!----------------------------------------------------------------------
#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'   :                                       FABM model
   !!----------------------------------------------------------------------
   !! trc_nam_fabm      : FABM initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE iom

   USE par_fabm
   USE trcsms_fabm
   USE st2D_fabm

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_rst_read_fabm   ! called by trcrst.F90 module
   PUBLIC   trc_rst_wri_fabm    ! called by trcrst.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL licence (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_rst_read_fabm
      INTEGER :: jn
      REAL(wp), POINTER, DIMENSION(:,:) :: pdata !, toto
      REAL(wp), POINTER :: pdata_s
      REAL(wp), POINTER, DIMENSION(:,:,:) :: pdata_i

      DO jn=1,jp_fabm_surface
         CALL iom_get( numrtr, jpdom_auto, 'fabm_st2Db'//TRIM(model%surface_state_variables(jn)%name), fabm_st2Db(:,:,jn) )
      END DO
      DO jn=1,jp_fabm_bottom
         CALL iom_get( numrtr, jpdom_auto, 'fabm_st2Db'//TRIM(model%bottom_state_variables(jn)%name), fabm_st2Db(:,:,jp_fabm_surface+jn) )
      END DO

      DO jn=1,jp_fabm_surface
         CALL iom_get( numrtr, jpdom_auto, 'fabm_st2Dn'//TRIM(model%surface_state_variables(jn)%name), fabm_st2Dn(:,:,jn) )
      END DO
      DO jn=1,jp_fabm_bottom
         CALL iom_get( numrtr, jpdom_auto, 'fabm_st2Dn'//TRIM(model%bottom_state_variables(jn)%name), fabm_st2Dn(:,:,jp_fabm_surface+jn) )
      END DO

      ! Mokrane: horizontal Variables 
      DO jn = 1, size(model%horizontal_diagnostic_variables)
        IF (model%horizontal_diagnostic_variables(jn)%part_of_state) THEN
                pdata => model%get_horizontal_diagnostic_data(jn)
                !write(numout,*) "La valeur de pdata for ", jn, " variable : " , TRIM(model%horizontal_diagnostic_variables(jn)%name), " est : ", pdata
                CALL iom_get( numrtr, jpdom_auto, TRIM(model%horizontal_diagnostic_variables(jn)%name), pdata )
                !toto => model%get_horizontal_diagnostic_data(jn)
                !write(numout,*) "La valeur de pdata dans trc_rst_read_fabm après pointage for ", jn, " variable : " , TRIM(model%horizontal_diagnostic_variables(jn)%name), " est : ", toto
        END IF
      END DO

      ! scalar variables:
      !DO jn = 1, size(model%scalar_diagnostic_variables)
         !IF (model%scalar_diagnostic_variables(jn)%part_of_state) THEN
          !       pdata_s => model%get_scalar_diagnostic_data(jn)
          !       CALL iom_get( numrtr, jpdom_auto, TRIM(model%scalar_diagnostic_variables(jn)%name), pdata_s )
          !       write(numout,*) "La valeur de pdata for ", jn, " variable : " , &
          !               TRIM(model%scalar_diagnostic_variables(jn)%name), " est : ", pdata_s
         !END IF
      !END DO

      ! interior variables
      DO jn = 1, size(model%interior_diagnostic_variables)
        IF (model%interior_diagnostic_variables(jn)%part_of_state) THEN
                pdata_i => model%get_interior_diagnostic_data(jn)
                CALL iom_get( numrtr, jpdom_auto, TRIM(model%interior_diagnostic_variables(jn)%name), pdata_i )
        END IF
      END DO



   END SUBROUTINE trc_rst_read_fabm

   SUBROUTINE trc_rst_wri_fabm(kt)
      INTEGER, INTENT( in ) ::   kt    ! ocean time-step index

      INTEGER :: jn
      REAL(wp), POINTER, DIMENSION(:,:) :: pdata
      REAL(wp), POINTER :: pdata_s
      REAL(wp), POINTER, DIMENSION(:,:,:) :: pdata_i
      REAL(wp):: pdata_local(3,3)

      DO jn=1,jp_fabm_surface
         CALL iom_rstput( kt, nitrst, numrtw, 'fabm_st2Db'//TRIM(model%surface_state_variables(jn)%name), fabm_st2Db(:,:,jn) )
      END DO
      DO jn=1,jp_fabm_bottom
         CALL iom_rstput( kt, nitrst, numrtw, 'fabm_st2Db'//TRIM(model%bottom_state_variables(jn)%name), fabm_st2Db(:,:,jp_fabm_surface+jn) )
      END DO

      DO jn=1,jp_fabm_surface
         CALL iom_rstput( kt, nitrst, numrtw, 'fabm_st2Dn'//TRIM(model%surface_state_variables(jn)%name), fabm_st2Dn(:,:,jn) )
      END DO
      DO jn=1,jp_fabm_bottom
         CALL iom_rstput( kt, nitrst, numrtw, 'fabm_st2Dn'//TRIM(model%bottom_state_variables(jn)%name), fabm_st2Dn(:,:,jp_fabm_surface+jn) )
      END DO

      ! Mokrane: Diagnostic variables

      DO jn = 1, size(model%horizontal_diagnostic_variables)
         IF (model%horizontal_diagnostic_variables(jn)%part_of_state) THEN
                pdata => model%get_horizontal_diagnostic_data(jn)

                CALL iom_rstput( kt, nitrst, numrtw,TRIM(model%horizontal_diagnostic_variables(jn)%name), pdata)
         END IF
      END DO

      ! SCALAR VARIABLES

      DO jn = 1, size(model%scalar_diagnostic_variables)
         IF (model%scalar_diagnostic_variables(jn)%part_of_state) THEN
                pdata_s => model%get_scalar_diagnostic_data(jn)
     !           !pdata_local = pdata
                CALL iom_rstput( kt, nitrst, numrtw,TRIM(model%scalar_diagnostic_variables(jn)%name), pdata_s)
     !           write(numout,*) "WRITE : The value of pdata_s for ", jn, " variable : " , &
     !                   TRIM(model%scalar_diagnostic_variables(jn)%name)," est : ", pdata_s
         END IF
      END DO

      ! INTERIOR VARIABLES
      DO jn = 1, size(model%interior_diagnostic_variables)
        IF (model%interior_diagnostic_variables(jn)%part_of_state) THEN
                pdata_i => model%get_interior_diagnostic_data(jn)
                CALL iom_rstput( kt, nitrst, numrtw,TRIM(model%interior_diagnostic_variables(jn)%name), pdata_i)
        END IF
      END DO

   END SUBROUTINE trc_rst_wri_fabm

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                             No FABM
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rst_read_fabm
   END  SUBROUTINE trc_rst_read_fabm

   SUBROUTINE trc_rst_wri_fabm(kt)
      INTEGER, INTENT( in ) ::   kt    ! ocean time-step index
   END SUBROUTINE trc_rst_wri_fabm
#endif

   !!======================================================================
END MODULE trcrst_fabm
