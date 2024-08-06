MODULE trcini_fabm
   !!======================================================================
   !!                         ***  MODULE trcini_fabm  ***
   !! TOP :   initialisation of the FABM tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                               FABM tracers
   !!----------------------------------------------------------------------
   !! trc_ini_fabm   : FABM model initialisation
   !!----------------------------------------------------------------------
   USE in_out_manager ! I/O manager
   USE iom            ! I/O module
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE par_fabm
   USE st2D_fabm
   USE trcsms_fabm
   USE fabm, ONLY: fabm_create_model, type_fabm_variable
   USE fabm_driver
   USE inputs_fabm, ONLY: initialize_inputs,link_inputs, &
     type_input_variable,type_input_data,type_river_data, &
     first_input_data,first_river_data
     
#if defined key_git_version
   USE fabm_version,ONLY: fabm_commit_id=>git_commit_id, &
                          fabm_branch_name=>git_branch_name
   USE fabm_types,ONLY: type_version,first_module_version
#endif


   IMPLICIT NONE
   PRIVATE

#if defined key_git_version
#include "gitversion.h90"
   CHARACTER(len=*),parameter :: git_commit_id = _NEMO_COMMIT_ID_
   CHARACTER(len=*),parameter :: git_branch_name = _NEMO_BRANCH_
#endif

   PUBLIC   trc_ini_fabm   ! called by trcini.F90 module
   PUBLIC   nemo_fabm_configure

   TYPE,extends(type_base_driver) :: type_nemo_fabm_driver
   contains
      procedure :: fatal_error => nemo_fabm_driver_fatal_error
      procedure :: log_message => nemo_fabm_driver_log_message
   end type

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL licence     (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE nemo_fabm_configure()
      INTEGER :: jn
      INTEGER, PARAMETER :: xml_unit = 1979

! TODO  rem.inp
      TYPE (type_input_data),POINTER :: input_data
      TYPE (type_river_data),POINTER :: river_data
      CLASS (type_input_variable),POINTER :: input_pointer

      ! Allow FABM to parse fabm.yaml. This ensures numbers of variables are known.
      model => fabm_create_model()

      jp_fabm = size(model%interior_state_variables)
      jp_fabm_bottom = size(model%bottom_state_variables)
      jp_fabm_surface = size(model%surface_state_variables)
!      jp_fabm0 = jptra + 1
!      jp_fabm1 = jptra + jp_fabm
      jp_fabm_m1=0
!      jptra = jptra + jp_fabm
      jp_dia2d = jp_dia2d + size(model%horizontal_diagnostic_variables)
      jp_dia3d = jp_dia3d + size(model%interior_diagnostic_variables)
      jp_diabio = jp_diabio + jp_fabm

! TODO  remove inputs
      ! Read inputs (river and additional 2D forcing) from fabm_input.nml
      ! This must be done before writing field_def_fabm.xml, as that file
      ! also describes the additional input variables.
      call initialize_inputs

      IF (lwp) THEN
         ! write field_def_fabm.xml on lead process,  TODO i  not lon lead pocess then ... call sleep(for a second or two)  ?
         OPEN(UNIT=xml_unit,FILE='field_def_fabm.xml',ACTION='WRITE',STATUS='REPLACE')

         WRITE (xml_unit,1000) '<field_definition level="1" prec="4" operation="average" enabled=".TRUE." default_value="1.e20" >'

         WRITE (xml_unit,1000) ' <field_group id="ptrc_T" grid_ref="grid_T_3D">'
         DO jn=1,jp_fabm
            CALL write_variable_xml(xml_unit,model%interior_state_variables(jn))
#if defined key_trdtrc
            CALL write_trends_xml(xml_unit,model%interior_state_variables(jn))
#endif
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) ' <field_group id="sf_T" grid_ref="grid_T_2D">'
         DO jn=1,jp_fabm_surface
            CALL write_variable_xml(xml_unit,model%surface_state_variables(jn))
         END DO
         DO jn=1,jp_fabm_bottom
            CALL write_variable_xml(xml_unit,model%bottom_state_variables(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) ' <field_group id="diad_T" grid_ref="grid_T_2D">'
         DO jn=1,size(model%interior_diagnostic_variables)
            CALL write_variable_xml(xml_unit,model%interior_diagnostic_variables(jn),3)
         END DO
         DO jn=1,size(model%horizontal_diagnostic_variables)
            CALL write_variable_xml(xml_unit,model%horizontal_diagnostic_variables(jn),2)
         END DO
         DO jn=1,size(model%interior_state_variables)
            WRITE (xml_unit,'(A)') '  <field id="'//TRIM(model%interior_state_variables(jn)%name)//'_VINT" long_name="depth-integrated ' &
                    //TRIM(model%interior_state_variables(jn)%long_name)//'" unit="'//TRIM(model%interior_state_variables(jn)%units)//'*m" default_value="0.0" />'
         END DO
         DO jn=1,size(model%interior_diagnostic_variables)
            WRITE (xml_unit,'(A)') '  <field id="'//TRIM(model%interior_diagnostic_variables(jn)%name)//'_VINT" long_name="depth-integrated &
                    '//TRIM(model%interior_diagnostic_variables(jn)%long_name)//'" unit="'//TRIM(model%interior_diagnostic_variables(jn)%units)//'*m" default_value="0.0" />'
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) ' <field_group id="fabm_scalar" grid_ref="grid_0">'

         DO jn=1,size(model%conserved_quantities)
            CALL write_variable_xml(xml_unit,model%conserved_quantities(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'
!-----     Mokrane
         WRITE (xml_unit,1000) ' <field_group id="fabm_input" grid_ref="grid_T_2D">'
          input_data => first_input_data
          DO WHILE (ASSOCIATED(input_data))
            input_pointer => input_data
            CALL write_input_xml(xml_unit,input_pointer)
             input_data => input_data%next
          END DO
         river_data => first_river_data
         DO WHILE (ASSOCIATED(river_data))
           input_pointer => river_data
           CALL write_input_xml(xml_unit,input_pointer,3)
            river_data => river_data%next
          END DO
          WRITE (xml_unit,1000) ' </field_group>'
!------- Mokrane
         WRITE (xml_unit,1000) '</field_definition>'

         CLOSE(xml_unit)
      END IF
      IF( lk_mpp )   CALL mppsync !Ensure field_def_fabm is ready.

1000 FORMAT (A)

   END SUBROUTINE nemo_fabm_configure

   SUBROUTINE write_variable_xml(xml_unit,variable,flag_grid_ref)
      INTEGER,INTENT(IN) :: xml_unit
      INTEGER,INTENT(IN),OPTIONAL :: flag_grid_ref
      CLASS (type_fabm_variable),INTENT(IN) :: variable

      CHARACTER(LEN=20) :: missing_value,string_dimensions
      INTEGER :: number_dimensions

      ! Check variable dimension for grid_ref specificaiton.
      ! Default is to not specify the grid_ref in the field definition.
      IF (present(flag_grid_ref)) THEN
          number_dimensions=flag_grid_ref
      ELSE
          number_dimensions=-1 !default, don't specify grid_ref
      ENDIF

      WRITE (missing_value,'(E10.3)') variable%missing_value
      WRITE (string_dimensions,'(I1)') number_dimensions
      SELECT CASE (number_dimensions)
      CASE (3)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
      CASE (2)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_2D"/>'
      CASE (0)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="1point"/>'
      CASE (-1)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" />'
      CASE default
         IF(lwp) WRITE(numout,*) ' trc_ini_fabm: Failing to initialise output of variable '//TRIM(variable%name)//': Output of '//TRIM(ADJUSTL(string_dimensions))//'-dimensional variables not supported!!!'
      END SELECT

   END SUBROUTINE write_variable_xml

   SUBROUTINE write_trends_xml(xml_unit,variable,flag_grid_ref)
      INTEGER,INTENT(IN) :: xml_unit
      INTEGER,INTENT(IN),OPTIONAL :: flag_grid_ref
      CLASS (type_fabm_variable),INTENT(IN) :: variable

      INTEGER :: number_dimensions,i
      CHARACTER(LEN=20) :: missing_value,string_dimensions
#if defined key_tracer_budget
      CHARACTER(LEN=3),DIMENSION(10),PARAMETER :: trd_tags = (/ &
        'LDF','BBL','FOR','ZDF','DMP','SMS','ATF', &
        'RDB','RDN','VMV' /)
      CHARACTER(LEN=3),DIMENSION(3),PARAMETER :: trd_e3t_tags = (/ &
        'XAD','YAD','ZAD' /)
#else
      CHARACTER(LEN=3),DIMENSION(13),PARAMETER :: trd_tags = (/ &
        'XAD','YAD','ZAD','LDF','BBL','FOR','ZDF','DMP','SMS','ATF', &
        'RDB','RDN','VMV' /)
#endif

      ! Check variable dimension for grid_ref specificaiton.
      ! Default is to not specify the grid_ref in the field definition.
      IF (present(flag_grid_ref)) THEN
          number_dimensions=flag_grid_ref
      ELSE
          number_dimensions=-1 !default, don't specify grid_ref
      ENDIF

      WRITE (missing_value,'(E10.3)') -2.E20
      WRITE (string_dimensions,'(I1)') number_dimensions
      SELECT CASE (number_dimensions)
      CASE (3)
        DO i=1,size(trd_tags)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(trd_tags(i))//'_'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//' '//TRIM(trd_tags(i))//' trend" &
                 unit="'//TRIM(variable%units)//'/s" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
        END DO
#if defined key_tracer_budget
        DO i=1,size(trd_e3t_tags)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(trd_e3t_tags(i))//'_'//TRIM(variable%name)//'_e3t" long_name="'//TRIM(variable%long_name)//' cell depth integrated '//TRIM(trd_e3t_tags(i))//' trend" &
                 unit="'//TRIM(variable%units)//'*m/s" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
        END DO
        WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'_e3t" long_name="'//TRIM(variable%long_name)//' cell depth integrated" &
                unit="'//TRIM(variable%units)//'*m" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
#endif
      CASE (-1)
        DO i=1,size(trd_tags)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(trd_tags(i))//'_'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//' '//TRIM(trd_tags(i))//' trend" &
                 unit="'//TRIM(variable%units)//'/s" default_value="'//TRIM(ADJUSTL(missing_value))//'" />'
        END DO
#if defined key_tracer_budget
        DO i=1,size(trd_e3t_tags)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(trd_e3t_tags(i))//'_'//TRIM(variable%name)//'_e3t" long_name="'//TRIM(variable%long_name)//' cell depth integrated '//TRIM(trd_e3t_tags(i))//' trend" &
                 unit="'//TRIM(variable%units)//'*m/s" default_value="'//TRIM(ADJUSTL(missing_value))//'" />'
        END DO
        WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'_e3t" long_name="'//TRIM(variable%long_name)//' cell depth integrated" unit="'//TRIM(variable%units)//'*m" default_value="'//TRIM(ADJUSTL(missing_value))//'" />'
#endif
      CASE default
         IF(lwp) WRITE(numout,*) ' trc_ini_fabm: Failing to initialise trends of variable '//TRIM(variable%name)//': Output of '//TRIM(ADJUSTL(string_dimensions))//'-dimensional trends not supported!!!'
      END SELECT

   END SUBROUTINE write_trends_xml

!--------- Mokrane --------

   SUBROUTINE write_input_xml(xml_unit,variable,flag_grid_ref)
      INTEGER,INTENT(IN) :: xml_unit
      INTEGER,INTENT(IN),OPTIONAL :: flag_grid_ref
      CLASS(type_input_variable),POINTER,INTENT(IN) :: variable

      INTEGER :: number_dimensions,i
      CHARACTER(LEN=20) :: missing_value,string_dimensions

      ! Check variable dimension for grid_ref specificaiton.
      ! Default is to not specify the grid_ref in the field definition.
      IF (present(flag_grid_ref)) THEN
          number_dimensions=flag_grid_ref
      ELSE
          number_dimensions=-1 !default, don't specify grid_ref
      ENDIF
 
       WRITE (missing_value,'(E10.3)') -2.E20
       WRITE (string_dimensions,'(I1)') number_dimensions
       SELECT CASE (number_dimensions)
       CASE (3)
         WRITE (xml_unit,'(A)') '  <field id="'//'INP_'//TRIM(variable%sf(1)%clvar)//'" long_name="'//TRIM(variable%sf(1)%clvar)//' input" unit="" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
       CASE (-1)
         WRITE (xml_unit,'(A)') '  <field id="'//'INP_'//TRIM(variable%sf(1)%clvar)//'" long_name="'//TRIM(variable%sf(1)%clvar)//' input" unit="" default_value="'//TRIM(ADJUSTL(missing_value))//'" />'
       CASE default
          IF(lwp) WRITE(numout,*) ' trc_ini_fabm: Failing to initialise input diagnostic of variable '//TRIM(variable%sf(1)%clvar)//': Output of '//TRIM(ADJUSTL(string_dimensions))//'-dimensional diagnostic not supported!!!'
       END SELECT
 
   END SUBROUTINE write_input_xml

!-----------Mokrane--------------

   SUBROUTINE trc_ini_fabm( Kmm )
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_fabm  ***
      !!
      !! ** Purpose :   initialization for FABM model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)  ::  Kmm  ! time level indices
#if defined key_git_version
      TYPE (type_version),POINTER :: version
#endif
      INTEGER :: jn
      INTEGER :: inum      ! local integers
      REAL(wp), DIMENSION(jpi,jpj) ::  fld_inp_2d   ! 2D workspace
      CHARACTER(len=40) ::   cl_file_sed, cl_file_surf  ! local names

      cl_file_sed    = 'fabm_input_sed.nc'
      cl_file_surf   = 'fabm_input_suf.nc'

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_fabm: initialisation of FABM model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
#if defined key_git_version
      IF(lwp) WRITE(numout,*) ' NEMO version:   ',git_commit_id,' (',git_branch_name,' branch)'
      IF(lwp) WRITE(numout,*) ' FABM version:   ',fabm_commit_id,' (',fabm_branch_name,' branch)'
      version => first_module_version
      do while (associated(version))
         IF(lwp) WRITE(numout,*)  ' '//trim(version%module_name)//' version:   ',trim(version%version_string)
         version => version%next
      end do
#endif

      ! Allocate FABM arrays
      ! Constant initial values from the fabm namelist are set to the arrays here
      IF(trc_sms_fabm_alloc(Kmm) /= 0) CALL ctl_stop( 'STOP', 'trc_ini_fabm: unable to allocate FABM arrays' )

      ! Log mapping of FABM states:
      IF (lwp) THEN
         IF (jp_fabm.gt.0) WRITE(numout,*) " FABM tracers:"
         DO jn=1,jp_fabm
            WRITE(numout,*) "   State",jn,":",trim(model%interior_state_variables(jn)%name), &
               " (",trim(model%interior_state_variables(jn)%long_name), &
               ") [",trim(model%interior_state_variables(jn)%units),"]"
         ENDDO
         IF (jp_fabm_surface.gt.0) WRITE(numout,*) "FABM seasurface states:"
         DO jn=1,jp_fabm_surface
            WRITE(numout,*) "   State",jn,":",trim(model%surface_state_variables(jn)%name), &
               " (",trim(model%surface_state_variables(jn)%long_name), &
               ") [",trim(model%surface_state_variables(jn)%units),"]"
         ENDDO
         IF (jp_fabm_bottom.gt.0) WRITE(numout,*) "FABM seafloor states:"
         DO jn=1,jp_fabm_bottom
            WRITE(numout,*) "   State",jn,":",trim(model%bottom_state_variables(jn)%name), &
               " (",trim(model%bottom_state_variables(jn)%long_name), &
               ") [",trim(model%bottom_state_variables(jn)%units),"]"
         ENDDO
      ENDIF

    IF( .NOT. ln_rsttr ) tr(:,:,:,jp_fabm0:jp_fabm1,Kmm) = 1.     
    DO jn=1,jp_fabm
        tr(:,:,:,jn,Kmm)=model%interior_state_variables(jn)%initial_value * tmask
    ENDDO

! initial conditions for the bottom and surface state variables
   IF(ln_trcdta ) THEN
       IF(lwp) WRITE(numout,*) ' initialization of surface and bottom from files', TRIM(cl_file_surf),TRIM(cl_file_sed),' nsed,nbot',jp_fabm_surface,jp_fabm_bottom
       IF (jp_fabm_surface .gt. 0) THEN
          DO jn=1,jp_fabm_surface
           CALL iom_open ( TRIM(cl_file_surf), inum )
           CALL iom_get  ( inum, jpdom_global, TRIM(model%surface_state_variables(jn)%name), fld_inp_2d, 1 )
           CALL iom_close( inum )
           fabm_st2Dn(:,:,jn) = fld_inp_2d * tmask(:,:,1)
          END DO
       END IF

       IF (jp_fabm_bottom .gt. 0) THEN
          DO jn=1,jp_fabm_bottom
            IF(lwp) WRITE(numout,*) ' reading   from file:',  TRIM(cl_file_surf), '   parameter: ',   TRIM(model%bottom_state_variables(jn)%name)
           CALL iom_open ( TRIM(cl_file_sed), inum )
           CALL iom_get  ( inum, jpdom_global, TRIM(model%bottom_state_variables(jn)%name), fld_inp_2d, 1 )
           CALL iom_close( inum )
           fabm_st2Dn(:,:,jp_fabm_surface+jn) = fld_inp_2d * tmask(:,:,1)
          END DO
       END IF
   ELSE
       IF (jp_fabm_surface .gt. 0) THEN
          DO jn=1,jp_fabm_surface
           fabm_st2Dn(:,:,jn) = model%surface_state_variables(jn)%initial_value * tmask(:,:,1)
          END DO
       END IF
       IF (jp_fabm_bottom .gt. 0) THEN
          DO jn=1,jp_fabm_bottom
            fabm_st2Dn(:,:,jp_fabm_surface+jn) = model%bottom_state_variables(jn)%initial_value * tmask(:,:,1)
          END DO
       END IF
   END IF
   fabm_st2Db = fabm_st2Dn
   END SUBROUTINE trc_ini_fabm

   SUBROUTINE nemo_fabm_driver_fatal_error(self, location, message)
      CLASS (type_nemo_fabm_driver), INTENT(INOUT) :: self
      CHARACTER(len=*),              INTENT(IN)    :: location, message

      CALL ctl_stop('STOP', TRIM(location)//': '//TRIM(message))
      STOP
   END SUBROUTINE

   SUBROUTINE nemo_fabm_driver_log_message(self, message)
      CLASS (type_nemo_fabm_driver), INTENT(INOUT) :: self
      CHARACTER(len=*),              INTENT(IN)    :: message

      IF(lwp) WRITE (numout,*) TRIM(message)
   END SUBROUTINE

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No FABM model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE nemo_fabm_configure
   END SUBROUTINE nemo_fabm_configure

   SUBROUTINE trc_ini_fabm( Kmm )           ! Empty routine
   END SUBROUTINE trc_ini_fabm
#endif

   !!======================================================================
END MODULE trcini_fabm
