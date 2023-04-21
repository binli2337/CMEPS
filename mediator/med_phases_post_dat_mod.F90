module med_phases_post_dat_mod

  !-----------------------------------------------------------------------------
  ! Mediator phase for post dat calculations, maps dat->ice, dat->lnd, dat->ocn
  ! and dat->wav
  !-----------------------------------------------------------------------------

  implicit none
  private

  public :: med_phases_post_dat

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_post_dat(gcomp, rc)

    !---------------------------------------
    ! map dat to ocn and dat to ice and dat to land
    !---------------------------------------

    use NUOPC_Mediator        , only : NUOPC_MediatorGet
    use ESMF                  , only : ESMF_Clock, ESMF_ClockIsCreated
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_FieldBundleGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_internalstate_mod , only : InternalState, maintask, logunit
    use med_phases_history_mod, only : med_phases_history_write_comp
    use med_map_mod           , only : med_map_field_packed
    use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
    use med_utils_mod         , only : chkerr    => med_utils_ChkErr
    use med_internalstate_mod , only : compocn, compdat, compatm, compice, complnd, compwav
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Clock)    :: dClock
    character(len=*), parameter :: subname='(med_phases_post_dat)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    ! Get the internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! map dat to ocn
    if (is_local%wrap%med_coupling_active(compdat,compocn)) then
       call t_startf('MED:'//trim(subname)//' map_dat2ocn')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compdat,compdat), &
            FBDst=is_local%wrap%FBImp(compdat,compocn), &
            FBFracSrc=is_local%wrap%FBFrac(compdat), &
            field_normOne=is_local%wrap%field_normOne(compdat,compocn,:), &
            packed_data=is_local%wrap%packed_data(compdat,compocn,:), &
            routehandles=is_local%wrap%RH(compdat,compocn,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_dat2ocn')
    end if
    ! map dat->ice
    if (is_local%wrap%med_coupling_active(compdat,compice)) then
       call t_startf('MED:'//trim(subname)//' map_dat2ice')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compdat,compdat), &
            FBDst=is_local%wrap%FBImp(compdat,compice), &
            FBFracSrc=is_local%wrap%FBFrac(compdat), &
            field_normOne=is_local%wrap%field_normOne(compdat,compice,:), &
            packed_data=is_local%wrap%packed_data(compdat,compice,:), &
            routehandles=is_local%wrap%RH(compdat,compice,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_dat2ice')
    end if
    ! map dat->wav
    if (is_local%wrap%med_coupling_active(compdat,compwav)) then
       call t_startf('MED:'//trim(subname)//' map_dat2wav')
       call med_map_field_packed( &
            FBSrc=is_local%wrap%FBImp(compdat,compdat), &
            FBDst=is_local%wrap%FBImp(compdat,compwav), &
            FBFracSrc=is_local%wrap%FBFrac(compdat), &
            field_normOne=is_local%wrap%field_normOne(compdat,compwav,:), &
            packed_data=is_local%wrap%packed_data(compdat,compwav,:), &
            routehandles=is_local%wrap%RH(compdat,compwav,:), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf('MED:'//trim(subname)//' map_dat2wav')
    end if

    ! Write dat inst, avg or aux if requested in mediator attributes
    call NUOPC_MediatorGet(gcomp, driverClock=dClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_ClockIsCreated(dclock)) then
       call med_phases_history_write_comp(gcomp, compdat, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_post_dat

end module med_phases_post_dat_mod
