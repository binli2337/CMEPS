module esmFldsExchange_hafs_mom6_mod

  !---------------------------------------------------------------------
  ! This is a mediator specific routine that determines ALL possible
  ! fields exchanged between components and their associated routing,
  ! mapping and merging
  !---------------------------------------------------------------------

  implicit none
  public

  public :: esmFldsExchange_hafs_mom6

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine esmFldsExchange_hafs_mom6(gcomp, phase, rc)

    use ESMF
    use NUOPC
    use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
    use med_utils_mod         , only : chkerr => med_utils_chkerr
    use med_methods_mod       , only : fldchk => med_methods_FB_FldChk
    use med_internalstate_mod , only : InternalState
    use med_internalstate_mod , only : maintask, logunit
    use med_internalstate_mod , only : compmed, compatm, compdat, compocn, comprof, compwav, ncomps
    use med_internalstate_mod , only : mapbilnr
    use med_internalstate_mod , only : mapfcopy
    use med_internalstate_mod , only : mapfillv_bilnr
    use med_internalstate_mod , only : coupling_mode, mapnames
    use esmFlds               , only : med_fldList_type
    use esmFlds               , only : addfld_to => med_fldList_addfld_to
    use esmFlds               , only : addmrg_to => med_fldList_addmrg_to
    use esmFlds               , only : addfld_from => med_fldList_addfld_from
    use esmFlds               , only : addmap_from => med_fldList_addmap_from
    use esmFlds               , only : addfld_aoflux => med_fldList_addfld_aoflux
    use esmFlds               , only : addmap_aoflux => med_fldList_addmap_aoflux

    use med_internalstate_mod , only : InternalState, maintask, logunit

    ! input/output parameters:
    type(ESMF_GridComp)              :: gcomp
    character(len=*) , intent(in)    :: phase
    integer          , intent(inout) :: rc

    ! local variables:
    type(InternalState) :: is_local
    integer             :: i, n, maptype, maptype_dat
    character(len=CX)   :: msgString
    character(len=CL)   :: cvalue
    character(len=CS)   :: fldname
    character(len=CS), allocatable :: flds(:), wflds(:), oflds(:), aflds(:), dflds(:)
    character(len=*) , parameter   :: subname='(esmFldsExchange_hafs_mom6)'
    !--------------------------------------

    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Set maptype according to coupling_mode
    if (trim(coupling_mode) == 'hafs_mom6') then
      maptype = mapfillv_bilnr
      maptype_dat = mapfcopy
    end if
    write(msgString,'(A,i6,A)') trim(subname)//': maptype is ',maptype,', '//mapnames(maptype)
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    write(msgString,'(A,i6,A)') trim(subname)//': maptype_dat is ',maptype_dat,', '//mapnames(maptype_dat)
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    !=====================================================================
    ! scalar information
    !=====================================================================

    if (phase == 'advertise') then
       call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,ncomps
          call addfld_from(n, trim(cvalue))
          call addfld_to(n, trim(cvalue))
       end do
    end if

    !=====================================================================
    ! Mediator fields
    !=====================================================================

    ! masks from components
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compocn)) call addfld_from(compocn, 'So_omask')
    else
       if (fldchk(is_local%wrap%FBImp(compocn,compocn),'So_omask',rc=rc)) then
          call ESMF_LogWrite(trim(subname)//": Field connected "//"So_omask", &
          ESMF_LOGMSG_INFO)
       else
          call ESMF_LogSetError(ESMF_FAILURE, &
             msg=trim(subname)//": Field is not connected "//"So_omask", &
             line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    end if

    !=====================================================================
    ! FIELDS TO ATMOSPHERE (compatm)
    !=====================================================================

    ! to atm: unmerged surface temperatures from ocn
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compatm)) then
          call addfld_from(compocn, 'So_t')
          call addfld_to(compatm, 'So_t')
       end if
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'So_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          call addmap_from(compocn, 'So_t', compatm, maptype, 'none', 'unset')
          call addmrg_to(compatm, 'So_t', mrg_from=compocn, mrg_fld='So_t', mrg_type='copy')
       end if
    end if

    ! to atm: surface roughness length from wav
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compwav) .and. is_local%wrap%comp_present(compatm)) then
          call addfld_from(compwav, 'Sw_z0')
          call addfld_to(compatm, 'Sw_z0')
       end if
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Sw_z0', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav,compwav), 'Sw_z0', rc=rc)) then
          call addmap_from(compwav, 'Sw_z0', compatm, maptype, 'none', 'unset')
          call addmrg_to(compatm, 'Sw_z0', mrg_from=compwav, mrg_fld='Sw_z0', mrg_type='copy')
       end if
    end if

    !=====================================================================
    ! FIELDS TO OCEAN (compocn)
    !=====================================================================
    ! to ocn: from sw from atm and sw net from ice (custom merge in med_phases_prep_ocn)
    ! - downward direct  near-infrared ("n" or "i") incident solar radiation
    ! - downward diffuse near-infrared ("n" or "i") incident solar radiation
    ! - downward direct visible ("v") incident solar radiation
    ! - downward diffuse visible ("v") incident solar radiation


    ! to ocn: sw from atm and sw from cdeps (custom merge in med_phases_prep_ocn)
    ! - downward direct  near-infrared ("n" or "i") incident solar radiation
    ! - downward diffuse near-infrared ("n" or "i") incident solar radiation
    ! - downward direct visible ("v") incident solar radiation
    ! - downward diffuse visible ("v") incident solar radiation
    allocate(aflds(4))
    allocate(dflds(4))
    allocate(oflds(4))
    aflds = (/'Faxa_swndr'    , 'Faxa_swndf'    , 'Faxa_swvdr'    , 'Faxa_swvdf'/)
    dflds = (/'Faxd_swndr'    , 'Faxd_swndf'    , 'Faxd_swvdr'    , 'Faxd_swvdf'/)
    oflds = (/'Foxx_swnet_idr', 'Foxx_swnet_idf', 'Foxx_swnet_vdr', 'Foxx_swnet_vdf'/)
    do n = 1,size(oflds)
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
             call addfld_from(compatm, trim(aflds(n)))
             call addfld_to(compocn, trim(oflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(oflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(aflds(n)), rc=rc)) then
             call addmap_from(compatm, trim(aflds(n)), compocn, maptype, 'none', 'unset')
          end if
       end if
    end do

    do n = 1,size(oflds)
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compdat) .and. is_local%wrap%comp_present(compocn)) then
             call addfld_from(compdat, trim(dflds(n)))
             call addfld_to(compocn, trim(oflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(oflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compdat,compdat), trim(dflds(n)), rc=rc)) then
             call addmap_from(compdat, trim(dflds(n)), compocn, maptype_dat, 'none', 'unset')
          end if
       end if
    end do
    deallocate(oflds)
    deallocate(aflds)
    deallocate(dflds)

    ! to ocn: rain, net longwave, sensible heat and latent heat via auto merge
    allocate(oflds(4))
    allocate(dflds(4))
    allocate(aflds(4))
    aflds = (/'Faxa_rain', 'Faxa_lwnet', 'Faxa_sen', 'Faxa_lat'/)
    dflds = (/'Faxd_rain', 'Faxd_lwnet', 'Faxd_sen', 'Faxd_lat'/)
    oflds = (/'Faxa_rain', 'Foxx_lwnet', 'Foxx_sen', 'Foxx_evap'/)
    do n = 1,size(oflds)
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn) &
             .and. is_local%wrap%comp_present(compdat)) then
             call addfld_from(compatm, trim(aflds(n)))
             call addfld_from(compdat, trim(dflds(n)))
             call addfld_to(compocn, trim(oflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(oflds(n)), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compdat,compdat), trim(dflds(n)), rc=rc) .and. &
             fldchk(is_local%wrap%FBImp(compatm,compatm), trim(aflds(n)), rc=rc)) then
             call addmap_from(compdat, trim(dflds(n)), compocn, maptype_dat, 'none', 'unset')
             call addmap_from(compatm, trim(aflds(n)), compocn, maptype, 'none', 'unset')
          end if
       end if
    end do
    deallocate(oflds)
    deallocate(dflds)
    deallocate(aflds)

    ! to ocn: merge surface stress and sea level pressure Sa_pslv (custom merge calculation in med_phases_prep_ocn)
    allocate(aflds(3))
    allocate(dflds(3))
    allocate(oflds(3))
    aflds = (/'Faxa_taux', 'Faxa_tauy', 'Sa_pslv'/)
    dflds = (/'Faxd_taux', 'Faxd_tauy', 'Sd_pslv'/)
    oflds = (/'Foxx_taux', 'Foxx_tauy', 'Sa_pslv'/)
    do n = 1,size(oflds)
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn) &
          .and. is_local%wrap%comp_present(compdat)) then
             call addfld_from(compatm, trim(aflds(n)))
             call addfld_from(compdat, trim(dflds(n)))
             call addfld_to(compocn, trim(oflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)     , trim(oflds(n)), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compdat,compdat), trim(dflds(n)), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm,compatm), trim(aflds(n)), rc=rc)) then
             call addmap_from(compdat, trim(dflds(n)), compocn, maptype_dat, 'none', 'unset')
             call addmap_from(compatm, trim(aflds(n)), compocn, maptype, 'none', 'unset')
          end if
       end if
    end do
    deallocate(aflds)
    deallocate(dflds)
    deallocate(oflds)

    !=====================================================================
    ! FIELDS TO WAV (compwav)
    !=====================================================================

    ! to wav - 10m winds from atm and datm (custom merge calculation in
    ! med_phases_prep_wav)
    allocate(aflds(2))
    allocate(dflds(2))
    allocate(oflds(2))
    aflds = (/'Sa_u10m', 'Sa_v10m'/)
    dflds = (/'Sd_u10m', 'Sd_u10m'/)
    wflds = (/'Sa_u10m', 'Sa_u10m'/)
    do n = 1,size(flds)
       !fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compwav) &
          .and. is_local%wrap%comp_present(compdat)) then
             call addfld_from(compatm, trim(aflds(n)))
             call addfld_from(compdat, trim(dflds(n)))
             call addfld_to(compwav, trim(wflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compwav)     , trim(wflds(n)), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compdat,compdat), trim(dflds(n)), rc=rc) .and. &
          fldchk(is_local%wrap%FBImp(compatm,compatm), trim(aflds(n)), rc=rc)) then
             call addmap_from(compdat, trim(dflds(n)), compwav, maptype_dat, 'none', 'unset')
             call addmap_from(compatm, trim(aflds(n)), compwav, maptype, 'none', 'unset')
            ! call addmrg_to(compwav, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(aflds)
    deallocate(dflds)
    deallocate(wflds)

  end subroutine esmFldsExchange_hafs_mom6

end module esmFldsExchange_hafs_mom6_mod
