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
    use med_internalstate_mod , only : mastertask, logunit
    use med_internalstate_mod , only : compmed, compatm, compocn, compice, comprof, compwav, ncomps
    use med_internalstate_mod , only : compdat
    use med_internalstate_mod , only : mapbilnr, mapconsf, mapconsd, mappatch
    use med_internalstate_mod , only : mapfcopy, mapnstod, mapnstod_consd, mapnstod_consf
    use med_internalstate_mod , only : mapconsf_aofrac, mapfillv_bilnr, mapbilnr_nstod
    use med_internalstate_mod , only : coupling_mode, mapnames
    use esmFlds               , only : med_fldList_type
    use esmFlds               , only : addfld => med_fldList_AddFld
    use esmFlds               , only : addmap => med_fldList_AddMap
    use esmFlds               , only : addmrg => med_fldList_AddMrg
    use esmflds               , only : fldListTo, fldListFr, fldListMed_aoflux, fldListMed_ocnalb
    use med_internalstate_mod , only : InternalState, mastertask, logunit

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
    character(len=CS), allocatable :: flds(:), oflds(:), aflds(:), dflds(:)
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
      maptype_dat = mapbilnr
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
          call addfld(fldListFr(n)%flds, trim(cvalue))
          call addfld(fldListTo(n)%flds, trim(cvalue))
       end do
    end if

    !=====================================================================
    ! Mediator fields
    !=====================================================================

    ! masks from components
    if (phase == 'advertise') then
       !if (is_local%wrap%comp_present(compice)) call addfld(fldListFr(compice)%flds, 'Si_imask')
       if (is_local%wrap%comp_present(compocn)) call addfld(fldListFr(compocn)%flds, 'So_omask')
!BL
!       if (is_local%wrap%comp_present(compocn)) call addfld(fldListTo(compatm)%flds, 'So_ofrac')
!BL
    else
       if ( fldchk(is_local%wrap%FBexp(compice)        , trim(fldname), rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), trim(fldname), rc=rc)) then
          call addmap(fldListFr(compocn)%flds, 'So_omask', compice,  mapfcopy, 'unset', 'unset')
       end if
    end if

    !=====================================================================
    ! FIELDS TO ATMOSPHERE (compatm)
    !=====================================================================

    ! to atm: fractions (computed in med_phases_prep_atm)
    if (phase == 'advertise') then
       ! ofrac used by atm
       if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compatm)) then
          call addfld(fldListFr(compatm)%flds, 'Sa_ofrac')
       end if
    end if

    ! to atm: unmerged from ice
    ! - zonal surface stress, meridional surface stress
    ! - surface latent heat flux,
    ! - surface sensible heat flux
    ! - surface upward longwave heat flux
    ! - evaporation water flux from water
    ! - mean ice volume per unit area
    ! - mean snow volume per unit area
    ! - surface temperatures

    ! to atm: unmerged surface temperatures from ocn
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compatm)) then
          call addfld(fldListFr(compocn)%flds, 'So_t')
          call addfld(fldListTo(compatm)%flds, 'So_t')
       end if
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'So_t', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compocn,compocn), 'So_t', rc=rc)) then
          !call addmap(fldListFr(compocn)%flds, 'So_t', compatm, maptype, 'one', 'unset')
          call addmap(fldListFr(compocn)%flds, 'So_t', compatm, maptype, 'none', 'unset')
          call addmrg(fldListTo(compatm)%flds, 'So_t', mrg_from=compocn, mrg_fld='So_t', mrg_type='copy')
       end if
    end if

    ! to atm: unmerged from mediator, merge will be done under FV3/CCPP composite step
    ! - zonal surface stress, meridional surface stress
    ! - surface latent heat flux,
    ! - surface sensible heat flux
    ! - surface upward longwave heat flux
    ! - evaporation water flux from water, not in the list do we need to send it to atm?

    ! to atm: surface roughness length from wav
    if (phase == 'advertise') then
       if (is_local%wrap%comp_present(compwav) .and. is_local%wrap%comp_present(compatm)) then
          call addfld(fldListFr(compwav)%flds, 'Sw_z0')
          call addfld(fldListTo(compatm)%flds, 'Sw_z0')
       end if
    else
       if ( fldchk(is_local%wrap%FBexp(compatm)        , 'Sw_z0', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compwav,compwav), 'Sw_z0', rc=rc)) then
          call addmap(fldListFr(compwav)%flds, 'Sw_z0', compatm, maptype, 'one', 'unset')
          call addmrg(fldListTo(compatm)%flds, 'Sw_z0', mrg_from=compwav, mrg_fld='Sw_z0', mrg_type='copy')
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

    ! to ocn: sea level pressure from atm
! Sa_pslv may be used to compute pressure gradient
!    if (phase == 'advertise') then
!       if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn)) then
!          call addfld(fldListFr(compatm)%flds, 'Sa_pslv')
!          call addfld(fldListFr(compdat)%flds, 'Sd_pslv')
!          call addfld(fldListTo(compocn)%flds, 'Sa_pslv')
!       end if
!    else
!       if ( fldchk(is_local%wrap%FBexp(compocn)        , 'Sa_pslv', rc=rc) .and. &
!            fldchk(is_local%wrap%FBImp(compdat,compdat), 'Sd_pslv', rc=rc) .and. &
!            fldchk(is_local%wrap%FBImp(compatm,compatm), 'Sa_pslv', rc=rc)) then
!          call addmap(fldListFr(compdat)%flds, 'Sd_pslv', compocn, maptype_dat, 'one', 'unset')
!          call addmap(fldListFr(compatm)%flds, 'Sa_pslv', compocn, maptype, 'one', 'unset')
!       end if
!    end if

    ! to ocn: from sw from atm and sw from cdeps (custom merge in med_phases_prep_ocn)
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
             call addfld(fldListFr(compatm)%flds, trim(aflds(n)))
             call addfld(fldListTo(compocn)%flds, trim(oflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(oflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(aflds(n)), rc=rc)) then
             !call addmap(fldListFr(compatm)%flds, trim(aflds(n)), compocn, maptype, 'one', 'unset')
             call addmap(fldListFr(compatm)%flds, trim(aflds(n)), compocn, maptype, 'none', 'unset')
          end if
       end if
    end do

    do n = 1,size(oflds)
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compdat) .and. is_local%wrap%comp_present(compocn)) then
             call addfld(fldListFr(compdat)%flds, trim(dflds(n)))
             call addfld(fldListTo(compocn)%flds, trim(oflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(oflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compdat,compdat), trim(dflds(n)), rc=rc)) then
             !call addmap(fldListFr(compdat)%flds, trim(dflds(n)), compocn, maptype_dat, 'one', 'unset')
             call addmap(fldListFr(compdat)%flds, trim(dflds(n)), compocn, maptype_dat, 'none', 'unset')
          end if
       end if
    end do
    deallocate(oflds)
    deallocate(aflds)
    deallocate(dflds)

    ! to ocn: rain and snow via auto merge
    allocate(oflds(4))
    allocate(dflds(4))
    allocate(aflds(4))
    aflds = (/'Faxa_rain', 'Faxa_lwnet', 'Faxa_sen', 'Faxa_lat'/)
    dflds = (/'Faxd_rain', 'Faxd_lwnet', 'Faxd_sen', 'Faxd_lat'/)
    oflds = (/'Foxx_rain', 'Foxx_lwnet', 'Foxx_sen', 'Faxa_evap'/)
    do n = 1,size(oflds)
!       fldname = trim(oflds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn) &
             .and. is_local%wrap%comp_present(compdat)) then
             call addfld(fldListFr(compatm)%flds, trim(aflds(n)))
             call addfld(fldListFr(compdat)%flds, trim(dflds(n)))
             call addfld(fldListTo(compocn)%flds, trim(oflds(n)))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(oflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compdat,compdat), trim(dflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(aflds(n)), rc=rc)) then
             !call addmap(fldListFr(compdat)%flds, trim(dflds(n)), compocn, maptype_dat, 'one', 'unset')
             !call addmap(fldListFr(compatm)%flds, trim(aflds(n)), compocn, maptype, 'one', 'unset')
             call addmap(fldListFr(compdat)%flds, trim(dflds(n)), compocn, maptype_dat, 'none', 'unset')
             call addmap(fldListFr(compatm)%flds, trim(aflds(n)), compocn, maptype, 'none', 'unset')
          end if
       end if
    end do
    deallocate(oflds)
    deallocate(dflds)
    deallocate(aflds)

    if (trim(coupling_mode) == 'hafs_mom6' ) then
       ! to ocn: merge surface stress (custom merge calculation in med_phases_prep_ocn)
       allocate(aflds(2))
       allocate(dflds(2))
       allocate(oflds(2))
       aflds = (/'Faxa_taux', 'Faxa_tauy'/)
       dflds = (/'Faxd_taux', 'Faxd_tauy'/)
       oflds = (/'Foxx_taux', 'Foxx_tauy'/)
       do n = 1,size(oflds)
          if (phase == 'advertise') then
             if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compocn) &
             .and. is_local%wrap%comp_present(compdat)) then
                   call addfld(fldListFr(compatm)%flds, trim(aflds(n)))
                   call addfld(fldListFr(compdat)%flds, trim(dflds(n)))
                   call addfld(fldListTo(compocn)%flds, trim(oflds(n)))
             end if
          else
             if ( fldchk(is_local%wrap%FBexp(compocn)     , trim(oflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compdat,compdat), trim(dflds(n)), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(aflds(n)), rc=rc)) then
                !call addmap(fldListFr(compdat)%flds, trim(dflds(n)), compocn, maptype_dat, 'one', 'unset')
                !call addmap(fldListFr(compatm)%flds, trim(aflds(n)), compocn, maptype, 'one', 'unset')
                call addmap(fldListFr(compdat)%flds, trim(dflds(n)), compocn, maptype_dat, 'none', 'unset')
                call addmap(fldListFr(compatm)%flds, trim(aflds(n)), compocn, maptype, 'none', 'unset')
             end if
          end if
       end do
       deallocate(aflds)
       deallocate(dflds)
       deallocate(oflds)

    endif

    ! to ocn: salt flux from ice
    allocate(flds(3))
    flds = (/'Fioi_meltw', 'Fioi_melth', 'Fioi_salt '/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compocn)) then
             call addfld(fldListFr(compice)%flds, trim(fldname))
             call addfld(fldListTo(compocn)%flds, trim(fldname))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice,compice), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compice)%flds, trim(fldname), compocn,  mapfcopy, 'unset', 'unset')
             call addmrg(fldListTo(compocn)%flds, trim(fldname), &
                  mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy_with_weights', mrg_fracname='ifrac')
          end if
       end if
    end do
    deallocate(flds)

    ! to ocn: partitioned stokes drift from wav
    allocate(flds(6))
    flds = (/'Sw_ustokes1', 'Sw_ustokes2', 'Sw_ustokes3', &
             'Sw_vstokes1', 'Sw_vstokes2', 'Sw_vstokes3'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compwav) .and. is_local%wrap%comp_present(compocn)) then
             call addfld(fldListFr(compwav)%flds, trim(fldname))
             call addfld(fldListTo(compocn)%flds, trim(fldname))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compocn)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compwav,compwav), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compwav)%flds, trim(fldname), compocn, mapbilnr_nstod, 'one', 'unset')
             call addmrg(fldListTo(compocn)%flds, trim(fldname), mrg_from=compwav, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
    end do
    deallocate(flds)

    !=====================================================================
    ! FIELDS TO WAV (compwav)
    !=====================================================================

    ! to wav - 10m winds and bottom temperature from atm
    allocate(flds(3))
    flds = (/'Sa_u10m', 'Sa_v10m', 'Sa_tbot'/)
    do n = 1,size(flds)
       fldname = trim(flds(n))
       if (phase == 'advertise') then
          if (is_local%wrap%comp_present(compatm) .and. is_local%wrap%comp_present(compwav)) then
             call addfld(fldListFr(compatm)%flds, trim(fldname))
             call addfld(fldListTo(compwav)%flds, trim(fldname))
          end if
       else
          if ( fldchk(is_local%wrap%FBexp(compwav)        , trim(fldname), rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compatm), trim(fldname), rc=rc)) then
             call addmap(fldListFr(compatm)%flds, trim(fldname), compwav, mapnstod_consf, 'one', 'unset')
             call addmrg(fldListTo(compwav)%flds, trim(fldname), mrg_from=compatm, mrg_fld=trim(fldname), mrg_type='copy')
          end if
       end if
     end do
     deallocate(flds)

     ! to wav: sea ice fraction, thickness and floe diameter
!     allocate(flds(3))
!     flds = (/'Si_ifrac   ', 'Si_floediam', 'Si_thick   '/)
!     do n = 1,size(flds)
!        fldname = trim(flds(n))
!        if (phase == 'advertise') then
!           if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compwav)) then
!              call addfld(fldListFr(compice)%flds, trim(fldname))
!              call addfld(fldListTo(compwav)%flds, trim(fldname))
!           end if
!        else
!           if ( fldchk(is_local%wrap%FBexp(compwav)        , trim(fldname), rc=rc) .and. &
!                fldchk(is_local%wrap%FBImp(compice,compice), trim(fldname), rc=rc)) then
!               call addmap(fldListFr(compice)%flds, trim(fldname), compwav, mapbilnr_nstod, 'one', 'unset')
!               call addmrg(fldListTo(compwav)%flds, trim(fldname), mrg_from=compice, mrg_fld=trim(fldname), mrg_type='copy')
!           end if
!        end if
!      end do
!      deallocate(flds)

     ! to wav: zonal sea water velocity from ocn
     ! to wav: meridional sea water velocity from ocn
     ! to wav: surface temperature from ocn
     allocate(flds(3))
     flds = (/'So_u', 'So_v', 'So_t'/)
     do n = 1,size(flds)
        fldname = trim(flds(n))
        if (phase == 'advertise') then
           if (is_local%wrap%comp_present(compocn) .and. is_local%wrap%comp_present(compwav)) then
              call addfld(fldListFr(compocn)%flds, trim(fldname))
              call addfld(fldListTo(compwav)%flds, trim(fldname))
           end if
        else
           if ( fldchk(is_local%wrap%FBexp(compwav)        , trim(fldname), rc=rc) .and. &
                fldchk(is_local%wrap%FBImp(compocn,compocn), trim(fldname), rc=rc)) then
              call addmap(fldListFr(compocn)%flds, trim(fldname), compwav, mapbilnr_nstod , 'one', 'unset')
              call addmrg(fldListTo(compwav)%flds, trim(fldname), mrg_from=compocn, mrg_fld=trim(fldname), mrg_type='copy')
           end if
        end if
     end do
     deallocate(flds)

  end subroutine esmFldsExchange_hafs_mom6

end module esmFldsExchange_hafs_mom6_mod
