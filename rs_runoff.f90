module array
  REAL(8), dimension(:), allocatable, save::Cl, Cp, Sl, Sp, San, Sar, Di, Fmub, &
                                             Ftub, Ffm, Free1, Car, Bw, DmC, Dzbp, &
                                             Cave, Q, Zbmin, Z, Ux, Dzbf_b, Dzbf_g, Free2, Free3, Rpr, Reti, gari_length, &
                                             Pro_rachi, I_rachi_ini, I_rachi_tot, Qdep, Ldep_b, Ldep_g, &
                                             Bdep, Briv1, Briv2, Hei_gari, Pslo, B_area, W_depth, S_rate, &
                                             I_gari_ini, I_gari_tot, Bare_r, Width_b, Depth_10, &
                                             Times_ft, Dep_ini_b, Dep_ini_g, Granite_r, Diam_b, Ini_time_ft, Tot_time_ft, Pdep, &
                                             Pro_gari, Qdep_g, Qdep_b, Bdep_g, Bdep_b, Qsup_g, Qsup_b, &
                                             Npool, Vsed, Vmax, Vout, Vsup, Vout_tot, Qsup_b_tot, Qsup_g_tot, Other_r, &
                                             Use_res, Tse_res, strahler, Farea1, &
                                             Farea2, Farea3, Farea4, Zb0, Zdif, Hrsum, Hraverage, Hrmax, Amedas, Ame, rain, &
                                             Pro_rachi_total, past_pro_rachi_total, Depositable, &
                                             Sulen, Sualp, Suvol, Suelv, Suthecr, Suwcr, Sutheini, Suwini, Suthet, Sudthet, Suwt, &
                                             Sufs, Sufsmin, Surs, Suvs, Sutcol, SuproHB, SuRpr, &
                                             Pro_accum, Sup_accum, zb_lowest, Vwc, Vw, Ub1, Ub2, Tlag, Vwe, r_srock, r_hrock, &
                                             Rprsu, R_as, R_bs, R_cc, Vx, Hmax, Suh, Suq, Suinfin, Cdeb, &
                                             Water_accum, Cc, H, Zb, Et, Dzb1, Emb, Dm, Bwcom, Zbcom, Sup_vol, &
                                             Nbcom, Etcom, Embcom, Dzb1com, Dzbp1com, Dmcom, &
                                             Da, Db, Rka, Rkb, Poros2, srcol, rprmax, & !20171220suzuki全てparaから移動
                                             Surate, colrate, Suarea, susr, Sraverage, susrmax, Cs, Cscol, D, Susfq, &
                                             susfqmax, hanran, Zbcomini, &!!!!201801suzuki
                                             Suwid, Psar, Psl, Psp, Pqr, Qbtot, Qstot, Acsec, Suwidsum, &
                                             Qmc, Qfp, Umc, Ufp, Parate, &
                                             Qmax, & !最大流量
                                             Qtot, & !累積流量
                                             Qbmax, & !最大掃流砂量
                                             Qsmax, & !最大流量
                                             Zbmax, & !最大河床位
                                             Caprmax, & !最大断面平均濃度
                                             wlmax, &
                                             hhmax
  REAL(8), dimension(:, :), allocatable, save::Hr, Qr, Qbx, FmC, Qsu, &
                                                Ca, Capr, Qs, Prec, Prec2, &
                                                Fmubb, Ftubb, Ffmm_b, Ffmm_g, Srate, Pro_rate_gari, Pro_rate_rachi, &
                                                ti_rachi, ti_gari, &
                                                Thaw_time_aut, Thaw_time_spr, Dmc_pool, Vsed2, &
                                                Vout2, Vsup2, Vmax2, Ffmm_b1, Ffmm_b2, Cd, Qd, F0, dpro, Tlagrec, &
                                                Qsup_brec, Qsup_grec, &
                                                nftsurf, depth10, Prec3, Fm, Ft, Ff_b, Fmcom, Ftcom, &
                                                Sr, colHr, colQr, Qrsum, &!!!!201801suzuki
                                                Suhr, Suqr, Suqrcol, Suhrcol, &!!!斜面要素モデル用
                                                Caprmc, Caprfp, Qtotc, Qsufp, Qsmc, Qsfp, &
                                                Dzbp1, Wcon
  REAL(8), dimension(:, :, :), allocatable, save::Qsrec, Qbrec, Fd, Fdcom
  INTEGER, dimension(:), allocatable :: Sucat, SuMeshID, SuRainID, SuRnum, &
                                        SuCol, SuHbasin, Time_now, Time_past, Iav, &
                                        Jud1, Jud2, Jud3, Jud4, Jud5, Nb, Istate, colj, jud6, jud7, jud8, jud9, &
                                        Iin1, Iin2, Iin3, Iout1, Iright, Ileft
  integer, dimension(:), allocatable, save :: JDamtype ! 0ならばダムなし、1ならば不透過、2ならば透過
  real(8), dimension(:), allocatable, save :: Vdam, & !砂防堰堤の貯砂容量
                                              Vdep, &   !砂防堰堤の堆砂量
                                              Hdam   !砂防堰堤の高さ
end module

module para
  Real(8), save :: Fscin, Cmns, Cmn, Dt, Q0, &
                   Regeme, Dm1, Em, Ed, T, Data_time, Out_time, &
                   Hsini, Hcini, Qmin, Tp, Pd, Ps, &
                   Rate_rain, c_sf, S, Poros, G, C_karman, Ara, &
                   Hmin, T_end, zbdown, coltot, SuBls, pi, &
                   Slope_water_ini, Tot_rain, &
                   Dep_water, Cmnf, da0, db0, rka0, rkb0, dzbmin
  Integer, save :: Nxc, Nyc, Nys, Nxs, Nxc2, Nr, Nl, &
                   J1, J2, J3, J4, J5, J6, J7, Np, Nbl, Ntrate, Str_type, &
                   Nsu, Time, Time_max, J8, J10, J11, J12, J13, J14, Nspr, Tspr
  real(8), save :: epsiron, lpool, vpool, lpave, &
                   vpave, Qsup_b_tot_sum, Qbx_sum, Fscol, Qsuptot
end module

Program simhis
  use array
  use para
  Implicit none
  Real(8) Tdata, Tout
  character(1000) directory

  !OPEN FILES
  !input_data
  ! write(*,*) "input the name of directory"
  ! read(*,*) directory
  directory = "."
  Tdata = 0.d0
  Tout = 0.d0

  call open_files_input(directory)

  WRITE (*, *) 'Setting Initial Conditions...'

  if (str_type == 2) then
    READ (1102, *) Nxc
  else
    READ (1108, *) Nxc
  end if
  Np = 1
  READ (1006, *) Nr        !number_of_rainfall data
  READ (1006, *) Tp        !time_interval of rainfall data
  rate_rain = 1.0
  Nys = 10
  Fscin = 0.d0
  READ (1001, *) Cmns      !斜面のマニングの粗度係数
  READ (1001, *) Cmn       !流路のマニングの粗度係数
  READ (1001, *) Da0        !A層の厚さ(m)
  READ (1001, *) Db0        !B層の厚さ(m)
  READ (1001, *) Rka0       !A層の浸透係数(m/s)
  READ (1001, *) Rkb0       !B層の浸透係数(m/s)
  Cmnf = Cmn
  READ (1001, *) Hsini     !初期斜面水深
  READ (1001, *) Hcini     !初期流路水深
  READ (1001, *) Dt        !時間ステップ
  READ (1001, *) Q0        !下流端における流量(川幅決定用)
  READ (1001, *) Regeme    !レジーム定数
  READ (1001, *) Dm1       !一様砂使用時の平均粒径
  READ (1001, *) Data_time
  READ (1001, *) out_time
  READ (1001, *) dzbmin !

  Em = 0.2
  Ed = 0.2
  Nl = 10
  Nbl = 4
  T = 0.0
  Fscol = 1.0
  Subls = 10
  Nxc2 = 1
  ! READ (1001, *) J2        !計算の種類、流れ計算のみ=1,河床変動計算=2,step_pool_model=3
  J2 = 2
  ! READ (1001, *) J3        !河床材料の種類、一様砂=1、混合砂=2
  J3 = 1
  Np = 1
  J4 = 1
  J5 = 1
  J6 = 2
  J10 = 2
  J11 = 1
  J12 = 1
  J14 = 1
  ! READ (1104, *) Nsu       !the number of slope units
  Nsu=1

  Nspr = 1
  Tspr = int(Tp)
  J13 = 1

  Nxs = Nxc*2
  C_sf = 0.795d0             !静止摩擦係数
  S = 1.65d0                 !土砂の水中比重
  Poros = 0.4d0              !空隙率(河床）
!    Poros2=0.45d0           !空隙率(斜面）
  G = 9.81d0                 !重力加速度
  C_Karman = 0.4d0           !カルマン定数
  Ara = 6.d0                 !平均流に対する対数則定数
  Hmin = 0.0               !最小水深?
  T_end = Tp*Nr - 1        !計算終了時刻
  Fscin = Fscin/1000.d0/3600.d0
  Time_max = int(T_end/Dt - T)
  Time = 0
  epsiron = 0.41d0
  pi = 4.0d0*datan(1.0d0)

  allocate (Cl(Nxc), Cp(Nxc), Sl(Nxs), Sp(Nxs), San(Nxs), Sar(Nxs), &
            Di(Np), Fmub(Np), Ftub(Np), Free1(3000), Iin1(Nxc), Iin2(Nxc), Iin3(Nxc), &
            Ffm(Np), Iright(Nxc), Ileft(Nxc), Car(Nxc), Bw(Nxc), DmC(Nxc), Dzbp(Np), &
            Cave(Nxc), Q(Nxc), Zbmin(Nxc), Z(Nxc), Ux(Nxc), Iout1(Nxc), &
            Dzbf_b(Nxc), Dzbf_g(Nxc), Free2(3000), Free3(3000), Rpr(Nxs), Reti(Nr), &
            Qmax(Nxc), Qtot(Nxc), Hei_gari(Nxs), Pslo(Nxs), Pdep(Nxs), &
            B_area(Nxs), W_depth(Nxs), S_rate(Nxs), I_gari_ini(Nxs), &
            I_gari_tot(Nxs), Diam_b(Nxs), Ini_time_ft(Nxs), Tot_time_ft(Nxs), &
            Pro_gari(Nxs), Npool(Nxc), Vsed(Nxc), Vmax(Nxc), Vout(Nxc), &
            Vsup(Nxc), Vout_tot(Nxc), Use_res(Nxc), Tse_res(Nxc), Farea1(Nxc), &
            Farea2(Nxc), Farea3(Nxc), Farea4(Nxc), Zb0(Nxc), Zdif(Nxc), Hrsum(Nxs), &
            Hraverage(Nxs), Hrmax(Nxs), Amedas(Nxs), Ame(50), rain(100), &
            Sulen(Nsu), Sualp(Nsu), Suvol(Nsu), Suelv(Nsu), &
            Suthecr(Nsu), Suwcr(Nsu), Sutheini(Nsu), Suwini(Nsu), Suthet(Nsu), &
            Sudthet(Nsu), Suwt(Nsu), Sufs(Nsu), Sufsmin(Nsu), Surs(Nsu), &
            Suvs(Nsu), Sutcol(Nsu), SuproHB(Nxs), SuRpr(Nsu), Pro_accum(Nxs), &
            Sup_accum(Nxs), zb_lowest(Nxc), Vwc(Nxc), Vw(Nxc), Ub1(Nxc), Ub2(Nxc), &
            Tlag(Nxc), Vwe(Nxc), Rprsu(Nsu), R_as(Nxs), R_bs(Nxs), R_cc(Nxc))
  allocate (Sucat(Nsu), SuMeshID(Nsu), SuRainID(Nr), SuRnum(Nsu), &
            SuCol(Nsu), SuHbasin(Nsu))
  allocate (H(Nxc), Zb(Nxc), Hr(Nxs, Nys), Qr(Nxs, Nys), Qbx(Np, Nxc), &
            Et(Nxc), Prec(Nr, Nxs), &
            Prec2(Nr, Nxs), Fmubb(Nxc, Np), &
            Ftubb(Nxc, Np), Ffmm_b(Nxc, Np), Ffmm_g(Nxc, Np), Nb(Nxc), FmC(Np, Nxc), &
            Emb(Nxc), Dzbp1(Np, Nxc), &
            Qsu(Np, Nxc), &
            Ca(Np, Nxc), Capr(Np, Nxc), Qs(Np, Nxc), Dm(Nxc), &
            Dzb1(Nxc), Srate(Nr, Nxc), &
            Thaw_time_spr(Nxs, Nr), Ffmm_b1(Nxc, Np), Ffmm_b2(Nxc, Np), &
            Cd(Nxc, Np), Qd(Nxc, Np), &
            Prec3(Nr, Nsu))
  allocate (Fm(Np, Nxc), Ft(Np, Nxc), Ff_b(Np, Nxc))
  allocate (Fd(Np, Nl, Nxc))
  allocate (strahler(Nxc))
  allocate (gari_length(Nxs), Pro_rachi(Nxs), I_rachi_ini(Nxs), &
            I_rachi_tot(Nxs))
  allocate (Pro_rate_gari(Nxc, Nr), Pro_rate_rachi(Nxc, Nr), &
            ti_rachi(Nxc, Nr), ti_gari(Nxc, Nr))
  allocate (Qdep(Nxc), Ldep_b(Nxs), Ldep_g(Nxs), Bdep(Nxc), &
            Briv1(Nxc), Briv2(Nxc))
  allocate (Bare_r(Nxs), Width_b(Nxs), Depth_10(Nxs), Times_ft(Nxs), &
            Dep_ini_b(Nxs), Dep_ini_g(Nxs), Granite_r(Nxs), Other_r(Nxs))
  allocate (Qdep_g(Nxs), Qdep_b(Nxs), Bdep_g(Nxs), Bdep_b(Nxs), &
            Qsup_g(Nxs), Qsup_b(Nxs), Qsup_b_tot(Nxs), Qsup_g_tot(Nxs))
  allocate (Dmc_pool(Nxc, 10000), Vsed2(Nxc, 10000), Vout2(Nxc, 10000), &
            Vsup2(Nxc, 10000), Vmax2(Nxc, 10000))
  allocate (Pro_rachi_total(Nxs), past_pro_rachi_total(Nxs), Depositable(Nxs))
  allocate (Time_now(Nxc), Time_past(Nxc))
  allocate (r_srock(Nxs), r_hrock(Nxs), nftsurf(Nxs, Nr), depth10(Nxs, Nr))
  allocate (Vx(Nxc))
  allocate (Hmax(Nxc), Iav(Nsu))
  allocate (Suh(Nsu), Suq(Nsu), Suinfin(Nsu))
  allocate (Jud1(Nxc), Jud2(Nxc), Jud3(Nxc), Jud4(Nxc), Jud5(Nxc), jud6(Nsu), &
            jud8(Nsu), jud7(Nsu), jud9(Nsu), Cdeb(Nxc), Water_accum(Nxc) &
            , Cc(Nxc))
  allocate (Bwcom(Nxc), Zbcom(Nxc), Istate(Nxc), Sup_vol(Nxc), Nbcom(Nxc), &
            Etcom(Nxc), Embcom(Nxc), Dzb1com(Nxc), Dzbp1com(Np), &
            Fdcom(Nxc, Np, Nl), Fmcom(Nxc, Np), Ftcom(Nxc, Np), Dmcom(Nxc))
  allocate (Rka(Nxs), Rkb(Nxs), Da(Nxs), Db(Nxs), D(Nxs), Poros2(Nxs), &
            Sr(Nxs, Nys), susr(Nsu), Sraverage(Nxs), susrmax(Nsu))
  !20171220鈴木
  allocate (Suarea(Nsu), colrate(Nxs), Surate(Nsu), colHr(Nxs, Nys), &
            colQr(Nxs, Nys), colj(Nxs), Qrsum(Nxs, Nys), srcol(Nsu), Cs(Nsu), &
            Cscol(Nsu))
  allocate (susfq(Nsu), susfqmax(Nsu), rprmax(Nsu), hanran(Nxc), &
            Zbcomini(Nxc), Acsec(Nxc))
  allocate (Suwid(Nsu), Psar(Nxs), Psl(Nxs), Psp(Nxs), Pqr(Nxs), Suhr(Nsu, Nys), &
            Suqr(Nsu, Nys), Qbtot(Nxc), Qstot(Nxc), Suwidsum(Nxs))
  allocate (Qmc(Nxc), Qfp(Nxc), Umc(Nxc), Ufp(Nxc), Caprmc(Nxc, Np), &
            Caprfp(Nxc, Np), Qtotc(Nxc, Np), Qsufp(Nxc, Np), Qsmc(Nxc, Np), &
            Qsfp(Nxc, Np), Suhrcol(Nsu, Nys), Suqrcol(Nsu, Nys))
  allocate (Parate(Nxs))
  allocate (Qbmax(Nxc), Qsmax(Nxc), Zbmax(Nxc), Caprmax(Nxc), wlmax(Nxc), hhmax(Nxc))
  allocate (Wcon(Np, Nxc))
  call open_files_output(directory)
  call Initial_data_ucs
  WRITE (*, *) '計算開始!'
  do
    call rainfall_ucs
    call Channelflow
    call bedload
    call beddeformation
    if (T + Dt*0.01d0 < T_end) then !T=終了時刻 以外
      if ((Tdata + Dt*0.01d0 > data_time) .or. (T <= Dt*0.01d0)) then !T=終了時刻 以外
        call Results
        Tdata = Dt
      ELSE
        Tdata = Tdata + dt
      end if
      if ((Tout + Dt*0.01d0 > out_time) .or. (T <= Dt*0.01d0)) then
        call Screen
        Tout = Dt
      ELSE
        Tout = Tout + dt
      end if
      T = T + dt
      Time = Time + 1
      cycle
    ELSE !T=終了時刻
      call results
      call screen
      stop 'END of Calculation!!!!!'
    end if
  end do

contains

  subroutine beddeformation
    use array
    use para
    implicit none
    Integer I, L
    Real(8) Dzb, Qbin1, Qbin2, Qbin3, Qbout, Qsin1, Qsin2, Qsin3, Qsout, Dasp, Qsup

    !---- Calculation of bed elevation ----*
    Do I = 1, Nxc
      Dzb = 0.d0
      If (Z(I) > Zbmin(I) + Dm1*Hmin) then
        Dzb1(I) = 0.d0
        do L = 1, Np
          if (Iin1(I) > 0) then
            Qbin1 = Qbx(L, Iin1(I))
          ELSE
            Qbin1 = 0.d0
          end if
          if (Iin2(I) > 0) then
            Qbin2 = Qbx(L, Iin2(I))
          ELSE
            Qbin2 = 0.d0
          end if
          if (Iin3(I) > 0) then
            Qbin3 = Qbx(L, Iin3(I))
          ELSE
            Qbin3 = 0.d0
          end if
          Qbout = Qbx(L, I)
          if (Iin1(I) > 0) then
            Qsin1 = Qs(L, Iin1(I))
          ELSE
            Qsin1 = 0.d0
          end if
          if (Iin2(I) > 0) then
            Qsin2 = Qs(L, Iin2(I))
          ELSE
            Qsin2 = 0.d0
          end if
          if (Iin3(I) > 0) then
            Qsin3 = Qs(L, Iin3(I))
          ELSE
            Qsin3 = 0.d0
          end if
          Qsout = Qs(L, I)
          Qsup = 0.d0
          Qsuptot = 0.d0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DAsp = Dt/(1.d0 - Poros)*(1.d0/(Cl(I))*(Qbin1 + Qbin2 + Qbin3 - Qbout + Qsup))
          if (H(I) > 0.d0) then
            Dzbp1(L, I) = DAsp/Bw(I)    !L粒径階の河床変動量
            Dzb1(I) = Dzb1(I) + Dzbp1(L, I)
          ELSE
            Dzb1(I) = 0.d0
          end if
        end do
        if (H(I) > 0.d0) then
          Zb(I) = Zb(I) + Dzb1(I)
          Zbcom(I) = Zbcom(I) + Dzb1com(I)
        end if
      end if
    end do

    do I = 1, Nxc
      Zbmin(I) = 10000000.d0
      Zbmin(I) = Dmin1(Zbmin(I), Zb(I))
    end do
    do I = 1, Nxc
      Z(I) = H(I) + Zb(I)
    end do

    do I = 1, Nxc
      !      write(*,*)Cp(I)
      if (J6 == 1) then
        !       if(I == Nxc2)then
        if (Iout1(I) == 0) then
          !                Zbdown=Zb(I,1)-Cp(I)*Cl(I)
          Cp(I) = (Zb(I) - Zbdown)/Cl(I)
        ELSEIF (J6 == 1) then
          Cp(I) = (Zb(I) - Zb(Iout1(I)))/Cl(I)
        end if
        !if(Cp(I) <= 0.00001d0) Cp(I)=0.00001d0
        if (Cp(I) <= 0.0001d0) Cp(I) = 0.0001d0 !20161007
      end if
    end do

    do I = 1, Nxc
      If (Z(I) <= Zbmin(I) + Dm1*Hmin) then
        Z(I) = Zbmin(I) + Dm1*Hmin
        Q(I) = 0.d0
      end if
      H(I) = Z(I) - Zb(I)
      If (H(I) <= Dm1*Hmin) H(I) = Dm1*Hmin
    end do

  end subroutine beddeformation

  subroutine bedload
    use array
    use para
    implicit none
    integer I, L
    Real(8) U, Theta, C_kc, Ex, Ras, Aw, Se, Us, Ts, U_flux, Cf, Use, Tse, Tsc, Tscm, Usc, Uscm

    do I = 1, Nxc
      DmC(I) = Dm1
      do L = 1, Np
        Fmc(L, I) = 1.d0
      end do
    end do

    do I = 1, Nxc
      if (H(I) > 0.d0) then
        Dm(I) = DmC(I)
        do L = 1, Np
          Fm(L, I) = Fmc(L, I)
        end do
      end if
    end do

    !---- Calculation of sediment discharge ----*
    do I = 1, Nxc
      U = Ux(I)
      Theta = Datan(Cp(I))
      C_Kc = Dmax1(1.d0 + 1.d0/C_sf*((1.d0/S + 1.d0)*Dtan(Theta)), 0.5D0)
      Ex = 1.d0
      do L = 1, Np
        If (Z(I) > Zbmin(I) + Dm1*Hmin) then
          Ras = H(I)
          Aw = H(I)*Bw(I)
          Se = Dabs(Cmn**2.d0*Q(I)**2.d0/(Ras**(4.d0/3.d0)*Aw**2.d0))
          Us = Dsqrt(G*Ras*Se)
          Ts = Us**2.d0/(S*G*Di(L))
          Ux(I) = Q(I)/(Bw(I)*H(I))
          !~                    Use=Us_ave(Ux(I),Ras,DmC(I)*(1.+2.*Ts),C_Karman,Ara,Cmn)
          U_flux = Dabs(Ux(I))
          Cf = Dabs(Cmn**2.d0*9.8d0/Ras**(1.d0/3.d0))**0.5d0
          If (Ras > DmC(I)*(1.d0 + 2.d0*Ts)) then
            Use = U_flux/(Ara + 1.d0/C_Karman*Dlog(Ras/(DmC(I)*(1.d0 + 2.d0*Ts))))
          Else
            Use = U_flux/Ara
          End if
          Tse = Use**2.d0/(S*G*Di(L))
          DmC(I) = DmC(I)*100.d0
          If (DmC(I) >= 0.303d0) then
            Uscm = (80.9d0*DmC(I))**0.5d0
          Elseif (DmC(I) >= 0.118d0) then
            Uscm = (134.6d0*DmC(I)**(31.d0/22.d0))**0.5d0
          Elseif (DmC(I) >= 0.0565d0) then
            Uscm = (55.0d0*DmC(I))**0.5d0
          Elseif (DmC(I) >= 0.0065d0) then
            Uscm = (8.41d0*DmC(I)**(11.d0/32.d0))**0.5d0
          Elseif (DmC(I) < 0.0065) then
            Uscm = (226.0d0*DmC(I))**0.5d0
          else
            Uscm=0.d0
            write(*,*) "err"
            read(*,*)
          End if
          DmC(I) = DmC(I)/100.d0
          Uscm = Uscm/100.d0
          Tscm = Uscm**2.d0/(S*G*DmC(I))
          If ((Di(L)/DmC(I)) <= 0.4d0) then
            Usc = (0.85d0*Tscm*S*G*DmC(I))**0.5d0
          Else
            Usc = (1.64d0*Tscm*S*G*Di(L))**0.5d0/(Dlog10(19.d0*Di(L)/DmC(I)))
          End if
          Tsc = Usc**2.d0/(S*G*Di(L))
          if (Usc >= Us) then
            Qbx(L, I) = 0.d0
          Else
            Qbx(L, I) = 17.d0*(S*G*Di(L)**3.d0)**0.5d0*Tse**1.5d0*(1.d0 - Tsc/Ts)*(1.d0 - Usc/Us)*Bw(I)*Ex
            if (Zb(I) <= zb_lowest(I)) then
              Qbx(L, I) = 0.d0
            end if
          end if
        ELSE
          Qbx(L, I) = 0.d0
        end if
      end do
    end do

  end subroutine bedload

  subroutine Channelflow
    use array
    use para
    implicit none
    Real(8) Qin1, Qin2, Qin3, Qout, Qr2
    Integer I

    do I = 1, Nxc
      Q(I) = 1./Cmn*Bw(I)*Cp(I)**0.5d0*H(I)**(5.d0/3.d0)
    end do

    do I = 1, Nxc
      if (Iin1(I) > 0) then
        Qin1 = Q(Iin1(I))
      ELSE
        Qin1 = 0.0d0
      end if
      if (Iin2(I) > 0) then
        Qin2 = Q(Iin2(I))
      ELSE
        Qin2 = 0.d0
      end if
      if (Iin3(I) > 0) then
        Qin3 = Q(Iin3(I))
      ELSE
        Qin3 = 0.d0
      end if
      Qout = Q(I)
      Qr2 = Qr(Iright(I), Nys) + Qr(Ileft(I), Nys)
      Free1(I) = Z(I) + Dt*(1.d0/(Bw(I)*Cl(I))*(Qin1 + Qin2 + Qin3 - Qout) + 1./Bw(I)*Qr2)
    end do

    do I = 1, Nxc
      Z(I) = Free1(I)
      If (Z(I) <= Zbmin(I) + Dm1*Hmin) then
        Z(I) = Zbmin(I) + Dm1*Hmin
        Q(I) = 0.d0
      end if
      H(I) = Z(I) - Zb(I)
      If (H(I) <= Dm1*Hmin) then
        H(I) = Dm1*Hmin
      end if
    end do
  end subroutine Channelflow

  subroutine Initial_data_ucs
    use array
    use para
    implicit none
    Real(8) Car0, Prectemp, zb02, dummy, carpro
    Integer I, I1, I2, II, IIII, K, M, N, iost
    character(30) dummychar

    !watershed_condition
    Zdif = 0.d0
    READ (1102, *)
    DO I = 1, Nxc
      READ (1102, *) dummy, Iout1(I), Iin1(I), Iin2(I), Cl(I), Zb(I), Zb02
      if (Iout1(I) < 0) then
        Iout1(I) = 0
      end if
      if (J10 == 2) Bwcom(I) = 0.d0
      !Iin3(I)=Iin3(I)/2
      Ileft(I) = 2*I - 1
      Iright(I) = 2*I
      !        Zb(I,1)=Zb(I,1)+Zdif(I)
      Zb0(I) = Zb(I)
      if (Cl(I) < 1.d0) Cl(I) = 1.d0
      Cp(I) = (Zb(I) - (Zb02 + Zdif(I)))/Cl(I)
    end do

    read (1103, *)
    DO
      READ (1103, *, iostat=iost) II, Sar(II), Sp(II)
      if (iost < 0) exit
    end do
    DO I = 1, Nxc
      I1 = 2*I - 1
      I2 = 2*I
      if (Sar(I1) < 1.d0) Sar(I1) = 1.d0
      if (Sar(I2) < 1.d0) Sar(I2) = 1.d0
      if (Sp(I1) < 0.001d0) Sp(I1) = 0.001d0
      if (Sp(I2) < 0.001d0) Sp(I2) = 0.001d0
      Sl(I1) = Sar(I1)/(Cl(I))/Dcos(Sp(I1)*2*3.14159265358979323846D0/360.d0)
      Sl(I2) = Sar(I2)/(Cl(I))/Dcos(Sp(I2)*2*3.14159265358979323846D0/360.d0)

      if (Sl(I1) < 1.) Sl(I1) = 1.
      if (Sl(I2) < 1.) Sl(I2) = 1.
    end do
    !calclate catchmnent area
    DO I = 1, Nxc
      Car(I) = 0.d0
    end do
    DO
      Carpro = 1.d0
      DO I = 1, Nxc
        Iiii = 2*I
        if ((Iin1(I) == 0) .and. (Iin2(I) == 0) .and. (Iin3(I) == 0)) then
          Car(I) = Sar(Iiii) + Sar(Iiii - 1)
        ELSEIF (Iin1(I) > 0 .and. Iin2(I) > 0 .and. Iin3(I) > 0) then
          if (Car(Iin1(I))*Car(Iin2(I))*Car(Iin3(I)) /= 0.) then
            Car(I) = Car(Iin1(I)) + Car(Iin2(I)) + Car(Iin3(I)) + Sar(Iiii) + Sar(Iiii - 1)
          end if
        ELSEIF (Iin1(I) > 0 .and. Iin2(I) > 0) then
          if (Car(Iin1(I))*Car(Iin2(I)) /= 0.) then
            Car(I) = Car(Iin1(I)) + Car(Iin2(I)) + Sar(Iiii) + Sar(Iiii - 1)
          end if
        ELSEIF (Iin1(I) > 0 .and. Iin2(I) == 0) then
          if (Car(Iin1(I)) /= 0.) then
            Car(I) = Car(Iin1(I)) + Sar(Iiii) + Sar(Iiii - 1)
          end if
        end if
        if (Car(I) > 0.1d0) then
          Carpro = Carpro*1.d0
        else
          Carpro = Carpro*0.d0
        end if
        !Carpro=Carpro*Car(I)
      end do
      if (carpro > 0.1d0) exit
      !        write(*,*) carpro
    end do
    Car0 = maxval(Car)

    ! READ (1002, *)
    DO I = 1, Nxc

      if (Bw(I) == 0.d0) Bw(I) = Regeme*Dsqrt(Q0*Car(I)/Car0)
      if (Hmax(I) == 0.d0) Hmax(I) = Bw(I)/2.d0
      Zbcom(I) = Zb(I) + Hmax(I)
      Zbcomini(I) = Zbcom(I)
      Zbmin(I) = Zb(I)
      zb_lowest(I) = Zb(I) + dzbmin
      Zb(I) = Zb(I)
    end do

    Do I = 1, Nxc
      if (Iout1(I) == 0) then
        Zbdown = Zb(I) - Cp(I)*Cl(I)
        Cp(I) = (Zb(I) - Zbdown)/Cl(I)
        if (Cp(I) <= 0.00001d0) Cp(I) = 0.00001d0
      ELSE
        Cp(I) = (Zb(I) - Zb(Iout1(I)))/Cl(I)
        if (Cp(I) <= 0.00001d0) Cp(I) = 0.00001d0
      end if
    end do

    DO I = 1, Nxc
      H(I) = Hcini*(Car(I)/Car0)**0.6d0 !20170310
      Z(I) = H(I) + Zb(I)
      Q(I) = 1./Cmn*Bw(I)*Cp(I)**0.5d0*H(I)**(5.d0/3.d0)  !20170523'
      Emb(I) = Em
      Nb(I) = Nbl
      Et(I) = Em
      Embcom(I) = Em
      Nbcom(I) = Nbl
      Etcom(I) = Em
    end do

    !Unit slopes
    ! read (1008, *)
    do M = 1, Nxs
      ! read (1008, *) dummy, Rka(M), Rkb(M), Da(M), Db(M), Poros2(M), Parate(M)
      Rka(M) = rka0
      Rkb(M) = rkb0
      Da(M) = da0
      Db(M) = db0
      Poros2(M) = 0.56d0
      Parate(M) = 0.d0
      D(M) = Da(M) + Db(M)
    end do

    DO M = 1, Nxs
      DO N = 1, Nys
        if (N == 1) then
          Hr(M, N) = 0.d0
        ELSEIF (N == 2) then
          Hr(M, N) = Hsini*0.5
        ELSE
          Hr(M, N) = Hsini
        end if
        Qr(M, N) = 0.d0
      end do
    end do

    !rainfall
    READ (1006, *)
    DO II = 1, Nr
      READ (1006, *) dummychar, prectemp
      DO M = 1, Nxs
        Prec(II, M) = rate_rain*prectemp/1000.d0/3600.d0
      end do
    end do

    !grainsize_distribution
    DO I = 1, Nxc
      Di(1) = Dm1
      Fm(1, I) = 1.d0
      Ft(1, I) = 1.d0
      Ff_b(1, I) = 1.d0
      DO K = 1, Nl
        Fd(1, K, I) = 1.d0
      end do
    end do

    DO I = 1, Nxc
      write (1209, *) I, Bw(I)
      write (1208, *) I, Car(I)
    end do

    write (*, *) "END READING ALL FILES!!!"

  end subroutine Initial_data_ucs

  subroutine open_files_input(directory)
    use array
    use para
    Implicit none
    character(1000), intent(in) :: directory
    character(1000) filetemp, filename
    integer iost

    filetemp = "/conditions/parameters.dat"
    filename = "./"//trim(directory)//trim(filetemp)
    OPEN (1001, FILE=filename, STATUS='old', action='read') !初期条件もろもろ

    filetemp = "/conditions/rainfall.dat"
    filename = "./"//trim(directory)//trim(filetemp)
    OPEN (1006, FILE=filename, STATUS='old', action='read') !

    !river_data
    filetemp = "/watershed_data/unitchannels_02.dat"
    filename = "./"//trim(directory)//trim(filetemp)
    OPEN (1102, FILE=filename, STATUS='old', action='read', iostat=iost)
    str_type = 2

    filetemp = "/watershed_data/unit_slopes.dat"
    filename = "./"//trim(directory)//trim(filetemp)
    OPEN (1103, FILE=filename, STATUS='old', action='read', iostat=iost)


  end subroutine open_files_input

  subroutine open_files_output(directory)
    use array
    use para
    Implicit none
    character(1000), intent(in) :: directory
    character(1000) filetemp, filename

    !results
    !流量
    filetemp = "/results/Q.dat"
    filename = "./"//trim(directory)//trim(filetemp)
    OPEN (1201, FILE=filename, STATUS='UNKNOWN')
    !掃流砂量
    filetemp = "/results/Qb.dat"
    filename = "./"//trim(directory)//trim(filetemp)
    OPEN (1202, FILE=filename, STATUS='UNKNOWN')
    !河床位
    filetemp = "/results/Zb.dat"
    filename = "./"//trim(directory)//trim(filetemp)
    OPEN (1205, FILE=filename, STATUS='UNKNOWN')
    !集水面積
    filetemp = "/results/Car.dat"
    filename = "./"//trim(directory)//trim(filetemp)
    OPEN (1208, FILE=filename, STATUS='UNKNOWN')
    !川幅
    filetemp = "/results/Bw.dat"
    filename = "./"//trim(directory)//trim(filetemp)
    OPEN (1209, FILE=filename, STATUS='UNKNOWN')
    !河床変動量
    filetemp = "/results/Dzb.dat"
    filename = "./"//trim(directory)//trim(filetemp)
    OPEN (1220, FILE=filename, STATUS='UNKNOWN')
    !水深
    filetemp = "/results/H.dat"
    filename = "./"//trim(directory)//trim(filetemp)
    OPEN (1221, FILE=filename, STATUS='UNKNOWN')
  end subroutine open_files_output

  subroutine rainfall_ucs
    use array
    use para
    implicit none
    Real(8) Dy, Alp, Theta
    Integer II, M, N

    II = int((T/Tp - 0.0000001d0)) + 1  !steprain割り切れるとき（正時)は前の時刻のデータを使用する．
    do M = 1, Nxs
      Rpr(M) = Prec(II, M)    !steprain
    end do

    do M = 1, Nxs
      Dy = Sl(M)/(Nys - 1)
      Qr(M, 1) = 0.
      Theta = Sp(M)*2.d0*3.14159265358979323846D0/360.d0
      do N = 2, Nys
        if (Hr(M, N) < Db(M)) then
          Qr(M, N) = Rkb(M)*Hr(M, N)*Dsin(Theta)
          Alp = Poros2(M)
        ELSEIF (Hr(M, N) < Da(M) + Db(M)) then
          Qr(M, N) = Rka(M)*(Hr(M, N) - Db(M))*Dsin(Theta) + Rkb(M)*Db(M)*Dsin(Theta)
          Alp = Poros2(M)
        ELSE
          Qr(M, N) = Rka(M)*Da(M)*Dsin(Theta) + Rkb(M)*Db(M)* &
                     Dsin(Theta) + 1./Cmns*Dsqrt(Dsin(Theta))*(Hr(M, N) - (Da(M) + Db(M)))**(5./3.)
          Alp = 1.
        end if
        Hr(M, N) = Hr(M, N) + Dt/Alp*(-(Qr(M, N) - Qr(M, N - 1))/Dy + (Rpr(M) - Fscin)*Dcos(Theta))
        if (Hr(M, N) <= 0.d0) Hr(M, N) = 0.d0
      end do
    end do
    Qrsum(1:Nxs, 1:Nys) = Qr(1:Nxs, 1:Nys)

    do M = 1, Nxs
      Hrsum(M) = 0.d0
      do N = 2, Nys
        if (Hr(M, N) > Da(M) + Db(M)) then
          Hrsum(M) = Hrsum(M) + Da(M) + Db(M)
        ELSE
          Hrsum(M) = Hrsum(M) + Hr(M, N)
        end if
      end do
      Hraverage(M) = Hrsum(M)/(Nys - 1)
      if (Hraverage(M) > Hrmax(M)) Hrmax(M) = Hraverage(M)
    end do

  end subroutine rainfall_ucs

  subroutine results
    use array
    use para
    implicit none
    integer I, L

    write (*, *) 'writing files... T=', T

    !Q
    if (T <= Dt*0.01) Write (1201, '(A15, 3000(I15))') 'T(s)', (I, I=1, Nxc)
    Write (1201, '(F15.3, 3000(E15.7))') T, (Q(I), I=1, Nxc)

    !Qb
    do I = 1, Nxc
      Free1(I) = 0.d0
      do L = 1, Np
        Free1(I) = Free1(I) + Qbx(L, I)
      end do
    end do
    if (T <= Dt*0.01) Write (1202, '(A15, 3000(I15))') 'T(s)', (I, I=1, Nxc)
    Write (1202, '(F15.3, 3000(E15.7))') T, (Free1(I), I=1, Nxc)

    !Zb
    if (T <= Dt*0.01) Write (1205, '(A15, 3000(I15))') 'T(s)', (I, I=1, Nxc)
    Write (1205, '(F15.3, 3000(E15.7))') T, (Zb(I), I=1, Nxc)

    !Dzb
    do I = 1, Nxc
      Free1(I) = Zb(I) - Zb0(I)
    end do
    if (T <= Dt*0.01) Write (1220, '(A15, 3000(I15))') 'T(s)', (I, I=1, Nxc)
    Write (1220, '(F15.3, 3000(E15.7))') T, (Free1(I), I=1, Nxc)

    !H
    if (T <= Dt*0.01) Write (1221, '(A15, 3000(I15))') 'T(s)', (I, I=1, Nxc)
    Write (1221, '(F15.3, 3000(E15.7))') T, (H(I), I=1, Nxc)
  end subroutine results

  subroutine screen
    use array
    use para
    implicit none
    real(8) Rpr11

    Write (*, *)
    Write (*, *)
    Write (*, *)
    Write (*, *)
    Write (*, *)
    Write (*, *)
    Write (*, *) 'Sediment Runoff model'
    Write (*, '(A17, F15.5, A6, A17, F15.7, A6)') 'calculation time=', T, ' s    '
    Rpr11 = Rpr(Nxs)*1000.d0*3600.d0
    Write (*, '(A17, F15.5, A6, A17, F15.7, A6)') 'END TIME=', T_end, ' s    ', 'Rainfall', Rpr11, ' mm/h '
    Write (*, '(A17, F15.5, A6, A17, F15.7, A6)') 'dt=', Dt, ' s    '
    Write (*, *) 'TIME=', TIME, 'TIME_MAX=', TIME_MAX
    Write (*, '(A17, F15.5, A6, A17, F15.7, A6)') 'Manning(slope)=', Cmns, ' s-m  '
    Write (*, '(A17, F15.5, A6, A17, F15.7, A6)') 'Manning(river)=', Cmn, ' s-m  ', 'Em', Em, ' m    '
    Write (*, '(A17, F15.5, A6, A17, F15.7, A6)') 'Q=', Q(Nxc2), ' m3/s '
  end subroutine screen

End Program simhis
