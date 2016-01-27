
MODULE module_mpas_atm_dyn
!
   IMPLICIT NONE
!
CONTAINS
!
$$$ ! replace with atm_compute_dyn_tend_work()
  SUBROUTINE wsm62D(t, q                                          &   
                   ,qci, qrs, den, p, delz                        &
                   ,delt,g, cpd, cpv, rd, rv, t0c                 &
                   ,ep1, ep2, qmin                                &
                   ,XLS, XLV0, XLF0, den0, denr                   &
                   ,cliq,cice,psat                                &
                   ,lat                                           &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte                     &
                   ,snow,snowncv                                  &
                   ,graupel,graupelncv                            &
                                                                  )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!
!  This code is a 6-class GRAUPEL phase microphyiscs scheme (WSM6) of the 
!  Single-Moment MicroPhyiscs (WSMMP). The WSMMP assumes that ice nuclei
!  number concentration is a function of temperature, and seperate assumption
!  is developed, in which ice crystal number concentration is a function
!  of ice amount. A theoretical background of the ice-microphysics and related
!  processes in the WSMMPs are described in Hong et al. (2004).
!  All production terms in the WSM6 scheme are described in Hong and Lim (2006).
!  All units are in m.k.s. and source/sink terms in kgkg-1s-1.
!
!  WSM6 cloud scheme
!
!  Coded by Song-You Hong and Jeong-Ock Jade Lim (Yonsei Univ.)
!           Summer 2003
!
!  Implemented by Song-You Hong (Yonsei Univ.) and Jimy Dudhia (NCAR)
!           Summer 2004
!
!  History :  semi-lagrangian scheme sedimentation(JH), and clean up
!             Hong, August 2009
!
!  Reference) Hong, Dudhia, Chen (HDC, 2004) Mon. Wea. Rev.
!             Hong and Lim (HL, 2006) J. Korean Meteor. Soc.
!             Dudhia, Hong and Lim (DHL, 2008) J. Meteor. Soc. Japan
!             Lin, Farley, Orville (LFO, 1983) J. Appl. Meteor.
!             Rutledge, Hobbs (RH83, 1983) J. Atmos. Sci.
!             Rutledge, Hobbs (RH84, 1984) J. Atmos. Sci.
!             Juang and Hong (JH, 2010) Mon. Wea. Rev.
!
  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte,  &
                                      lat
  REAL, DIMENSION( its:ite , kts:kte ),                           &
        INTENT(INOUT) ::                                          &
                                                               t
  REAL, DIMENSION( its:ite , kts:kte, 2 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qci
  REAL, DIMENSION( its:ite , kts:kte, 3 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qrs
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                               q
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(IN   ) ::                                          &
                                                             den, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                             cpd, &
                                                             cpv, &
                                                             t0c, &
                                                            den0, &
                                                              rd, &
                                                              rv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( ims:ime ),                                     &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr
  REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv
  REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                 graupel, &
                                                      graupelncv
! LOCAL VAR
  REAL, DIMENSION( its:ite , kts:kte , 3) ::                      &
                                                              rh, &
                                                              qs, &
                                                          rslope, &
                                                         rslope2, &
                                                         rslope3, &
                                                         rslopeb, &
                                                         qrs_tmp, & 
                                                            falk, &
                                                            fall, &
                                                           work1
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                           fallc, &
                                                           falkc, &
                                                          work1c, &
                                                          work2c, &
                                                           workr, &
                                                           worka
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                         den_tmp, &
                                                        delz_tmp
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                           pigen, &
                                                           pidep, &
                                                           pcond, &
                                                           prevp, &
                                                           psevp, &
                                                           pgevp, &
                                                           psdep, &
                                                           pgdep, &
                                                           praut, &
                                                           psaut, &
                                                           pgaut, &
                                                           piacr, &
                                                           pracw, &
                                                           praci, &
                                                           pracs, &
                                                           psacw, &
                                                           psaci, &
                                                           psacr, &
                                                           pgacw, &
                                                           pgaci, &
                                                           pgacr, &
                                                           pgacs, &
                                                           paacw, &
                                                           psmlt, &
                                                           pgmlt, &
                                                           pseml, &
                                                           pgeml
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                            qsum, &
                                                              xl, &
                                                             cpm, &
                                                           work2, &
                                                          denfac, &
                                                             xni, &
                                                         denqrs1, &
                                                         denqrs2, &
                                                         denqrs3, &
                                                          denqci, & 
                                                          delta2, &
                                                          delta3, &
                                                           n0sfac
  REAL, DIMENSION( its:ite ) ::                          delqrs1, &
                                                         delqrs2, &
                                                         delqrs3, &
                                                           delqi  
  REAL, DIMENSION( its:ite ) ::                        tstepsnow, &
                                                      tstepgraup
  INTEGER, DIMENSION( its:ite ) ::                         mstep, &
                                                           numdt
  LOGICAL, DIMENSION( its:ite ) ::                        flgcld
  REAL  ::                                                        &
            cpmcal, xlcal, diffus,                                &
            viscos, xka, venfac, conden, diffac,                  &
            x, y, z, a, b, c, d, e,                               &
            qdt, holdrr, holdrs, holdrg, supcol, supcolt, pvt,    &
            coeres, supsat, dtcld, xmi, eacrs, satdt,             &
            qimax, diameter, xni0, roqi0,                         &
            fallsum, fallsum_qsi, fallsum_qg,                     &
            vt2i,vt2r,vt2s,vt2g,acrfac,egs,egi,                   &
            xlwork2, factor, source, value,                       &
            xlf, pfrzdtc, pfrzdtr, supice, alpha2
  REAL  :: vt2ave
  REAL  :: holdc, holdci
  INTEGER :: i, j, k, mstepmax,                                   &
            iprt, latd, lond, loop, loops, ifsat, n, idim, kdim
  INTEGER :: itest,ktest
! Temporaries used for inlining fpvs function
  REAL  :: dldti, xb, xai, tr, xbi, xa, hvap, cvap, hsub, dldt, ttp
! variables for optimization
  REAL, DIMENSION( its:ite ) ::                             tvec1
  REAL                       ::                              temp
!
!=================================================================
!   compute internal functions
!
      cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
      xlcal(x) = xlv0-xlv1*(x-t0c)
!----------------------------------------------------------------
!     diffus: diffusion coefficient of the water vapor
!     viscos: kinematic viscosity(m2s-1)
!     Optimizatin : A**B => exp(log(A)*(B))
!
      diffus(x,y) = 8.794e-5 * exp(log(x)*(1.81)) / y        ! 8.794e-5*x**1.81/y
      viscos(x,y) = 1.496e-6 * (x*sqrt(x)) /(x+120.)/y  ! 1.496e-6*x**1.5/(x+120.)/y
      xka(x,y) = 1.414e3*viscos(x,y)*y
      diffac(a,b,c,d,e) = d*a*a/(xka(c,d)*rv*c*c)+1./(e*diffus(c,b))
      venfac(a,b,c) = exp(log((viscos(b,c)/diffus(b,a)))*((.3333333)))         &
                     /sqrt(viscos(b,c))*sqrt(sqrt(den0/c))
      conden(a,b,c,d,e) = (max(b,qmin)-c)/(1.+d*d/(rv*e)*c/(a*a))
!
#ifdef ALIGN_OK
!DIR$ ASSUME_ALIGNED t:64,qci:64,qrs:64,q:64,den:64,p:64,delz:64,rain:64,rainncv:64,sr:64
!DIR$ ASSUME_ALIGNED snow:64,snowncv:64,graupel:64,graupelncv:64,mstep:64,flgcld:64
!DIR$ ASSUME_ALIGNED worka:64,pgacs:64,n0sfac:64,work2:64,psmlt:64,rslope:64
!DIR$ ASSUME_ALIGNED rslope2:64,rslope3:64,rslopeb:64,cpm:64,pgmlt:64
!DIR$ ASSUME_ALIGNED tstepsnow:64,tstepgraup:64
#endif
!
      idim = ite-its+1
      kdim = kte-kts+1
itest=979
ktest=1
!
!----------------------------------------------------------------
!     padding 0 for negative values generated by dynamics
!
      do k = kts, kte
        do i = its, ite
          qci(i,k,1) = max(qci(i,k,1),0.0)
          qrs(i,k,1) = max(qrs(i,k,1),0.0)
          qci(i,k,2) = max(qci(i,k,2),0.0)
          qrs(i,k,2) = max(qrs(i,k,2),0.0)
          qrs(i,k,3) = max(qrs(i,k,3),0.0)
        enddo
      enddo
!
!----------------------------------------------------------------
!     latent heat for phase changes and heat capacity. neglect the
!     changes during microphysical process calculation
!     emanuel(1994)
!
      do k = kts, kte
        do i = its, ite
          cpm(i,k) = cpmcal(q(i,k))
          xl(i,k) = xlcal(t(i,k))
        enddo
      enddo
      do k = kts, kte
        do i = its, ite
          delz_tmp(i,k) = delz(i,k)
          den_tmp(i,k) = den(i,k)
        enddo
      enddo
!
!----------------------------------------------------------------
!    initialize the surface rain, snow, graupel
!
      do i = its, ite
        rainncv(i) = 0.
        if(PRESENT (snowncv) .AND. PRESENT (snow)) snowncv(i,lat) = 0.
        if(PRESENT (graupelncv) .AND. PRESENT (graupel)) graupelncv(i,lat) = 0.
        sr(i) = 0.
! new local array to catch step snow and graupel
        tstepsnow(i) = 0.
        tstepgraup(i) = 0.
      enddo
!
!----------------------------------------------------------------
!     compute the minor time steps.
!
      loops = max(nint(delt/dtcldcr),1)
      dtcld = delt/loops
      if(delt.le.dtcldcr) dtcld = delt
!
      do loop = 1,loops
!
!----------------------------------------------------------------
!     initialize the large scale variables
!
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
      do i = its, ite
        mstep(i) = 1
        flgcld(i) = .true.
      enddo

      do k = kts, kte
        do i = its, ite
          denfac(i,k) = sqrt(den0/den(i,k))
        enddo
      enddo
!
! Inline expansion for fpvs
!         qs(i,k,1) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!         qs(i,k,2) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      do k = kts, kte
        do i = its, ite
          tr=ttp/t(i,k)
          qs(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          qs(i,k,1) = min(qs(i,k,1),0.99*p(i,k))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          rh(i,k,1) = max(q(i,k) / qs(i,k,1),qmin)
          tr=ttp/t(i,k)
          if(t(i,k).lt.ttp) then
            qs(i,k,2)=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
          else
            qs(i,k,2)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          endif
          qs(i,k,2) = min(qs(i,k,2),0.99*p(i,k))
          qs(i,k,2) = ep2 * qs(i,k,2) / (p(i,k) - qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
          rh(i,k,2) = max(q(i,k) / qs(i,k,2),qmin)
        enddo
      enddo
!
!----------------------------------------------------------------
!     initialize the variables for microphysical physics
!
!
      do k = kts, kte
        do i = its, ite
          prevp(i,k) = 0.
          psdep(i,k) = 0.
          pgdep(i,k) = 0.
          praut(i,k) = 0.
          psaut(i,k) = 0.
          pgaut(i,k) = 0.
          pracw(i,k) = 0.
          praci(i,k) = 0.
          piacr(i,k) = 0.
          psaci(i,k) = 0.
          psacw(i,k) = 0.
          pracs(i,k) = 0.
          psacr(i,k) = 0.
          pgacw(i,k) = 0.
          paacw(i,k) = 0.
          pgaci(i,k) = 0.
          pgacr(i,k) = 0.
          pgacs(i,k) = 0.
          pigen(i,k) = 0.
          pidep(i,k) = 0.
          pcond(i,k) = 0.
          psmlt(i,k) = 0.
          pgmlt(i,k) = 0.
          pseml(i,k) = 0.
          pgeml(i,k) = 0.
          psevp(i,k) = 0.
          pgevp(i,k) = 0.
          falk(i,k,1) = 0.
          falk(i,k,2) = 0.
          falk(i,k,3) = 0.
          fall(i,k,1) = 0.
          fall(i,k,2) = 0.
          fall(i,k,3) = 0.
          fallc(i,k) = 0.
          falkc(i,k) = 0.
          xni(i,k) = 1.e3
        enddo
      enddo
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
      do k = kts, kte
        do i = its, ite
          temp = (den(i,k)*max(qci(i,k,2),qmin))
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
        enddo
      enddo
!
!----------------------------------------------------------------
!     compute the fallout term:
!     first, vertical terminal velosity for minor loops
!----------------------------------------------------------------
      do k = kts, kte
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
        enddo
      enddo
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, & 
                     work1,its,ite,kts,kte)
!
      do k = kte, kts, -1
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
          workr(i,k) = work1(i,k,1)
          qsum(i,k) = max( (qrs(i,k,2)+qrs(i,k,3)), 1.E-15)
          IF ( qsum(i,k) .gt. 1.e-15 ) THEN
            worka(i,k) = (work1(i,k,2)*qrs(i,k,2) + work1(i,k,3)*qrs(i,k,3)) &
                      /qsum(i,k)
          ELSE
            worka(i,k) = 0.
          ENDIF
          denqrs1(i,k) = den(i,k)*qrs(i,k,1)
          denqrs2(i,k) = den(i,k)*qrs(i,k,2)
          denqrs3(i,k) = den(i,k)*qrs(i,k,3)
          if(qrs(i,k,1).le.0.0) workr(i,k) = 0.0
        enddo
      enddo
      call _NISLFV_RAIN_PLM_(its,ite,kts,kte,den_tmp,denfac,t,delz_tmp,workr,denqrs1,  &
                           delqrs1,dtcld,1,1)
      call _NISLFV_RAIN_PLM6_(its,ite,kts,kte,den_tmp,denfac,t,delz_tmp,worka,         & 
                           denqrs2,denqrs3,delqrs2,delqrs3,dtcld,1,1)
      do k = kts, kte
        do i = its, ite
          qrs(i,k,1) = max(denqrs1(i,k)/den(i,k),0.)
          qrs(i,k,2) = max(denqrs2(i,k)/den(i,k),0.)
          qrs(i,k,3) = max(denqrs3(i,k)/den(i,k),0.)
          fall(i,k,1) = denqrs1(i,k)*workr(i,k)/delz(i,k)
          fall(i,k,2) = denqrs2(i,k)*worka(i,k)/delz(i,k)
          fall(i,k,3) = denqrs3(i,k)*worka(i,k)/delz(i,k)
        enddo
      enddo
      do i = its, ite
        fall(i,1,1) = delqrs1(i)/delz(i,1)/dtcld
        fall(i,1,2) = delqrs2(i)/delz(i,1)/dtcld
        fall(i,1,3) = delqrs3(i)/delz(i,1)/dtcld
      enddo
      do k = kts, kte
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
        enddo
      enddo
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                     work1,its,ite,kts,kte)
!
      do k = kte, kts, -1 
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
          supcol = t0c-t(i,k)
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(t(i,k).gt.t0c) then
!---------------------------------------------------------------
! psmlt: melting of snow [HL A33] [RH83 A25]
!       (T>T0: S->R)
!---------------------------------------------------------------
            xlf = xlf0
            work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
            if(qrs(i,k,2).gt.0.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psmlt(i,k) = xka(t(i,k),den(i,k))/xlf*(t0c-t(i,k))*pi/2.       &
                         *n0sfac(i,k)*(precs1*rslope2(i,k,2)                 &
                         +precs2*work2(i,k)*coeres)
              psmlt(i,k) = min(max(psmlt(i,k)*dtcld/mstep(i),                &
                          -qrs(i,k,2)/mstep(i)),0.)
              qrs(i,k,2) = qrs(i,k,2) + psmlt(i,k)
              qrs(i,k,1) = qrs(i,k,1) - psmlt(i,k)
              t(i,k) = t(i,k) + xlf/cpm(i,k)*psmlt(i,k)
            endif
!---------------------------------------------------------------
! pgmlt: melting of graupel [HL A23]  [LFO 47]
!       (T>T0: G->R)
!---------------------------------------------------------------
            if(qrs(i,k,3).gt.0.) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgmlt(i,k) = xka(t(i,k),den(i,k))/xlf                          &
                         *(t0c-t(i,k))*(precg1*rslope2(i,k,3)                &
                         +precg2*work2(i,k)*coeres)
              pgmlt(i,k) = min(max(pgmlt(i,k)*dtcld/mstep(i),                &
                          -qrs(i,k,3)/mstep(i)),0.)                          
              qrs(i,k,3) = qrs(i,k,3) + pgmlt(i,k)
              qrs(i,k,1) = qrs(i,k,1) - pgmlt(i,k)
              t(i,k) = t(i,k) + xlf/cpm(i,k)*pgmlt(i,k)
            endif
          endif
        enddo
      enddo
!---------------------------------------------------------------
! Vice [ms-1] : fallout of ice crystal [HDC 5a]
!---------------------------------------------------------------
      do k = kte, kts, -1
        do i = its, ite
          if(qci(i,k,2).le.0.) then
            work1c(i,k) = 0.
          else
            xmi = den(i,k)*qci(i,k,2)/xni(i,k)
            diameter  = max(min(dicon * sqrt(xmi),dimax), 1.e-25)
            work1c(i,k) = 1.49e4*exp(log(diameter)*(1.31))
          endif
        enddo
      enddo
!
!  forward semi-laglangian scheme (JH), PCM (piecewise constant),  (linear)
!
      do k = kte, kts, -1
        do i = its, ite
          denqci(i,k) = den(i,k)*qci(i,k,2)
        enddo
      enddo
      call _NISLFV_RAIN_PLM_(its,ite,kts,kte,den_tmp,denfac,t,delz_tmp,work1c,denqci,  &
                           delqi,dtcld,1,0)
      do k = kts, kte
        do i = its, ite
          qci(i,k,2) = max(denqci(i,k)/den(i,k),0.)
        enddo
      enddo
      do i = its, ite
        fallc(i,1) = delqi(i)/delz(i,1)/dtcld
      enddo
!
!----------------------------------------------------------------
!      rain (unit is mm/sec;kgm-2s-1: /1000*delt ===> m)==> mm for wrf
!
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
      do i = its, ite
        fallsum = fall(i,kts,1)+fall(i,kts,2)+fall(i,kts,3)+fallc(i,kts)
        fallsum_qsi = fall(i,kts,2)+fallc(i,kts)
        fallsum_qg = fall(i,kts,3)
        if(fallsum.gt.0.) then
          rainncv(i) = fallsum*delz(i,kts)/denr*dtcld*1000. + rainncv(i)
          rain(i) = fallsum*delz(i,kts)/denr*dtcld*1000. + rain(i)
        endif
        if(fallsum_qsi.gt.0.) then
          tstepsnow(i)   = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.            &
                           +tstepsnow(i)
        IF ( PRESENT (snowncv) .AND. PRESENT (snow)) THEN
          snowncv(i,lat) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.            & 
                           +snowncv(i,lat)
          snow(i,lat) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000. + snow(i,lat)
        ENDIF
        endif
        if(fallsum_qg.gt.0.) then
          tstepgraup(i)  = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.            &
                           +tstepgraup(i)
        IF ( PRESENT (graupelncv) .AND. PRESENT (graupel)) THEN
          graupelncv(i,lat) = fallsum_qg*delz(i,kts)/denr*dtcld*1000.          &   
                              + graupelncv(i,lat)
          graupel(i,lat) = fallsum_qg*delz(i,kts)/denr*dtcld*1000. + graupel(i,lat)
        ENDIF
        endif
!       if(fallsum.gt.0.)sr(i)=(snowncv(i,lat) + graupelncv(i,lat))/(rainncv(i)+1.e-12)
        if(fallsum.gt.0.)sr(i)=(tstepsnow(i) + tstepgraup(i))/(rainncv(i)+1.e-12)
      enddo
!
!---------------------------------------------------------------
! pimlt: instantaneous melting of cloud ice [HL A47] [RH83 A28]
!       (T>T0: I->C)
!---------------------------------------------------------------
      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)
          xlf = xls-xl(i,k)
          if(supcol.lt.0.) xlf = xlf0
          if(supcol.lt.0.and.qci(i,k,2).gt.0.) then
            qci(i,k,1) = qci(i,k,1) + qci(i,k,2)
            t(i,k) = t(i,k) - xlf/cpm(i,k)*qci(i,k,2)
            qci(i,k,2) = 0.
          endif
!---------------------------------------------------------------
! pihmf: homogeneous freezing of cloud water below -40c [HL A45]
!        (T<-40C: C->I)
!---------------------------------------------------------------
          if(supcol.gt.40..and.qci(i,k,1).gt.0.) then
            qci(i,k,2) = qci(i,k,2) + qci(i,k,1)
            t(i,k) = t(i,k) + xlf/cpm(i,k)*qci(i,k,1)
            qci(i,k,1) = 0.
          endif
!---------------------------------------------------------------
! pihtf: heterogeneous freezing of cloud water [HL A44]
!        (T0>T>-40C: C->I)
!---------------------------------------------------------------
          if(supcol.gt.0..and.qci(i,k,1).gt.qmin) then
!           pfrzdtc = min(pfrz1*(exp(pfrz2*supcol)-1.)                         &
!              *den(i,k)/denr/xncr*qci(i,k,1)**2*dtcld,qci(i,k,1))
            supcolt=min(supcol,50.)
            pfrzdtc = min(pfrz1*(exp(pfrz2*supcolt)-1.)                        &
            *den(i,k)/denr/xncr*qci(i,k,1)*qci(i,k,1)*dtcld,qci(i,k,1))
            qci(i,k,2) = qci(i,k,2) + pfrzdtc
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtc
            qci(i,k,1) = qci(i,k,1)-pfrzdtc
          endif
!---------------------------------------------------------------
! pgfrz: freezing of rain water [HL A20] [LFO 45]
!        (T<T0, R->G)
!---------------------------------------------------------------
          if(supcol.gt.0..and.qrs(i,k,1).gt.0.) then
!           pfrzdtr = min(20.*pi**2*pfrz1*n0r*denr/den(i,k)                    &
!                 *(exp(pfrz2*supcol)-1.)*rslope3(i,k,1)**2                    &
!                 *rslope(i,k,1)*dtcld,qrs(i,k,1))
            temp = rslope3(i,k,1)
            temp = temp*temp*rslope(i,k,1)
            supcolt=min(supcol,50.)
            pfrzdtr = min(20.*(pi*pi)*pfrz1*n0r*denr/den(i,k)                  &
                  *(exp(pfrz2*supcolt)-1.)*temp*dtcld,                         &
                  qrs(i,k,1))
            qrs(i,k,3) = qrs(i,k,3) + pfrzdtr
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtr
            qrs(i,k,1) = qrs(i,k,1)-pfrzdtr
          endif
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     update the slope parameters for microphysics computation
!
      do k = kts, kte
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
        enddo
      enddo
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                     work1,its,ite,kts,kte)
!------------------------------------------------------------------
!     work1:  the thermodynamic term in the denominator associated with
!             heat conduction and vapor diffusion
!             (ry88, y93, h85)
!     work2: parameter associated with the ventilation effects(y93)
!
      do k = kts, kte
        do i = its, ite
          work1(i,k,1) = diffac(xl(i,k),p(i,k),t(i,k),den(i,k),qs(i,k,1))
          work1(i,k,2) = diffac(xls,p(i,k),t(i,k),den(i,k),qs(i,k,2))
          work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
        enddo
      enddo
!
!===============================================================
!
! warm rain processes
!
! - follows the processes in RH83 and LFO except for autoconcersion
!
!===============================================================
!
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
          supsat = max(q(i,k),qmin)-qs(i,k,1)
          satdt = supsat/dtcld
!---------------------------------------------------------------
! praut: auto conversion rate from cloud to rain [HDC 16]
!        (C->R)
!---------------------------------------------------------------
          if(qci(i,k,1).gt.qc0) then
            praut(i,k) = qck1*qci(i,k,1)**(7./3.)
            praut(i,k) = min(praut(i,k),qci(i,k,1)/dtcld)
          endif
!---------------------------------------------------------------
! pracw: accretion of cloud water by rain [HL A40] [LFO 51]
!        (C->R)
!---------------------------------------------------------------
          if(qrs(i,k,1).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            pracw(i,k) = min(pacrr*rslope3(i,k,1)*rslopeb(i,k,1)               &
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!---------------------------------------------------------------
! prevp: evaporation/condensation rate of rain [HDC 14]
!        (V->R or R->V)
!---------------------------------------------------------------
          if(qrs(i,k,1).gt.0.) then
            coeres = rslope2(i,k,1)*sqrt(rslope(i,k,1)*rslopeb(i,k,1))
            prevp(i,k) = (rh(i,k,1)-1.)*(precr1*rslope2(i,k,1)                 &
                         +precr2*work2(i,k)*coeres)/work1(i,k,1)
            if(prevp(i,k).lt.0.) then
              prevp(i,k) = max(prevp(i,k),-qrs(i,k,1)/dtcld)
              prevp(i,k) = max(prevp(i,k),satdt/2)
            else
              prevp(i,k) = min(prevp(i,k),satdt/2)
            endif
          endif
        enddo
      enddo
!
!===============================================================
!
! cold rain processes
!
! - follows the revised ice microphysics processes in HDC
! - the processes same as in RH83 and RH84  and LFO behave
!   following ice crystal hapits defined in HDC, inclduing
!   intercept parameter for snow (n0s), ice crystal number
!   concentration (ni), ice nuclei number concentration
!   (n0i), ice diameter (d)
!
!===============================================================
!
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
          supcol = t0c-t(i,k)
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          supsat = max(q(i,k),qmin)-qs(i,k,2)
          satdt = supsat/dtcld
          ifsat = 0
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
!         xni(i,k) = min(max(5.38e7*(den(i,k)                                  &
!                      *max(qci(i,k,2),qmin))**0.75,1.e3),1.e6)
          temp = (den(i,k)*max(qci(i,k,2),qmin))
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
          eacrs = exp(0.07*(-supcol))
!
          xmi = den(i,k)*qci(i,k,2)/xni(i,k)
          diameter  = min(dicon * sqrt(xmi),dimax)
          vt2i = 1.49e4*diameter**1.31
          vt2r=pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt2s=pvts*rslopeb(i,k,2)*denfac(i,k)
          vt2g=pvtg*rslopeb(i,k,3)*denfac(i,k)
          qsum(i,k) = max( (qrs(i,k,2)+qrs(i,k,3)), 1.E-15)
          if(qsum(i,k) .gt. 1.e-15) then
          vt2ave=(vt2s*qrs(i,k,2)+vt2g*qrs(i,k,3))/(qsum(i,k))
          else
          vt2ave=0.
          endif
          if(supcol.gt.0.and.qci(i,k,2).gt.qmin) then
            if(qrs(i,k,1).gt.qcrmin) then
!-------------------------------------------------------------
! praci: Accretion of cloud ice by rain [HL A15] [LFO 25]
!        (T<T0: I->R)
!-------------------------------------------------------------
              acrfac = 2.*rslope3(i,k,1)+2.*diameter*rslope2(i,k,1)            &
                      +diameter**2*rslope(i,k,1)
              praci(i,k) = pi*qci(i,k,2)*n0r*abs(vt2r-vt2i)*acrfac/4.
              praci(i,k) = min(praci(i,k),qci(i,k,2)/dtcld)
!-------------------------------------------------------------
! piacr: Accretion of rain by cloud ice [HL A19] [LFO 26]
!        (T<T0: R->S or R->G)
!-------------------------------------------------------------
              piacr(i,k) = pi**2*avtr*n0r*denr*xni(i,k)*denfac(i,k)            &
                          *g6pbr*rslope3(i,k,1)*rslope3(i,k,1)                 &
                          *rslopeb(i,k,1)/24./den(i,k)
              piacr(i,k) = min(piacr(i,k),qrs(i,k,1)/dtcld)
            endif
!-------------------------------------------------------------
! psaci: Accretion of cloud ice by snow [HDC 10]
!        (T<T0: I->S)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.qcrmin) then
              acrfac = 2.*rslope3(i,k,2)+2.*diameter*rslope2(i,k,2)            &
                      +diameter**2*rslope(i,k,2)
              psaci(i,k) = pi*qci(i,k,2)*eacrs*n0s*n0sfac(i,k)                 &
                          *abs(vt2ave-vt2i)*acrfac/4.
              psaci(i,k) = min(psaci(i,k),qci(i,k,2)/dtcld)
            endif
!-------------------------------------------------------------
! pgaci: Accretion of cloud ice by graupel [HL A17] [LFO 41]
!        (T<T0: I->G)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.qcrmin) then
              egi = exp(0.07*(-supcol))
              acrfac = 2.*rslope3(i,k,3)+2.*diameter*rslope2(i,k,3)            &
                      +diameter**2*rslope(i,k,3)
              pgaci(i,k) = pi*egi*qci(i,k,2)*n0g*abs(vt2ave-vt2i)*acrfac/4.
              pgaci(i,k) = min(pgaci(i,k),qci(i,k,2)/dtcld)
            endif
          endif
!-------------------------------------------------------------
! psacw: Accretion of cloud water by snow  [HL A7] [LFO 24]
!        (T<T0: C->S, and T>=T0: C->R)
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            psacw(i,k) = min(pacrc*n0sfac(i,k)*rslope3(i,k,2)*rslopeb(i,k,2)   &    
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pgacw: Accretion of cloud water by graupel [HL A6] [LFO 40]
!        (T<T0: C->G, and T>=T0: C->R)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            pgacw(i,k) = min(pacrg*rslope3(i,k,3)*rslopeb(i,k,3)               &
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! paacw: Accretion of cloud water by averaged snow/graupel 
!        (T<T0: C->G or S, and T>=T0: C->R) 
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qrs(i,k,3).gt.qcrmin) then
            paacw(i,k) = (qrs(i,k,2)*psacw(i,k)+qrs(i,k,3)*pgacw(i,k))         & 
                        /(qsum(i,k))
           endif      
!-------------------------------------------------------------
! pracs: Accretion of snow by rain [HL A11] [LFO 27]
!         (T<T0: S->G)
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qrs(i,k,1).gt.qcrmin) then
            if(supcol.gt.0) then
              acrfac = 5.*rslope3(i,k,2)*rslope3(i,k,2)*rslope(i,k,1)          &
                      +2.*rslope3(i,k,2)*rslope2(i,k,2)*rslope2(i,k,1)         &
                      +.5*rslope2(i,k,2)*rslope2(i,k,2)*rslope3(i,k,1)
              pracs(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2r-vt2ave)          &
                          *(dens/den(i,k))*acrfac
              pracs(i,k) = min(pracs(i,k),qrs(i,k,2)/dtcld)
            endif
!-------------------------------------------------------------
! psacr: Accretion of rain by snow [HL A10] [LFO 28]
!         (T<T0:R->S or R->G) (T>=T0: enhance melting of snow)
!-------------------------------------------------------------
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,2)            &
                    +2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,2)           &
                    +.5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,2)
            psacr(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2ave-vt2r)            &
                        *(denr/den(i,k))*acrfac
            psacr(i,k) = min(psacr(i,k),qrs(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pgacr: Accretion of rain by graupel [HL A12] [LFO 42]
!         (T<T0: R->G) (T>=T0: enhance melting of graupel)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qrs(i,k,1).gt.qcrmin) then
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,3)            &
                    +2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,3)           &
                    +.5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,3)
            pgacr(i,k) = pi**2*n0r*n0g*abs(vt2ave-vt2r)*(denr/den(i,k))        &
                        *acrfac
            pgacr(i,k) = min(pgacr(i,k),qrs(i,k,1)/dtcld)
          endif
!
!-------------------------------------------------------------
! pgacs: Accretion of snow by graupel [HL A13] [LFO 29]
!        (S->G): This process is eliminated in V3.0 with the 
!        new combined snow/graupel fall speeds
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qrs(i,k,2).gt.qcrmin) then
            pgacs(i,k) = 0.
          endif
          if(supcol.le.0) then
            xlf = xlf0
!-------------------------------------------------------------
! pseml: Enhanced melting of snow by accretion of water [HL A34]
!        (T>=T0: S->R)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0.)                                               &
              pseml(i,k) = min(max(cliq*supcol*(paacw(i,k)+psacr(i,k))         &
                          /xlf,-qrs(i,k,2)/dtcld),0.)
!-------------------------------------------------------------
! pgeml: Enhanced melting of graupel by accretion of water [HL A24] [RH84 A21-A22]
!        (T>=T0: G->R)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0.)                                               &
              pgeml(i,k) = min(max(cliq*supcol*(paacw(i,k)+pgacr(i,k))         &
                          /xlf,-qrs(i,k,3)/dtcld),0.)
          endif
          if(supcol.gt.0) then
!-------------------------------------------------------------
! pidep: Deposition/Sublimation rate of ice [HDC 9]
!       (T<T0: V->I or I->V)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.and.ifsat.ne.1) then
              pidep(i,k) = 4.*diameter*xni(i,k)*(rh(i,k,2)-1.)/work1(i,k,2)
              supice = satdt-prevp(i,k)
              if(pidep(i,k).lt.0.) then
                pidep(i,k) = max(max(pidep(i,k),satdt/2),supice)
                pidep(i,k) = max(pidep(i,k),-qci(i,k,2)/dtcld)
              else
                pidep(i,k) = min(min(pidep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)).ge.abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! psdep: deposition/sublimation rate of snow [HDC 14]
!        (T<T0: V->S or S->V)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psdep(i,k) = (rh(i,k,2)-1.)*n0sfac(i,k)*(precs1*rslope2(i,k,2)   &    
                           + precs2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)
              if(psdep(i,k).lt.0.) then
                psdep(i,k) = max(psdep(i,k),-qrs(i,k,2)/dtcld)
                psdep(i,k) = max(max(psdep(i,k),satdt/2),supice)
              else
                psdep(i,k) = min(min(psdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)).ge.abs(satdt))          &
                ifsat = 1
            endif
!-------------------------------------------------------------
! pgdep: deposition/sublimation rate of graupel [HL A21] [LFO 46]
!        (T<T0: V->G or G->V)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgdep(i,k) = (rh(i,k,2)-1.)*(precg1*rslope2(i,k,3)               &
                              +precg2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)
              if(pgdep(i,k).lt.0.) then
                pgdep(i,k) = max(pgdep(i,k),-qrs(i,k,3)/dtcld)
                pgdep(i,k) = max(max(pgdep(i,k),satdt/2),supice)
              else
                pgdep(i,k) = min(min(pgdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)+pgdep(i,k)).ge.          &
                abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! pigen: generation(nucleation) of ice from vapor [HL 50] [HDC 7-8]
!       (T<T0: V->I)
!-------------------------------------------------------------
            if(supsat.gt.0.and.ifsat.ne.1) then
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)-pgdep(i,k)
              xni0 = 1.e3*exp(0.1*supcol)
              roqi0 = 4.92e-11*xni0**1.33
              pigen(i,k) = max(0.,(roqi0/den(i,k)-max(qci(i,k,2),0.))/dtcld)
              pigen(i,k) = min(min(pigen(i,k),satdt),supice)
            endif
!
!-------------------------------------------------------------
! psaut: conversion(aggregation) of ice to snow [HDC 12]
!        (T<T0: I->S)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.) then
              qimax = roqimax/den(i,k)
              psaut(i,k) = max(0.,(qci(i,k,2)-qimax)/dtcld)
            endif
!
!-------------------------------------------------------------
! pgaut: conversion(aggregation) of snow to graupel [HL A4] [LFO 37]
!        (T<T0: S->G)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0.) then
              alpha2 = 1.e-3*exp(0.09*(-supcol))
              pgaut(i,k) = min(max(0.,alpha2*(qrs(i,k,2)-qs0)),qrs(i,k,2)/dtcld)
            endif
          endif
!
!-------------------------------------------------------------
! psevp: Evaporation of melting snow [HL A35] [RH83 A27]
!       (T>=T0: S->V)
!-------------------------------------------------------------
          if(supcol.lt.0.) then
            if(qrs(i,k,2).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psevp(i,k) = (rh(i,k,1)-1.)*n0sfac(i,k)*(precs1                  &
                           *rslope2(i,k,2)+precs2*work2(i,k)                   &
                           *coeres)/work1(i,k,1)
              psevp(i,k) = min(max(psevp(i,k),-qrs(i,k,2)/dtcld),0.)
            endif
!-------------------------------------------------------------
! pgevp: Evaporation of melting graupel [HL A25] [RH84 A19]
!       (T>=T0: G->V)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgevp(i,k) = (rh(i,k,1)-1.)*(precg1*rslope2(i,k,3)               &
                         +precg2*work2(i,k)*coeres)/work1(i,k,1)
              pgevp(i,k) = min(max(pgevp(i,k),-qrs(i,k,3)/dtcld),0.)
            endif
          endif
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     check mass conservation of generation terms and feedback to the
!     large scale
!
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
          delta2(i,k)=0.
          delta3(i,k)=0.
          if(qrs(i,k,1).lt.1.e-4.and.qrs(i,k,2).lt.1.e-4) delta2(i,k)=1.
          if(qrs(i,k,1).lt.1.e-4) delta3(i,k)=1.
        enddo
      enddo
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
!
!     cloud water
!
          value = max(qmin,qci(i,k,1))
          source = (praut(i,k)+pracw(i,k)+paacw(i,k)+paacw(i,k))*dtcld
          if(t(i,k).le.t0c) then
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
            endif
          else
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
            endif
          endif
        enddo
      enddo
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
!
!     cloud ice
!
          if(t(i,k).le.t0c) then
            value = max(qmin,qci(i,k,2))
            source = (psaut(i,k)-pigen(i,k)-pidep(i,k)+praci(i,k)+psaci(i,k)   &
                     +pgaci(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              psaut(i,k) = psaut(i,k)*factor
              pigen(i,k) = pigen(i,k)*factor
              pidep(i,k) = pidep(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
            endif
          endif
        enddo
      enddo
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
!
!     rain
!
          value = max(qmin,qrs(i,k,1))
          if(t(i,k).le.t0c) then
            source = (-praut(i,k)-prevp(i,k)-pracw(i,k)+piacr(i,k)+psacr(i,k)  &
                      +pgacr(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
            endif
          else
            source = (-paacw(i,k)-praut(i,k)+pseml(i,k)+pgeml(i,k)-pracw(i,k)  &
                      -paacw(i,k)-prevp(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif
          endif
        enddo
      enddo
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
!
!     snow
!
          if(t(i,k).le.t0c) then
            value = max(qmin,qrs(i,k,2))
            source = -(psdep(i,k)+psaut(i,k)-pgaut(i,k)+paacw(i,k)+piacr(i,k)  &
                     *delta3(i,k)+praci(i,k)*delta3(i,k)                       &
                     -pracs(i,k)*(1.-delta2(i,k))                              &
                     +psacr(i,k)*delta2(i,k)+psaci(i,k)-pgacs(i,k) )*dtcld
            if (source.gt.value) then
              factor = value/source
              psdep(i,k) = psdep(i,k)*factor
              psaut(i,k) = psaut(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif
          else
            value = max(qcrmin,qrs(i,k,2))
            source=(pgacs(i,k)-pseml(i,k)-psevp(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              psevp(i,k) = psevp(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
            endif
          endif
        enddo
      enddo
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
!
!     graupel
!
          if(t(i,k).le.t0c) then
            value = max(qmin,qrs(i,k,3))
            source = -(pgdep(i,k)+pgaut(i,k)                                   &
                     +piacr(i,k)*(1.-delta3(i,k))+praci(i,k)*(1.-delta3(i,k))  &
                     +psacr(i,k)*(1.-delta2(i,k))+pracs(i,k)*(1.-delta2(i,k))  &
                     +pgaci(i,k)+paacw(i,k)+pgacr(i,k)+pgacs(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgdep(i,k) = pgdep(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif
          else
            value = max(qcrmin,qrs(i,k,3))
            source=-(pgacs(i,k)+pgevp(i,k)+pgeml(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              pgevp(i,k) = pgevp(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif
          endif
        enddo
      enddo
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
!
!     update
!
          if(t(i,k).le.t0c) then
            work2(i,k)=-(prevp(i,k)+psdep(i,k)+pgdep(i,k)+pigen(i,k)+pidep(i,k))
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)                 &
                           +paacw(i,k)+paacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)                 &
                           +prevp(i,k)-piacr(i,k)-pgacr(i,k)                   &
                           -psacr(i,k))*dtcld,0.)
            qci(i,k,2) = max(qci(i,k,2)-(psaut(i,k)+praci(i,k)                 &
                           +psaci(i,k)+pgaci(i,k)-pigen(i,k)-pidep(i,k))       &
                           *dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psdep(i,k)+psaut(i,k)+paacw(i,k)      &
                           -pgaut(i,k)+piacr(i,k)*delta3(i,k)                  &
                           +praci(i,k)*delta3(i,k)+psaci(i,k)-pgacs(i,k)       &
                           -pracs(i,k)*(1.-delta2(i,k))+psacr(i,k)*delta2(i,k))&
                           *dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgdep(i,k)+pgaut(i,k)                 &
                           +piacr(i,k)*(1.-delta3(i,k))                        &
                           +praci(i,k)*(1.-delta3(i,k))                        &
                           +psacr(i,k)*(1.-delta2(i,k))                        &
                           +pracs(i,k)*(1.-delta2(i,k))+pgaci(i,k)+paacw(i,k)  &
                           +pgacr(i,k)+pgacs(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xls*(psdep(i,k)+pgdep(i,k)+pidep(i,k)+pigen(i,k))       &
                      -xl(i,k)*prevp(i,k)-xlf*(piacr(i,k)+paacw(i,k)           &
                      +paacw(i,k)+pgacr(i,k)+psacr(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          else
            work2(i,k)=-(prevp(i,k)+psevp(i,k)+pgevp(i,k))
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)                 &
                    +paacw(i,k)+paacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)                 &
                    +prevp(i,k)+paacw(i,k)+paacw(i,k)-pseml(i,k)               &
                    -pgeml(i,k))*dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psevp(i,k)-pgacs(i,k)                 &
                    +pseml(i,k))*dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgacs(i,k)+pgevp(i,k)                 &
                    +pgeml(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xl(i,k)*(prevp(i,k)+psevp(i,k)+pgevp(i,k))              &
                      -xlf*(pseml(i,k)+pgeml(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          endif
        enddo
      enddo
!
! Inline expansion for fpvs
!         qs(i,k,1) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!         qs(i,k,2) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      do k = kts, kte
        do i = its, ite
          tr=ttp/t(i,k)
          qs(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          qs(i,k,1) = min(qs(i,k,1),0.99*p(i,k))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          tr=ttp/t(i,k)
          if(t(i,k).lt.ttp) then
            qs(i,k,2)=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
          else
            qs(i,k,2)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          endif
          qs(i,k,2) = min(qs(i,k,2),0.99*p(i,k))
          qs(i,k,2) = ep2 * qs(i,k,2) / (p(i,k) - qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
        enddo
      enddo
!
!----------------------------------------------------------------
!  pcond: condensational/evaporational rate of cloud water [HL A46] [RH83 A6]
!     if there exists additional water vapor condensated/if
!     evaporation of cloud water is not enough to remove subsaturation
!
      do k = kts, kte
        do i = its, ite
          work1(i,k,1) = conden(t(i,k),q(i,k),qs(i,k,1),xl(i,k),cpm(i,k))
          work2(i,k) = qci(i,k,1)+work1(i,k,1)
          pcond(i,k) = min(max(work1(i,k,1)/dtcld,0.),max(q(i,k),0.)/dtcld)
          if(qci(i,k,1).gt.0..and.work1(i,k,1).lt.0.)                          &
            pcond(i,k) = max(work1(i,k,1),-qci(i,k,1))/dtcld
          q(i,k) = q(i,k)-pcond(i,k)*dtcld
          qci(i,k,1) = max(qci(i,k,1)+pcond(i,k)*dtcld,0.)
          t(i,k) = t(i,k)+pcond(i,k)*xl(i,k)/cpm(i,k)*dtcld
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     padding for small values
!
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
          if(qci(i,k,1).le.qmin) qci(i,k,1) = 0.0
          if(qci(i,k,2).le.qmin) qci(i,k,2) = 0.0
        enddo
      enddo
      enddo                  ! big loops
  END SUBROUTINE wsm62d

  SUBROUTINE firstTouch(fEdge, dvEdge, dcEdge, invDcEdge, invDvEdge,           &
    invAreaCell, invAreaTriangle, meshScalingDel2, meshScalingDel4,            &
    weightsOnEdge, zgrid, rho_edge, rho_zz, ru, u, v, tend_u, divergence,      &
    vorticity, ke, pv_edge, theta_m, rw, tend_rho, rt_diabatic_tend,           &
    tend_theta, tend_w, w, cqw, rb, rr, pp, pressure_b, zz, zxu, cqu,          &
    h_divergence, kdiff, edgesOnCell_sign, edgesOnVertex_sign, rw_save,        &
    ru_save, theta_m_save, exner, pp_fwd, rr_save, scalars, tend_u_euler,      &
    tend_w_euler, tend_theta_euler, deriv_two, cellsOnEdge, verticesOnEdge,    &
    edgesOnCell, edgesOnEdge, cellsOnCell, edgesOnVertex, nEdgesOnCell,        &
    nEdgesOnEdge, latCell, latEdge, angleEdge, u_init, advCellsForEdge,        &
    nAdvCellsForEdge, adv_coefs, adv_coefs_3rd, rdzu, rdzw, fzm, fzp, qv_init, &
    t_init, pzp, pzm, ur_cell, vr_cell, defc_a, defc_b, tend_w_pgf,            &
    tend_w_buoy)

  IMPLICIT NONE
!
! Mimic memory access patterns of full MPAS atm dynamics while setting arrays 
! to zero.  
!
  real (kind=RKIND), dimension(:) :: fEdge
  real (kind=RKIND), dimension(:) :: dvEdge
  real (kind=RKIND), dimension(:) :: dcEdge
  real (kind=RKIND), dimension(:) :: invDcEdge
  real (kind=RKIND), dimension(:) :: invDvEdge
  real (kind=RKIND), dimension(:) :: invAreaCell
  real (kind=RKIND), dimension(:) :: invAreaTriangle
  real (kind=RKIND), dimension(:) :: meshScalingDel2
  real (kind=RKIND), dimension(:) :: meshScalingDel4
  real (kind=RKIND), dimension(:,:) :: weightsOnEdge
  real (kind=RKIND), dimension(:,:) :: zgrid
  real (kind=RKIND), dimension(:,:) :: rho_edge
  real (kind=RKIND), dimension(:,:) :: rho_zz
  real (kind=RKIND), dimension(:,:) :: ru
  real (kind=RKIND), dimension(:,:) :: u
  real (kind=RKIND), dimension(:,:) :: v
  real (kind=RKIND), dimension(:,:) :: tend_u
  real (kind=RKIND), dimension(:,:) :: divergence
  real (kind=RKIND), dimension(:,:) :: vorticity
  real (kind=RKIND), dimension(:,:) :: ke
  real (kind=RKIND), dimension(:,:) :: pv_edge
  real (kind=RKIND), dimension(:,:) :: theta_m
  real (kind=RKIND), dimension(:,:) :: rw
  real (kind=RKIND), dimension(:,:) :: tend_rho
  real (kind=RKIND), dimension(:,:) :: rt_diabatic_tend
  real (kind=RKIND), dimension(:,:) :: tend_theta
  real (kind=RKIND), dimension(:,:) :: tend_w
  real (kind=RKIND), dimension(:,:) :: w
  real (kind=RKIND), dimension(:,:) :: cqw
  real (kind=RKIND), dimension(:,:) :: rb
  real (kind=RKIND), dimension(:,:) :: rr
  real (kind=RKIND), dimension(:,:) :: pp
  real (kind=RKIND), dimension(:,:) :: pressure_b
  real (kind=RKIND), dimension(:,:) :: zz
  real (kind=RKIND), dimension(:,:) :: zxu
  real (kind=RKIND), dimension(:,:) :: cqu
  real (kind=RKIND), dimension(:,:) :: h_divergence
  real (kind=RKIND), dimension(:,:) :: kdiff
  real (kind=RKIND), dimension(:,:) :: edgesOnCell_sign
  real (kind=RKIND), dimension(:,:) :: edgesOnVertex_sign
  real (kind=RKIND), dimension(:,:) :: rw_save
  real (kind=RKIND), dimension(:,:) :: ru_save
  real (kind=RKIND), dimension(:,:) :: theta_m_save
  real (kind=RKIND), dimension(:,:) :: exner
  real (kind=RKIND), dimension(:,:) :: pp_fwd
  real (kind=RKIND), dimension(:,:) :: rr_save
  real (kind=RKIND), dimension(:,:,:) :: scalars
  real (kind=RKIND), dimension(:,:) :: tend_u_euler
  real (kind=RKIND), dimension(:,:) :: tend_w_euler
  real (kind=RKIND), dimension(:,:) :: tend_theta_euler
  real (kind=RKIND), dimension(:,:,:) :: deriv_two
  integer, dimension(:,:) :: cellsOnEdge
  integer, dimension(:,:) :: verticesOnEdge
  integer, dimension(:,:) :: edgesOnCell
  integer, dimension(:,:) :: edgesOnEdge
  integer, dimension(:,:) :: cellsOnCell
  integer, dimension(:,:) :: edgesOnVertex
  integer, dimension(:) :: nEdgesOnCell
  integer, dimension(:) :: nEdgesOnEdge
  real (kind=RKIND), dimension(:) :: latCell
  real (kind=RKIND), dimension(:) :: latEdge
  real (kind=RKIND), dimension(:) :: angleEdge
  real (kind=RKIND), dimension(:) :: u_init

  integer, dimension(:,:) :: advCellsForEdge
  integer, dimension(:) :: nAdvCellsForEdge
  real (kind=RKIND), dimension(:,:) :: adv_coefs
  real (kind=RKIND), dimension(:,:) :: adv_coefs_3rd

  real (kind=RKIND), dimension(:) :: rdzu
  real (kind=RKIND), dimension(:) :: rdzw
  real (kind=RKIND), dimension(:) :: fzm
  real (kind=RKIND), dimension(:) :: fzp
  real (kind=RKIND), dimension(:) :: qv_init
  real (kind=RKIND), dimension(:,:) :: t_init

  real (kind=RKIND), dimension(:,:) :: pzp
  real (kind=RKIND), dimension(:,:) :: pzm

  real (kind=RKIND), dimension(:,:) :: ur_cell
  real (kind=RKIND), dimension(:,:) :: vr_cell

  real (kind=RKIND), dimension(:,:) :: defc_a
  real (kind=RKIND), dimension(:,:) :: defc_b

  real (kind=RKIND), dimension(:,:) :: tend_w_pgf
  real (kind=RKIND), dimension(:,:) :: tend_w_buoy

! $$$ TODO:  make this really do first-touch after OpenMP is re-introduced
#ifdef _OPENMP
  print *,'ERROR:  must update firstTouch routine for OpenMP'
  stop
#endif
  fEdge = 0.
  dvEdge = 0.
  dcEdge = 0.
  invDcEdge = 0.
  invDvEdge = 0.
  invAreaCell = 0.
  invAreaTriangle = 0.
  meshScalingDel2 = 0.
  meshScalingDel4 = 0.
  weightsOnEdge = 0.
  zgrid = 0.
  rho_edge = 0.
  rho_zz = 0.
  ru = 0.
  u = 0.
  v = 0.
  tend_u = 0.
  divergence = 0.
  vorticity = 0.
  ke = 0.
  pv_edge = 0.
  theta_m = 0.
  rw = 0.
  tend_rho = 0.
  rt_diabatic_tend = 0.
  tend_theta = 0.
  tend_w = 0.
  w = 0.
  cqw = 0.
  rb = 0.
  rr = 0.
  pp = 0.
  pressure_b = 0.
  zz = 0.
  zxu = 0.
  cqu = 0.
  h_divergence = 0.
  kdiff = 0.
  edgesOnCell_sign = 0.
  edgesOnVertex_sign = 0.
  rw_save = 0.
  ru_save = 0.
  theta_m_save = 0.
  exner = 0.
  pp_fwd = 0.
  rr_save = 0.
  scalars = 0.
  tend_u_euler = 0.
  tend_w_euler = 0.
  tend_theta_euler = 0.
  deriv_two = 0.
  cellsOnEdge = 0
  verticesOnEdge = 0
  edgesOnCell = 0
  edgesOnEdge = 0
  cellsOnCell = 0
  edgesOnVertex = 0
  nEdgesOnCell = 0
  nEdgesOnEdge = 0
  latCell = 0.
  latEdge = 0.
  angleEdge = 0.
  u_init = 0.
  advCellsForEdge = 0
  nAdvCellsForEdge = 0
  adv_coefs = 0.
  adv_coefs_3rd = 0.
  rdzu = 0.
  rdzw = 0.
  fzm = 0.
  fzp = 0.
  qv_init = 0.
  t_init = 0.
  pzp = 0.
  pzm = 0.
  ur_cell = 0.
  vr_cell = 0.
  defc_a = 0.
  defc_b = 0.
  tend_w_pgf = 0.
  tend_w_buoy = 0.

  END SUBROUTINE firstTouch

END MODULE module_mpas_atm_dyn




  PROGRAM atm_compute_dyn_tend_work_kernel

  USE module_mpas_atm_dyn

  IMPLICIT NONE

! timers
#include <gptl.inc>
  INTEGER :: ret
  REAL*8  :: totaltime
  INTEGER :: handle = 0
$$$ ! get rid of chunk.h and literalk.h in favor of "use mpas_atm_dimensions"

  integer :: nCells, nEdges, nVertices, nVertLevels, nCellsSolve, &
             nEdgesSolve, vertexDegree, maxEdges, maxEdges2, &
             num_scalars, moist_start, moist_end
$$$ ! set up RKIND and StrKIND
  real (kind=RKIND), allocatable, dimension(:) :: fEdge
  real (kind=RKIND), allocatable, dimension(:) :: dvEdge
  real (kind=RKIND), allocatable, dimension(:) :: dcEdge
  real (kind=RKIND), allocatable, dimension(:) :: invDcEdge
  real (kind=RKIND), allocatable, dimension(:) :: invDvEdge
  real (kind=RKIND), allocatable, dimension(:) :: invAreaCell
  real (kind=RKIND), allocatable, dimension(:) :: invAreaTriangle
  real (kind=RKIND), allocatable, dimension(:) :: meshScalingDel2
  real (kind=RKIND), allocatable, dimension(:) :: meshScalingDel4
  real (kind=RKIND), allocatable, dimension(:,:) :: weightsOnEdge
  real (kind=RKIND), allocatable, dimension(:,:) :: zgrid
  real (kind=RKIND), allocatable, dimension(:,:) :: rho_edge
  real (kind=RKIND), allocatable, dimension(:,:) :: rho_zz
  real (kind=RKIND), allocatable, dimension(:,:) :: ru
  real (kind=RKIND), allocatable, dimension(:,:) :: u
  real (kind=RKIND), allocatable, dimension(:,:) :: v
  real (kind=RKIND), allocatable, dimension(:,:) :: tend_u
  real (kind=RKIND), allocatable, dimension(:,:) :: divergence
  real (kind=RKIND), allocatable, dimension(:,:) :: vorticity
  real (kind=RKIND), allocatable, dimension(:,:) :: ke
  real (kind=RKIND), allocatable, dimension(:,:) :: pv_edge
  real (kind=RKIND), allocatable, dimension(:,:) :: theta_m
  real (kind=RKIND), allocatable, dimension(:,:) :: rw
  real (kind=RKIND), allocatable, dimension(:,:) :: tend_rho
  real (kind=RKIND), allocatable, dimension(:,:) :: rt_diabatic_tend
  real (kind=RKIND), allocatable, dimension(:,:) :: tend_theta
  real (kind=RKIND), allocatable, dimension(:,:) :: tend_w
  real (kind=RKIND), allocatable, dimension(:,:) :: w
  real (kind=RKIND), allocatable, dimension(:,:) :: cqw
  real (kind=RKIND), allocatable, dimension(:,:) :: rb
  real (kind=RKIND), allocatable, dimension(:,:) :: rr
  real (kind=RKIND), allocatable, dimension(:,:) :: pp
  real (kind=RKIND), allocatable, dimension(:,:) :: pressure_b
  real (kind=RKIND), allocatable, dimension(:,:) :: zz
  real (kind=RKIND), allocatable, dimension(:,:) :: zxu
  real (kind=RKIND), allocatable, dimension(:,:) :: cqu
  real (kind=RKIND), allocatable, dimension(:,:) :: h_divergence
  real (kind=RKIND), allocatable, dimension(:,:) :: kdiff
  real (kind=RKIND), allocatable, dimension(:,:) :: edgesOnCell_sign
  real (kind=RKIND), allocatable, dimension(:,:) :: edgesOnVertex_sign
  real (kind=RKIND), allocatable, dimension(:,:) :: rw_save
  real (kind=RKIND), allocatable, dimension(:,:) :: ru_save
  real (kind=RKIND), allocatable, dimension(:,:) :: theta_m_save
  real (kind=RKIND), allocatable, dimension(:,:) :: exner
  real (kind=RKIND), allocatable, dimension(:,:) :: pp_fwd
  real (kind=RKIND), allocatable, dimension(:,:) :: rr_save
  real (kind=RKIND), allocatable, dimension(:,:,:) :: scalars
  real (kind=RKIND), allocatable, dimension(:,:) :: tend_u_euler
  real (kind=RKIND), allocatable, dimension(:,:) :: tend_w_euler
  real (kind=RKIND), allocatable, dimension(:,:) :: tend_theta_euler
  real (kind=RKIND), allocatable, dimension(:,:,:) :: deriv_two
  integer, allocatable, dimension(:,:) :: cellsOnEdge
  integer, allocatable, dimension(:,:) :: verticesOnEdge
  integer, allocatable, dimension(:,:) :: edgesOnCell
  integer, allocatable, dimension(:,:) :: edgesOnEdge
  integer, allocatable, dimension(:,:) :: cellsOnCell
  integer, allocatable, dimension(:,:) :: edgesOnVertex
  integer, allocatable, dimension(:) :: nEdgesOnCell
  integer, allocatable, dimension(:) :: nEdgesOnEdge
  real (kind=RKIND), allocatable, dimension(:) :: latCell
  real (kind=RKIND), allocatable, dimension(:) :: latEdge
  real (kind=RKIND), allocatable, dimension(:) :: angleEdge
  real (kind=RKIND), allocatable, dimension(:) :: u_init

  integer, allocatable, dimension(:,:) :: advCellsForEdge
  integer, allocatable, dimension(:) :: nAdvCellsForEdge
  real (kind=RKIND), allocatable, dimension(:,:) :: adv_coefs
  real (kind=RKIND), allocatable, dimension(:,:) :: adv_coefs_3rd

  real (kind=RKIND), allocatable, dimension(:) :: rdzu
  real (kind=RKIND), allocatable, dimension(:) :: rdzw
  real (kind=RKIND), allocatable, dimension(:) :: fzm
  real (kind=RKIND), allocatable, dimension(:) :: fzp
  real (kind=RKIND), allocatable, dimension(:) :: qv_init
  real (kind=RKIND), allocatable, dimension(:,:) :: t_init

  real (kind=RKIND), allocatable, dimension(:,:) :: pzp
  real (kind=RKIND), allocatable, dimension(:,:) :: pzm

  real (kind=RKIND) :: cf1, cf2, cf3
  real (kind=RKIND) :: r_earth
  real (kind=RKIND), allocatable, dimension(:,:) :: ur_cell
  real (kind=RKIND), allocatable, dimension(:,:) :: vr_cell

  real (kind=RKIND), allocatable, dimension(:,:) :: defc_a
  real (kind=RKIND), allocatable, dimension(:,:) :: defc_b

  real (kind=RKIND), allocatable, dimension(:,:) :: tend_w_pgf
  real (kind=RKIND), allocatable, dimension(:,:) :: tend_w_buoy

  real (kind=RKIND) :: coef_3rd_order, c_s, smdiv
  logical :: config_mix_full
  character (len=StrKIND) :: config_horiz_mixing
  real (kind=RKIND) :: config_del4u_div_factor
  real (kind=RKIND) :: config_h_theta_eddy_visc4
  real (kind=RKIND) :: config_h_mom_eddy_visc4
  real (kind=RKIND) :: config_visc4_2dsmag
  real (kind=RKIND) :: config_len_disp

  integer, intent(in) :: rk_step
  real (kind=RKIND), intent(in) :: dt

  integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
  integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd

! LOCAL VAR
$$$ ! clean up
  INTEGER :: i,j,k
  INTEGER :: ios, kunit
  CHARACTER(LEN=64) :: fn  ! file name

! EXTERNAL FUNCTIONS
#ifdef _OPENMP
  integer, external :: omp_get_max_threads
#endif

  call gptlprocess_namelist ('GPTLnamelist', 77, ret)
  ret = gptlinitialize ()
  ret = gptlstart('Total')

  ! TODO:  Read only INTENT(IN) variables and write only INTENT(OUT) 
  ! TODO:  variables. For now, "best is enemy of good enough".  

  ! read input data
  PRINT *,'atm_compute_dyn_tend_work():  read input state'
  fn = 'atm_compute_dyn_tend_work_in.dat'
  kunit=31
  open (kunit,file=trim(fn),form="unformatted",action='read', &
        iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR: failed to open input file ',trim(fn), &
               ' . stopping'
    stop
  endif

  ! read array dimensions
  read(kunit) nCells
  read(kunit) nEdges
  read(kunit) nVertices
  read(kunit) nVertLevels
  read(kunit) nCellsSolve
  read(kunit) nEdgesSolve
  read(kunit) vertexDegree
  read(kunit) maxEdges
  read(kunit) maxEdges2
  read(kunit) num_scalars
  read(kunit) moist_start
  read(kunit) moist_end

  ! allocate arrays
  allocate(fEdge(nEdges+1))
  allocate(dvEdge(nEdges+1)) 
  allocate(dcEdge(nEdges+1))
  allocate(invDcEdge(nEdges+1))
  allocate(invDvEdge(nEdges+1))
  allocate(invAreaCell(nCells+1))
  allocate(invAreaTriangle(nVertices+1))
  allocate(meshScalingDel2(nEdges+1))
  allocate(meshScalingDel4(nEdges+1))
  allocate(weightsOnEdge(maxEdges2,nEdges+1))
  allocate(zgrid(nVertLevels+1,nCells+1))
  allocate(rho_edge(nVertLevels,nEdges+1))
  allocate(rho_zz(nVertLevels,nCells+1))
  allocate(ru(nVertLevels,nEdges+1))
  allocate(u(nVertLevels,nEdges+1))
  allocate(v(nVertLevels,nEdges+1))
  allocate(tend_u(nVertLevels,nEdges+1))
  allocate(divergence(nVertLevels,nCells+1))
  allocate(vorticity(nVertLevels,nVertices+1))
  allocate(ke(nVertLevels,nCells+1))
  allocate(pv_edge(nVertLevels,nEdges+1))
  allocate(theta_m(nVertLevels,nCells+1))
  allocate(rw(nVertLevels+1,nCells+1))
  allocate(tend_rho(nVertLevels,nCells+1))
  allocate(rt_diabatic_tend(nVertLevels,nCells+1))
  allocate(tend_theta(nVertLevels,nCells+1))
  allocate(tend_w(nVertLevels+1,nCells+1))
  allocate(w(nVertLevels+1,nCells+1))
  allocate(cqw(nVertLevels,nCells+1))
  allocate(rb(nVertLevels,nCells+1))
  allocate(rr(nVertLevels,nCells+1))
  allocate(pp(nVertLevels,nCells+1))
  allocate(pressure_b(nVertLevels,nCells+1))
  allocate(zz(nVertLevels,nCells+1))
  allocate(zxu(nVertLevels,nEdges+1))
  allocate(cqu(nVertLevels,nEdges+1))
  allocate(h_divergence(nVertLevels,nCells+1))
  allocate(kdiff(nVertLevels,nCells+1))
  allocate(edgesOnCell_sign(maxEdges,nCells+1))
  allocate(edgesOnVertex_sign(vertexDegree,nVertices+1))
  allocate(rw_save(nVertLevels+1,nCells+1))
  allocate(ru_save(nVertLevels,nEdges+1))
  allocate(theta_m_save(nVertLevels,nCells+1))
  allocate(exner(nVertLevels,nCells+1))
  allocate(pp_fwd(nVertLevels,nCells+1))
  allocate(rr_save(nVertLevels,nCells+1))
  allocate(scalars(num_scalars,nVertLevels,nCells+1))
  allocate(tend_u_euler(nVertLevels,nEdges+1))
  allocate(tend_w_euler(nVertLevels+1,nCells+1))
  allocate(tend_theta_euler(nVertLevels,nCells+1))
  allocate(deriv_two(15,2,nEdges+1))
  allocate(cellsOnEdge(2,nEdges+1))
  allocate(verticesOnEdge(2,nEdges+1))
  allocate(edgesOnCell(maxEdges,nCells+1))
  allocate(edgesOnEdge(maxEdges2,nEdges+1))
  allocate(cellsOnCell(maxEdges,nCells+1))
  allocate(edgesOnVertex(vertexDegree,nVertices+1))
  allocate(nEdgesOnCell(nCells+1))
  allocate(nEdgesOnEdge(nEdges+1))
  allocate(latCell(nCells+1))
  allocate(latEdge(nEdges+1))
  allocate(angleEdge(nEdges+1))
  allocate(u_init(nVertLevels))

  allocate(advCellsForEdge(15,nEdges+1))
  allocate(nAdvCellsForEdge(nEdges+1))
  allocate(adv_coefs(15,nEdges+1))
  allocate(adv_coefs_3rd(15,nEdges+1))

  allocate(rdzu(nVertLevels))
  allocate(rdzw(nVertLevels))
  allocate(fzm(nVertLevels))
  allocate(fzp(nVertLevels))
  allocate(qv_init(nVertLevels))
  allocate(t_init(nVertLevels,nCells+1))

  allocate(pzp(nVertLevels,nCells+1))
  allocate(pzm(nVertLevels,nCells+1))

  allocate(ur_cell(nVertLevels,nCells+1))
  allocate(vr_cell(nVertLevels,nCells+1))

  allocate(defc_a(maxEdges,nCells+1))
  allocate(defc_b(maxEdges,nCells+1))

  allocate(tend_w_pgf(nVertLevels+1,nCells+1))
  allocate(tend_w_buoy(nVertLevels+1,nCells+1))

  ! read arrays
  read(kunit) fEdge
  read(kunit) dvEdge
  read(kunit) dcEdge
  read(kunit) invDcEdge
  read(kunit) invDvEdge
  read(kunit) invAreaCell
  read(kunit) invAreaTriangle
  read(kunit) meshScalingDel2
  read(kunit) meshScalingDel4
  read(kunit) weightsOnEdge
  read(kunit) zgrid
  read(kunit) rho_edge
  read(kunit) rho_zz
  read(kunit) ru
  read(kunit) u
  read(kunit) v
  read(kunit) tend_u
  read(kunit) divergence
  read(kunit) vorticity
  read(kunit) ke
  read(kunit) pv_edge
  read(kunit) theta_m
  read(kunit) rw
  read(kunit) tend_rho
  read(kunit) rt_diabatic_tend
  read(kunit) tend_theta
  read(kunit) tend_w
  read(kunit) w
  read(kunit) cqw
  read(kunit) rb
  read(kunit) rr
  read(kunit) pp
  read(kunit) pressure_b
  read(kunit) zz
  read(kunit) zxu
  read(kunit) cqu
  read(kunit) h_divergence
  read(kunit) kdiff
  read(kunit) edgesOnCell_sign
  read(kunit) edgesOnVertex_sign
  read(kunit) rw_save
  read(kunit) ru_save
  read(kunit) theta_m_save
  read(kunit) exner
  read(kunit) pp_fwd
  read(kunit) rr_save
  read(kunit) scalars
  read(kunit) tend_u_euler
  read(kunit) tend_w_euler
  read(kunit) tend_theta_euler
  read(kunit) deriv_two
  read(kunit) cellsOnEdge
  read(kunit) verticesOnEdge
  read(kunit) edgesOnCell
  read(kunit) edgesOnEdge
  read(kunit) cellsOnCell
  read(kunit) edgesOnVertex
  read(kunit) nEdgesOnCell
  read(kunit) nEdgesOnEdge
  read(kunit) latCell
  read(kunit) latEdge
  read(kunit) angleEdge
  read(kunit) u_init
  read(kunit) advCellsForEdge
  read(kunit) nAdvCellsForEdge
  read(kunit) adv_coefs
  read(kunit) adv_coefs_3rd
  read(kunit) rdzu
  read(kunit) rdzw
  read(kunit) fzm
  read(kunit) fzp
  read(kunit) qv_init
  read(kunit) t_init 
  read(kunit) pzp
  read(kunit) pzm
  read(kunit) cf1
  read(kunit) cf2
  read(kunit) cf3
  read(kunit) r_earth
  read(kunit) ur_cell
  read(kunit) vr_cell
  read(kunit) defc_a
  read(kunit) defc_b
  read(kunit) tend_w_pgf
  read(kunit) tend_w_buoy
  read(kunit) coef_3rd_order
  read(kunit) c_s
  read(kunit) smdiv
  read(kunit) config_mix_full
  read(kunit) config_horiz_mixing
  read(kunit) config_del4u_div_factor
  read(kunit) config_h_theta_eddy_visc4
  read(kunit) config_h_mom_eddy_visc4
  read(kunit) config_visc4_2dsmag
  read(kunit) config_len_disp
  read(kunit) rk_step
  read(kunit) dt
  read(kunit) cellStart
  read(kunit) cellEnd
  read(kunit) vertexStart
  read(kunit) vertexEnd
  read(kunit) edgeStart
  read(kunit) edgeEnd
  read(kunit) cellSolveStart
  read(kunit) cellSolveEnd
  read(kunit) vertexSolveStart
  read(kunit) vertexSolveEnd
  read(kunit) edgeSolveStart
  read(kunit) edgeSolveEnd
  close(kunit)

  ! first-touch for OpenMP
  call firstTouch(fEdge, dvEdge, dcEdge, invDcEdge, invDvEdge,                 &
    invAreaCell, invAreaTriangle, meshScalingDel2, meshScalingDel4,            &
    weightsOnEdge, zgrid, rho_edge, rho_zz, ru, u, v, tend_u, divergence,      &
    vorticity, ke, pv_edge, theta_m, rw, tend_rho, rt_diabatic_tend,           &
    tend_theta, tend_w, w, cqw, rb, rr, pp, pressure_b, zz, zxu, cqu,          &
    h_divergence, kdiff, edgesOnCell_sign, edgesOnVertex_sign, rw_save,        &
    ru_save, theta_m_save, exner, pp_fwd, rr_save, scalars, tend_u_euler,      &
    tend_w_euler, tend_theta_euler, deriv_two, cellsOnEdge, verticesOnEdge,    &
    edgesOnCell, edgesOnEdge, cellsOnCell, edgesOnVertex, nEdgesOnCell,        &
    nEdgesOnEdge, latCell, latEdge, angleEdge, u_init, advCellsForEdge,        &
    nAdvCellsForEdge, adv_coefs, adv_coefs_3rd, rdzu, rdzw, fzm, fzp, qv_init, &
    t_init, pzp, pzm, ur_cell, vr_cell, defc_a, defc_b, tend_w_pgf,            &
    tend_w_buoy)

$$$ ! here...  swap out WSM6 for atm_compute_dyn_tend_work()

  ! minimize timer overhead inside OpenMP loop
  ret = gptlinit_handle ('WSM62D', handle)
  ! call WSM6
  ret = gptlstart('WSM62D+OpenMP')
!$OMP PARALLEL DO &
!$OMP PRIVATE ( j,ret ) &
!$OMP SCHEDULE(runtime)
  do j = jjts,jjte
    ret = gptlstart_handle ('WSM62D', handle)
    CALL wsm62D(t(iits,kts,j), q(iims,kms,j)                  &
               ,qci(iits,kts,1,j), qrs(iits,kts,1,j)          &
               ,den(iims,kms,j)                               &
               ,p(iims,kms,j), delz(iims,kms,j)               &
               ,delt,g, cpd, cpv, rd, rv, t0c                 &
               ,ep1, ep2, qmin                                &
               ,XLS, XLV0, XLF0, den0, denr                   &
               ,cliq,cice,psat                                &
               ,j                                             &
               ,rain(iims,j),rainncv(iims,j)                  &
               ,sr(iims,j)                                    &
               ,iids,iide, jjds,jjde, kds,kde                 &
               ,iims,iime, jjms,jjme, kms,kme                 &
               ,iits,iite, jjts,jjte, kts,kte                 &
               ,snow,snowncv                                  &
               ,graupel,graupelncv                            &
                                                              )
    ret = gptlstop_handle ('WSM62D', handle)
  enddo  ! j loop
!$OMP END PARALLEL DO
   ret = gptlstop('WSM62D+OpenMP')

  ! write output data
  PRINT *,'atm_compute_dyn_tend_work():  write output state'
  fn = 'atm_compute_dyn_tend_work_kernel_out.dat'
  kunit=31
  open (kunit,file=trim(fn),form="unformatted",action='write', &
        iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR: failed to open output file ',trim(fn), &
               ' . stopping'
    stop
  endif
  write(kunit) nCells
  write(kunit) nEdges
  write(kunit) nVertices
  write(kunit) nVertLevels
  write(kunit) nCellsSolve
  write(kunit) nEdgesSolve
  write(kunit) vertexDegree
  write(kunit) maxEdges
  write(kunit) maxEdges2
  write(kunit) num_scalars
  write(kunit) moist_start
  write(kunit) moist_end

  write(kunit) fEdge
  write(kunit) dvEdge
  write(kunit) dcEdge
  write(kunit) invDcEdge
  write(kunit) invDvEdge
  write(kunit) invAreaCell
  write(kunit) invAreaTriangle
  write(kunit) meshScalingDel2
  write(kunit) meshScalingDel4
  write(kunit) weightsOnEdge
  write(kunit) zgrid
  write(kunit) rho_edge
  write(kunit) rho_zz
  write(kunit) ru
  write(kunit) u
  write(kunit) v
  write(kunit) tend_u
  write(kunit) divergence
  write(kunit) vorticity
  write(kunit) ke
  write(kunit) pv_edge
  write(kunit) theta_m
  write(kunit) rw
  write(kunit) tend_rho
  write(kunit) rt_diabatic_tend
  write(kunit) tend_theta
  write(kunit) tend_w
  write(kunit) w
  write(kunit) cqw
  write(kunit) rb
  write(kunit) rr
  write(kunit) pp
  write(kunit) pressure_b
  write(kunit) zz
  write(kunit) zxu
  write(kunit) cqu
  write(kunit) h_divergence
  write(kunit) kdiff
  write(kunit) edgesOnCell_sign
  write(kunit) edgesOnVertex_sign
  write(kunit) rw_save
  write(kunit) ru_save
  write(kunit) theta_m_save
  write(kunit) exner
  write(kunit) pp_fwd
  write(kunit) rr_save
  write(kunit) scalars
  write(kunit) tend_u_euler
  write(kunit) tend_w_euler
  write(kunit) tend_theta_euler
  write(kunit) deriv_two
  write(kunit) cellsOnEdge
  write(kunit) verticesOnEdge
  write(kunit) edgesOnCell
  write(kunit) edgesOnEdge
  write(kunit) cellsOnCell
  write(kunit) edgesOnVertex
  write(kunit) nEdgesOnCell
  write(kunit) nEdgesOnEdge
  write(kunit) latCell
  write(kunit) latEdge
  write(kunit) angleEdge
  write(kunit) u_init
  write(kunit) advCellsForEdge
  write(kunit) nAdvCellsForEdge
  write(kunit) adv_coefs
  write(kunit) adv_coefs_3rd
  write(kunit) rdzu
  write(kunit) rdzw
  write(kunit) fzm
  write(kunit) fzp
  write(kunit) qv_init
  write(kunit) t_init 
  write(kunit) pzp
  write(kunit) pzm
  write(kunit) cf1
  write(kunit) cf2
  write(kunit) cf3
  write(kunit) r_earth
  write(kunit) ur_cell
  write(kunit) vr_cell
  write(kunit) defc_a
  write(kunit) defc_b
  write(kunit) tend_w_pgf
  write(kunit) tend_w_buoy
  write(kunit) coef_3rd_order
  write(kunit) c_s
  write(kunit) smdiv
  write(kunit) config_mix_full
  write(kunit) config_horiz_mixing
  write(kunit) config_del4u_div_factor
  write(kunit) config_h_theta_eddy_visc4
  write(kunit) config_h_mom_eddy_visc4
  write(kunit) config_visc4_2dsmag
  write(kunit) config_len_disp
  write(kunit) rk_step
  write(kunit) dt
  write(kunit) cellStart
  write(kunit) cellEnd
  write(kunit) vertexStart
  write(kunit) vertexEnd
  write(kunit) edgeStart
  write(kunit) edgeEnd
  write(kunit) cellSolveStart
  write(kunit) cellSolveEnd
  write(kunit) vertexSolveStart
  write(kunit) vertexSolveEnd
  write(kunit) edgeSolveStart
  write(kunit) edgeSolveEnd
  close(kunit)

  ! deallocate arrays
  deallocate(fEdge(nEdges+1))
  deallocate(dvEdge(nEdges+1)) 
  deallocate(dcEdge(nEdges+1))
  deallocate(invDcEdge(nEdges+1))
  deallocate(invDvEdge(nEdges+1))
  deallocate(invAreaCell(nCells+1))
  deallocate(invAreaTriangle(nVertices+1))
  deallocate(meshScalingDel2(nEdges+1))
  deallocate(meshScalingDel4(nEdges+1))
  deallocate(weightsOnEdge(maxEdges2,nEdges+1))
  deallocate(zgrid(nVertLevels+1,nCells+1))
  deallocate(rho_edge(nVertLevels,nEdges+1))
  deallocate(rho_zz(nVertLevels,nCells+1))
  deallocate(ru(nVertLevels,nEdges+1))
  deallocate(u(nVertLevels,nEdges+1))
  deallocate(v(nVertLevels,nEdges+1))
  deallocate(tend_u(nVertLevels,nEdges+1))
  deallocate(divergence(nVertLevels,nCells+1))
  deallocate(vorticity(nVertLevels,nVertices+1))
  deallocate(ke(nVertLevels,nCells+1))
  deallocate(pv_edge(nVertLevels,nEdges+1))
  deallocate(theta_m(nVertLevels,nCells+1))
  deallocate(rw(nVertLevels+1,nCells+1))
  deallocate(tend_rho(nVertLevels,nCells+1))
  deallocate(rt_diabatic_tend(nVertLevels,nCells+1))
  deallocate(tend_theta(nVertLevels,nCells+1))
  deallocate(tend_w(nVertLevels+1,nCells+1))
  deallocate(w(nVertLevels+1,nCells+1))
  deallocate(cqw(nVertLevels,nCells+1))
  deallocate(rb(nVertLevels,nCells+1))
  deallocate(rr(nVertLevels,nCells+1))
  deallocate(pp(nVertLevels,nCells+1))
  deallocate(pressure_b(nVertLevels,nCells+1))
  deallocate(zz(nVertLevels,nCells+1))
  deallocate(zxu(nVertLevels,nEdges+1))
  deallocate(cqu(nVertLevels,nEdges+1))
  deallocate(h_divergence(nVertLevels,nCells+1))
  deallocate(kdiff(nVertLevels,nCells+1))
  deallocate(edgesOnCell_sign(maxEdges,nCells+1))
  deallocate(edgesOnVertex_sign(vertexDegree,nVertices+1))
  deallocate(rw_save(nVertLevels+1,nCells+1))
  deallocate(ru_save(nVertLevels,nEdges+1))
  deallocate(theta_m_save(nVertLevels,nCells+1))
  deallocate(exner(nVertLevels,nCells+1))
  deallocate(pp_fwd(nVertLevels,nCells+1))
  deallocate(rr_save(nVertLevels,nCells+1))
  deallocate(scalars(num_scalars,nVertLevels,nCells+1))
  deallocate(tend_u_euler(nVertLevels,nEdges+1))
  deallocate(tend_w_euler(nVertLevels+1,nCells+1))
  deallocate(tend_theta_euler(nVertLevels,nCells+1))
  deallocate(deriv_two(15,2,nEdges+1))
  deallocate(cellsOnEdge(2,nEdges+1))
  deallocate(verticesOnEdge(2,nEdges+1))
  deallocate(edgesOnCell(maxEdges,nCells+1))
  deallocate(edgesOnEdge(maxEdges2,nEdges+1))
  deallocate(cellsOnCell(maxEdges,nCells+1))
  deallocate(edgesOnVertex(vertexDegree,nVertices+1))
  deallocate(nEdgesOnCell(nCells+1))
  deallocate(nEdgesOnEdge(nEdges+1))
  deallocate(latCell(nCells+1))
  deallocate(latEdge(nEdges+1))
  deallocate(angleEdge(nEdges+1))
  deallocate(u_init(nVertLevels))

  deallocate(advCellsForEdge(15,nEdges+1))
  deallocate(nAdvCellsForEdge(nEdges+1))
  deallocate(adv_coefs(15,nEdges+1))
  deallocate(adv_coefs_3rd(15,nEdges+1))

  deallocate(rdzu(nVertLevels))
  deallocate(rdzw(nVertLevels))
  deallocate(fzm(nVertLevels))
  deallocate(fzp(nVertLevels))
  deallocate(qv_init(nVertLevels))
  deallocate(t_init(nVertLevels,nCells+1))

  deallocate(pzp(nVertLevels,nCells+1))
  deallocate(pzm(nVertLevels,nCells+1))

  deallocate(ur_cell(nVertLevels,nCells+1))
  deallocate(vr_cell(nVertLevels,nCells+1))

  deallocate(defc_a(maxEdges,nCells+1))
  deallocate(defc_b(maxEdges,nCells+1))

  deallocate(tend_w_pgf(nVertLevels+1,nCells+1))
  deallocate(tend_w_buoy(nVertLevels+1,nCells+1))

  ret = gptlstop('Total')
  ret = gptlget_wallclock ('Total', 0, totaltime)  ! The "0" is thread number
  print*,''
  print*,'Total time =' , totaltime

  ! print timing info
  ret = gptlpr (0)
  ret = gptlpr_summary (0)

  END PROGRAM atm_compute_dyn_tend_work_kernel

