    module explicit_lib
    use properties
    implicit none
    
    contains
    
!-----------------------------------------------------------------------
!****************************** RK-Step *******************************
!-----------------------------------------------------------------------
    subroutine rk_step(stage, time, fintime, err_prev)
    integer, intent(inout) :: stage
    real(dp), intent(in) :: time, fintime
    real(dp), intent(inout) :: err_prev
    real(dp) :: error_ni(nr, nz_loc), error_ne(nr, nz_loc), &
                error_nt(nr, nz_loc), normerr_ni, normerr_ne, &
                normerr_nt, err_curr, scaling_factor, &
                abs_tol = 1e-4, rel_tol = 1e-4
    
    ! merson 4("5") adaptive time-stepping
    if (stage == 1) then
        ni_mi = ni_org + ki(1,:,:)*deltime/3.0_dp
        ne_mi = ne_org + ke(1,:,:)*deltime/3.0_dp
        nt_mi = nt_org + kt(1,:,:)*deltime/3.0_dp
        
    else if (stage == 2) then
        
        ni_mi = ni_org + deltime * ( &
                     ki(1,:,:) / 6.0_dp + ki(2,:,:)/6.0_dp )
        ne_mi = ne_org + deltime * ( &
                     ke(1,:,:) / 6.0_dp + ke(2,:,:)/6.0_dp )
        nt_mi = nt_org + deltime * ( &
                     kt(1,:,:) / 6.0_dp + kt(2,:,:)/6.0_dp )
        
    else if (stage == 3) then
        ni_mi = ni_org + deltime * ( &
                     ki(1,:,:) / 8.0_dp + ki(3,:,:) * 3.0_dp/8.0_dp )
        ne_mi = ne_org + deltime * ( &
                     ke(1,:,:) / 8.0_dp + ke(3,:,:) * 3.0_dp/8.0_dp )
        nt_mi = nt_org + deltime * ( &
                     kt(1,:,:) / 8.0_dp + kt(3,:,:) * 3.0_dp/8.0_dp )
            
    else if (stage == 4) then
        ni_mi = ni_org + deltime * ( ki(1,:,:) / 2.0_dp &
                     - ki(3,:,:) * 3.0_dp / 2.0_dp + ki(4,:,:) * 2.0_dp )
        ne_mi = ne_org + deltime * ( ke(1,:,:) / 2.0_dp &
                     - ke(3,:,:) * 3.0_dp / 2.0_dp + ke(4,:,:) * 2.0_dp )
        nt_mi = nt_org + deltime * ( kt(1,:,:) / 2.0_dp &
                     - kt(3,:,:) * 3.0_dp / 2.0_dp + kt(4,:,:) * 2.0_dp )
    else
        ni_pl = ni_org + deltime * ( ki(1,:,:) / 6.0_dp &
                    + ki(4,:,:) * 2.0_dp / 3.0_dp + ki(5,:,:) / 6.0_dp )
        ne_pl = ne_org + deltime * ( ke(1,:,:) / 6.0_dp &
                    + ke(4,:,:) * 2.0_dp / 3.0_dp + ke(5,:,:) / 6.0_dp )
        nt_pl = nt_org + deltime * ( kt(1,:,:) / 6.0_dp &
                    + kt(4,:,:) * 2.0_dp / 3.0_dp + kt(5,:,:) / 6.0_dp )
        
        
        error_ni = abs(deltime *( ki(1,:,:)*2.0_dp/30.0_dp &
                    - ki(3,:,:)*3.0_dp/10.0_dp + ki(4,:,:)*4.0_dp/15.0_dp &
                    - ki(5,:,:)/30.0_dp ))
        error_ne = abs(deltime *( ke(1,:,:)*2.0_dp/30.0_dp &
                    - ke(3,:,:)*3.0_dp/10.0_dp + ke(4,:,:)*4.0_dp/15.0_dp &
                    - ke(5,:,:)/30.0_dp ))
        error_nt = abs(deltime *( kt(1,:,:)*2.0_dp/30.0_dp &
                    - kt(3,:,:)*3.0_dp/10.0_dp + kt(4,:,:)*4.0_dp/15.0_dp &
                    - kt(5,:,:)/30.0_dp ))
        
        normerr_ni = maxval(error_ni/(abs_tol+rel_tol*abs(ni_org)))
        normerr_ne = maxval(error_ne/(abs_tol+rel_tol*abs(ne_org)))
        normerr_nt = maxval(error_nt/(abs_tol+rel_tol*abs(nt_org)))
        
        err_curr = (normerr_ni**2. + normerr_ne**2. + normerr_nt**2.)**0.5
        
        scaling_factor = 0.8_dp * err_curr**(-0.7_dp / 4._dp) &
                                * err_prev**( 0.4_dp / 4._dp)
        scaling_factor = min(2.5,max(0.3,scaling_factor))
        
        deltime = scaling_factor*deltime
        deltime = max(min(fintime - time, deltime),1.d-16)
        
        if (deltime == 1.d-16) then
            write(*,*) 'minimum timestep reached; finishing simulation'
            write(*,'(es10.2)') err_curr
            call petscfinalize(ierr)
            stop
        end if
        
        if (err_curr .le. 1.0_dp) then
            err_prev = err_curr
        else
            stage = 0
        end if
    end if
    
    end subroutine

!-----------------------------------------------------------------------
!****************************** Get Flux ******************************
!-----------------------------------------------------------------------
    subroutine get_flux(flux, dphi, dh, q, mu, D, n_pl, n_mi)
    real(dp), intent(inout) :: flux
    integer, intent(in) :: q
    real(dp), intent(in) :: dphi, dh, mu, D, n_pl, n_mi
    real(dp) :: tol
    
    tol = 1e-12_dp
    
    if (abs(dphi) < tol) then
         
         flux = D * (n_pl - n_mi) / dh
         
    else if (dphi < 0) then
    
         flux = -q * mu * (dphi / dh) * (n_pl - n_mi * exp(mu * dphi / D)) / &
                (1.0_dp - exp(mu * dphi / D))
    else
    
         flux = -q * mui * (dphi / dh) * (n_pl * exp(-mu * dphi / D) &
                -n_mi) / (exp(-mu * dphi / D) - 1.0_dp)
    
    end if
    end subroutine get_flux

!-----------------------------------------------------------------------
!***************************** Continuity *****************************
!-----------------------------------------------------------------------
    subroutine continuity(i, j, side_r, side_z, stage)
    use properties
    integer, intent(in):: i, j, side_r, side_z, stage
    real(dp):: dz_pl, dz_mi, dr_pl, dr_mi, &
               dphi, Ez, Er, Te, ve, anode_coef, flux_coef, &
               fluxi_pl, fluxi_mi, fluxi_r, fluxi_z, &
               fluxe_pl, fluxe_mi, fluxe_r, fluxe_z, &
               fluxt_pl, fluxt_mi, fluxt_r, fluxt_z, &
               mue_pl, mue_mi, De_pl, De_mi, &
               mut_pl, mut_mi, Dt_pl, Dt_mi, &
               dfluxi_dz = 0, dfluxe_dz = 0, dfluxt_dz, &
               dfluxi_dr = 0, dfluxe_dr = 0, dfluxt_dr, &
               term_sie, term_st1, term_st2, term_st3
    
    Te = get_Te( nt_mi(i,j), ne_mi(i,j) )

    ! z-dir center
    if (nz > 1) then
        if (side_z == 0) then
            ! grid spacing
            dz_pl = z(j+1) - z(j)
            dz_mi = z(j)   - z(j-1)

            ! rates and coefficients
            mue_pl = 0.5_dp * ( mue(i,j) + mue(i,j+1) )
            mue_mi = 0.5_dp * ( mue(i,j) + mue(i,j-1) )
            De_pl  = 0.5_dp * ( De( i,j) + De( i,j+1) )
            De_mi  = 0.5_dp * ( De( i,j) + De( i,j-1) )
            
            mut_pl = 0.5_dp * ( mut(i,j) + mut(i,j+1) )
            mut_mi = 0.5_dp * ( mut(i,j) + mut(i,j-1) )
            Dt_pl  = 0.5_dp * ( Dt( i,j) + Dt( i,j+1) )
            Dt_mi  = 0.5_dp * ( Dt( i,j) + Dt( i,j-1) )

            ! evaluate fluxes (in z-direction)
            dphi = phi(i,j+1) - phi(i,j)
            Ez   = -dphi / dz_pl
            
            call get_flux( fluxi_pl, dphi, dz_pl, 1, mui, Di, &
                           ni_mi(i,j), ni_mi(i,j+1) )
            call get_flux( fluxe_pl, dphi, dz_pl, -1, mue_pl, De_pl, &
                           ne_mi(i,j), ne_mi(i,j+1) )
            call get_flux( fluxt_pl, dphi, dz_pl, -1, mut_pl, Dt_pl, &
                           nt_mi(i,j), nt_mi(i,j+1) )

            dphi = phi(i,j) - phi(i,j-1)
            Ez   = (Ez - dphi / dz_mi) / 2.0_dp
            
            call get_flux( fluxi_mi, dphi, dz_mi, 1, mui, Di, &
                           ni_mi(i,j-1), ni_mi(i,j) )
            call get_flux( fluxe_mi, dphi, dz_mi, -1, mue_mi, De_mi, &
                           ne_mi(i,j-1), ne_mi(i,j) )
            call get_flux( fluxt_mi, dphi, dz_mi, -1, mut_mi, Dt_mi, &
                           nt_mi(i,j-1), nt_mi(i,j) )

            fluxi_z = 0.5_dp * (fluxi_pl + fluxi_mi)
            fluxe_z = 0.5_dp * (fluxe_pl + fluxe_mi)
            fluxt_z = 0.5_dp * (fluxt_pl + fluxt_mi)

            ! z-direction flux gradients
            dfluxi_dz = 2.0_dp * (fluxi_pl - fluxi_mi) / (dz_pl + dz_mi)
            dfluxe_dz = 2.0_dp * (fluxe_pl - fluxe_mi) / (dz_pl + dz_mi)
            dfluxt_dz = 2.0_dp * (fluxt_pl - fluxt_mi) / (dz_pl + dz_mi)

        ! z-dir right electrode
        else if (side_z > 0) then
        
            ! grid spacing
            dz_mi = z(j) - z(j-1)
            dphi  = phi(i,j) - phi(i,j-1)
            Ez    = -dphi / dz_mi

            ! rates and coefficients
            mue_mi = 0.5_dp * ( mue(i,j) + mue(i,j-1) )
            De_mi  = 0.5_dp * ( De( i,j) + De( i,j-1) )
            
            mut_mi = 0.5_dp * ( mut(i,j) + mut(i,j-1) )
            Dt_mi  = 0.5_dp * ( Dt( i,j) + Dt( i,j-1) )

            ! evaluate fluxes (in x-direction)
            call get_flux( fluxi_mi, dphi, dz_mi,  1, mui, Di, &
                           ni_mi(i,j-1), ni_mi(i,j) )
            call get_flux( fluxe_mi, dphi, dz_mi, -1, mue_mi, De_mi, &
                           ne_mi(i,j-1), ne_mi(i,j) )
            call get_flux( fluxt_mi, dphi, dz_mi, -1, mut_mi, Dt_mi, &
                           nt_mi(i,j-1), nt_mi(i,j) )

            ! Electrode
            if (side_z == 2) then
            
                ! dot product: n.gamma
                ! n = 1
                if (Ez > 0) then
                ! defines alpha in flux bc
                    flux_coef = 1.0_dp
                else
                    flux_coef = 0.0_dp
                end if

                ! check if node is part of cathode
                if ((phiL > phiR) .and. (side_z == 2)) then
                    anode_coef = 1.0_dp
                else
                    anode_coef = 0.0_dp
                end if

                ve = sqrt( (8.0_dp * echarg * Te) / (pi * me)) * tau0 / l0
                fluxi_z = 0.25_dp * vi * ni_mi(i,j) &
                          + flux_coef * mui * ni_mi(i,j) * Ez
                fluxe_z = 0.25_dp * ve * ne_mi(i,j) &
                          - anode_coef * gamma * fluxi_z
                fluxt_z = 1.0_dp/3 * ve * nt_mi(i,j) &
                          - 4.0_dp/3 * anode_coef * gamma * Te * fluxi_z

                dfluxi_dz = 2.0_dp * (fluxi_z - fluxi_mi) / dz_mi
                dfluxe_dz = 2.0_dp * (fluxe_z - fluxe_mi) / dz_mi
                dfluxt_dz = 2.0_dp * (fluxt_z - fluxt_mi) / dz_mi

            ! Vacuum
            else
                fluxi_z = fluxi_mi
                fluxe_z = fluxe_mi
                fluxt_z = fluxt_mi
                
                dfluxi_dz = 0.0_dp
                dfluxe_dz = 0.0_dp
                dfluxt_dz = 0.0_dp
            end if

        ! z-dir left
        else if (side_z < 0) then
            ! grid spacing
            dz_pl =  z(j+1) - z(j)
            dphi  =  phi(i,j+1) - phi(i,j)
            Ez    = -dphi / dz_pl

            ! rates and coefficients
            mue_pl = 0.5_dp * ( mue(i,j) + mue(i,j+1) )
            de_pl  = 0.5_dp * ( De( i,j) + De( i,j+1) )
            
            mut_pl = 0.5_dp * ( mut(i,j) + mut(i,j+1) )
            Dt_pl  = 0.5_dp * ( Dt( i,j) + Dt( i,j+1) )
            
            ! evaluate fluxes (in x-direction)
            call get_flux( fluxi_pl, dphi, dz_pl, 1, mui, Di, &
                           ni_mi(i,j), ni_mi(i,j+1) )
            call get_flux( fluxe_pl, dphi, dz_pl, -1, mue_pl, De_pl, &
                           ne_mi(i,j), ne_mi(i,j+1) )
            call get_flux( fluxt_pl, dphi, dz_pl, -1, mut_pl, Dt_pl, &
                           nt_mi(i,j), nt_mi(i,j+1) )

            ! Electrode
            if (side_z == -2) then
                ! dot product: n.gamma
                ! n = -1
                if (-Ez > 0) then
                    ! defines alpha in flux bc
                    flux_coef = 1.0_dp
                else
                    flux_coef = 0.0_dp
                end if

                if ( (phiL < phiR) .and. (side_z == -2) ) then
                    anode_coef = 1.0_dp
                else
                    anode_coef = 0.0_dp
                end if

                ve = sqrt( (8.0_dp * echarg * Te) / (pi * me)) * tau0 / l0
                
                fluxi_z = -0.25_dp * vi * ni_mi(i,j) &
                          + flux_coef * mui * ni_mi(i,j) * Ez
                fluxe_z = -0.25_dp * ve * ne_mi(i,j) &
                          - anode_coef * gamma * fluxi_z
                fluxt_z = -1.0_dp/3 * ve * nt_mi(i,j) &
                          - 4.0_dp/3 * anode_coef * gamma * Te * fluxi_z

                dfluxi_dz = 2.0_dp * (fluxi_pl - fluxi_z) / dz_pl
                dfluxe_dz = 2.0_dp * (fluxe_pl - fluxe_z) / dz_pl
                dfluxt_dz = 2.0_dp * (fluxt_pl - fluxt_z) / dz_pl

            ! Vacuum
            else
                fluxi_z = fluxi_pl
                fluxe_z = fluxe_pl
                fluxt_z = fluxt_pl
                
                dfluxi_dz = 0.0_dp
                dfluxe_dz = 0.0_dp
                dfluxt_dz = 0.0_dp
            end if
        end if
    end if

    ! r-dir center
    if (nr > 1) then
        if (side_r == 0) then
            ! grid spacing
            dr_pl = r(i+1) - r(i)
            dr_mi = r(i)   - r(i-1)

            ! rates and coefficients
            mue_pl = 0.5_dp * ( mue(i,j) + mue(i+1,j) )
            mue_mi = 0.5_dp * ( mue(i,j) + mue(i-1,j) )
            De_pl  = 0.5_dp * ( De( i,j) + De( i+1,j) )
            De_mi  = 0.5_dp * ( De( i,j) + De( i-1,j) )
            
            mut_pl = 0.5_dp * ( mut(i,j) + mut(i+1,j) )
            mut_mi = 0.5_dp * ( mut(i,j) + mut(i-1,j) )
            Dt_pl  = 0.5_dp * ( Dt( i,j) + Dt( i+1,j) )
            Dt_mi  = 0.5_dp * ( Dt( i,j) + Dt( i-1,j) )

            ! evaluate fluxes
            dphi = phi(i+1,j) - phi(i,j)
            Er   = -dphi / dr_pl
            
            call get_flux( fluxi_pl, dphi, dr_pl, 1, mui, Di, &
                           ni_mi(i,j), ni_mi(i+1,j) )
            call get_flux( fluxe_pl, dphi, dr_pl, -1, mue_pl, De_pl, &
                           ne_mi(i,j), ne_mi(i+1,j) )
            call get_flux( fluxt_pl, dphi, dr_pl, -1, mut_pl, Dt_pl, &
                           nt_mi(i,j), nt_mi(i+1,j) )

            dphi = phi(i,j) - phi(i-1,j)
            Er   = (Er - dphi / dr_mi) / 2.0_dp
            
            call get_flux( fluxi_mi, dphi, dr_mi, 1, mui, Di, &
                           ni_mi(i-1,j), ni_mi(i,j) )
            call get_flux( fluxe_mi, dphi, dr_mi, -1, mue_mi, De_mi, &
                           ne_mi(i-1,j), ne_mi(i,j) )
            call get_flux( fluxt_mi, dphi, dr_mi, -1, mut_mi, Dt_mi, &
                           nt_mi(i-1,j), nt_mi(i,j) )

            fluxi_r = 0.5_dp * (fluxi_pl + fluxi_mi)
            fluxe_r = 0.5_dp * (fluxe_pl + fluxe_mi)
            fluxt_r = 0.5_dp * (fluxt_pl + fluxt_mi)

            ! flux gradients
            dfluxi_dr = 2.0_dp * (fluxi_pl - fluxi_mi) / (dr_pl + dr_mi) &
                        + fluxi_r / r(i)
            dfluxe_dr = 2.0_dp * (fluxe_pl - fluxe_mi) / (dr_pl + dr_mi) &
                        + fluxe_r / r(i)
            dfluxt_dr = 2.0_dp * (fluxt_pl - fluxt_mi) / (dr_pl + dr_mi) &
                        + fluxt_r / r(i)

        ! r-dir right sdie
        else if (side_r > 0) then
        
            ! grid spacing
            dr_mi =  r(i) - r(i-1)
            dphi  =  phi(i,j) - phi(i-1,j)
            Er    = -dphi / dr_mi

            ! rates and coefficients
            mue_mi = 0.5_dp * ( mue(i,j) + mue(i-1,j) )
            De_mi  = 0.5_dp * ( De( i,j) + De( i-1,j) )
            
            mut_mi = 0.5_dp * ( mut(i,j) + mut(i-1,j) )
            Dt_mi  = 0.5_dp * ( Dt( i,j) + Dt( i-1,j) )

            ! evaluate fluxes
            call get_flux( fluxi_mi, dphi, dr_mi, 1, mui, Di, &
                           ni_mi(i-1,j), ni_mi(i,j) )
            call get_flux( fluxe_mi, dphi, dr_mi, -1, mue_mi, De_mi, &
                           ne_mi(i-1,j), ne_mi(i,j) )
            call get_flux( fluxt_mi, dphi, dr_mi, -1, mut_mi, Dt_mi, &
                           nt_mi(i-1,j), nt_mi(i,j) )

            ! Electrode
            if (side_r == 2) then
            
                ! dot product: n.gamma
                ! n = 1
                if (Er > 0) then
                ! defines alpha in flux bc
                    flux_coef = 1.0_dp
                else
                    flux_coef = 0.0_dp
                end if

                ! check if node is part of cathode
                if ((phiL > phiR) .and. (side_r == 2)) then
                    anode_coef = 1.0_dp
                else
                    anode_coef = 0.0_dp
                end if

                ve = sqrt( (8.0_dp * echarg * Te) / (pi * me)) * tau0 / l0
                fluxi_r = 0.25_dp * vi * ni_mi(i,j) &
                          + flux_coef * mui * ni_mi(i,j) * Er
                fluxe_r = 0.25_dp * ve * ne_mi(i,j) &
                          - anode_coef * gamma * fluxi_r
                fluxt_r = 1.0_dp/3 * ve * nt_mi(i,j) &
                          - 4.0_dp/3 * anode_coef * gamma * Te * fluxi_r

                dfluxi_dr = 2.0_dp * (fluxi_r - fluxi_mi) / dr_mi &
                            + fluxi_r / r(i)
                dfluxe_dr = 2.0_dp * (fluxe_r - fluxe_mi) / dr_mi &
                            + fluxe_r / r(i)
                dfluxt_dr = 2.0_dp * (fluxt_r - fluxt_mi) / dr_mi &
                            + fluxt_r / r(i)

            ! Vacuum
            else
                fluxi_r = 0.0_dp
                fluxe_r = 0.0_dp
                fluxt_r = 0.0_dp
                
                dfluxi_dr = 2.0_dp * (fluxi_r - fluxi_mi) / dr_mi &
                            + fluxi_r / r(i)
                dfluxe_dr = 2.0_dp * (fluxe_r - fluxe_mi) / dr_mi &
                            + fluxe_r / r(i)
                dfluxt_dr = 2.0_dp * (fluxt_r - fluxt_mi) / dr_mi &
                            + fluxt_r / r(i)
            end if

        ! r-dir left
        else if (side_r < 0) then
            ! grid spacing
            dr_pl =  r(i+1) - r(i)
            dphi  =  phi(i+1,j) - phi(i,j)
            Er    = -dphi / dr_pl

            ! rates and coefficients
            mue_pl = 0.5_dp * ( mue(i,j) + mue(i+1,j) )
            de_pl  = 0.5_dp * ( De( i,j) + De( i+1,j) )
            
            mut_pl = 0.5_dp * ( mut(i,j) + mut(i+1,j) )
            Dt_pl  = 0.5_dp * ( Dt( i,j) + Dt( i+1,j) )
            
            ! evaluate fluxes (in x-direction)
            call get_flux( fluxi_pl, dphi, dr_pl, 1, mui, Di, &
                           ni_mi(i,j), ni_mi(i+1,j) )
            call get_flux( fluxe_pl, dphi, dr_pl, -1, mue_pl, De_pl, &
                           ne_mi(i,j), ne_mi(i+1,j) )
            call get_flux( fluxt_pl, dphi, dr_pl, -1, mut_pl, Dt_pl, &
                           nt_mi(i,j), nt_mi(i+1,j) )

            ! Electrode
            if (side_r == -2) then
                ! dot product: n.gamma
                ! n = -1
                if (-Er > 0) then
                    ! defines alpha in flux bc
                    flux_coef = 1.0_dp
                else
                    flux_coef = 0.0_dp
                end if

                if ( (phiL < phiR) .and. (side_r == -2) ) then
                    anode_coef = 1.0_dp
                else
                    anode_coef = 0.0_dp
                end if

                ve = sqrt( (8.0_dp * echarg * Te) / (pi * me)) * tau0 / l0
                
                fluxi_r = -0.25_dp * vi * ni_mi(i,j) &
                          + flux_coef * mui * ni_mi(i,j) * Er
                fluxe_r = -0.25_dp * ve * ne_mi(i,j) &
                          - anode_coef * gamma * fluxi_r
                fluxt_r = -1.0_dp/3 * ve * nt_mi(i,j) &
                          -4.0_dp/3 * anode_coef * gamma * Te * fluxi_r

                dfluxi_dr = (fluxi_pl - fluxi_r) / dr_pl
                dfluxe_dr = (fluxe_pl - fluxe_r) / dr_pl
                dfluxt_dr = (fluxt_pl - fluxt_r) / dr_pl

            ! Vacuum
            else
                fluxi_r = 0.0_dp
                fluxe_r = 0.0_dp
                fluxt_r = 0.0_dp
                
                dfluxi_dr = 4.0_dp * fluxi_pl / dr_pl
                dfluxe_dr = 4.0_dp * fluxe_pl / dr_pl
                dfluxt_dr = 4.0_dp * fluxt_pl / dr_pl
            end if
        end if
    end if

    ! evaluate source terms
    term_sie = k_ir(i,j) * ninf/n0 * ne_mi(i,j) - beta * ni_mi(i,j) * ne_mi(i,j)
   
    ! evaluate expression
    ki(stage,i,j) = -dfluxi_dz - dfluxi_dr + term_sie
    ke(stage,i,j) = -dfluxe_dz - dfluxe_dr + term_sie
    
    ! -e flux_e . E
    term_st1 = -echarg * phi0 / (3.0_dp/2 * kb*T0) &
               * ( fluxe_z * Ez + fluxe_r * Er)

    ! -me/mg nu_e (Te - Tg)
    term_st2 = ne_mi(i,j) * nu(i,j) * me/mi * (Te - Ti/T0)

    ! reactions
    term_st3 = h_ir * k_ir(i,j) * ninf/n0 * ne_mi(i,j) &
               + h_ex * k_ex(i,j) * ninf/n0 * ne_mi(i,j)

    ! evaluate expression
    kt(stage,i,j) = -dfluxt_dz - dfluxt_dr + term_st1 - term_st2 - term_st3
     

    if (ki(stage,i,j) /= ki(stage,i,j)) then
        write(*,*) 'i,j:', i, j
        write(*,*) 'ki(stage,i,j)',ki(stage,i,j)
        write(*,*) 'side_z', side_z, 'side_r', side_r
        call PetscFinalize(ierr)
        stop
    end if
    if (ke(stage,i,j) /= ke(stage,i,j)) then
        write(*,*) 'i,j:', i, j
        write(*,*) 'ke(stage,i,j)',ke(stage,i,j)
        write(*,*) 'side_z', side_z, 'side_r', side_r
        call PetscFinalize(ierr)
        stop
    end if
    if (kt(stage,i,j) /= kt(stage,i,j)) then
        write(*,*) 'i,j:', i, j
        write(*,*) 'ke(stage,i,j)',kt(stage,i,j)
        write(*,*) 'side_z', side_z, 'side_r', side_r
        call PetscFinalize(ierr)
        stop
    end if
    
    end subroutine
    end module
