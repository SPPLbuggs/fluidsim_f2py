    module laplace_lib
    use properties
    use explicit_lib
    implicit none

    contains
!-----------------------------------------------------------------------
!******************************* Laplace ******************************
!-----------------------------------------------------------------------
    subroutine laplace(i, j, side_r, side_z, b_temp)
    
    integer, intent(in):: i, j, side_r, side_z
    real(dp), intent(out):: b_temp
    real(dp):: dz_pl = 0, dz_mi = 0, dr_pl= 0, dr_mi= 0, &
               dphi_dz = 0, dphi_dr = 0, dphi, &
               dfluxi_dz = 0, dfluxe_dz = 0, dfluxi_dr = 0, dfluxe_dr = 0, &
               fluxi_mi, fluxi_pl, fluxe_mi, fluxe_pl, &
               mue_pl, mue_mi, De_pl, De_mi, term_s = 0

    if (nz > 1) then
        ! Z-dir Center
        if (side_z .eq. 0) then
            dz_pl = z(i+1) - z(i)
            dz_mi = z(i)   - z(i-1)

            ! term: d^2(phi)/dx^2 (x-term)
            dphi_dz = 2.0_dp * phi(i+1,j) / ( dz_pl * (dz_pl + dz_mi) ) &
                     -2.0_dp * phi(i,j)   / ( dz_pl * dz_mi ) &
                     +2.0_dp * phi(i-1,j) / ( dz_mi * (dz_pl + dz_mi) )
                
            ! Rates and coefficients
            mue_pl = 0.5_dp * ( mue(i,j) + mue(i+1,j) )
            mue_mi = 0.5_dp * ( mue(i,j) + mue(i-1,j) )
            De_pl  = 0.5_dp * ( De( i,j) + De( i+1,j) )
            De_mi  = 0.5_dp * ( De( i,j) + De( i-1,j) )    
                
            dphi = phi(i+1,j) - phi(i,j)        
            call get_flux( fluxi_pl, dphi, dz_pl, 1, mui, Di, &
                           ni_mi(i,j), ni_mi(i+1,j) )
            call get_flux( fluxe_pl, dphi, dz_pl, -1, mue_pl, De_pl, &
                           ne_mi(i,j), ne_mi(i+1,j) )
                           
            dphi = phi(i,j) - phi(i-1,j)
            call get_flux( fluxi_mi, dphi, dz_mi, 1, mui, Di, &
                           ni_mi(i-1,j), ni_mi(i,j) )
            call get_flux( fluxe_mi, dphi, dz_mi, -1, mue_mi, De_mi, &
                           ne_mi(i-1,j), ne_mi(i,j) )
    
            ! X-direction flux gradients
            dfluxe_dz = 2.0_dp * (fluxe_pl - fluxe_mi) / (dz_pl + dz_mi)
            dfluxi_dz = 2.0_dp * (fluxi_pl - fluxi_mi) / (dz_pl + dz_mi)
    
        ! X-dir left (vacuum)
        else if (side_z < 0) then
            dz_pl = z(i+1) - z(i)

            ! BC is E_x = 0
            dphi_dz = 2.0_dp * (phi(i+1,j) - phi(i,j)) / dz_pl**2.0_dp
    
        ! X-dir right (vacuum) is fixed phi = 0
        else if (side_z > 0) then
            dz_mi = z(i) - z(i-1)
            
            dphi_dz = 2.0_dp * (phi(i-1,j) - phi(i,j)) / dz_mi**2.0_dp
        
        end if
    end if
    
    if (nr > 1) then
        ! r-dir Center
        if (side_r .eq. 0) then
            dr_pl = r(j+1) - r(j)
            dr_mi = r(j)   - r(j-1)

            ! term: d^2(phi)/dr^2 (r-term)
            dphi_dr = 2.0_dp * phi(i,j+1) / ( dr_pl * (dr_pl + dr_mi) ) &
                     -2.0_dp * phi(i,j)   / ( dr_pl * dr_mi ) &
                     +2.0_dp * phi(i,j-1) / ( dr_mi * (dr_pl + dr_mi) ) &
				     + (phi(i,j+1) - phi(i,j-1)) / ( r(j) * (dr_pl + dr_mi) )
                
            ! Rates and coefficients
            mue_pl = 0.5_dp * ( mue(i,j) + mue(i,j+1) )
            mue_mi = 0.5_dp * ( mue(i,j) + mue(i,j-1) )
            De_pl  = 0.5_dp * ( De( i,j) + De( i,j+1) )
            De_mi  = 0.5_dp * ( De( i,j) + De( i,j-1) )   
                
            ! r-direction flux
            dphi = phi(i,j+1) - phi(i,j)        
            call get_flux( fluxi_pl, dphi, dr_pl, 1, mui, Di, &
                           ni_mi(i,j), ni_mi(i,j+1) )
            call get_flux( fluxe_pl, dphi, dr_pl, -1, mue_pl, De_pl, &
                           ne_mi(i,j), ne_mi(i,j+1) )
                           
            dphi = phi(i,j) - phi(i,j-1)
            call get_flux( fluxi_mi, dphi, dr_mi, 1, mui, Di, &
                           ni_mi(i-1,j), ni_mi(i,j) )
            call get_flux( fluxe_mi, dphi, dr_mi, -1, mue_mi, De_mi, &
                           ne_mi(i-1,j), ne_mi(i,j) )
    
            ! X-direction flux gradients
            dfluxe_dr = 2.0_dp * (fluxe_pl - fluxe_mi) / (dr_pl + dr_mi) &
                       +0.5_dp * (fluxe_pl + fluxe_mi) / r(j)
            dfluxi_dr = 2.0_dp * (fluxi_pl - fluxi_mi) / (dr_pl + dr_mi) &
                       +0.5_dp * (fluxi_pl + fluxi_mi) / r(j)
    
        ! r-dir left (vacuum)
        else if (side_r < 0) then
            dr_pl= r(j+1) - r(j)

            ! BC is E_r = 0
            dphi_dr = 2.0_dp * ( phi(i,j+1) - phi(i,j) ) / dr_pl**2.0_dp
                         
            ! rates and coefficients
            mue_pl = 0.5d0 * ( mue(i,j) + mue(i,j+1) )
            De_pl  = 0.5d0 * ( De( i,j) + De( i,j+1) )
            
            ! evaluate fluxes (in x-direction)
            dphi = phi(i,j+1) - phi(i,j)        
            call get_flux( fluxi_pl, dphi, dr_pl, 1, mui, Di, &
                           ni_mi(i,j), ni_mi(i,j+1) )
            call get_flux( fluxe_pl, dphi, dr_pl, -1, mue_pl, De_pl, &
                           ne_mi(i,j), ne_mi(i,j+1) )
                
            dfluxe_dr = 4.0_dp * fluxe_pl / dr_pl
            dfluxi_dr = 4.0_dp * fluxi_pl / dr_pl
    
        ! r-dir right (vacuum) is fixed phi = 0
        else if (side_r > 0) then
            dr_mi = r(j) - r(j-1)
            
            dphi_dr = 2.0_dp * (phi(i,j-1) - phi(i,j)) / dr_mi**2.0_dp
        
        end if
    end if
    
    ! Source term: e/epsilon*(n_i-n_e + dt*(del . flux_e - del . flux_i))
    term_s = L0**2 * n0 * echarg / (eps * phi0) * ( ni_mi(i,j) - ne_mi(i,j) &
             + deltime * ( dfluxe_dz + dfluxe_dr - dfluxi_dz - dfluxi_dr ) )

    b_temp = (dphi_dz + dphi_dr + term_s) * max(dz_pl,dz_mi,dr_pl,dr_mi)
    

    if (term_s /= term_s) then
        write(*,*)
        write(*,*) rank, timestep
        write(*,*) i,j, glob_node(i,j), side_z, side_r
        write(*,111) term_s
        write(*,111) ni_mi(i,j), ne_mi(i,j)
        write(*,111) dfluxe_dz, dfluxe_dr, dfluxi_dz, dfluxi_dr
        write(*,111) fluxe_pl, fluxe_mi
        write(*,111) mue_pl, De_pl, mue_mi, De_mi
        write(*,111) fluxi_mi, dphi, dz_mi, mui, Di, &
                           ni_mi(i-1,j), ni_mi(i,j)
        read(*,*) term_s
    end if

    if (b_temp /= b_temp) then
        write(*,*)
        write(*,*) rank, timestep
        write(*,*) i,j, glob_node(i,j), side_z, side_r
        write(*,111) b_temp
        write(*,111) dphi_dz, dphi_dr, term_s, max(dz_pl,dz_mi,dr_pl,dr_mi)
        read(*,*) term_s
    end if
111 format(10es14.6)
    end subroutine
    
!-----------------------------------------------------------------------
!****************************** Jacobian ******************************
!-----------------------------------------------------------------------

    subroutine jacobian(i_loc, j_loc, side_r, side_z, &
                    b_temp, cols, A_temp)
    integer, intent(in):: i_loc, j_loc, side_r, side_z
    integer, intent(inout):: cols(5)
    real(dp), intent(in):: b_temp
    real(dp), intent(inout):: A_temp(1,5)
    logical:: zero_perturb
    real(dp):: perturb, b_pert, temp
    integer:: I, J, K, width, k_start, k_stop, cols_idx
    integer, dimension(5,2):: stencil

    ! initialize
    cols_idx = 0
    temp = 0
    Perturb = 1e-4_dp
    stencil = 0
    width = 0
    
    k_start = 1
    k_stop  = 5
    if (nz .eq. 1) k_start = 3
    if (nr .eq. 1) k_stop  = 3
    
    do K = k_start, k_stop
        if ((K .eq. 1) .and. (side_z .ge. 0)) then
            if (glob_node(i_loc-1,j_loc) > 0) then
                width = width + 1
                stencil(width,1) = -1
                stencil(width,2) =  0
            end if
        else if ((K .eq. 2) .and. (side_z .le. 0)) then
            if (glob_node(i_loc+1,j_loc) > 0) then
                width = width + 1
                stencil(width,1) = 1
                stencil(width,2) = 0
            end if
        else if (K .eq. 3) then
            width = width + 1
            stencil(width,1) = 0
            stencil(width,2) = 0
        else if ((K .eq. 4) .and. (side_r .le. 0)) then
            if (glob_node(i_loc,j_loc+1) > 0) then
                width = width + 1
                stencil(width,1) = 0
                stencil(width,2) = 1
            end if
        else if ((K .eq. 5) .and. (side_r .ge. 0)) then
            if (glob_node(i_loc,j_loc-1) > 0) then
                width = width + 1
                stencil(width,1) =  0
                stencil(width,2) = -1
            end if
        end if
    end do
    
    DO k=1, width
        I = i_loc + stencil(k,1)
        J = j_loc + stencil(k,2)
        
        zero_perturb = .false.
        cols_idx = cols_idx + 1
        temp = phi(I,J)
        if (Abs(phi(I,J)) > 0) then
            phi(I,J) = phi(I,J) + &
                phi(I,J)*perturb
        else
            zero_perturb = .true.
            phi(I,J) = perturb
        end if
        call laplace(i_loc, j_loc, side_r, side_z, b_pert)
        if (.not. zero_perturb) then
            phi(I,J) = temp
        else
            phi(I,J) = 1
        end if
        
        cols(cols_idx) = glob_node(i_loc + stencil(k,1), &
                                     j_loc + stencil(k,2)) - 1
        A_temp(1,cols_idx) = (b_pert - b_temp)/(phi(I,J)*perturb)
        
        if (Zero_Perturb) then
            phi(I,J) = temp
        end if
    end do
    end subroutine jacobian
    end module laplace_lib

