    module mod
    use properties
    use laplace_lib
    use petsc_lib
    use explicit_lib
    implicit none
    
    !f2py integer:: neqn, nz, nr, nz_loc
    !f2py real(8):: phiL, phiR
    !f2py integer, allocatable :: type_z(:,:), type_r(:,:), glob_idx(:,:)
    !f2py integer, allocatable :: loc_idx(:,:), glob_node(:,:)
    !f2py real(8), allocatable :: z(:), r(:), phi(:,:)
    !f2py real(8), allocatable :: ni_mi(:,:), ni_pl(:,:), ni_org(:,:)
    !f2py real(8), allocatable :: ne_mi(:,:), ne_pl(:,:), ne_org(:,:)
    !f2py real(8), allocatable :: nt_mi(:,:), nt_pl(:,:), nt_org(:,:)
    
    contains
!-----------------------------------------------------------------------
!***************************** Initialize *****************************
!-----------------------------------------------------------------------
    subroutine initialize(rel_err, abs_err, max_iter)
    integer, intent(in) :: max_iter
    real(dp), intent(in) :: rel_err, abs_err
    integer:: i_glob, j_glob, i_loc, j_loc, node, cols(5)
    real(dp):: A_temp(1,5), b_temp, soln
    
    
    call petsc_initialize
    
    allocate( ki(5, nr, nz_loc), ke(5, nr, nz_loc), kt(5, nr, nz_loc), &
              nu(nr, nz_loc), mue(nr, nz_loc), mut(nr, nz_loc), &
              De(nr, nz_loc), Dt(nr,nz_loc), k_ir(nr, nz_loc), &
              k_ex(nr, nz_loc), k_sc(nr, nz_loc), k_si(nr, nz_loc) )
    
    ni_mi  = ni_pl
    ni_org = ni_pl
    ne_mi  = ne_pl
    ne_org = ne_pl
    nt_mi  = nt_pl
    nt_org = nt_pl
    
    call update_coef
    
    ! update boundary conditions
    !call update_bc
    
    ! assemble A and b
    do node = Istart, Iend-1
        i_loc  = loc_idx(node+1, 1)
        j_loc  = loc_idx(node+1, 2)
        i_glob = glob_idx(node+1,1)
        j_glob = glob_idx(node+1,2)
               
        b_temp =  0.0_dp
        A_temp =  0.0_dp
        cols   = -5
        
        call laplace(i_loc, j_loc, type_z(i_glob,j_glob), &
                     type_r(i_glob,j_glob), b_temp)
        call jacobian(i_loc, j_loc, i_glob, j_glob, &
                      type_z(i_glob,j_glob), type_r(i_glob,j_glob), &
                      b_temp, cols, A_temp)
        
        call MatSetValues(A, 1, node, 5, cols, A_temp, insert_values, ierr)
        call VecSetValues(b, 1, node, -b_temp, insert_values, ierr)
    end do
    
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    
    call vecassemblybegin(b, ierr)
    call vecassemblyend(b, ierr)

    ! create linear solver
    call KSPCreate(comm,ksp,ierr)
    call KSPSetOperators(ksp,A,A,ierr)
    call KSPSetFromOptions(ksp,ierr)
    !call KSPSetTolerances(ksp,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL, & 
    !                      PETSC_DEFAULT_REAL,max_iter,ierr)
    
    end subroutine
!-----------------------------------------------------------------------
!******************************** View ********************************
!-----------------------------------------------------------------------
    subroutine view
    
    call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    end subroutine

!-----------------------------------------------------------------------
!******************************** Run *********************************
!-----------------------------------------------------------------------
    subroutine run()
    real(dp):: fintime = 1
    integer  :: timestep = 0, i_loc, j_loc, i_glob, j_glob, node, stage, recv_status
    real(dp) :: time = 0, print_time = 1e-2, b_temp, &
                err_prev = 0, soln = 0
    
    do
        timestep = timestep + 1
        
        ! solve implicit system
        call KSPSolve(ksp,b,x,ierr)
        
        ! update loc implicit variables
        do node = Istart, Iend - 1
            call VecGetValues(x, 1, node, soln, ierr)
            i_loc = loc_idx(node + 1, 1)
            j_loc = loc_idx(node + 1, 2)
            phi(i_loc, j_loc) = phi(i_loc, j_loc) + soln
        end do
        
        ! communicate loc boundary nodes
        if (rank .ne. 0) then
            call MPI_Send(phi(:,2), nr, MPI_REAL8, &
                          rank-1, 16, comm, ierr)
        end if
        if (rank .ne. tasks-1) then
            call MPI_Send(phi(:,nz_loc-1), nr, MPI_REAL8, &
                          rank+1, 16, comm, ierr)
        end if
        
        if (rank .ne. tasks-1) then
            call MPI_Recv(phi(:,nz_loc), nr, MPI_REAL8, &
                          rank-1, 16, comm, recv_status, ierr)
        end if
        if (rank .ne. 0) then
            call MPI_Send(phi(:,1), nr, MPI_REAL8, &
                          rank+1, 16, comm, recv_status, ierr)
        end if
        
        ! solve explicit system
        stage = 0
        do
            stage = stage + 1
            if (stage > 1) exit
            
            do node = Istart, Iend - 1
                i_loc  = loc_idx( node + 1, 1)
                j_loc  = loc_idx( node + 1, 2)
                i_glob = glob_idx(node + 1, 1)
                j_glob = glob_idx(node + 1, 2)
                call continuity(i_loc, j_loc, type_r(i_glob,j_glob), &
                                type_z(i_glob,j_glob), stage)
            end do
            
            !call rk_step(stage, time, fintime, err_prev)
        end do
        
        do node = Istart, Iend - 1
            i_loc  = loc_idx( node + 1, 1)
            j_loc  = loc_idx( node + 1, 2)
            ni_mi = ni_org + deltime * ki(1,i_loc,j_loc)
            ne_mi = ne_org + deltime * ke(1,i_loc,j_loc)
            nt_mi = nt_org + deltime * kt(1,i_loc,j_loc)
        end do
        
        ! communicate loc boundary nodes
        if (rank .ne. 0) then
            call MPI_Send(ni_pl(:,2), nr, MPI_REAL8, &
                          rank-1, 16, comm, ierr)
            call MPI_Send(ne_pl(:,2), nr, MPI_REAL8, &
                          rank-1, 16, comm, ierr)
            call MPI_Send(nt_pl(:,2), nr, MPI_REAL8, &
                          rank-1, 16, comm, ierr)
        end if
        if (rank .ne. tasks-1) then
            call MPI_Send(ni_pl(:,nz_loc-1), nr, MPI_REAL8, &
                          rank+1, 16, comm, ierr)
            call MPI_Send(ne_pl(:,nz_loc-1), nr, MPI_REAL8, &
                          rank+1, 16, comm, ierr)
            call MPI_Send(nt_pl(:,nz_loc-1), nr, MPI_REAL8, &
                          rank+1, 16, comm, ierr)
        end if
        
        if (rank .ne. tasks-1) then
            call MPI_Recv(ni_pl(:,nz_loc), nr, MPI_REAL8, &
                          rank-1, 16, comm, recv_status, ierr)
            call MPI_Recv(ne_pl(:,nz_loc), nr, MPI_REAL8, &
                          rank-1, 16, comm, recv_status, ierr)
            call MPI_Recv(nt_pl(:,nz_loc), nr, MPI_REAL8, &
                          rank-1, 16, comm, recv_status, ierr)
        end if
        if (rank .ne. 0) then
            call MPI_Send(ni_pl(:,1), nr, MPI_REAL8, &
                          rank+1, 16, comm, recv_status, ierr)
            call MPI_Send(ne_pl(:,1), nr, MPI_REAL8, &
                          rank+1, 16, comm, recv_status, ierr)
            call MPI_Send(nt_pl(:,1), nr, MPI_REAL8, &
                          rank+1, 16, comm, recv_status, ierr)
        end if
        
        time = time + deltime
        
        ni_mi  = ni_pl
        ni_org = ni_pl
        ne_mi  = ne_pl
        ne_org = ne_pl
        nt_mi  = nt_pl
        nt_org = nt_pl
        
        call update_coef
        
        if ( (timestep - int(timestep/1000)*1000 == 0) .and. (rank == 0) ) &
            write(6,90) timestep/1000, deltime, time, fintime
90 format('Step # ', i5, 'e3 dT: ', e9.2, ' Time: ', e9.2, ' Final Time: ', f7.2)
        
        if (time >= fintime) exit
        
        ! update boundary nodes
        
        ! assemble b for next timestep
        do node = Istart, Iend-1
            i_loc  = loc_idx(node+1, 1)
            j_loc  = loc_idx(node+1, 2)
            i_glob = glob_idx(node+1,1)
            j_glob = glob_idx(node+1,2)
        
            b_temp = 0.0_dp
        
            call laplace(i_loc, j_loc, type_z(i_glob,j_glob), &
                         type_r(i_glob,j_glob), b_temp)
        
            call VecSetValues(b, 1, node, -b_temp, insert_values, ierr)
        end do

        call vecassemblybegin(b, ierr)
        call vecassemblyend(b, ierr)
    end do
    call petsc_finalize
    end subroutine
    
!-----------------------------------------------------------------------
!**************************** Check Error *****************************
!-----------------------------------------------------------------------
    
    subroutine checkerr

    call MatMult(A,x,u,ierr)
    call VecAYPX(u,none,b,ierr)
    call VecNorm(u,NORM_2,norm,ierr)
    call KSPGetIterationNumber(ksp,its,ierr)
    if (rank .eq. 0) then
        if (norm .gt. 1.e-12) then
            write(6,100) norm,its
        else
            write(6,110) its
        endif
    endif
100 format('Norm of error ',e11.4,' iterations ',i5)
110 format('Norm of error < 1.e-12 iterations ',i5)
    end subroutine   

    end module
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
