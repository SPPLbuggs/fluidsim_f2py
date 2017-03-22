    module mod
    use properties
    use laplace_lib
    use petsc_lib
    use explicit_lib
    implicit none
    
    !f2py integer:: neqn, nr, nz, nr_loc
    !f2py real(8):: phiL, phiR
    !f2py integer, allocatable :: type_z(:,:), type_r(:,:), glob_node(:,:)
    !f2py real(8), allocatable :: z(:), r(:), phi(:,:)
    !f2py real(8), allocatable :: ni_pl(:,:), ne_pl(:,:), nt_pl(:,:)
    !f2py real(8), allocatable :: phi_save(:,:,:), ne_save(:,:,:)
    
    
    contains
!-----------------------------------------------------------------------
!***************************** Initialize *****************************
!-----------------------------------------------------------------------

    subroutine initialize(rel_tol, abs_tol)
    !f2py real(8), intent(in):: rel_tol, abs_tol
    real(dp), intent(in) :: rel_tol, abs_tol
    integer:: i, j, rows, cols(5)
    real(dp):: A_temp(1,5), b_temp
    
    call petsc_initialize
    
    allocate( ki(5, nz, nr_loc), ke(5, nz, nr_loc), kt(5, nz, nr_loc), &
              nu(nz, nr_loc), mue(nz, nr_loc), mut(nz, nr_loc), &
              De(nz, nr_loc), Dt(nz, nr_loc), k_ir(nz, nr_loc), &
              k_ex(nz, nr_loc), k_sc(nz, nr_loc), k_si(nz, nr_loc) )
    
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
    do j = Jstart, Jend
        do i = 1, nz
            
            b_temp =  0.0_dp
            A_temp =  0.0_dp
            cols   = -1
            rows   = glob_node(i,j) - 1
            
            if (rows < 0) cycle
               
            call laplace(i, j, type_r(i,j), type_z(i,j), b_temp)
            call jacobian(i, j, type_r(i,j), type_z(i,j), b_temp, cols, A_temp)
            
            call MatSetValues(A, 1, rows, 5, cols, A_temp, insert_values, ierr)
            call VecSetValues(b, 1, rows, -b_temp, insert_values, ierr)
        end do
    end do

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    
    call vecassemblybegin(b, ierr)
    call vecassemblyend(b, ierr)
    
    ! create linear solver
    call KSPCreate(comm,ksp,ierr)
    call KSPSetOperators(ksp,A,A,ierr)
    call KSPSetFromOptions(ksp,ierr)
    call KSPSetTolerances(ksp,rel_tol,abs_tol, & 
                          PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
    
    end subroutine
!-----------------------------------------------------------------------
!******************************** View ********************************
!-----------------------------------------------------------------------
    subroutine view
    integer :: wait
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)
    read(*,*) wait
    end subroutine

!-----------------------------------------------------------------------
!******************************** Run *********************************
!-----------------------------------------------------------------------
    subroutine run(fintime,savetime,nsave)
    !f2py integer, intent(in):: nsave
    !f2py real(8), intent(in):: fintime, savetime(nsave)
    integer, intent(in) :: nsave
    real(dp), intent(in):: fintime, savetime(nsave)
    integer  :: i, j, node, stage, recv_status, Nprint = 1000, ksave = 1
    real(dp) :: time = 0, stime, itime, etime, ftime, b_temp, &
                err_prev = 1, soln = 0,  &
                avg_itime = 0, avg_etime = 0, avg_stime = 0
    
    do
        call cpu_time(stime)
        timestep = timestep + 1
        
        ! solve implicit system
        !call view
        call KSPSolve(ksp,b,x,ierr)
        call KSPGetIterationNumber(ksp,its,ierr)
        
        call cpu_time(itime)
        ! update loc implicit variables
        do j = Jstart, Jend
            do i = 1, nz
                node = glob_node(i,j) - 1

                if (node .ge. 0) then
                    call VecGetValues(x, 1, node, soln, ierr)
                    phi(i, j) = phi(i, j) + soln
                end if
            end do
        end do
        
        ! communicate loc boundary nodes               
        if (rank .ne. 0) &
            call MPI_Send(phi(:,2), nz, MPI_REAL8, &
                          rank-1, 15, comm, ierr)
        
        if (rank .ne. tasks -1) &
            call MPI_Send(phi(:,nr_loc-1), nz, MPI_REAL8, &
                          rank+1, 16, comm, ierr)
        
        if (rank .ne. tasks-1) &
            call MPI_Recv(phi(:,nr_loc), nz, MPI_REAL8, &
                          rank+1, 15, comm, recv_status, ierr)
        
        if (rank .ne. 0) &
            call MPI_Recv(phi(:,1), nz, MPI_REAL8, &
                          rank-1, 16, comm, recv_status, ierr)
        
        call MPI_Barrier( comm, ierr )

        ! solve explicit system
        stage = 0
        ki = 0
        kt = 0
        ke = 0
        do
            stage = stage + 1
            if (stage > 5) exit
            
            call update_coef
            
            do j = Jstart, Jend
                do i = 1, nz
                    call continuity(i, j, type_r(i,j), type_z(i,j), stage)
                end do
            end do
            
            call rk_step(stage, time, fintime, err_prev)
        end do
        
        ! communicate loc boundary nodes
        if (rank .ne. 0) then
            call MPI_Send(ni_pl(:,2), nz, MPI_REAL8, &
                          rank-1, 15, comm, ierr)
            call MPI_Send(ne_pl(:,2), nz, MPI_REAL8, &
                          rank-1, 15, comm, ierr)
            call MPI_Send(nt_pl(:,2), nz, MPI_REAL8, &
                          rank-1, 15, comm, ierr)
        end if
        if (rank .ne. tasks -1) then
            call MPI_Send(ni_pl(:,nr_loc-1), nz, MPI_REAL8, &
                          rank+1, 16, comm, ierr)
            call MPI_Send(ne_pl(:,nr_loc-1), nz, MPI_REAL8, &
                          rank+1, 16, comm, ierr)
            call MPI_Send(nt_pl(:,nr_loc-1), nz, MPI_REAL8, &
                          rank+1, 16, comm, ierr)
        end if
        if (rank .ne. tasks-1) then
            call MPI_Recv(ni_pl(:,nr_loc), nz, MPI_REAL8, &
                          rank+1, 15, comm, recv_status, ierr)
            call MPI_Recv(ne_pl(:,nr_loc), nz, MPI_REAL8, &
                          rank+1, 15, comm, recv_status, ierr)
            call MPI_Recv(nt_pl(:,nr_loc), nz, MPI_REAL8, &
                          rank+1, 15, comm, recv_status, ierr)
        end if
        if (rank .ne. 0) then
            call MPI_Recv(ni_pl(:,1), nz, MPI_REAL8, &
                          rank-1, 16, comm, recv_status, ierr)
            call MPI_Recv(ne_pl(:,1), nz, MPI_REAL8, &
                          rank-1, 16, comm, recv_status, ierr)
            call MPI_Recv(nt_pl(:,1), nz, MPI_REAL8, &
                          rank-1, 16, comm, recv_status, ierr)
        end if
        
        call MPI_Barrier( comm, ierr )
        
        call cpu_time(etime)
        
        time = time + deltime
        
        ni_mi  = ni_pl
        ni_org = ni_pl
        ne_mi  = ne_pl
        ne_org = ne_pl
        nt_mi  = nt_pl
        nt_org = nt_pl
        
        if (time >= savetime(ksave)) then
            phi_save(:,:,ksave) = phi
            ne_save(:,:,ksave) = ne_pl
            ksave = ksave+1
        end if
        
        call cpu_time(ftime)
        if (timestep == 1) then
            avg_stime = ftime - stime
            avg_itime = itime - stime
            avg_etime = etime - itime
        else
            avg_stime = avg_stime - 0.01_dp * (avg_stime - ftime + stime)
            avg_itime = avg_itime - 0.01_dp * (avg_itime - itime + stime)
            avg_etime = avg_etime - 0.01_dp * (avg_etime - etime + itime)
        end if
        
        if ( (timestep - int(timestep/Nprint)*Nprint == 0) .and. (rank == 0) ) then
            write(6,90) timestep/Nprint, int(log10(real(Nprint))), deltime, time, fintime
            call KSPGetIterationNumber(ksp,its,ierr)
            write(6,91) avg_stime, avg_itime, avg_etime, its
            write(*,*)
        end if
90 format('Step # ', i4,'e',i1 '  dT:', es9.2, '  Time:', es9.2, '  Final Time:', f4.1)
91 format('  time/s: ', f5.3, '  iTime: ', f5.3, ' eTime: ', f5.3, '  ksp:', i4)

        if (time >= fintime) exit
        
        ! update boundary nodes
        
        ! assemble b for next timestep
        call VecSet(b,0,ierr)
        do j = Jstart, Jend
            do i = 1, nz
                b_temp =  0.0_dp
                node   = glob_node(i,j) - 1
                
                call laplace(i, j, type_r(i,j), type_z(i,j), b_temp)
                call VecSetValues(b, 1, node, -b_temp, insert_values, ierr)
            end do
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
