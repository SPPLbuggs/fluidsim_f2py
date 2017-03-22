    module petsc_lib
    use properties
    implicit none

    contains 
    subroutine petsc_initialize
    
    integer:: i, j, loc_neqn

    call PetscInitialize(petsc_null_character, ierr)
    comm = PETSC_COMM_WORLD
        
    none = -1.0
    
    call MPI_Comm_rank(comm,rank,ierr)
    call MPI_Comm_size(comm,tasks,ierr)
        
    ! create parallel matrix
    call MatCreate(comm,A,ierr)
    
    Jstart = 2
    Jend   = nr_loc-1
    if (rank == 0) Jstart = 1
    if (rank == tasks-1) Jend = nr_loc

    loc_neqn = 0
    do j = Jstart, Jend
        do i = 1, nz
            if (glob_node(i, j) > 0) loc_neqn = loc_neqn+1
        end do
    end do
                
    call MatSetSizes(A,loc_neqn,loc_neqn,neqn,neqn,ierr)
    call MatSetUp(A,ierr)
        
    call MatSetFromOptions(A,ierr)
    call MatSeqAIJSetPreallocation(A, 5, petsc_null_integer, ierr)
    call MatSetOption(A, mat_ignore_zero_entries, petsc_true, ierr)

    ! find parallel partitioning range
    call MatGetOwnershipRange(A,Istart,Iend,ierr)

    ! create parallel vectors
    call VecCreateMPI(comm,loc_neqn,neqn,u,ierr)
    call VecCreateMPI(comm,loc_neqn,neqn,b,ierr)
    call VecCreateMPI(comm,loc_neqn,neqn,x,ierr)
    call VecSetFromOptions(u,ierr)
    call VecSetFromOptions(b,ierr)
    call VecSetFromOptions(x,ierr)
    call VecSetOption(b, vec_ignore_negative_indices, petsc_true, ierr)
    
    end subroutine

    subroutine petsc_finalize
    
    call KSPDestroy(ksp,ierr)
    call VecDestroy(u,ierr)
    call VecDestroy(x,ierr)
    call VecDestroy(b,ierr)
    call MatDestroy(A,ierr)
    call PetscFinalize(ierr)

    end subroutine
    end module
