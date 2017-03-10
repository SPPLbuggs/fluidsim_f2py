    module properties
    implicit none
    integer, parameter:: dp=selected_real_kind(15)
    
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>

    MPI_Comm comm
    PetscMPIInt rank,tasks
    Mat A
    Vec b,u,x
    KSP ksp
    PC pc
    PetscInt Istart, Iend, its
    PetscReal  norm
    PetscErrorCode ierr
    PetscScalar none
    
    integer  :: neqn, nz, nr, nz_loc
    real(dp) :: phiL, phiR, deltime = 1e-5
    integer, allocatable :: type_z(:,:), type_r(:,:), glob_idx(:,:), &
                            loc_idx(:,:), glob_node(:,:)
    real(dp), allocatable :: z(:), r(:), phi(:,:), ki(:,:,:), ke(:,:,:), kt(:,:,:), &
                             ni_mi(:,:), ni_pl(:,:), ni_org(:,:), &
                             ne_mi(:,:), ne_pl(:,:), ne_org(:,:), &
                             nt_mi(:,:), nt_pl(:,:), nt_org(:,:)
    
    ! non-dimensional parameters
    real(dp), parameter:: L0   = 1e-3_dp
    real(dp), parameter:: Phi0 = 1e3_dp
    real(dp), parameter:: Tau0 = 1e-6_dp
    real(dp), parameter:: n0 = 1e16_dp
    real(dp), parameter:: T0 = 11604.52_dp

    ! fundamental constants
    real(dp), parameter:: eps = 8.85418782e-12_dp
    real(dp), parameter:: echarg = 1.60217646e-19_dp
    real(dp), parameter:: Pi = 4.0_dp * atan(1.0_dp)
    real(dp), parameter:: kb = 1.3806503e-23_dp

    ! case properties
    real(dp), parameter:: Tg = 300.0_dp
    real(dp), parameter:: p = 3.0_dp
    real(dp), parameter:: ninf = p * 101325.0_dp / 760.0_dp / kb / Tg

    ! electron properties
    real(dp), parameter:: me = 9.10938188e-31_dp

    ! positive ion properties
    real(dp), parameter:: mi  = 1.67262178e-27_dp * 39.948_dp
    real(dp), parameter:: mui = 1.45e3_dp / p * 1e-4_dp * phi0 * tau0 / l0**2
    real(dp), parameter:: Ti  = Tg
    real(dp), parameter:: vi  = Sqrt( (8.0_dp * kb*Ti) / (pi * mi) ) * tau0 / l0
    real(dp), parameter:: Di  = mui * (kb*Ti / echarg) / phi0

    ! metastable argon properties
    real(dp), parameter:: mm = mi
    real(dp), parameter:: Dm = 2.42e18_dp / ninf / 1e2_dp * tau0 / l0**2
    real(dp), parameter:: kr  = 2e-7_dp / 1.0e6_dp * n0 * tau0
    real(dp), parameter:: kmp = 6.2e-10_dp / 1.0e6_dp * n0 * tau0
    real(dp), parameter:: k2q = 3e-15_dp / 1.0e6_dp * n0 * tau0
    real(dp), parameter:: k3q = 1.1e-31_dp / 1.0e12_dp * n0**2 * tau0

    ! reactions
    real(dp), parameter:: H_ir  =  15.8_dp  / 1.5_dp
    real(dp), parameter:: H_ex  =  11.56_dp / 1.5_dp
    real(dp), parameter:: H_si  =   4.14_dp / 1.5_dp
    real(dp), parameter:: H_sc  = -11.56_dp / 1.5_dp
    real(dp), parameter:: beta  =   1e-13_dp * Tau0 * n0
    real(dp), parameter:: gamma =   0.075_dp
  
    ! Updatable variables
    real(dp), allocatable, dimension(:, :) :: mue, mut
    real(dp), allocatable, dimension(:, :) :: De, Dt
    real(dp), allocatable, dimension(:, :) :: nu
    real(dp), allocatable, dimension(:, :) :: k_ir, k_ex, k_sc, k_si

    contains
    
    subroutine update_coef
    integer:: node, i, j
    real(dp) :: Te
    do node = Istart, Iend-1
        i  = loc_idx(node+1, 1)
        j  = loc_idx(node+1, 2)
        Te = get_Te(nt_mi(i,j), ne_mi(i,j))
        
        mue(i,j) = get_mue(Te)
        mut(i,j) = get_mut(Te)
        De(i,j)  = get_De(Te)
        Dt(i,j)  = get_Dt(Te)
        nu(i,j)   = get_nu(Te)
        k_ir(i,j) = get_k_ir(Te)
        k_sc(i,j) = get_k_sc(Te)
        k_si(i,j) = get_k_si(Te)
        k_ex(i,j) = get_k_ex(Te)
    end do
    end subroutine
    
    subroutine update_bc
    integer :: i_loc, j_loc, i_glob, j_glob, node
    !assuming DC bc right now.
    ! phiR = 
    ! phiL = 
    do j_glob = 1, nz
        do i_glob = 1, nr
            node = glob_node(i_glob, j_glob)
            if ( (node >= Istart) .and. (node < Iend) ) then
                if (type_z(i_glob, j_glob) == -2) then
                    i_loc = loc_idx(node+1,1)
                    j_loc = loc_idx(node+1,2)
                    phi(i_loc, j_loc) = phiR
                end if
                if (type_z(i_glob, j_glob) == 2) then
                    i_loc = loc_idx(node+1,1)
                    j_loc = loc_idx(node+1,2)
                    phi(i_loc, j_loc) = phiR
                end if
             end if
         end do
     end do
     end subroutine
    
    function get_Te(nt,ne)
    real(dp):: get_Te
    real(dp), intent(in) :: nt, ne
    
    if ((nt >= 0.0_dp) .and. (ne >= 0.0_dp)) then
        get_Te = nt / ne
    else
        get_Te = 1e-8_dp
    end if
    
    return
    end function
    
    function get_mue(T)
    real(dp):: get_mue
    real(dp), intent(in):: T
    real(dp):: x,a,b,c,d,e,f,g,h,k
    a = 55.0126564455
    b =  0.685594575551
    c = -0.383637563328
    d = -0.0340209762821
    e =  0.0276748887423
    f =  0.00242420301108
    g = -0.0012183881312
    h = -5.6471677382e-05
    k =  2.11510269874e-05
    x = log(T*3./2.)
    if (x > log(200.)) x = log(200.)
    if (x < log(1.d-2)) x = log(1.d-2)
    get_mue = exp(a + b*x + c*x**2. + d*x**3. + e*x**4. &
        + f*x**5. + g*x**6. + h*x**7. + k*x**8.) / ninf * tau0 * phi0 / l0**2
    return
    end function get_mue
  
    function get_mut(T)
    real(dp):: get_mut
    real(dp), intent(in):: T
    real(dp):: x,a,b,c,d,e,f,g,h,k
    a = 55.7894575628
    b =  0.428129671316
    c = -0.429648313902
    d = -0.00193575524377
    e =  0.034430796495
    f =  0.000678700242152
    g = -0.00153955670552
    h = -2.38690441212d-05
    k =  2.58672391609d-05
    x = log(T*3./2.)
    if (x > log(200.)) x = log(200.)
    if (x < log(1.d-2)) x = log(1.d-2)
    get_mut = exp(a + b*x + c*x**2. + d*x**3. + e*x**4. &
        + f*x**5. + g*x**6. + h*x**7. + k*x**8.) / ninf * tau0 * phi0 / l0**2
    return
    end function get_mut
  
    function get_De(T)
    real(dp):: get_De
    real(dp), intent(in):: T
    real(dp):: x,a,b,c,d,e,f,g,h,k
    a = 54.6077441165
    b =  1.68552023011
    c = -0.383788369738
    d = -0.0340244943028
    e =  0.0276980470142
    f =  0.00242512124501
    g = -0.00121973056843
    h = -5.65037595561e-05
    k =  2.1177126055e-05
    x = log(T*3./2.)
    if (x > log(200.)) x = log(200.)
    if (x < log(1.d-2)) x = log(1.d-2)
    get_De = exp(a + b*x + c*x**2. + d*x**3. + e*x**4. &
        + f*x**5. + g*x**6. + h*x**7. + k*x**8.) / ninf * tau0 / l0**2
    return
    end function get_De

    function get_Dt(T)
    integer :: get_Dt
    real(dp), intent(in):: T
    real(dp):: x,a,b,c,d,e,f,g,h,k
    a = 55.3840545867
    b =  1.42801298837
    c = -0.429629341572
    d = -0.00191349710732
    e =  0.0344269595677
    f =  0.000677020293459
    g = -0.001539248252
    h = -2.3830286892e-05
    k =  2.58597001685e-05
    x = log(T*3./2.)
    if (x > log(200.)) x = log(200.)
    if (x < log(1.d-2)) x = log(1.d-2)
    get_Dt = exp(a + b*x + c*x**2. + d*x**3. + e*x**4. &
        + f*x**5. + g*x**6. + h*x**7. + k*x**8.) / ninf * tau0 / l0**2
    return
    end function get_Dt

    function get_k_ex(T)
    real(dp):: get_k_ex
    real(dp), intent(in):: T
    real(dp):: x,a,b,c,d,e,f,g,h,k
    a = -50.3785102239 
    b =  19.0129183764
    c = -11.7950315424
    d =   7.41674013553
    e =  -3.84148086698
    f =   1.2962229976
    g =  -0.259359346989
    h =   0.0279182131315
    k =  -0.00124438710099
    x = log(T*3./2.)
    if (x > log(200.)) x = log(200.)
    if (x < log(1.d-2)) x = log(1.d-2)
    get_k_ex = exp(a + b*x + c*x**2. + d*x**3. + e*x**4. &
        + f*x**5. + g*x**6. + h*x**7. + k*x**8.) * tau0 * n0
    return
    end function get_k_ex

    function get_k_ir(T)
    real(dp):: get_k_ir
    real(dp), intent(in):: T
    real(dp):: x,a,b,c,d,e,f,g,h,k
    a = -56.2478139553
    b =  26.6123052468
    c = -15.9868576469
    d =   7.97316507041
    e =  -3.18287109994
    f =   0.890666459851
    g =  -0.156771317788
    h =   0.0153819279555
    k =  -0.000638729430911
    x = log(T*3./2.)
    if (x > log(200.)) x = log(200.)
    if (x < log(1.d-2)) x = log(1.d-2)
    get_k_ir = exp(a + b*x + c*x**2. + d*x**3. + e*x**4. &
        + f*x**5. + g*x**6. + h*x**7. + k*x**8.) * tau0 * n0
    return
    end function get_k_ir

    function get_nu(T)
    real(dp):: get_nu
    real(dp), intent(in):: T
    real(dp):: x,a,b,c,d,e,f,g,h,k
    a = -32.3199262825
    b =   1.69388044254
    c =   0.0842378277404
    d =  -0.164556371047
    e =   0.00861307304011
    f =   0.00855257716132
    g =  -0.000983504760785
    h =  -0.000160952834008
    k =   2.37965210684e-05
    x = log(T*3./2.)
    if (x > log(200.)) x = log(200.)
    if (x < log(1.d-2)) x = log(1.d-2)
    get_nu = exp(a + b*x + c*x**2. + d*x**3. + e*x**4. &
        + f*x**5. + g*x**6. + h*x**7 + k*x**8.) * ninf * tau0
    return
    end function get_nu

    function get_k_sc(T)
    real(dp):: get_k_sc
    real(dp), intent(in):: T
    real(dp):: x,a,b,c,d,e,f
    a = -21.4827864151
    b =   0.457356923276
    c =  -0.555439231606
    d =   1.27257798891
    e =  -0.67840685073
    f =   0.10591014464
    x = log(T*3./2.)
    if (x > log(16.)) x = log(16.)
    if (x < log(1.d-2)) x = log(1.d-2)
    get_k_sc = exp(a + b*x + c*x**2. + d*x**3. + e*x**4. &
        + f*x**5. ) / 1.0d6 * n0 * tau0
    return
    end function get_k_sc

    function get_k_si(T)
    real(dp):: get_k_si
    real(dp), intent(in):: T
    real(dp):: x,a,b,c,d,e
    a = -43.1347385848
    b =  43.9905424566
    c = -28.1169537586
    d =   8.28853856817
    e =  -0.931626144207
    x = log(T*3./2.)
    if (x > log(16.)) x = log(16.)
    if (x < log(1.d-2)) x = log(1.d-2)
    get_k_si = exp(a + b*x + c*x**2. + d*x**3. + e*x**4.) &
        / 1.0d6 * n0 * tau0
    return
    end function get_k_si

    end module
