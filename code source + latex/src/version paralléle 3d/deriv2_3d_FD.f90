program deriv2_3d
    use mpi
    use decomp_2d
    implicit none

    integer :: Nx, Ny, Nz
    integer :: i, j, k, ierr
    integer :: p_row, p_col, rank, size
    integer :: local_xstart(3), local_xend(3)
    real(kind=8), allocatable :: fi(:, :, :), lap_fi(:, :, :)
    real(kind=8) :: dx, dy, dz
    
    ! Taille globale du domaine
    Nx = 128
    Ny = 128
    Nz = 128
    p_row = 4
    p_col = 4

    dx = 1.0 / Nx  ! Supposons un domaine normalisé [0,1]
    dy = 1.0 / Ny 
    dz = 1.0 / Nz 

    ! Initialisation MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! Initialisation du domaine decomp_2d
    call decomp_2d_init(Nx, Ny, Nz, p_row, p_col)

    ! Obtenir les indices de début et de fin du sous-domaine local
    call decomp_2d_get(local_xstart, local_xend)

    ! Allocation des tableaux locaux
    allocate(fi(local_xstart(1):local_xend(1), local_xstart(2):local_xend(2), local_xstart(3):local_xend(3)))
    allocate(lap_fi(local_xstart(1):local_xend(1), local_xstart(2):local_xend(2), local_xstart(3):local_xend(3)))

    ! Initialisation du champ fi
    do k = local_xstart(3), local_xend(3)
        do j = local_xstart(2), local_xend(2)
            do i = local_xstart(1), local_xend(1)
                fi(i, j, k) = sin(real(i) * dx) * sin(real(j) * dy) * sin(real(k) * dz)
            end do
        end do
    end do

    ! Calcul de la dérivée seconde en x avec différences finies centrées (ordre 2)
    ! Utilisation de la méthode de différences finies centrées en 3D
    do k = local_xstart(3), local_xend(3)
        do j = local_xstart(2), local_xend(2)
            do i = local_xstart(1)+1, local_xend(1)-1
                lap_fi(i, j, k) = (fi(i+1, j, k) - 2.0 * fi(i, j, k) + fi(i-1, j, k)) / (dx**2)
            end do
        end do
    end do

    ! Finalisation de la décomposition de domaine et MPI
    call decomp_2d_finalize()
    call MPI_Finalize(ierr)

end program deriv2_3d

