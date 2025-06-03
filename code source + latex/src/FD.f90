program poisson_2D
    implicit none
    integer, parameter :: N = 100      ! Taille de la grille
    integer :: i, j, iter
    integer, parameter :: max_iter = 10000
    real(8), parameter :: L = 1.0      ! Domaine [0,L]x[0,L]
    real(8), parameter :: tol = 1.0E-6 ! Critère d'arrêt
    real(8) :: dx, dy, h2, error
    real(8), dimension(N, N) :: phi, phi_new, S

    ! Définition du pas spatial
    dx = L / (N-1)
    dy = dx  ! On suppose un maillage carré
    h2 = dx**2

    ! Initialisation de phi et de S
    phi = 0.0  ! Conditions initiales : phi nul partout
    S = 1.0    ! Terme source constant pour ce test

    ! Boucle de Gauss-Seidel avec conditions périodiques en x uniquement
    do iter = 1, max_iter
        error = 0.0

        ! Mise à jour des points intérieurs (i=2:N-1, j=2:N-1)
        do i = 2, N-1
            do j = 2, N-1
                phi_new(i, j) = (S(i, j) * h2 + phi(i+1, j) + phi(i-1, j) + phi(i, j+1) + phi(i, j-1)) / 4.0
                error = max(error, abs(phi_new(i, j) - phi(i, j)))
            end do
        end do

        ! Conditions périodiques en x (i=1 et i=N) pour j=2:N-1
        do j = 2, N-1
            ! Bord i=1 (utilise i=N pour la périodicité)
            phi_new(1, j) = (S(1, j) * h2 + phi(2, j) + phi(N, j) + phi(1, j+1) + phi(1, j-1)) / 4.0
            ! Bord i=N (utilise i=1 pour la périodicité)
            phi_new(N, j) = (S(N, j) * h2 + phi(N-1, j) + phi(1, j) + phi(N, j+1) + phi(N, j-1)) / 4.0
            error = max(error, abs(phi_new(1, j) - phi(1, j)))
            error = max(error, abs(phi_new(N, j) - phi(N, j)))
        end do

        ! Mettre à jour phi (les bords j=1 et j=N restent à 0)
        phi = phi_new

        ! Vérification de la convergence
        if (error < tol) then
            print *, "Convergence atteinte après", iter, "itérations"
            exit
        end if
    end do

    ! Sauvegarde des résultats
    open(10, file="solution.dat", status="unknown")
    do i = 1, N
        do j = 1, N
            write(10,*) i*dx, j*dy, phi(i,j)
        end do
        write(10,*)
    end do
    close(10)

    print *, "Résultat sauvegardé dans solution.dat"
end program poisson_2D