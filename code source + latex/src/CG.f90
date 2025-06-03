program solve_poisson
  implicit none
  integer, parameter :: nx = 1000, ny = 1000
  real(8) :: dx, dy
  real(8), dimension(nx, ny) :: phi, S
  real(8) :: tol
  integer :: max_iter, i, j

  ! Initialiser dx et dy
  dx = 1.0d0 / (nx - 1)
  dy = 1.0d0 / (ny - 1)

  ! Initialiser le problème
  max_iter = 1000
  tol = 1.0E-6

  ! Initialiser la source S et la solution phi (exemple avec une source aléatoire)
  S = 1.0
  phi = 0.0

  ! Résoudre l'équation de Poisson
  call gradient_conjugue(phi, S, nx, ny, dx, dy, max_iter, tol)

  ! Afficher une partie du résultat
  print *, 'Solution de phi à la position (1, 1) : ', phi(1, 1)

contains

  ! Sous-routine pour calculer le Laplacien
  subroutine laplacien(phi, lap, nx, ny, dx, dy)
    implicit none
    real(8), dimension(:,:), intent(in) :: phi
    real(8), dimension(:,:), intent(out) :: lap
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: dx, dy
    integer :: i, j

    ! Initialiser laplacien à zéro
    lap = 0.0d0

    ! Calculer la matrice du laplacien avec conditions périodiques
    do j = 1, ny
      do i = 1, nx
        lap(i,j) = (phi(mod(i, nx)+1, j) - 2.0d0*phi(i,j) + phi(mod(i-2, nx)+1, j)) / (dx**2) + &
                   (phi(i, mod(j, ny)+1) - 2.0d0*phi(i,j) + phi(i, mod(j-2, ny)+1)) / (dy**2)
      end do
    end do
  end subroutine laplacien

  ! Sous-routine pour enregistrer la solution phi dans un fichier binaire (.dat)
  subroutine save_phi(phi, nx, ny, filename)
    implicit none
    real(8), dimension(:,:), intent(in) :: phi
    integer, intent(in) :: nx, ny
    character(len=*), intent(in) :: filename
    integer :: iunit

    ! Ouvrir le fichier en mode écriture binaire (unformatted)
    open(unit=iunit, file=filename, form='unformatted', access='stream')

    ! Écrire les dimensions dans le fichier (nx et ny)
    write(iunit) nx, ny

    ! Écrire les données de phi dans le fichier
    write(iunit) phi

    ! Fermer le fichier
    close(iunit)
  end subroutine save_phi

  ! Sous-routine pour calculer le produit scalaire entre deux matrices
  function produit_scalaire(A, B, nx, ny) result(res)
    implicit none
    real(8), dimension(:,:), intent(in) :: A, B
    integer, intent(in) :: nx, ny
    real(8) :: res
    integer :: i, j

    res = 0.0d0
    do j = 1, ny
      do i = 1, nx
        res = res + A(i,j) * B(i,j)
      end do
    end do
  end function produit_scalaire

  ! Sous-routine pour la méthode du gradient conjugué
  subroutine gradient_conjugue(phi, S, nx, ny, dx, dy, max_iter, tol)
    implicit none
    real(8), dimension(:,:), intent(inout) :: phi
    real(8), dimension(:,:), intent(in) :: S
    integer, intent(in) :: nx, ny, max_iter
    real(8), intent(in) :: dx, dy, tol
    real(8), dimension(:,:), allocatable :: r, p, Ap
    real(8) :: alpha, beta, r_old, r_new
    integer :: iter, i, j

    ! Initialisation des variables
    allocate(r(nx, ny), p(nx, ny), Ap(nx, ny))
    ! Appel à laplacien et assignation à r
    call laplacien(phi, r, nx, ny, dx, dy)  ! Remplacer r = S - laplacien(phi) par un appel à laplacien
    r = S - r  ! Calculer le résidu

    p = r
    r_old = produit_scalaire(r, r, nx, ny)

    ! Boucle principale du gradient conjugué
    do iter = 1, max_iter
      ! Calcul de Ap = A * p (appliqué au Laplacien)
      call laplacien(p, Ap, nx, ny, dx, dy)
      
      ! Calcul de alpha
      alpha = r_old / produit_scalaire(p, Ap, nx, ny)
      
      ! Mise à jour de phi
      phi = phi + alpha * p
      
      ! Mise à jour du résidu r
      r = r - alpha * Ap
      r_new = produit_scalaire(r, r, nx, ny)
      
      ! Test de convergence
      if (sqrt(r_new) < tol) exit
      
      ! Calcul de beta
      beta = r_new / r_old
      p = r + beta * p
      r_old = r_new
    end do
  ! Sauvegarder les résultats dans un fichier binaire
  call save_phi(phi, nx, ny, 'solution_cg.dat')    
    ! Libération de la mémoire allouée
    deallocate(r, p, Ap)
  end subroutine gradient_conjugue

end program solve_poisson
