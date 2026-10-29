!---------------------------------------------------------------------
! Program: klimontovich
! Equation: df/dt + v*df/dx + (qE/m)*df/dv = 0
! Domain:   1D space (periodic), 1D velocity
! Scheme:   Upwind (x) + Central (v), explicit Euler in time
!---------------------------------------------------------------------
program klimontovich
  implicit none

  !--- Physical constants ---
  real(8), parameter :: q  = 1.602176634d-19   ! charge (C)
  real(8), parameter :: m  = 1.676221923d-27   ! mass (kg)
  real(8), parameter :: pi = 3.14159265358979323846d0

  !--- Grid parameters (input) ---
  integer :: nx, nv, nsteps
  real(8) :: L, v_max, dx, dv, dt, E_ext
  real(8) :: CFL

  !--- IC parameters ---
  real(8) :: sigma_x, v_th

  !--- Arrays (with ghost cells) ---
  real(8), allocatable :: f(:,:), f_new(:,:)
  real(8), allocatable :: x_arr(:), v_arr(:)

  !--- Diagnostics ---
  real(8) :: N_particles, P_momentum, U_energy
  integer :: i, j, istep

  print *, "------------------------------------------------------------------"
  print *, "  KLIMONTOVICH EQUATION  (1D-1V, prescribed E-field)"
  print *, "  df/dt + v*df/dx + (qE/m)*df/dv = 0"
  print *, "------------------------------------------------------------------"
  print *, ""

  !---input from user ---
  print *, "Enter grid size in x (Nx):"
  read(*,*) nx
  print *, "Enter grid size in v (Nv):"
  read(*,*) nv
  print *, "Enter domain length L (m):"
  read(*,*) L
  print *, "Enter max velocity v_max (m/s):"
  read(*,*) v_max
  print *, "Enter external E-field (V/m):"
  read(*,*) E_ext
  print *, "Enter number of time steps:"
  read(*,*) nsteps
  print *, "Enter bump width sigma_x (m):"
  read(*,*) sigma_x
  print *, "Enter thermal velocity v_th (m/s):"
  read(*,*) v_th

  dx = L / nx
  dv = 2.0d0 * v_max / nv

  ! CFL condition
  CFL = dx / (2.0d0 * v_max)
  dt  = 0.5d0 * CFL
  print *, ""
  print *, "dx = ", dx, " m"
  print *, "dv = ", dv, " m/s"
  print *, "dt = ", dt, " s  (CFL-limited)"
  print *, ""

  !--- Allocate nx for ghost cells ---
    allocate(f(0:nx+1, 0:nv+1), f_new(0:nx+1, 0:nv+1))
  allocate(x_arr(0:nx+1), v_arr(0:nv+1))   

  !--- Build grids (interior: 1..nx, 1..nv) ---
  do i = 1, nx
    x_arr(i) = (i - 1) * dx
  end do
  do j = 1, nv
    v_arr(j) = -v_max + (j - 1) * dv
  end do

  !--- Initial condition: Gaussian bump ---
  f = 0.0d0
  do j = 1, nv
    do i = 1, nx
      f(i,j) = exp(-((x_arr(i) - L/2.0d0)**2) / (sigma_x**2)) &
             * exp(-(v_arr(j)**2) / (v_th**2))
    end do
  end do

  !--- Time loop ---
  do istep = 1, nsteps

    f_new = 0.0d0

    do j = 1, nv
      do i = 1, nx
        if (v_arr(j) >= 0.0d0) then
          f_new(i,j) = f(i,j) &
            - v_arr(j) * dt * (f(i+1,j) - f(i,j)) / dx &
            + (q * E_ext / m) * dt * (f(i,j+1) - f(i,j-1)) / (2.0d0 * dv)
        else
          f_new(i,j) = f(i,j) &
            - v_arr(j) * dt * (f(i,j) - f(i-1,j)) / dx &
            + (q * E_ext / m) * dt * (f(i,j+1) - f(i,j-1)) / (2.0d0 * dv)
        end if
      end do
    end do

    ! Periodic boundary in x (ghost cells)
    f_new(0,1:nv)    = f_new(nx,1:nv)
    f_new(nx+1,1:nv) = f_new(1,1:nv)
    ! Zero-gradient boundary in v (ghost cells)
    f_new(1:nx,0)    = f_new(1:nx,1)
    f_new(1:nx,nv+1) = f_new(1:nx,nv)

    f = f_new

    !--- Diagnostics every 100 steps ---
    if (mod(istep, 100) == 0 .or. istep == 1) then
      N_particles = 0.0d0
      P_momentum  = 0.0d0
      U_energy    = 0.0d0
      do j = 1, nv
        do i = 1, nx
          N_particles = N_particles + f(i,j) * dx * dv
          P_momentum  = P_momentum  + f(i,j) * v_arr(j) * dx * dv
          U_energy    = U_energy    + f(i,j) * 0.5d0 * v_arr(j)**2 * dx * dv
        end do
      end do
      print *, "Step ", istep, "  N = ", N_particles, &
               "  <v> = ", P_momentum/N_particles, &
               "  U = ", U_energy
    end if

  end do

  print *, ""
  print *, "---------------------------------------------------------------"
  print *, "  FINAL DIAGNOSTICS"
  print *, "---------------------------------------------------------------"
  print *, "Total particles (integral of f): ", N_particles
  print *, "Mean velocity <v> (m/s):          ", P_momentum / N_particles
  print *, "Kinetic energy U (arb. units):    ", U_energy
  print *, ""
  print *, "Note: <v> should shift by (qE/m)*t over the simulation."
  print *, "Expected shift = (q*E/m)*t_total = ", (q*E_ext/m) * nsteps * dt, " m/s"

  deallocate(f, f_new, x_arr, v_arr)
end program klimontovich   
