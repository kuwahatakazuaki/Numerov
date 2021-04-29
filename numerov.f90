!--------------------
! - Solves the Schrodinger equation in one dimension with Numerov methods
! - You need "input.dat" and "pes.inp" file
! - Last modified on 2018/10/24
! - Written by Kazuaki Kuwahata
!---------------------
module globals
  implicit none
  integer, save :: nx
  real(8), save :: xmin, xmax, mass,kine, h, h12
  real(8), save, allocatable :: f(:)
  real(8), parameter :: amu2au = 1.6605655D+4/9.1093897D0 ! atomic unit mass to atomic unit
  real(8), parameter :: au2cm = 2.19466D+5 ! atomic unit to wave number
contains

end module globals

program main
use globals
implicit none
integer :: i, ni, itr, itr_max
real(8), parameter :: eps = 1.0d-5
real(8) :: eI, de, bis1, bis2, kine1, kine2, Penergy, Zenergy
real(8), allocatable :: psi(:)
integer :: Ulog
! Penergy: Potential energy, Zenergy: Zero-point energy

! xmin, xmax is the range of x, mass is reduced mass
! nx is the number of steps. nx must be 'EVEN'
! eI is the initial kinetic energy, de is the delta of kinetic energy
! itr_max is the maximum of the iteration
call ReadInp(xmin,xmax,mass,nx,eI,de,itr_max)
! Reading potential energy from "pes.inp"
allocate(f(0:nx))
call ReadPotential
if ( eI == 0.0d0 ) then
  eI = minval(f) + de
  print *, "eI is ", eI
end if
h = (xmax-xmin)/dble(nx)
h12 = h*h/12.0d0
kine = eI
allocate(psi(0:nx))
  psi(0) = 0.0d0
  psi(1) = 1.d-4
open(23, file='norm.out', status='replace');close(23)


call numerov(psi)
bis1 = psi(nx); bis2 = bis1
! open(22, file='log.out', status='replace')
open(NewUnit=Ulog, file='log.out', status='replace')
  write(Ulog,'("# Starting step search")')
  call printlog(1,kine,bis1)

  do itr = 1, itr_max
    kine = kine + de
    bis1 = bis2
    call numerov(psi)
    bis2 = psi(nx)
    call printlog(itr+1,kine,bis2)
    if ( bis1 * bis2 <= 0) exit
  end do
  if ( itr >= itr_max ) then; goto 922; end if

  write(Ulog,'("# Starting bisection search")')
  kine2 = kine
  kine1 = kine - de
  do itr = 1, itr_max
    kine = (kine1+kine2)/2.0d0
    call numerov(psi)
    bis2 = psi(nx)
    call printlog(itr,kine,bis2)
    if (abs(bis2) <= eps) exit
    if (bis2*bis1 <= 0.0d0) then
      kine2 = kine
    else
      kine1 = kine
      bis1 = bis2
    end if
  end do
  if ( itr >= itr_max ) then; goto 923; end if
  write(Ulog,'("# Searching is finished",/)')
  call potential_energy(psi,Penergy,Zenergy)
  write(Ulog,'(4A20)') "Kinetic", "Potential","Zero-point", "Zero+Potential"
  write(Ulog,'(4f20.7)') kine, Penergy, Zenergy, Zenergy+Penergy
close(Ulog)

call print_psi(psi)

print *, "Normal termination Data is saved in wave.out"
stop
922 stop "Erro first serching"
923 stop "Erro bisection serching"
contains

  subroutine printlog(itr,kine,bis)
    integer, intent(in) :: itr
    real(8), intent(in) :: kine, bis
    write(Ulog,'("No.", I3,"   Kinetic = ",f13.8, "   Bound = ",f13.8)') itr, kine,bis
  end subroutine printlog

end program main

subroutine ReadInp(xmin,xmax,mass,nx,eI,de,itr_max)
  implicit none
  real(8) :: xmax,xmin,mass,eI,de
  integer :: nx, itr_max
  character(len=20) keyword
  open(20, file='input.dat', status='old', err=920)
    read(20,*) keyword
      if ( keyword == "xmin" ) then; read(20,*) xmin; else; goto 901; end if
    read(20,*) keyword
      if ( keyword == "xmax" ) then; read(20,*) xmax; else; goto 901; end if
    read(20,*) keyword
      if ( keyword == "mass" ) then; read(20,*) mass; else; goto 901; end if
    read(20,*) keyword
      if ( keyword == "nx" ) then; read(20,*) nx; else; goto 901; end if
    read(20,*) keyword
      if ( keyword == "eI" ) then; read(20,*) eI; else; goto 901; end if
    read(20,*) keyword
      if ( keyword == "de" ) then; read(20,*) de; else; goto 901; end if
    read(20,*) keyword
      if ( keyword == "itr_max" ) then; read(20,*) itr_max; else; goto 901; end if
  close(20)
  return
  901 stop "Check the imput parameter!!"
  920 stop "There is no input.dat"
end subroutine ReadInp

subroutine numerov(psi)
  use globals
  implicit none
  integer :: i, n
  real(8), intent(inout) :: psi(0:nx)
  real(8), external :: k2

  do i = 0, nx-2, 1
    psi(i+2) = 2.0*(1.0-5.0*h12*k2(i+1))*psi(i+1) - (1.0+h12*k2(i))*psi(i)
    psi(i+2) = psi(i+2) / (1.0+h12*k2(i+2))
  end do
  call normalize(psi)
end subroutine numerov

real(8) function k2(i)
  use globals
  implicit none
  integer, intent(in) :: i
  k2 = 2.0*amu2au*mass*(kine - f(i))
end function k2

! --- Reading from external potential ---
subroutine ReadPotential
  use globals
  implicit none
  integer :: i
  real(8) :: dummy
  open(31, file='pes.inp', status='old', err=931)
    do i = 0, nx
      read(31,*) dummy, f(i)
    end do
  close(31)
  return
  931 stop "Thre is no pes.inp"
end subroutine ReadPotential

subroutine normalize(psi)
  use globals
  implicit none
  integer :: i
  real(8), intent(inout) :: psi(0:nx)
  real(8) :: norm
  norm = psi(0)**2 + psi(nx)**2

  open(23, file='norm.out', position='append')
    write(23,*) norm
  close(23)

  do i = 1, nx-3, 2
    norm = norm + 4.0d0*psi(i)**2 + 2.0d0*psi(i+1)**2
  end do
  norm = norm + 4.0d0*psi(nx-1)**2
  norm = 1.0d0/sqrt(norm*h/3.0d0)
  psi(:) = psi(:) * norm
end subroutine normalize

subroutine print_psi(psi)
  use globals
  implicit none
  real(8),intent(in) :: psi(0:nx)
  real(8) :: x
  integer :: i
  open(21,file='wave.out', status='replace')
  do i = 0, nx
    x = xmin + h*dble(i) ! xmin + h*dble(i)
    write(21,*) x, psi(i), f(i)
  end do
  close(21)
end subroutine print_psi

subroutine potential_energy(psi,Penergy,Zenergy)
  use globals
  implicit none
  integer :: i
  real(8), intent(in) :: psi(0:nx)
  real(8), intent(out) :: Penergy, Zenergy
  real(8) :: minE
  Penergy = 0.0d0
  minE = f(0)
  do i = 0, nx, 1
    Penergy = Penergy + f(i)*psi(i)*psi(i)*h
    if (minE > f(i)) minE = f(i)
  end do
  Zenergy = kine - minE
end subroutine potential_energy



!! Making a potential
!real(8) function f(x)
!  implicit none
!  real(8) :: x
!! ---- harmonic potential ----
!  f = 0.1*x*x
!
!!!----- box type potential ------
!!  if ( -0.5d0 < x .and. x < 0.5d0 ) then
!!    f = 0.0d0
!!  else
!!    f = 5.0d0
!!  end if
!end function f

