!--------------------
! - Solves the Schrodinger equation in one dimension with Numerov methods
! - You need "input.dat" and "pes.inp" file
! - Last modified on 2021/05/03
! - Written by Kazuaki Kuwahata
!---------------------


program main
implicit none
real(8), parameter :: eps = 1.0d-5
real(8), parameter :: amu2au = 1.6605655D+4/9.1093897D0 ! atomic unit mass to atomic unit
real(8), parameter :: au2cm = 2.19466D+5 ! atomic unit to wave number
real(8), parameter :: ang2au = 1.d0/0.52918d0
real(8), parameter :: au2ang = 0.52918d0
integer :: i, nx, itr, itr_max
real(8) :: xmin, xmax, mass,kine, h, h12
real(8) :: eI, de, bis1, bis2, kine1, kine2, Penergy, Zenergy
real(8), allocatable :: psi(:)
real(8), allocatable :: f(:)
integer :: Ulog, Unor

call ReadInp
call ReadPotential

h = (xmax-xmin)/dble(nx) * ang2au
! h = (xmax-xmin)/dble(nx)
h12 = h*h/12.0d0
kine = eI
allocate(psi(0:nx))
  psi(0) = 0.0d0
  psi(1) = 1.d-4
open(newunit=Unor, file='norm.out', status='replace');close(Unor)

call calc_numerov
bis1 = psi(nx); bis2 = bis1
open(NewUnit=Ulog, file='log.out', status='replace')
  write(Ulog,'("# Starting step search")')
  call printlog(1,bis1)

  do itr = 1, itr_max
    kine = kine + de
    bis1 = bis2
    call calc_numerov
    bis2 = psi(nx)
    call printlog(itr+1,bis2)
    if ( bis1*bis2 <= 0 ) exit
  end do
  if ( itr >= itr_max ) goto 922

  write(Ulog,'("# Starting bisection search")')
  kine2 = kine
  kine1 = kine - de
  do itr = 1, itr_max
    kine = 0.5 * (kine1 + kine2)
    call calc_numerov
    bis2 = psi(nx)
    call printlog(itr,bis2)
    if (abs(bis2) <= eps) exit
    if (bis2*bis1 <= 0.0d0) then
      kine2 = kine
    else
      kine1 = kine
      bis1 = bis2
    end if
  end do
  if ( itr >= itr_max ) goto 923

  write(Ulog,'("# Searching is finished",/)')
  call calc_potential
  write(Ulog,'(4A20)') "Kinetic", "Potential","Zero-point", "Zero+Potential"
  write(Ulog,'(4f20.7)') kine, Penergy, Zenergy, Zenergy+Penergy
close(Ulog)

call print_psi

stop 'Normal termination Data is saved in "wave.out"'
922 stop "Erro first serching"
923 stop "Erro bisection serching"
contains

subroutine printlog(step,bis)
  integer, intent(in) :: step
  real(8), intent(in) :: bis
  write(Ulog,'("No.", I3,"   Kinetic = ",f13.8, "   Bound = ",f13.8)') step, kine, bis
end subroutine printlog

subroutine print_psi
  integer :: Upsi
  real(8) :: x
  open(NewUnit=Upsi, file='wave.out', status='replace')
    do i = 0, nx
      x = xmin + h*dble(i)*au2ang
      write(Upsi,*) x, psi(i)*sqrt(ang2au), f(i)
    end do
  close(Upsi)
end subroutine print_psi

subroutine ReadInp
  character(len=20) keyword
  open(20, file='input.dat', status='old', err=920)
    do
      read(20,*,end=101) keyword
      if ( keyword == "mass"      ) then
        read(20,*) mass
      else if ( keyword == "eI"   ) then
        read(20,*) ei
      else if ( keyword == "de"   ) then
        read(20,*) de
      else if ( keyword == "itr_max" ) then
        read(20,*) itr_max
      else if ( keyword == "EOF" ) then
        exit
!      else if ( keyword == "xmin" ) then
!        read(20,*) xmin
!      else if ( keyword == "xmax" ) then
!        read(20,*) xmax
!      else if ( keyword == "nx"   ) then
!        read(20,*) nx
      else
        goto 901
      end if
    end do
101 continue
  close(20)

  return
  920 continue
    print *, "There is no input.dat"
    open(30, file='input.dat', status='new')
      write(30,*) ' mass    # atomic mass of particle'
      write(30,*) '1.007825 '
      write(30,*) ' eI      # Initial kinetic energy'
      write(30,*) '0.0      '
      write(30,*) ' de      # Energy step for search'
      write(30,*) '0.01     '
      write(30,*) ' itr_max # The number of trials for search'
      write(30,*) '200      '
      write(30,*) ' EOF     '
    close(30)
  stop "ERROR!!"
  901 stop "Check the imput parameter!!"
end subroutine ReadInp

subroutine ReadPotential
  real(8) :: dummy
  real(8), allocatable :: xcoor(:)
  integer :: Upes
  nx = 0
  open(NewUnit=Upes, file='pes.inp', status='old', err=931)
    do
      read(Upes,*,err=900,end=100) dummy
      nx = nx + 1
    end do
100 continue
    nx = nx -1
    rewind(Upes)

allocate(f(0:nx))
allocate(xcoor(0:nx))

    do i = 0, nx
      read(Upes,*) xcoor(i), f(i)
    end do
  close(Upes)

  if ( eI == 0.0d0 ) then
    eI = minval(f) + 0.5d0*de
!    eI = minval(f)
    print *, "eI is ", eI
  end if
!  print *, xcoor(0), xcoor(nx)
  xmin = xcoor(0)
  xmax = xcoor(nx)
  return
  931 stop "Thre is no pes.inp"
  900 stop 'Check the format of "pes.inp"'
end subroutine ReadPotential

subroutine calc_numerov
  real(8) :: norm
  do i = 0, nx-2, 1
    psi(i+2) = 2.0*(1.0-5.0*h12*k2(i+1))*psi(i+1) - (1.0+h12*k2(i))*psi(i)
    psi(i+2) = psi(i+2) / (1.0+h12*k2(i+2))
  end do

! Normalized with Simpson's rule
  norm = 0.0d0
  norm = psi(0)**2 + psi(nx)**2 + 4.0d0*psi(nx-1)**2
  do i = 1, nx-3, 2
    norm = norm + 4.0d0*psi(i)**2 + 2.0d0*psi(i+1)**2
  end do
!  norm = sqrt(norm*h/3.0d0)
  norm = norm*h/3.0d0
  open(NewUnit=Unor, file='norm.out', position='append')
    write(Unor,*) norm
  close(Unor)
  psi(:) = psi(:) / sqrt(norm)
end subroutine calc_numerov

real(8) function k2(i)
  integer, intent(in) :: i
  k2 = 2.0 * amu2au * mass * ( kine - f(i) )
end function k2

subroutine calc_potential
  Penergy = 0.0d0
  do i = 0, nx, 1
    Penergy = Penergy + f(i)*psi(i)*psi(i)*h
  end do
  Zenergy = kine - minval(f)
end subroutine

end program main


!      if ( keyword == "xmin" ) then; read(20,*) xmin; else; goto 901; end if
!    read(20,*) keyword
!      if ( keyword == "xmax" ) then; read(20,*) xmax; else; goto 901; end if
!    read(20,*) keyword
!      if ( keyword == "mass" ) then; read(20,*) mass; else; goto 901; end if
!    read(20,*) keyword
!      if ( keyword == "nx" ) then; read(20,*) nx; else; goto 901; end if
!    read(20,*) keyword
!      if ( keyword == "ei" ) then; read(20,*) ei; else; goto 901; end if
!    read(20,*) keyword
!      if ( keyword == "de" ) then; read(20,*) de; else; goto 901; end if
!    read(20,*) keyword
!      if ( keyword == "itr_max" ) then; read(20,*) itr_max; else; goto 901; end if

