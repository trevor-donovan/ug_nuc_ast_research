subroutine read_data(a)                        !reads in the array
    implicit none
    real, dimension(1:3, 1:100), intent(out) :: a
    integer :: lunvar
    open(newunit=lunvar, file='num_densities_Np28_Nn28.dat', &
    & status='old', action='read')
    read(lunvar,*) a 
    close(unit=lunvar)
    return
end subroutine read_data 

subroutine integral_info(subint, num) !width and number of subintervals
    implicit none
    real, intent(out) :: subint
    integer, intent(out) :: num
    subint = .4 !the simpson subinterval is twice the step in the data
    num = 101
    return
end subroutine integral_info

function epp(n) result(val) !skyrme energy per particle at x=1/2
implicit none       !all quantities are in base units of MeV, fm, and 
    real, intent(in) :: n  !number density, not to be confused with program's n
    real :: val
    real :: a
    real :: b
    !model: SLY4; parameters are intended to be obvious 
    !from Boyang's papers' expressions
    real :: t0 = -2488.913
    real :: t1 = 486.818
    real :: t2 = -546.395
    real :: t31 = 13777.0
    real :: t32 = 0.0
    real :: t33 = 0.0
    real :: t4 = 0.0
    real :: t5 = 0.0
    real :: x1 = -.3438
    real :: x2 = -1.0
    real :: x5 = 0.0
    real :: gam = 0.0
    real :: del = 0.0
    real :: sig1 = .1666666667
    real :: sig2 = 0.0
    real :: sig3 = 0.0
    a = t1*(x1+2)+t2*(x2+2)
    b = t2*(x2+.5)-t1*(x1+.5)
    !ideally more efficient, but is now left like analytic form for clarity
    val = 74.959668903*n**(2.0/3.0) + .375*t0*n + .0625*(t31*n**(sig1 + &
    & 1.0) + t32*n**(sig2 + 1.0) + t33*n**(sig3 + 1.0)) + &
    & .4521910195*(a+b)*n**(5.0/3.0) + &
    & .4521910195*(1.5*t4)*n**(del+(5.0/3.0)) + &
    & .4521910195*(t5*(2.0*x5+2.5))*n**(gam+(5.0/3.0)) !E(n,1/2)
end function epp

function integrand(zi) result(integrand_array_element) 
!produces an array of densities for the main program to use
!this does a lot of file IO? code may be edited to instead use an imported array
implicit none
    real :: epp
    integer, intent(in) :: zi                     !index of a particular z value
    real :: n0 = .1596                         !SLy4 saturation density, fm**-3
    real :: integrand_array_element
    real, dimension(1:3, 1:100) :: dens_data
    call read_data(dens_data)
    integrand_array_element = 2*(epp(dens_data(2,zi)+dens_data(3,zi))- &
    & epp(n0))*(dens_data(2,zi)+dens_data(3,zi))
!    integrand_array_element = zi*.2
end function integrand


program simpson_integral                    !integrates the integrand array
implicit none


integer :: n                                !number of subintervals
real :: integ = 0.0                         !iterating value of the integral
integer :: i                                !looping variable
real :: subint_width                        !width of subinterval
integer :: lunvar                   !file-opening numbers left to the compiler
real :: ig1, igmid, ig2         !integrand at 3 points in the iterating interval
real :: integrand                           !The integrand function

real :: integ_less_accurate = 0.0 !riemann sum

real:: epp

call integral_info(subint_width, n)

!simpson integral that avoids redundant calls to the function 'integrand'
ig1 = integrand(1)
igmid = integrand(2)
ig2 = integrand(3)
integ = integ + (ig1 + ig2 + 4.0 * igmid) * subint_width / 6.0
!write(*,*) integ
do i=2,(n-1)/2 
ig1 = ig2 
igmid = integrand(2*i)
ig2 = integrand(2*i+1)
integ = integ + (ig1 + ig2 + 4.0 * igmid) * subint_width / 6.0
!write(*,*) ig1, igmid, ig2, integ
enddo

do i=1,n !the riemann sum version of the above answer
integ_less_accurate = integ_less_accurate + .2*integrand(i)
enddo

write(*,*) 'The primary result:'
write(*,*) 'sigma_0:', integ, 'MeV/fm**2'
write(*,*) 'The Riemann sum version (generally within a few % of simpson):'
write(*,*) integ_less_accurate, 'MeV/fm**2'

stop 0
end program simpson_integral