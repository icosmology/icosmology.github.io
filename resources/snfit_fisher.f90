module precision
  integer, parameter :: dl = kind(1.d0)
end module precision

module cosmo
use precision
real(dl) :: om, ol, ok, H0=70.d0
real(dl),parameter :: om_fid=0.303d0, ol_fid=0.788d0, delta=0.05d0
real(dl) :: c =3.d5
real(dl) :: DH
real(dl), parameter :: om_min=0.d0, om_max=3.d0
real(dl), parameter :: ol_min=0.d0, ol_max=3.d0
integer,  parameter :: Np=100
end module cosmo


program test
use precision
use cosmo
implicit none

integer i, j, iz, iom, iol
integer, parameter :: n=580

real(dl)    :: z_dat(n), mu_dat(n), invcov(n,n)
real(dl)    :: mu_om_plus(n),mu_om_minus(n)
real(dl)    :: mu_ol_plus(n),mu_ol_minus(n)
real(dl)    :: FD(2,n), Fisher(2,2), cov22(2,2), corr22(2,2)  
real(dl)    :: diff(n), chi2  
character   :: dummyc
real(dl)    :: dummyr

real(dl)    :: z

real(dl),external    :: E2, Einv, DC, DA, mu 
logical :: flag=.false.

character*256 :: root
real(dl) :: center(2)

open(unit=50, file='sn_z_mu_dmu_plow_union2.1.txt')

do i=1, n
 read(50,*) dummyc, z_dat(i), mu_dat(i), dummyr, dummyr
 !write(*,*) z_dat(i), mu(i)
end do

close(50)

open(unit=50, file='sn_wmat_nosys_union2.1.txt')

do i=1, n
  read(50,*) invcov(i,:)
  !write(*,*) invcov(i,i)
end do

close(50)

open(unit=50, file='test_mu.dat')
do i=1, n
  write(50,*) z_dat(i), mu_dat(i), 1.d0/sqrt(invcov(i,i))
end do
close(50)

DH=c/H0


!--------------Fisher ------------------------

om = om_fid*(1.d0+delta)
ol = ol_fid
ok = 1.d0-om-ol

do i=1, n
  mu_om_plus(i) = mu(z_dat(i))
end do


om = om_fid*(1.d0-delta)
ol = ol_fid
ok = 1.d0-om-ol

do i=1, n
  mu_om_minus(i) = mu(z_dat(i))
end do


om = om_fid
ol = ol_fid*(1.d0+delta)
ok = 1.d0-om-ol

do i=1, n
  mu_ol_plus(i) = mu(z_dat(i))
end do


om = om_fid
ol = ol_fid*(1.d0-delta)
ok = 1.d0-om-ol

do i=1, n
  mu_ol_minus(i) = mu(z_dat(i))
end do

FD(1,:) = (mu_om_plus(:)-mu_om_minus(:))/(2.d0*om_fid*delta)
FD(2,:) = (mu_ol_plus(:)-mu_ol_minus(:))/(2.d0*ol_fid*delta)

Fisher = matmul(matmul(FD,invcov),transpose(FD))

write(*,*) 'The Fisher matrix is,'

do i=1, 2
 write(*,*) Fisher(i,:)
end do

write(*,*) '----------------------------'

call inverse(Fisher,cov22,2)

write(*,*) 'The parameter covariance matrix is,'

do i=1, 2
 write(*,*) cov22(i,:)
end do

write(*,*) '----------------------------'
do i=1, 2
do j=1, 2
  corr22(i,j) = cov22(i,j)/sqrt(cov22(i,i))/sqrt(cov22(j,j))
end do
end do

write(*,*) 'The parameter correlation matrix is,'

do i=1, 2
 write(*,*) corr22(i,:)
end do

write(*,*) '----------------------------'

write(*,*) 'The error of Omega_M is', sqrt(cov22(1,1))
write(*,*) 'The error of Omega_L is', sqrt(cov22(2,2))

write(*,*) '----------------------------'

root = 'Fisher'

center(1) = om_fid
center(2) = ol_fid

call MakeContour(root,center,2,cov22,1,2,1)
call MakeContour(root,center,2,cov22,1,2,2)
call MakeContour(root,center,2,cov22,1,2,3)

!stop 
!--------------Fisher ------------------------


open(unit=50, file='chi2.dat')

do iom=1, Np

 om = om_min + (om_max-om_min)*dble(iom-1)/dble(Np-1) 
 
do iol=1, Np

 ol = ol_min + (ol_max-ol_min)*dble(iol-1)/dble(Np-1)
 
 ok = 1.d0-om-ol
 
 call check(flag)
 
 if(flag) then 
 
  do iz=1, n 
   diff(iz)  = mu(z_dat(iz))-mu_dat(iz)   
  end do ! loop for z
 
  chi2 = dot_product(diff,matmul(invcov,diff))
 
 else
 
  chi2 = 1.d10
 
 end if 
 
 write(50,*) om, ol, chi2
 
 !om = om_min * (om_max/om_min)**(dble(i-1)/dble(Np-1))

end do ! loop for ol 
end do ! loop for om


!z=0.5d0
!write(*,'(5e15.6)') z, E2(z), Einv(z), DC(z), DA(z), mu(z)

close(50)


end program test


subroutine check(flag)
use precision
use cosmo
implicit none

real(dl), parameter :: z_min=0.d0, z_max=3.d0
integer :: iz
integer, parameter :: nz=500
real(dl) :: z
real(dl), external :: E2, DA
logical :: flag

do iz=1, nz

 z = z_min +(z_max-z_min)*dble(iz-1)/dble(nz-1)
 
 if(E2(z)<0 .or. E2(z)==0 .or. DA(z)<0) then !! Thanks Yonghao
   flag = .false.
   return
 else  
   flag = .true.
 end if
 
end do

end subroutine check


function E2(z)
use precision
use cosmo
implicit none

real(dl) :: z, E2

E2 = om*(1.d0+z)**3.d0+ol+(1.d0-om-ol)*(1+z)**2.d0

end function E2

function Einv(z)
use precision
use cosmo
implicit none

real(dl) :: z, Einv
real(dl), external :: E2

Einv = 1.d0/sqrt(E2(z))

end function Einv

function DC(z)
use precision
use cosmo
implicit none

real(dl), external ::  rombint, Einv
real(dl), parameter :: tol=1.d-4
real(dl) :: DC, z

DC = rombint(Einv,0.d0,z,tol)
DC = DC*DH

end function DC

function DA(z)
use precision
use cosmo
implicit none

real(dl) :: DA, z, x
real(dl),external :: DC

if(ok>0) then

x= sqrt(ok)*DC(z)/DH

if(x<200.d0) then

DA = DH/(1.d0+z)/sqrt(ok)*sinh(sqrt(ok)*DC(z)/DH)

else

DA = 1.d0

end if

else if(ok<0) then

DA = DH/(1.d0+z)/sqrt(abs(ok))*sin(sqrt(abs(ok))*DC(z)/DH)

else

DA = DC(z)/(1.d0+z)

end if


end function DA


function Ldis(z)
use precision
use cosmo
implicit none

real(dl) :: Ldis, z
real(dl), external :: DA

Ldis = DA(z)*(1.d0+z)**2

end function Ldis


function mu(z)
use precision
use cosmo
implicit none

real(dl) :: mu, z
real(dl), external :: Ldis

 mu=5.d0*log10(Ldis(z))+25.d0

end function mu


       function rombint(f,a,b,tol)
        use Precision
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        implicit none
        integer, parameter :: MAXITER=20
        integer, parameter :: MAXJ=5
        dimension g(MAXJ+1)
        real(dl) f
        external f
        real(dl) :: rombint
        real(dl), intent(in) :: a,b,tol
        integer :: nint, i, k, jmax, j
        real(dl) :: h, gmax, error, g, g0, g1, fourj
!

        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
            go to 40
!  Calculate next trapezoidal rule approximation to integral.
          g0=0._dl
            do 20 k=1,nint
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1._dl
            do 30 j=1,jmax
!  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1._dl-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)  then
          write(*,*) 'Warning: Rombint failed to converge; '
          write (*,*)'integral, error, tol:', rombint,error, tol
        end if
        
        end function rombint

		
		
subroutine MakeContour(root,center,Np,Matrix,p,q,nsigma)
use Precision 
implicit none
integer Np
real(dl) :: Matrix(Np,Np)
real(dl) :: center(Np)
integer,intent(in) ::p,q, nsigma
integer M
real(dl),dimension(:,:),allocatable :: cov
integer i,j,k,ij
integer, parameter :: N=100
real(dl) theta, pi
real(dl),dimension(:,:),allocatable::x,y,contour
real(dl),dimension(:),allocatable::sigx,sigy,covxy,r,a2,b2,a,b
character*256 :: root

M= q-p+1

pi=4.0d0*datan(1.0d0)
ij=0
allocate(x(M*(M-1)/2,N))
allocate(y(M*(M-1)/2,N))
allocate(sigx(M*(M-1)/2))
allocate(sigy(M*(M-1)/2))
allocate(covxy(M*(M-1)/2))
allocate(r(M*(M-1)/2))
allocate(a2(M*(M-1)/2))
allocate(b2(M*(M-1)/2))
allocate(a(M*(M-1)/2))
allocate(b(M*(M-1)/2))

allocate(contour(M*(M-1),N))
allocate(cov(M,M))

do i=1,M
do j=1,M
if (nsigma==1) cov(i,j)=Matrix(p+i-1,p+j-1)*(1.5d0**2)
if (nsigma==2) cov(i,j)=Matrix(p+i-1,p+j-1)*(2.447d0**2)
if (nsigma==3) cov(i,j)=Matrix(p+i-1,p+j-1)*(3.4d0**2)
if (nsigma==4) cov(i,j)=Matrix(p+i-1,p+j-1)*(4.29d0**2)
end do
end do


do i=1,M-1
do j=i+1,M

ij=ij+1

sigx(ij)=cov(i,i)
sigy(ij)=cov(j,j)
covxy(ij)=cov(i,j)

r(ij)=0.5d0*atan(-2.d0*covxy(ij)/(sigx(ij)-sigy(ij)))

a2(ij)=sigx(ij)*sin(r(ij))**2+2.d0*covxy(ij)*sin(r(ij))*cos(r(ij))+sigy(ij)*cos(r(ij))**2
b2(ij)=sigx(ij)*cos(r(ij))**2-2.d0*covxy(ij)*sin(r(ij))*cos(r(ij))+sigy(ij)*sin(r(ij))**2

a(ij)=sqrt(a2(ij))
b(ij)=sqrt(b2(ij))


do k=1, N
theta = 2.d0*pi/(real(N)-1.d0)*(real(k)-1.d0)

x(ij,k) =  b(ij)*sin(theta)*cos(r(ij)) + a(ij)*cos(theta)*sin(r(ij))
y(ij,k) =  a(ij)*cos(theta)*cos(r(ij)) - b(ij)*sin(theta)*sin(r(ij))

x(ij,k)=x(ij,k)+center(p+i-1)
y(ij,k)=y(ij,k)+center(p+j-1)
!write(*,'(4i5,5e18.9)') i,j,k,ij,sigx(ij),sigy(ij),covxy(ij),x(ij,k),y(ij,k)
end do

end do 
end do

if (nsigma==1) open(unit=50,file=trim(root)//'_1sigma.dat')
if (nsigma==2) open(unit=50,file=trim(root)//'_2sigma.dat')
if (nsigma==3) open(unit=50,file=trim(root)//'_3sigma.dat')
if (nsigma==4) open(unit=50,file=trim(root)//'_4sigma.dat')


do k=1,N
do i=1,M*(M-1)/2
contour(2*i-1,k)=x(i,k)
contour(2*i,k)=y(i,k)
end do
end do

do k=1,N
write(50,'(500e15.6)') contour(:,k)
end do



deallocate(x)
deallocate(y)
deallocate(sigx)
deallocate(sigy)
deallocate(covxy)
deallocate(contour)
deallocate(r)
deallocate(a2)
deallocate(b2)
deallocate(a)
deallocate(b)

close(50)


end subroutine MakeContour




subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse		
