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
