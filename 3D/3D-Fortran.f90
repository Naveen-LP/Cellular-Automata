!################### 3-Dimensional static recrystallization and grain growth ####################

program threeD
implicit none

integer :: nx,ny,nz,i,j,k,steps,l,no,x,iii,ii,jj,kk,x1,y1,z1,stepsg,kink,curv,x2,y2,z2,m,l1,l2
real :: r1,a,b,q,grainid,c,count1,frac
real,dimension(:,:,:),allocatable :: cell,cell1,mn,Hi,CRSS,disl
real*8,dimension(:,:,:),allocatable :: prob,nuclei,velr,dist
integer,dimension(:,:,:),allocatable :: grain,grain1,grain2
real, allocatable, dimension(:) :: r,r2
integer :: values(1:8)
integer, dimension(:), allocatable :: seed
real :: Hc, c0, Qa, Runiv, T, tstep, Sn, epsilonc, ac, bc, gammala, modulus,Qb,boltz
real*8 :: activ,activ1,activ2,activ3,bb,D0
CHARACTER(128) :: buffer,buffer1
integer :: strlen, rows, cols, strlen1, rows1, cols1,asd,asd1
real, dimension(:,:), allocatable :: xid,xid1

!Parameters

 c0 = 1e-2 !Scaling parameter
 Qa = 170e3 !Activation Energy (J/mol)
 Runiv = 8.314 !Universal Gas constant
 T = 1173 !Temperature (K)
 tstep = 0.2 !Length of time step (s)
 Sn = 1 !Volume in whcih nucleas can appear
 epsilonc = 0.4 !critical deformation necessary to trigger nucleation
 ac = 0.1e-8 !Experimental Parameter
 bc = 0.9e-7
 gammala = 0.2 !Low angle grain boundary energy (J/m^2)
 
 Qb = 140e3 !Activation energy to compute velocity for recrystallization
 boltz = 1.38 !Boltzmann constant (10^-23 will be taken care in the expression)
 bb = 2.48e-10 !burges vector
 D0 = 1.1e-6 !Diffusion coefficient

!Modulus G

modulus = 8.1e10*(1 - 0.91*((T - 300)/1810))

nx = 128 !No. of points in x direction
ny = 128 !No. of points in y direction
nz = 128 !No. of points in z direction

steps = 6 !No. of recrystallization steps

allocate(cell(nx+1,ny+1,nz+1),cell1(nx+1,ny+1,nz+1),grain(nx+4,ny+4,nz+4),grain1(nx+4,ny+4,nz+4),grain2(nx+4,ny+4,nz+4))
allocate(prob(nx,ny,nz),nuclei(nx,ny,nz),mn(nx,ny,nz),Hi(nx,ny,nz),CRSS(nx,ny,nz),disl(nx,ny,nz),dist(nx,ny,nz),velr(nx,ny,nz))

!Initializaing all arrays to zero

 cell = 0
 cell1 = 0
 grain = 0
 grain1 = 0
 count1 = 0

!Grain ID data from 18_naveen_new1.out

OPEN (1, file = '128_naveen_new1.out', status='old', action='read')

!Count the number of columns

read(1,'(a)') buffer !read first line WITH SPACES INCLUDED
REWIND(1) !Get back to the file beginning

strlen = len(buffer) !Find the REAL length of a string read
do while (buffer(strlen:strlen) == ' ') 
  strlen = strlen - 1 
enddo

cols=0 !Count the number of spaces in the first line
do ii=0,strlen
  if (buffer(ii:ii) == ' ') then
    cols=cols+1
  endif
enddo

cols = cols+1

!Count the number of rows

rows = 0 !Count the number of lines in a file
DO
  READ(1,*,iostat=ii)
  IF (ii/=0) EXIT
  rows = rows + 1
END DO

REWIND(1)

allocate(xid(rows,cols))

do ii=1,rows,1
  read(1,*) xid(ii,:)
enddo

l1 = 0

do k = 1,128
do j = 1,128
do i = 1,128

l1 = l1 + 1

grain(i,j,k) = xid(l1,7)

end do
end do
end do

 CLOSE (1)

print*, 'Successfully read the grainID data'

grain1 = grain

 open(unit=99,file='grain'//trim(str(1))//'.vtk') 
 write(99,fmt='(A26)')'# vtk DataFile Version 2.0'
 write(99,fmt='(A2)')'# '
 write(99,fmt='(A5)')'ASCII'
 write(99,fmt='(A25)')'DATASET STRUCTURED_POINTS'
 write(99,fmt='(A11,3I5)')'DIMENSIONS ',nx+1,ny+1,nz+1 !npts1+1,npts2+1,npts3+1
 write(99,fmt='(A18)')'ASPECT_RATIO 1 1 1'
 write(99,fmt='(A12)')'ORIGIN 0 0 0'
 write(99,fmt='(A10,I8)')'CELL_DATA ', (nx)*(ny)*(nz) !npts1*npts2*npts3
 write(99,fmt='(A19)')'SCALARS grain int'
 write(99,fmt='(A20)')'LOOKUP_TABLE DEFAULT'

do kk=1,nz
do jj=1,ny

   write(99,fmt='(128(2I5))')grain(1:nx,jj,kk)

end do
end do

!Energy data (Hi) from tex_rolling.out

OPEN (111, file = 'rolling_data.out', status='old', action='read')

rows1 = nx*ny*nz
cols1 = 1

allocate(xid1(rows1,cols1))

do ii=1,rows1
  read(111,*) xid1(ii,:)
enddo

l2 = 0

do k = 1,128
do j = 1,128
do i = 1,128

l2 = l2 + 1

 CRSS(i,j,k) = xid1(l2,1)

end do
end do
end do

 CLOSE (111)

print*, 'Successfully obtained the CRSS data'

 CRSS = CRSS*1e6

activ1 = -Qa/(Runiv*T)
activ = dble((exp(activ1)))

! Calculating the probability

Hc = (epsilonc/(ac*epsilonc + bc))*gammala !Critical energy required to trigger nucleation

do i = 1,nx
do j = 1,ny
do k = 1,nz

Hi(i,j,k) = (CRSS(i,j,k))**2/modulus
mn(i,j,k) = c0*(Hi(i,j,k) - Hc)
nuclei(i,j,k) = mn(i,j,k)*activ
prob(i,j,k) = nuclei(i,j,k)*Sn*tstep !Probability of nucleation for every point

end do
end do
end do

print*, 'Successfully computed the probabilities'

 frac = 0

allocate(r(nx*ny*nz),r2(nx*ny*nz))

!To get different random numbers every time the code is run

call date_and_time(values=values)

call random_seed(size=k)
allocate(seed(1:k))
seed(:) = values(8)
call random_seed(put=seed)
call random_number(r)

no = 0 !No. of grains

l = 0 

m = 2127 !Newly nucleating grains's GrainID

do i = 1,nx
do j = 1,ny
do k = 1,nz

 l = l + 1

 r2(l) = r(l)

if (r2(l) <= prob(i,j,k)) then

 cell(i,j,k) = 1

 m = m + 1

 grain(i,j,k) = m

 no = no + 1

end if

end do
end do
end do

print*, "No. of newly nucleated grains =", no

 open(unit=99,file='grain'//trim(str(2))//'.vtk') 
 write(99,fmt='(A26)')'# vtk DataFile Version 2.0'
 write(99,fmt='(A2)')'# '
 write(99,fmt='(A5)')'ASCII'
 write(99,fmt='(A25)')'DATASET STRUCTURED_POINTS'
 write(99,fmt='(A11,3I5)')'DIMENSIONS ',nx+1,ny+1,nz+1 !npts1+1,npts2+1,npts3+1
 write(99,fmt='(A18)')'ASPECT_RATIO 1 1 1'
 write(99,fmt='(A12)')'ORIGIN 0 0 0'
 write(99,fmt='(A10,I8)')'CELL_DATA ', (nx)*(ny)*(nz) !npts1*npts2*npts3
 write(99,fmt='(A19)')'SCALARS grain int'
 write(99,fmt='(A20)')'LOOKUP_TABLE DEFAULT'

do kk=1,nz
do jj=1,ny

   write(99,fmt='(128(2I5))')grain(1:nx,jj,kk)

end do
end do

!Corrected Moore Neighborhood

a = 0.8 !Probability for Non-diagonal elements
b = 0.25 !Probability for diagonal elements

!Calculation of velocity

activ2 = -Qb/(Runiv*T)
activ3 = dble((exp(activ2)))

do i = 1,nx
do j = 1,ny
do k = 1,nz

velr(i,j,k) = 0.5*D0*(bb**2)*(CRSS(i,j,k)**2)*activ3*(1e23)/(boltz*T*modulus) !Velocity of the recrystallization front

dist(i,j,k) = velr(i,j,k)*tstep*(1e06)

dist(i,j,k) = int(dist(i,j,k))

end do
end do
end do

!Recrystallization step

do x = 1,steps

do i = 1,nx
do j = 1,ny
do k = 1,nz

if (cell(i,j,k) == 1) then

 cell1(i,j,k) = 1
 grain1(i,j,k) = grain(i,j,k)

else if (cell(i-1,j,k) == 1) then

 q = rand()

 if (q <= a) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i-1,j,k)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i+1,j,k) == 1) then

 q = rand()

 if (q <= a) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i+1,j,k)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i,j-1,k) == 1) then

 q = rand()

 if (q <= a) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i,j-1,k)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i,j+1,k) == 1) then

 q = rand()

 if (q <= a) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i,j+1,k)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i+1,j+1,k) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i+1,j+1,k)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i+1,j-1,k) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i+1,j-1,k)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i-1,j-1,k) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i-1,j-1,k)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i-1,j+1,k) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i-1,j+1,k)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i-1,j,k+1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i-1,j,k+1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i+1,j,k+1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i+1,j,k+1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i,j-1,k+1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i,j-1,k+1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i,j+1,k+1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i,j+1,k+1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i+1,j+1,k+1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i+1,j+1,k+1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i+1,j-1,k+1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i+1,j-1,k+1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i-1,j-1,k+1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i-1,j-1,k+1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i-1,j+1,k+1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i-1,j+1,k+1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i,j,k+1) == 1) then

 q = rand()

 if (q <= a) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i,j,k+1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i,j,k-1) == 1) then

 q = rand()

 if (q <= a) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i,j,k-1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i-1,j,k-1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i-1,j,k-1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i+1,j,k-1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i+1,j,k-1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i,j-1,k-1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i,j-1,k-1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i,j+1,k-1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i,j+1,k-1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i+1,j+1,k-1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i+1,j+1,k-1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i+1,j-1,k-1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i+1,j-1,k-1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i-1,j-1,k-1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i-1,j-1,k-1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

else if (cell(i-1,j+1,k-1) == 1) then

 q = rand()

 if (q <= b) then

  cell1(i,j,k) = 1
  grain1(i,j,k) = grain(i-1,j+1,k-1)

  count1 = count1 + 1

 else

  cell1(i,j,k) = 0
  grain1(i,j,k) = grain1(i,j,k)

 end if

end if

grain(i,j,k) = grain1(i,j,k)

end do
end do
end do

 cell = cell1

frac = count1/(nx*ny*nz) !Fraction of cell recrystallized 

!Copying the grain data into a vtk file

 open(unit=99,file='grain'//trim(str(x+2))//'.vtk') 
 write(99,fmt='(A26)')'# vtk DataFile Version 2.0'
 write(99,fmt='(A2)')'# '
 write(99,fmt='(A5)')'ASCII'
 write(99,fmt='(A25)')'DATASET STRUCTURED_POINTS'
 write(99,fmt='(A11,3I5)')'DIMENSIONS ',nx+1,ny+1,nz+1 !npts1+1,npts2+1,npts3+1
 write(99,fmt='(A18)')'ASPECT_RATIO 1 1 1'
 write(99,fmt='(A12)')'ORIGIN 0 0 0'
 write(99,fmt='(A10,I8)')'CELL_DATA ', (nx)*(ny)*(nz) !npts1*npts2*npts3
 write(99,fmt='(A19)')'SCALARS grain int'
 write(99,fmt='(A20)')'LOOKUP_TABLE DEFAULT'

do kk=1,nz
do jj=1,ny

   write(99,fmt='(128(2I5))')grain(1:nx,jj,kk)

end do
end do

end do

print*, 'Recrystallization is over'

print*, 'Grain growth has started'

stepsg = 5 !No. of grain growth steps

grain2 = grain

!Periodic Boundary conditions

do i = 1,nx
do j = 1,ny

grain2(nx+1,i,j) = grain2(1,i,j)
grain2(i,ny+1,j) = grain2(i,1,j)
grain2(i,j,nz+1) = grain2(i,j,1)

grain2(nx+2,i,j) = grain2(2,i,j)
grain2(i,ny+2,j) = grain2(i,2,j)
grain2(i,j,nz+2) = grain2(i,j,2)

grain2(0,j,i) = grain2(nx,j,i)
grain2(i,0,j) = grain2(i,ny,j)
grain2(i,j,0) = grain2(i,j,nz)

grain2(-1,j,i) = grain2(nx-1,j,i)
grain2(i,-1,j) = grain2(i,ny-1,j)
grain2(i,j,-1) = grain2(i,j,nz-1)

end do
end do


do i = 1,nx

grain2(-1,-1,i) = grain2(nx-1,ny-1,i)
grain2(-1,i,-1) = grain2(nx-1,i,nz-1)
grain2(i,-1,-1) = grain2(i,ny-1,nz-1)

grain2(0,0,i) = grain2(nx,ny,i)
grain2(0,i,0) = grain2(nx,i,nz)
grain2(i,0,0) = grain2(i,ny,nz)

grain2(nx+1,ny+1,i) = grain2(1,1,i)
grain2(nx+1,i,nz+1) = grain2(1,i,1)
grain2(i,ny+1,nz+1) = grain2(i,1,1)

grain2(nx+2,ny+2,i) = grain2(2,2,i)
grain2(nx+2,i,nz+2) = grain2(2,i,2)
grain2(i,ny+2,nz+2) = grain2(i,2,2)

grain2(0,-1,i) = grain2(nx,nx-1,i)
grain2(0,i,-1) = grain2(nx,i,nx-1)
grain2(-1,0,i) = grain2(nx-1,nx,i)
grain2(-1,i,0) = grain2(nx-1,i,nx)
grain2(i,0,-1) = grain2(i,nx,nx-1)
grain2(i,-1,0) = grain2(i,nx-1,nx)

grain2(0,nx+1,i) = grain2(nx,1,i)
grain2(0,i,nx+1) = grain2(nx,i,1)
grain2(nx+1,0,i) = grain2(1,nx,i)
grain2(nx+1,i,0) = grain2(1,i,nx)
grain2(i,0,nx+1) = grain2(i,nx,1)
grain2(i,nx+1,0) = grain2(i,1,nx)

grain2(0,nx+2,i) = grain2(nx,2,i)
grain2(0,i,nx+2) = grain2(nx,i,2)
grain2(nx+2,0,i) = grain2(2,nx,i)
grain2(nx+2,i,0) = grain2(2,i,nx)
grain2(i,0,nx+2) = grain2(i,nx,2)
grain2(i,nx+2,0) = grain2(i,2,nx)

grain2(-1,nx+1,i) = grain2(nx-1,1,i)
grain2(-1,i,nx+1) = grain2(nx-1,i,1)
grain2(nx+1,-1,i) = grain2(1,nx-1,i)
grain2(nx+1,i,-1) = grain2(1,i,nx-1)
grain2(i,-1,nx+1) = grain2(i,nx-1,1)
grain2(i,nx+1,-1) = grain2(i,1,nx-1)

grain2(-1,nx+2,i) = grain2(nx-1,2,i)
grain2(-1,i,nx+2) = grain2(nx-1,i,2)
grain2(nx+2,-1,i) = grain2(2,nx-1,i)
grain2(nx+2,i,-1) = grain2(2,i,nx-1)
grain2(i,-1,nx+2) = grain2(i,nx-1,2)
grain2(i,nx+2,-1) = grain2(i,2,nx-1)

grain2(nx+1,nx+2,i) = grain2(1,2,i)
grain2(nx+1,i,nx+2) = grain2(1,i,2)
grain2(nx+2,nx+1,i) = grain2(2,1,i)
grain2(nx+2,i,nx+1) = grain2(2,i,1)
grain2(i,nx+1,nx+2) = grain2(i,1,2)
grain2(i,nx+2,nx+1) = grain2(i,2,1)

end do

grain2(0,0,0) = grain2(nx,nx,nx)
grain2(0,0,-1) = grain2(nx,nx,nx-1)
grain2(0,0,nx+1) = grain2(nx,nx,1)
grain2(0,0,nx+2) = grain2(nx,nx,2)
grain2(0,-1,0) = grain2(nx,nx-1,nx)
grain2(0,-1,-1) = grain2(nx,nx-1,nx-1)
grain2(0,-1,nx+1) = grain2(nx,nx-1,1)
grain2(0,-1,nx+2) = grain2(nx,nx-1,2)
grain2(0,nx+1,0) = grain2(nx,1,nx)
grain2(0,nx+1,-1) = grain2(nx,1,nx-1)
grain2(0,nx+1,nx+1) = grain2(nx,1,1)
grain2(0,nx+1,nx+2) = grain2(nx,1,2)
grain2(0,nx+2,0) = grain2(nx,2,nx)
grain2(0,nx+2,-1) = grain2(nx,2,nx-1)
grain2(0,nx+2,nx+1) = grain2(nx,2,1)
grain2(0,nx+2,nx+2) = grain2(nx,2,2)

grain2(-1,0,0) = grain2(nx-1,nx,nx)
grain2(-1,0,-1) = grain2(nx-1,nx,nx-1)
grain2(-1,0,nx+1) = grain2(nx-1,nx,1)
grain2(-1,0,nx+2) = grain2(nx-1,nx,2)
grain2(-1,-1,0) = grain2(nx-1,nx-1,nx)
grain2(-1,-1,-1) = grain2(nx-1,nx-1,nx-1)
grain2(-1,-1,nx+1) = grain2(nx-1,nx-1,1)
grain2(-1,-1,nx+2) = grain2(nx-1,nx-1,2)
grain2(-1,nx+1,0) = grain2(nx-1,1,nx)
grain2(-1,nx+1,-1) = grain2(nx-1,1,nx-1)
grain2(-1,nx+1,nx+1) = grain2(nx-1,1,1)
grain2(-1,nx+1,nx+2) = grain2(nx-1,1,2)
grain2(-1,nx+2,0) = grain2(nx-1,2,nx)
grain2(-1,nx+2,-1) = grain2(nx-1,2,nx-1)
grain2(-1,nx+2,nx+1) = grain2(nx-1,2,1)
grain2(-1,nx+2,nx+2) = grain2(nx-1,2,2)

grain2(nx+1,0,0) = grain2(1,nx,nx)
grain2(nx+1,0,-1) = grain2(1,nx,nx-1)
grain2(nx+1,0,nx+1) = grain2(1,nx,1)
grain2(nx+1,0,nx+2) = grain2(1,nx,2)
grain2(nx+1,-1,0) = grain2(1,nx-1,nx)
grain2(nx+1,-1,-1) = grain2(1,nx-1,nx-1)
grain2(nx+1,-1,nx+1) = grain2(1,nx-1,1)
grain2(nx+1,-1,nx+2) = grain2(1,nx-1,2)
grain2(nx+1,nx+1,0) = grain2(1,1,nx)
grain2(nx+1,nx+1,-1) = grain2(1,1,nx-1)
grain2(nx+1,nx+1,nx+1) = grain2(1,1,1)
grain2(nx+1,nx+1,nx+2) = grain2(1,1,2)
grain2(nx+1,nx+2,0) = grain2(1,2,nx)
grain2(nx+1,nx+2,-1) = grain2(1,2,nx-1)
grain2(nx+1,nx+2,nx+1) = grain2(1,2,1)
grain2(nx+1,nx+2,nx+2) = grain2(1,2,2)

grain2(nx+2,0,0) = grain2(2,nx,nx)
grain2(nx+2,0,-1) = grain2(2,nx,nx-1)
grain2(nx+2,0,nx+1) = grain2(2,nx,1)
grain2(nx+2,0,nx+2) = grain2(2,nx,2)
grain2(nx+2,-1,0) = grain2(2,nx-1,nx)
grain2(nx+2,-1,-1) = grain2(2,nx-1,nx-1)
grain2(nx+2,-1,nx+1) = grain2(2,nx-1,1)
grain2(nx+2,-1,nx+2) = grain2(2,nx-1,2)
grain2(nx+2,nx+1,0) = grain2(2,1,nx)
grain2(nx+2,nx+1,-1) = grain2(2,1,nx-1)
grain2(nx+2,nx+1,nx+1) = grain2(2,1,1)
grain2(nx+2,nx+1,nx+2) = grain2(2,1,2)
grain2(nx+2,nx+2,0) = grain2(2,2,nx)
grain2(nx+2,nx+2,-1) = grain2(2,2,nx-1)
grain2(nx+2,nx+2,nx+1) = grain2(2,2,1)
grain2(nx+2,nx+2,nx+2) = grain2(2,2,2)

!Grain growth step

do x = 1,stepsg

do i = 1,nx
do j = 1,ny
do k = 1,nz

do x2 = -1,1
do y2 = -1,1
do z2 = -1,1

if (grain2(i,j,k) .NE. grain2(i+x2,j+y2,k+z2)) then !Identifying the grain boundary

grainid = grain2(i,j,k)

!Kink template to find the curvature at every point

kink = 0

do x1 = -2,2
do y1 = -2,2
do z1 = -2,2

if (grain2(i+x1,j+y1,k+z1) == grainid) then

kink = kink + 1

end if

end do
end do
end do

 curv = 75 - kink

 c = 0.6

if (curv <= 0) then

 q = rand()

 if (q <= c) then

  grain2(i+x2,j+y2,k+z2) = grain2(i,j,k)

 end if

end if

end if

end do
end do
end do

end do
end do
end do

!Copying the grain data into a vtk file

 open(unit=90,file='grain'//trim(str(x+steps+2))//'.vtk')
 write(90,fmt='(A26)')'# vtk DataFile Version 2.0'
 write(90,fmt='(A2)')'# '
 write(90,fmt='(A5)')'ASCII'
 write(90,fmt='(A25)')'DATASET STRUCTURED_POINTS'
 write(90,fmt='(A11,3I5)')'DIMENSIONS ',nx+1,ny+1,nz+1 !npts1+1,npts2+1,npts3+1
 write(90,fmt='(A18)')'ASPECT_RATIO 1 1 1'
 write(90,fmt='(A12)')'ORIGIN 0 0 0'
 write(90,fmt='(A10,I8)')'CELL_DATA ', (nx)*(ny)*(nz) !npts1*npts2*npts3
 write(90,fmt='(A19)')'SCALARS grain int'
 write(90,fmt='(A20)')'LOOKUP_TABLE DEFAULT'

do kk=1,nz
do jj=1,ny

   write(90,fmt='(128(2I5))')grain2(1:nx,jj,kk)

end do
end do

end do

print*, 'Grain growth is over'

contains

!Function used to copy the data to a new vtk file, every iteration

character(len=20) function str(zz)

  integer, intent(in) :: zz
  write (str, *) zz
  str  = adjustl(str)

end function str

end program
