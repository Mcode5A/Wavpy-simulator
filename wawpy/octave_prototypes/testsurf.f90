program testsurf
implicit none
integer*4 :: kph,Lsurf, i, j
real*8 :: x0coord,y0coord,U10,Ome,lambdaS,swh,thetaDW,thetaDS
integer*8, parameter :: nr2 = 4096
real :: surf2(nr2,nr2),rx(nr2),ry(nr2)
common/surface2/rx,ry,surf2


x0coord=0.0
y0coord=0.0
kph=7
Lsurf=4096
U10=4.0
Ome=0.84
lambdaS=400.0
swh=1.0
thetaDW=0.0
thetaDS=50.0

rx(:)=0.0
ry(:)=0.0
surf2(:,:)=0.0

call seasurfaceFFTphase(x0coord,y0coord,kph,Lsurf,U10,Ome,lambdaS,swh,thetaDW,thetaDS)

! do i=1,nr2
! do j=1,nr2
! write(*,*) rx(i),ry(j),surf2(i,j)
! enddo
! enddo

end program testsurf

