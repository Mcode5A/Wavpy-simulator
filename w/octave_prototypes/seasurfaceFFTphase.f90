subroutine seasurfaceFFTphase(x0coord,y0coord,kph,Lsurf,U10,Ome,lambdaS,swh,thetaDW,thetaDS)!,rx,ry,surf2)!nr2
  implicit none
!-----------------------------
! Routine to obtain the sea surface
! elevation profile h(x,y) for
! a sea surface spectrum.
!---------------------------------
! USAGE:
!  seasurfaceFFTphase(xoffset,yoffset,extraphase,surface)
!
!  extraphase is an integer number to
!   alterate the random phase
!---------------------------------
! E. Cardellach 2002/2003
! MODIFICATIONS:
!       12/12/02: input argument
!       04/28/03: random phase alterable
!       03/01/07: FFT
!       03/01/07: spectrum in the code
!       20/06/07: as routine to be called in soscf
!       18/03/11: swell added, from independent spectrum
!---------------------------------
! About Inverse Wave Age, Ome=U10/cp (there are other definitions!):
! fully developed sea: 0.84
! mature sea: 1.5
! young seas: 2-3
! From: Makin et al., Drag of the Sea Surface,
!       Boundary Layer Meteorology, 73, pp:159-182, 1995
!----------------------------------
  integer*4, intent(in) :: kph,Lsurf!,nr2
  real*8, intent(in) :: x0coord,y0coord,U10,Ome,lambdaS,swh,thetaDW,thetaDS
  !REAL*8, DIMENSION(:), INTENT(INOUT) :: rx,ry
  !REAL*8, DIMENSION(:,:), INTENT(INOUT) :: surf2
  integer*4 :: i,j,ii,k,ix,iy,ixf,iyf,thekph,thekph2
  real*8, allocatable, dimension(:,:) :: psi,phextra,z
  real*8, allocatable, dimension(:) :: kx,ky
  real*8 deltak, lres, ave, sdev, var
  real*8 x0,y0
  complex sumd
  complex, allocatable, dimension(:,:) :: h,temp,aux
  real*8 pi
  integer, parameter :: nr2 = 4096
  real*8 :: surfs(nr2,nr2)
  real :: surfw(nr2,nr2)
  real :: surf2(nr2,nr2),rx(nr2),ry(nr2)
  real*8 Rearth
  complex m_i
  common/surface2/rx,ry,surf2
  common/surface/Rearth
  pi  = 3.1415926535897932384626432
  m_i = (0., 1.)
  write(*,*) 'sono in seasurface, it will take some time now...'
  !-------------------------------------
  allocate(psi(nr2,nr2))
  allocate(phextra(nr2,nr2))
  allocate(z(nr2,nr2))
  allocate(kx(nr2))
  allocate(ky(nr2))
  allocate(h(nr2,nr2))
  allocate(temp(nr2,nr2))
  allocate(aux(nr2,nr2))
  !-------------------------------------------------------------
  h(:,:)=(0., 0.)
  ! if this patch starts at x=L, the origin must be set to x0=-L
  x0=-x0coord 
  y0=-x0coord 

  surfw(:,:)=0.0d0
  surfs(:,:)=0.0d0
  !--------- GENERATING WIND SPECTRA
  psi(:,:)=0.0d0
  kx(:)=0.0d0
  ky(:)=0.0d0
  phextra(:,:)=0.0d0
  deltak=2.0*pi/Lsurf       ! 0.015d0 old fixed number
                            ! now, for Lsurf=4096-> 0.001533...
                            !      for Lsurf=8192-> 0.000766...
  if (U10 .lt. 0.1) goto 100
  write(*,*) 'in seasurface... U10 thetaDW Ome',U10,thetaDW,Ome
  call myspectrum(psi,kx,ky,deltak,nr2/2,U10,thetaDW,Ome)
  !write(*,*) 'sono in myspecrtum'
  thekph=kph
  do i=1,nr2
  do ii=1,nr2
     k=i+ii
     thekph=mod(thekph+k,16777216)
     thekph2=1+int(ran(thekph)*10.0)
     phextra(i,ii)=ran(thekph2)*2.0*pi
     ! now computing the kernel of the transform, 
     ! note that sign of exponentials is "-" for phase0 but
     ! it's "+" for offset effect: so it works like (x-x0)*kx 
     rx(i)=(i-1)*2.0*pi/(nr2*deltak)+x0coord
     ry(ii)=(ii-1)*2.0*pi/(nr2*deltak)+y0coord
     h(i,ii)=sqrt(psi(i,ii))*exp(-m_i*phextra(i,ii))*exp(-m_i*kx(i)*rx(i))*exp(-m_i*ky(ii)*ry(ii))
  enddo
  enddo

  !-------------------------------------------
  ! FFT TO GENERATE THE WIND SURFACE z(x,y)
  !-------------------------------------------
  call fourrow(h,-1)
  !write(*,*) 'sono in FFT'
  temp=transpose(h)
  call fourrow(temp,-1)
  h=transpose(temp)
  
  do ix=1,nr2
    do iy=1,nr2
        !z=Real{InvFT[amplitude exp(extraphase)exp(x0,y0)exp(k*x)]}
        !  where amplitude=sqrt(2*Psi*deltak**2)
        !        amplitude=deltak*sqrt(2*Psi)
        ! we have computed the InvFT with sqrt(Psi), so
        ! we need to add: deltak*sqrt(2)
        !z(ix,iy) = deltak*sqrt(2.0)* real(h(ix,iy))
        !rx(ix)=(ix-1)*2.0*pi/(nr2*deltak)+x0coord
        !ry(iy)=(iy-1)*2.0*pi/(nr2*deltak)+y0coord
        surfw(ix,iy)=deltak*sqrt(2.0)* real(h(ix,iy))
     enddo
  enddo
  ! now h contains z(x,y). In just two more inverse FFT we can
  ! get the partial derivatives with respect to x-axis and y-axis
  ! FT(dh/dx)=2Pi kx FT(h)  --> dh/dx = IFT( 2 Pi kx FT(h))
  ! same for y-component
  ! ... to look carefully at this, will we save comp.time
  ! computing the derivatives this way rather than later on
  ! over each integration point? I guess so...
100 continue
  if (swh .lt. 0.01) goto 200
  ! NOW SWELL SPECTRUM
  psi(:,:)=0.0
  lres=2.0*pi/(nr2*deltak)
  call swell_piersmos(nr2,lres,lres,lambdaS,thetaDS,psi)
  ! shifting spectrum, so it's given following FFTs sorting (positive freq. first...)

  temp(:,:)=(0.0, 0.0)
  h(:,:)=(0.0, 0.0)

  do i=1,nr2
  do j=1,nr2
  if (j .gt. int(nr2/2.0) ) then
  temp(i,j-int(nr2/2.0))=cmplx(psi(i,j))
  else
  temp(i,j+int(nr2/2.0))=cmplx(psi(i,j))
  endif
  enddo
  enddo
  do j=1,nr2
  do i=1,nr2
  if (i .gt. int(nr2/2.0) ) then
  h(i-int(nr2/2.0),j)=temp(i,j)
  IF (ISNAN(real(h(i-int(nr2/2.0),j)))) THEN
  h(i-int(nr2/2.0),j)=0.0
  ENDIF
  else
  h(i+int(nr2/2.0),j)=temp(i,j)
  IF (ISNAN(real(h(i+int(nr2/2.0),j)))) THEN
  h(i+int(nr2/2.0),j)=0.0
  ENDIF
  endif
  enddo
  enddo
  ! now h contains the spectrum, in a complex array (although imaginary parts are all 0)
  temp(:,:)=0.0


  ! NOW THE SWELL SURFACE
  psi(:,:)=0.0
  aux(:,:)=(0.0, 0.0)
  thekph=3
  call RANDOM_SEED
  CALL RANDOM_SEED(SIZE=thekph)
  CALL RANDOM_NUMBER(psi)
  psi(:,:)=psi(:,:)*2.0*pi-pi
  aux=cmplx(psi)
  ! forward 2D fft for f:
  !f=fft(r,-1) (IDL code, OPPOSITE sign convention!)
  call fourrow(aux,1)
  temp=transpose(aux)
  call fourrow(temp,1)
  aux=transpose(temp)
  temp(:,:)=0.0

  do i=1,nr2
  do j=1,nr2
    temp(i,j)=2.0*aux(i,j)/abs(aux(i,j))*h(i,j)
    aux(i,j)=temp(i,j)
    !f=2.*f/abs(f)*shift(sp,dim/2,dim/2) (IDL code)
  enddo
  enddo

  ! now generating the surface
  !s=double(fft(f,1)) (IDL code, FFT OPPOSITE sign convention!)
  temp(:,:)=0.0
  call fourrow(aux,-1)
  temp=transpose(aux)
  call fourrow(temp,-1)
  surfs=real(transpose(temp))

  ! renormalize to obtain user desired SWH
  call moment(surfs,nr2,nr2,ave,sdev,var)
  surfs=surfs/(4.0*sdev)*swh

200 continue
  ! NOW FINAL SURFACE

  surf2=surfw+surfs ! total surface is the sum of wind and swell...
  do ix=1,nr2
    do iy=1,nr2
        rx(ix)=(ix-1)*2.0*pi/(nr2*deltak)+x0coord
        ry(iy)=(iy-1)*2.0*pi/(nr2*deltak)+y0coord
    enddo
  enddo


  return 
contains
!--------MANY SUBROUTINES-----------------------------
!--------subroutine ran:
function ran(idum)
implicit none
integer, intent(inout) :: idum
real*8                :: ran
integer, parameter :: ia=16807, im=2147483647, iq=127773, ir=2836
real*8, save :: am
integer, save :: ix=-1, iy=-1,l
if (idum <= 0 .or. iy < 0) then
   am=nearest(1.0,-1.0)/im
   iy=ior(ieor(888889999,abs(idum)),1)
   ix=ieor(777755555,abs(idum))
   idum=abs(idum)+1
endif
ix=ieor(ix,ishft(ix,13))
ix=ieor(ix,ishft(ix,-17))
ix=ieor(ix,ishft(ix,5))
l=iy/iq
iy=ia*(iy-l*iq)-ir*l
if (iy < 0) iy=iy+im
ran=am*ior(iand(im,ieor(ix,iy)),1)
end function ran
!----------------------------------------------------
SUBROUTINE fourrow(data,isign)
  IMPLICIT NONE
  COMPLEX, DIMENSION(:,:), INTENT(INOUT) :: data
  INTEGER, INTENT(IN) :: isign
  INTEGER :: n,i,istep,j,m,mmax,n2
  REAL*8 :: theta
  COMPLEX, DIMENSION(size(data,1)) :: temp
  COMPLEX :: w,wp
  COMPLEX :: ws
  REAL*8, PARAMETER :: pi  = 3.1415926535897932384626432
  n=size(data,2)
  call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow')
  n2=n/2
  j=n2
  do i=1,n-2
     if (j > i) call swap(data(:,j+1),data(:,i+1))
     m=n2
     do
          if (m < 2 .or. j < m) exit
          j=j-m
          m=m/2
     end do
     j=j+m
  end do
  mmax=1
  do
    if (n <= mmax) exit
    istep=2*mmax
    theta=pi/(isign*mmax)
    wp=cmplx(-2.0d0*sin(0.5*theta)**2,sin(theta))
    w=cmplx(1.0d0,0.0d0)
    do m=1,mmax
      ws=w
      do i=m,n,istep
         j=i+mmax
         temp=ws*data(:,j)
         data(:,j)=data(:,i)-temp
         data(:,i)=data(:,i)+temp
      end do
      w=w*wp+w
    end do
    mmax=istep
  end do
END SUBROUTINE fourrow
!--------------------------------------------------
SUBROUTINE swap(a,b)
COMPLEX, DIMENSION(:), INTENT(INOUT) :: a,b
COMPLEX, DIMENSION(SIZE(a)) :: dum
 dum=a
 a=b
 b=dum
END SUBROUTINE swap
!--------------------------------------------
SUBROUTINE assert(n1,string)
CHARACTER(LEN=*), INTENT(IN) :: string
LOGICAL, INTENT(IN) :: n1
if (.not. n1) then
  write (*,*) 'nrerror: an assertion failed with this tag:', &
      string
  STOP 'program terminated by assert1'
end if
END SUBROUTINE assert
!-------------------------------------------

subroutine myspectrum(psi,kx,ky,deltak,n,U10,thetaD,Ome)
!--------------------------------------------
! For a given set of sea state conditions,
! this program generates the 
! OMNIDIRECTIONAL sea spectra, S(k),
! as suggested in [Elfouhaily et al., 1997] 
! Also the DIRECTIONAL spectrum using 'unified'
! 'angular spreading functions'.
! 
! ASSUMPTIONS:
!   -Deep waters (c=sqrt(g/k))
!--------------------------------------------
!
! psi(n,n)=spectra (IN empty/OUT)
! kx(n),ky(n)=arrays of kx and ky (IN empty/OUT)
! deltak=step in kx and ky (SAME)
! n= half-length of spectra (i.e.length of pos or/and negative)
! U10: wind speed at 10m above the surface [m/s]
! theta: angle between wind and waves [deg]
! omega: wave age > 0.83
!---------------------------------------------
! Estel Cardellach, 12/17/2002
!                 , 02/03/2007, now subroutine
!---------------------------------------------
  implicit none
  real*8, INTENT(IN) :: U10,thetaD,Ome,deltak
  integer, INTENT(IN) :: n ! n is half-length of spectra (pos/neg)
  real*8, dimension(:,:), INTENT(INOUT) :: psi
  real*8, dimension(:), INTENT(INOUT) :: kx,ky
  real*8, parameter :: g = 9.81     ! gravitational acceleration
  integer i,ii,j,ix,iy,length_u,length_t,length_o
  real*8 S,Bl,Bh,Del,Phi,theta
  real*8 k,c,cp,alphap,kp,Jp,gamma
  real*8 cm,alpham,km,ufric,sigma2
  real*8 xp, xm, yp, ym
  real*8 a0,ap,am,ph

  !---------- Get arguments:
  theta = thetaD * 3.1416 / 180.0
  
  !---------- Wind dependent parameters, k-independent  
  cp = U10/(Ome*1.0)
  kp = g/(cp**2)     !**ATTENTION: assuming Deep Waters**  
  km = 370.0
  !ufric = sqrt(1e-3*(1.10+0.04*U10))*U10  !** from web page reference AOMIP
  ufric = sqrt(1e-3*(0.81+0.065*U10))*U10  !** from Elfouhaily's code
  alphap = 6e-3 * sqrt(Ome)
  gamma = gammaf(Ome,theta)
  sigma2 = (0.08*(1.0+4.0*(Ome*cos(theta))**(-3)))**2
  cm = 0.23
  alpham = alphamf(ufric,cm)
  !for spreading functions...
  a0=log(2.0)/4.0
  ap=4.0
  am=0.13*ufric/cm
  !------------------ k=0
  psi(:,:)=0.0d0
  psi(1,1)=0.0d0
  !------------------LOOP IN kx
  do i=1,2*n
     if (i.le.(n+1)) then
       ix=i
       kx(ix)=(i-1)*deltak
     else
       j=i-n-1
       ix=2*n-(j-1)
       kx(ix)=-j*deltak
     endif
     !---------------LOOP IN ky
     do ii=1,2*n
        if ((ii.eq.1).and.(i.eq.1)) goto 5
        if (ii.le.(n+1)) then
          iy=ii
          ky(iy)=(ii-1)*deltak
        else
          j=ii-n-1
          iy=2*n-(j-1)
          ky(iy)=-j*deltak
        endif
        k=sqrt(kx(ix)**2+ky(iy)**2)
        ph=atan2(ky(iy),kx(ix))
        !--- now omnidirectional
        xp = sqrt(k/kp)
        xm = k/km
        c = sqrt(g/k*(1.0+xm*xm))
        yp = Ome/sqrt(10.0)*(xp-1.0)
        Bl = 0.50*alphap*cp/c*exp(-yp)
        ym = (xm - 1.0)**2/4.0
        Bh = 0.50*alpham*cm/c*exp(-ym)
        Jp = jonswap(k,kp,sigma2,gamma)
        S = (Bl+Bh)*Jp/k**3.0
        ! S is omnidirectional spectra
        !--- now directional factor
        Del = tanh(a0 + ap*(c/cp)**2.5 + am*(cm/c)**2.5)
        
        !--- TOTAL
        Phi = 1.0/(2.0*3.1416)*(1.0+Del*cos(2.0*ph))
        psi(ix,iy) = 1.0/k *S*Phi
        write(*,*) kx(ix),ky(iy),psi(ix,iy)
5       continue 
     end do
6 end do
end subroutine myspectrum

!--------------------------- 
  function gammaf(O,th)
  ! gamma function INPUT PARAMETERS: Ome and theta (rad)
    real*8, intent(in) :: O
    real*8, intent(in) :: th
    real*8             :: Oc
    real*8             :: gammaf
    Oc = O*cos(th)
    if ((Oc >= 0.84) .and. (Oc < 1.0)) gammaf = 1.7
    if ((Oc >= 1.0) .and. (Oc < 5.0)) gammaf = 1.7 + 6.0*log10(Oc)
  end function gammaf
!-----------------------
  function alphamf(uf,cm)
  ! alpham function INPUT PARAMETERS: ufric and cm
  ! cm in Elfouhaily's formula is constant, but other
  ! authors consider it differently...that's why this
  ! parameter is entered as in the function
    real*8, intent(in) :: uf
    real*8, intent(in) :: cm
    real*8             :: alphamf
    if (uf < cm) alphamf = 1e-2*(1.0+log(uf/cm))
    if (uf > cm) alphamf = 1e-2*(1.0+3.0*log(uf/cm))
  end function alphamf
!-------------------------
  function jonswap(x,kkp,s2,gmm)
    implicit none
    real*8, intent(in) :: x,kkp,s2,gmm
    real*8             :: z,y,w
    real*8             :: jonswap
    y = exp(-5.0/4.0 * (kkp/x) * (kkp/x))
    w = sqrt(x/kkp) - 1.0
    z = exp(-0.5*w*w/s2)
    jonswap = y*gmm**z
  end function jonswap

subroutine swell_piersmos(ndim,resolr,resola,peak,dir,sp)
!-----------------------------------------------------------------
! fully developed seas' spectrum, here used for swell
! based on IDL code provided by Chapron (IFREMER)
!------------------------------------------------------------------
! Spectre Houle
! Pierson-Moskowitz
!------------------------------------------------------------------
implicit none
integer*4, INTENT(IN) :: ndim
real*8, INTENT(IN) :: resolr, resola, peak, dir
real*8, dimension(:,:), INTENT(INOUT) :: sp
integer*4 :: i,j
real*8 :: kmax,g, dir2
real*8, allocatable, dimension(:) :: kr,ka
real*8, allocatable, dimension(:,:) :: a,k
real*8, parameter :: pi = 3.1415926535897932384626432
allocate(kr(ndim))
allocate(ka(ndim))
allocate(a(ndim,ndim))
allocate(k(ndim,ndim))

! correct for direction (as it is, dir=0 produce propagation along Y-axis)
dir2=dir+90.0

do i=1,ndim
kr(i)=((i-1)-ndim/2.0)/(ndim*resolr)*2.*pi
ka(i)=((i-1)-ndim/2.0)/(ndim*resola)*2.*pi
enddo

kr(int(ndim/2.0)) = 1e-10
ka(int(ndim/2.0)) = 1e-10

g=9.81
kmax=2.0*pi/peak

do i=1,ndim
do j=1,ndim
a(i,j)=atan2(ka(i),kr(j))+dir2*pi/180.0
k(i,j)=sqrt(kr(j)**2+ka(i)**2)
sp(i,j)=0.0
if ( (a(i,j).lt.pi/2.0) .and. (a(i,j).gt.-pi/2.0) ) then
sp(i,j) = 0.0081/4.0*k(i,j)**(-4) * exp( -2.0*kmax**2/k(i,j)**2 )/(3.0*pi)*(cos(a(i,j)))**4
endif
enddo
enddo

end subroutine swell_piersmos

subroutine moment(dat,n1,n2,ave,sdev,var)
IMPLICIT NONE
INTEGER*4, INTENT(IN) :: n1, n2
REAL*8, dimension(:,:), INTENT(IN) :: dat
REAL*8, INTENT(OUT) :: ave,sdev,var
!Given an array of data, this routine returns its mean ave, standard
!deviation sdev, variance var, skewness skew, and kurtosis curt.
INTEGER*4 :: n
REAL*8 :: ep
REAL*8, allocatable, DIMENSION(:,:) :: p,s
allocate(s(n1,n2))
allocate(p(n1,n2))
s(:,:)=0.0
p(:,:)=0.0

n=n1*n2
if (n <= 1) then
write(*,*) 'moment: n must be at least 2'
stop
endif

ave=sum(dat(:,:))/n  !First pass to get the mean.
s=dat-ave !Second pass to get the first (absolute), second, third, and
ep=sum(s(:,:)) ! fourth moments of the deviation from the mean.
p(:,:)=s(:,:)*s(:,:)
var=sum(p(:,:))

var=(var-ep**2/n)/(n-1) !Corrected two-pass formula.
sdev=sqrt(var)
END subroutine moment


end subroutine seasurfaceFFTphase

