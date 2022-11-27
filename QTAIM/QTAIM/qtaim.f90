pure real*8 function mybinomial(n, k)

  implicit none

  integer*4,intent(in) :: n
  integer*4,intent(in) :: k

  mybinomial=gamma(real(n+1))/(gamma(real(k+1))*gamma(real(n-k+1)))

end function mybinomial


pure real*8 function hermite(n, x)

  implicit none

  integer*4, intent(in) :: n
  real*8, intent(in) :: x

  if (n==0) then
     hermite = 1
  else if (n==1) then
     hermite = (2*x)
  else if (n==2) then
     hermite = (-2 + 4*x**2)
  else if (n==3) then
     hermite = (-12*x + 8*x**3)
  else if (n==4) then
     hermite = (12 - 48*x**2 + 16*x**4)
  else if (n==5) then
     hermite = 120*x - 160*x**3 + 32*x**5
  else if (n==6) then
     hermite = (-120 + 720*x**2 - 480*x**4 + 64*x**6)
  else if (n==7) then
     hermite = (-1680*x + 3360*x**3 - 1344*x**5 + 128*x**7)
  else if (n==8) then
     hermite = (1680 - 13440*x**2 + 13440*x**4 - 3584*x**6 + 256*x**8)
  else if (n==9) then
     hermite = (30240*x - 80640*x**3 + 48384*x**5 - 9216*x**7 + 512*x**9)
  else if (n==10) then
     hermite = (-30240 + 302400*x**2 - 403200*x**4 + 161280*x**6 - 23040*x**8 + 1024*x**10)
  end if
end function hermite


subroutine read_wfn_init(filenamelength, filename, nmonpnn)

  implicit none

  integer*4, intent(inout) :: filenamelength
  character(len=filenamelength), intent(inout) :: filename
  integer*4, intent(inout) :: nmonpnn(3)
  integer*4 :: nmo
  integer*4 :: np
  integer*4 :: nn

  character*80 :: comment
  integer :: status = 0

  open(unit=5000, file=filename)
  read(unit=5000, fmt=5001) comment
  read(unit=5000, fmt=5002) nmo, np, nn
  close(unit=5000)

  nmonpnn(1) = nmo
  nmonpnn(2) = np
  nmonpnn(3) = nn

5001 format(a80)
5002 format(8x,10x,i5,13x,1x,i6,11x,4x,i5,7x)

end subroutine read_wfn_init


subroutine read_wfn(filename_length,filename,nmo,np,nn,xn,zn,centre,type,a,o,oe,c,ev)

  implicit none

  integer*4, intent(inout) :: filename_length
  character(len=filename_length), intent(inout) :: filename
  integer*4, intent(inout) :: nmo
  integer*4, intent(inout) :: np
  integer*4, intent(inout) :: nn
  real*8, intent(inout) :: xn(nn,3)
  real*8, intent(inout) :: zn(nn)
  integer*4, intent(inout) :: centre(np)
  integer*4, intent(inout) :: type(np)
  real*8, intent(inout) :: a(np)
  real*8, intent(inout) :: o(nmo)
  real*8, intent(inout) :: oe(nmo)
  real*8, intent(inout) :: c(np,nmo)
  real*8, intent(inout) :: ev(2)

  character*80 :: comment
  integer :: status = 0

  integer*4 :: iskip
  real*8 :: rskip

  real*8 :: e
  real*8 :: v

  integer*4 :: n,i,j

  open(unit=5000, file=filename)
  read(unit=5000, fmt=5001) comment
  read(unit=5000, fmt=5002) iskip, iskip, iskip

  do n=1,nn
     read(unit=5000, fmt=5003) xn(n,1), xn(n,2), xn(n,3), zn(n)
  end do

  read(unit=5000, fmt=5004) centre
  read(unit=5000, fmt=5005) type
  read(unit=5000, fmt=5006) a

  do i=1,nmo
     read(unit=5000,fmt=5007) o(i), oe(i)
     read(unit=5000,fmt=5008) (c(j,i), j=1,np)
  end do
  read(unit=5000,fmt=5009)
  read(unit=5000,fmt=5010) e, v

  close(unit=5000)

  ev(1) = e
  ev(2) = v

5001 format(a80)
5002 format(8x,10x,i5,13x,1x,i6,11x,4x,i5,7x)
5003 format(24x,3f12.8,10x,f5.1)
5004 format(18x,2x,20i3)
5005 format(16x,4x,20i3)
5006 format(9x,1x,1p,5e14.7)
5007 format(36x,f12.8,14x,f13.8 )
5008 format(1p,5e16.8)
5009 format(8x)
5010 format(8x,9x,f20.10,18x,f13.8)

end subroutine read_wfn

