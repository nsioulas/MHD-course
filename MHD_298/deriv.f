**********************************************************

		 subroutine deriv(uun,du,nx,ny)

*  calcule du/dx(1d)

	implicit real*8(a-h,o-z)
	include 'bunix.inc'
	dimension uun(2,nx/2+1,8),du(2,nx/2+1,8)

	do 1 k=1,8
	do 1 i=1,nx/2+1
		 du(1,i,k)=-rkx(i)*uun(2,i,k)
		 du(2,i,k)= rkx(i)*uun(1,i,k)
1	continue
	call ffts(nx,8,du,-1)
	return
	end
