***************************
	subroutine expansio

	implicit real*8(a-h,o-z)
	include 'bunix.inc'

	rayon=uu0*temp+1.
	rainv=1./rayon
	rainv2=rainv*rainv
	bx0=bx00*rainv
	fnu=uu0 *rainv
	sinar=sina*rainv
	ck2=cosa**2+sina**2*rainv2

	if (ifiso.eq.0)  then
	   area=rayon**2
	   sinard=sinar
	   ck2d=ck2
	else
	   area=1.
	   sinard=sina
	   ck2d=1.
	end if
	return
	end
