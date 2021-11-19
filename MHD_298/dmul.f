**********************************************************

	subroutine dmul(fbx,fby,fbz,nx,ny)

* fbx=-dyw1 -dzw2
* fby= dxw1 -dzw3
* fbz= dxw2 +dyw3
* => dxfbx+dyfby+dzfbz = 0

	implicit real*8(a-h,o-z)
	include 'bunix.inc'
	dimension fbx(2,nx/2+1),fby(2,nx/2+1)
	dimension fbz(2,nx/2+1)

	do 3  i=1,nx/2
	 w1r=fbx(1,i)
	 w1i=fbx(2,i)
	 w2r=fby(1,i)
	 w2i=fby(2,i)
	 w3r=fbz(1,i)
	 w3i=fbz(2,i)
	fbx(1,i)=+sinar*rkx(i)*w1i
	fbx(2,i)=-sinar*rkx(i)*w1r
	fby(1,i)=	-cosa*rkx(i)*w1i
	fby(2,i)=	+cosa*rkx(i)*w1r
	fbz(1,i)=	-cosa*rkx(i)*w2i-sinar*rkx(i)*w3i
	fbz(2,i)=	+cosa*rkx(i)*w2r+sinar*rkx(i)*w3r

3     continue
	end
