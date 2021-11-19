************************************************

		 subroutine dissi(vxn,vyn,vzn,d2v,ten,d2te,nx,ny)

*	calcule la partie explicite de la dissipation (tnl):
*	snu*(laplacien(u) + untrs*graddivu) (untrs=1/3 ou 0)
*	djsij=-snu ( k2 ui +1/3 kikjuj)

	implicit real*8(a-h,o-z)
	include 'bunix.inc'
	parameter (untrs=1./3.,qtr=4./3.,str=7./3.)
	dimension vxn(2,nx/2+1),vyn(2,nx/2+1)
	dimension vzn(2,nx/2+1),d2v(2,nx/2+1,3)
	dimension ten(2,nx/2+1),d2te(2,nx/2+1)

	do 1 l=1,2
	do 1 i=1,nx/2
	 ak2=rkx2(i)*ck2d
	 divu=rkx(i)*(cosa*vxn(l,i)+sinard*vyn(l,i))
	 d2v(l,i,1)=-snu*(ak2*vxn(l,i)+untrs*cosa*rkx(i)*divu)
	 d2v(l,i,2)=-snu*(ak2*vyn(l,i)+untrs*sinard*rkx(i)*divu)
	 d2v(l,i,3)=-snu*(ak2*vzn(l,i))
	 d2te(l,i)=-sk*ak2*ten(l,i)
1	continue

	call ffts  (nx,3,d2v,-1)
	call ffts  (nx,1,d2te,-1)

	return
	end
