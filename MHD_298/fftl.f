
****************************
		 subroutine fftl(nx,nf,f,ft,isg)

	implicit real*8(a-h,o-z)
	include 'bunix.inc'
	dimension f(nx+2,nf),ft(nx+2,nf)


*	f = reel,	ft=complexe
*	aller (isg=1):
*	fft(f)/n->ft
*	retour (isg=-1):
*	invfft(ft)->f

	isgi=-isg
	if (isg.lt.0) go to 20

* a:			  isg=1: reel -> complexe

	 nft=nf
	 inc=1
	 jump=nx+2

	do 1 l=1,nf
	 do 1 i=1,nx
 1		 ft(i,l)=f(i,l)
*	  call rfftmlt(ft,wf,ex,ifax,inc,jump,nx,nft,isgi)
	 call FFT991(ft,wf,ex,ifax,inc,jump,nx,nft,isgi)
	return

* b:		isg=-1: c -> reel

20	continue

	 nft=nf
	 inc=1
	 jump=nx+2
	 do 5  l=1,nf
	  do 5  i=1,nx
5		f(i,l)=ft(i,l)
*		call rfftmlt(f,wf,ex,ifax,inc,jump,nx,nft,isgi)
	call FFT991(f,wf,ex,ifax,inc,jump,nx,nft,isgi)
	return
	end
