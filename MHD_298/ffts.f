
**********************************************************

		 subroutine ffts(nx,nf,f,isg)

*	aller (isg=1):
*	fft(f)/n->f
*	retour (isg=-1):
*	invfft(f)->f

	implicit real*8(a-h,o-z)
	include 'bunix.inc'
	dimension f(nx+2,nf)

	isgi=-isg
	if (isg.lt.0) go to 20

* a:			  isg=1: reel -> complexe

	 nft=nf
	 inc=1
	 jump=nx+2

*	  call rfftmlt(f,wf,ex,ifax,inc,jump,nx,nft,isgi)
	 call FFT991(f,wf,ex,ifax,inc,jump,nx,nft,isgi)
	return

* b:		isg=-1: c -> reel

20	continue

	 nft=nf
	 inc=1
	 jump=nx+2
*	call rfftmlt(f,wf,ex,ifax,inc,jump,nx,nft,isgi)
	call FFT991(f,wf,ex,ifax,inc,jump,nx,nft,isgi)
	return
	end
