
***********		 calculs	   **************

	 subroutine avanti0(uun,fnl,vnl,cin,nx,ny,ifadams,romin,umax)


	implicit real*8(a-h,o-z)

*	dissipation completement explicite dans tnl
*	(la partie diagonale -laplacien- peut etre implicite ici, cf fac)
*	on avance en complexe
*	expansion anisotrope:
*	adams-bashforth: pour bx on traite l'expansion de facon exacte:
* dbx/dt = f - (dl/dt)/l bx ou l(t)=rayon=1+uu0.t, devient:
* bxnew= (3/2)dt(l/lnew)f-(1/2)dt(lold/lnew)fold + (l/lnew)bx
* avec lold=l(t-dtold), lnew=l(t+dt), de meme pour vy et vz
	include 'bunix.inc'
	complex *16  uun (nx/2+1,8)
	complex *16  fnl(nx/2+1,8)
	complex *16  vnl(nx/2+1,8)
	complex *16  cin (nx/2+1,8)
	dimension aamu(8),fmu(8)

	do 1000 i=1,8
		 aamu(i)=amu(i)
		 fmu(i)=0.
1000    continue
	aamu(2)=0.
	aamu(3)=0.
	aamu(4)=0.
	aamu(8)=0.

	dth=dt/2.
	if (dt.eq.dtold) then
	  trois=3.*dth
	  un=dth
	else
	  trois=(1.+0.5*dt/dtold)*dt
	  un=trois-dt
	end if
	k1=1
	k2=8
	if (ifns2D.eq.1) k2=3
	if (ifbur.eq.1) k1=2

	if (ifadams.eq.1) then

	rnew=1.+uu0*(temp+dt   )
	rold=1.+uu0*(temp-dtold)

*	aamu viscosite (partie diagonale)

	do 1 k=k1,k2
	if (k.eq.5.or.k.eq.3.or.k.eq.4) then
	 f   =rayon/rnew
	 fold=rold/rnew
	else
	 f   =1.
	 fold=1.
	end if
	do 1 i=1,nx/2
	 aa=dth*aamu(k)*ck2*rkx2(i)
	 uun(i,k) =
     &   ((trois*f*fnl(i,k)-un*fold*vnl(i,k))
     &     + f*uun(i,k) *(1-aa)) /(1+aa)
1	continue

	do 2 k=k1,k2
	do 2 i=1,nx/2
2	 vnl(i,k) = fnl(i,k)

	else

*	euler 2eme ordre

	if(leul.eq.2) then
	 dtt=dt/2.
	 f=1./(1.+uu0*dt/2.)
	 fci=f
	else
	 dtt=dt
	 fci=1./(1.+uu0*dt)
	 f=fci*(1.+uu0*dt/2.)
	end if

	if (leul.eq.2) then
	 do 3 k=k1,k2
	 do 3 i=1,nx/2
3	 vnl(i,k) = fnl(i,k)
	end if

	do 4 k=k1,k2

	 if (k.eq.5.or.k.eq.3.or.k.eq.4) then
	  ff=f
	  ffci=fci
	 else
	  ff=1.
	  ffci=1.
	 end if
	do 4 i=1,nx/2
	 uun(i,k) = ( dtt * ff*fnl(i,k) + ffci* cin(i,k))
     & /(1.+dtt*aamu(k)*ck2*rkx2(i))
4	continue

	end if

	dtold=dt
	return
	end
