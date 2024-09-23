! ------------------------------------------------ !
! Spatial mode linear analysis of alluvial bars   !
! ------------------------------------------------ !

subroutine ccuberoot(x,y)
	implicit none
	double precision, intent(in) :: y
	double precision, intent(out) :: x

	if( y>0.d0 ) then
		x = (y)**(1.d0/3.d0)
	else
		x = -(-y)**(1.d0/3.d0)
	end if
end subroutine ccuberoot

subroutine solvecubiceq(x,a,b,c,d)
	implicit none
	double precision, intent(in) :: a, b, c, d
	complex(8), dimension(0:2), intent(out) :: x

	double precision :: aa, bb, cc, p, q, dd, theta, u, v
	double precision, parameter :: pi = 3.141592d0
	complex(8), parameter :: iu=(0.d0,1.d0)

	aa = b/a
	bb = c/a
	cc = d/a
	p = bb-aa*aa/3.d0
	q = 2.d0*aa*aa*aa/27.d0-aa*bb/3.d0+cc
	dd = q*q/4.d0+p*p*p/27.d0

	if( dd<0.d0 ) then
		theta = datan2(sqrt(-dd),-q*0.5d0)
		x(0) = 2.d0*sqrt(-p/3.d0)*cos(theta/3.d0)-aa/3.d0
        x(1) = 2.d0*sqrt(-p/3.d0)*cos((theta+2.d0*pi)/3.d0)-aa/3.d0
        x(2) = 2.d0*sqrt(-p/3.d0)*cos((theta+4.d0*pi)/3.d0)-aa/3.d0
    else
        call ccuberoot(u,-q*0.5d0+sqrt(dd))
		call ccuberoot(v,-q*0.5d0-sqrt(dd))
        x(0) = u+v-aa/3.d0
        x(1) = -0.5d0*(u+v)+sqrt(3.d0)*0.5d0*iu*(u-v)-aa/3.d0
        x(2) = -0.5d0*(u+v)-sqrt(3.d0)*0.5d0*iu*(u-v)-aa/3.d0
	end if

end subroutine solvecubiceq

subroutine solvequarticeq(x,a,b,c,d,e)
	implicit none
	double precision, intent(in) :: a, b, c, d, e
	complex(8), dimension(0:3), intent(out) :: x

	double precision :: aa, bb, cc, dd, p, q, r, t
	double precision, parameter :: epsilon = 0.000001d0
	complex(8) :: m, w, w1, w2
	complex(8), dimension(0:2) :: xx

	aa = b/a
    bb = c/a
    cc = d/a
    dd = e/a
    p = -6.d0*(aa/4.d0)**2.d0+bb
    q = 8.d0*(aa/4.d0)**3.d0-2.d0*bb*aa/4.d0+cc
    r = -3.d0*(aa/4.d0)**4.d0+bb*(aa/4.d0)**2.d0-cc*aa/4.d0+dd

    call solvecubiceq(xx,1.d0,-p,-4.d0*r,4.d0*p*r-q*q)
	
    t = real(xx(0))
	w = t-p
    m = cdsqrt(w)
	w1 = -t-p+2.d0*q/m
	w2 = -t-p-2.d0*q/m
    x(0) = (-m+cdsqrt(w1))*0.5d0-aa/4.d0
    x(1) = (-m-cdsqrt(w1))*0.5d0-aa/4.d0
    x(2) = (m+cdsqrt(w2))*0.5d0-aa/4.d0
    x(3) = (m-cdsqrt(w2))*0.5d0-aa/4.d0


end subroutine solvequarticeq


	! --- subroutine for cal. critical Shields number by Iwagaki's formula --- !

subroutine usc(ddd,tsci,spec,g,snu00)
	implicit none
	double precision,intent(in) :: ddd, spec, g, snu00
	double precision,intent(out) :: tsci
	double precision :: rst, usci

	rst = sqrt(spec*g)/snu00*ddd**1.5

	if( rst<=2.14 )	usci = 0.14d0*spec*g*ddd
	if( rst>2.14  )	usci = (0.1235d0*spec*g)**(0.78125)*snu00**(0.4375)*ddd**(0.34375)
	if( rst>54.2  )	usci = 0.034d0*spec*g*ddd
	if( rst>162.7 )	usci = (0.01505d0*spec*g)**(1.136364)*snu00**(-0.272727273)*ddd**(1.40909091)
	if( rst>671   )	usci = 0.05d0*spec*g*ddd

	usci = sqrt(usci)
	tsci = usci**2.d0/(spec*g*ddd)

end subroutine usc

program bar_instability

	implicit none
	integer :: ii, nb, nl, i, j
	double precision :: lam, pi, g, snu0, l, n
	double precision :: ib, dis, diam, spec, b0
	double precision :: h0, nm, u0
	double precision :: fr, cf, ts0, tsc, gam0, phi0, beta, mu_s
	double precision :: phi, phi10
    double precision :: a, b, c, d, e
	double precision :: a11,a12,a13,a14
	double precision :: a21,a22,a23,a24
	double precision :: a31,a32,a33,a34
	double precision :: a41,a42,a43,a44,a442
	complex(8) :: temp
    complex(8), dimension(0:3) :: lll
	complex(8),parameter :: iu=(0.,1.)
	
	pi = 3.14159265358979d0
	g = 9.81d0
	snu0 = 0.000001d0
	
	! --- read conditions --- !

	open(1,file="cond.txt",status="old")
	
		read(1,*) ib		! slope of river bed
		read(1,*) dis		! discharge (m3/s)
		read(1,*) diam		! sediment diameter (m)
		read(1,*) spec		! specific weight of sediment in fluid
		read(1,*) b0		! half of channel width (m)
		read(1,*) mu_s		! kinematic friction coefs. of sediment
	
	close(1)
	
	! --- calculate parameters --- !
	
	nm = (2.5d0*diam)**(1.d0/6.d0)/(7.66d0*dsqrt(g))	! Manning - strickler
	h0 = (nm*dis/(2.d0*b0*ib**0.5d0))**0.6d0			! Uniform flow depth
	u0 = dis/(2d0*b0*h0)								! uniform flow velocity
	beta = b0/h0										! aspect ratio
	fr = u0/dsqrt(g*h0)					! Froude number
	cf = g*nm**2.d0/h0**0.333333d0		! friction coefficient
	ts0 = cf*u0**2.d0/(spec*g*diam)		! Shields number
	call usc(diam,tsc,spec,g,snu0)		! Critical Shields number
	phi0 = ts0-tsc
	gam0 = dsqrt(tsc/ts0)/mu_s	! coefs. of bed slope effect in n direction (Hasegawa)
	
	phi = tsc/ts0
	phi10 = (3.d0+phi**0.5d0)/(1.d0-phi)
	
	
	write(*,*) "h0= ", h0
	write(*,*) "u0= ", u0
	write(*,*) "Fr= ", fr
	write(*,*) "beta= ", beta
	write(*,*) "T*0= ", ts0
	write(*,*) "T*c= ", tsc
	write(*,*) "gamma= ", gam0

	open(12,file='ForcedLam.dat',status='unknown')
	
	nb = 1600
	nl = 800
	
	do n = 1, nb
	
		beta = dble(n)*0.125d0*0.25d0
				
        a12 = -0.5d0*pi
        a21 = 2.d0*cf*beta
        a23 = -4.d0/3.d0*cf*beta
        a32 = cf*beta
        a33 = 0.5d0*pi/fr**2.d0
        a34 = 0.5d0*pi/fr**2.d0
        a41 = phi10
        a42 = -0.5d0*pi
        a43 = -phi10/6.d0
        a44 = gam0*(0.5d0*pi)**2.d0/beta

        a = (a43-a41)/fr**2.d0
        b = a44-(a44-a32*a43+a32*a41)/fr**2.d0
        c = a32*a44-a23*a44+a21*a44+a12*a34*a43-a34*a42+(a34*a42-a32*a44-a33*a42-a12*a34*a41+a12*a33*a41)/fr**2.d0
        d = a21*a32*a44-a12*a33*a44-a23*a32*a44+a12*a21*a34*a43+a23*a34*a42-a21*a34*a42-a12*a23*a34*a41
        e = -a12*a21*a33*a44

        call solvequarticeq(lll,a,b,c,d,e)

		do i=0,2
			do j=0,2-i
				if( real(lll(j))<real(lll(j+1)) ) then
					temp = lll(j)
					lll(j) = lll(j+1)
					lll(j+1) = temp
				end if
			end do				
		end do

        write(12,'(9e16.6)') beta, real(lll(0)), imag(lll(0)), real(lll(1)), imag(lll(1)), real(lll(2)), imag(lll(2)), real(lll(3)), imag(lll(3))
				
	end do

	close(12)

end program bar_instability

