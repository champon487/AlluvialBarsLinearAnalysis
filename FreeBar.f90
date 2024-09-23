! ------------------------------------------------ !
! Temporal mode linear analysis of alluvial bars   !
! ------------------------------------------------ !

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

subroutine LinearStabilityAnalysis(gro,cel,lam,beta,m,cf,fr,phi10,gam0)
    implicit none
    integer, intent(in) :: m
    double precision, intent(in) :: lam, beta, cf, fr, phi10, gam0
    double precision, intent(out) :: gro, cel

    double precision, parameter :: pi=3.14159265358979d0
	complex :: a11,a12,a13,a14
	complex :: a21,a22,a23,a24
	complex :: a31,a32,a33,a34
	complex :: a41,a42,a43,a441,a442
	complex,parameter :: iu=(0.,1.)

    a11 = iu*lam
    a12 = (-1)**m*pi*0.5d0*dble(m)
    a13 = iu*lam
    
    a21 = 2.d0*cf*fr**2.d0*beta+iu*lam*fr**2.d0
    a23 = -4.d0/3.d0*cf*fr**2.d0*beta+iu*lam
    a24 = iu*lam
    
    a32 = cf*fr**2.d0*beta+iu*lam*fr**2.d0
    a33 = (-1)**(m+1)*pi*0.5d0*dble(m)
    a34 = (-1)**(m+1)*pi*0.5d0*dble(m)
    
    a41 = phi10*lam*iu
    a42 = (-1)**m*pi*0.5d0*dble(m)
    a43 = -phi10/6.d0*lam*iu
    a441 = gam0*(dble(m)*pi*0.5d0)**2.d0/beta

    a442 = ( -a11*a23*(a32*a441-a34*a42)+a11*a24*(a32*a43-a33*a42)	&
                -a12*a21*(a33*a441-a34*a43)-a12*a23*a34*a41+a12*a24*a33*a41	&
                +a13*a21*(a32*a441-a34*a42)-a13*a24*a32*a41 )				&
                        /(a11*a23*a32+a12*a21*a33-a13*a21*a32)

    gro = real(a442)
    cel = -imag(a442)

end subroutine 

program bar_instability

	implicit none
	integer :: m, ii, nb, nl
	double precision :: lam, pi, g, snu0, l, n
	double precision :: ib, dis, diam, spec, b0
	double precision :: h0, nm, u0
	double precision :: fr, cf, ts0, tsc, gam0, phi0, beta, beta0, mu_s
	double precision :: phi, phi10
    double precision :: gro, cel
	
	g = 9.81d0
	snu0 = 0.000001d0
	
	! --- read conditions --- !

	open(1,file="cond.txt",status="old")
	
		read(1,*) ib		! slope of river bed
		read(1,*) dis		! discharge (m3/s)
		read(1,*) diam		! sediment diameter (m)
		read(1,*) spec		! specific weight of sediment in fluid
		read(1,*) b0		! half of channel width (m)
		read(1,*) mu_s		! friction coefs. of sediment
		read(1,*) m 		! mode of sand bar
	
	close(1)
	
	! --- calculate parameters --- !
	
	nm = (2.5d0*diam)**(1.d0/6.d0)/(7.66d0*dsqrt(g))	! Manning's n by Manning - strickler
	h0 = (nm*dis/(2.d0*b0*ib**0.5d0))**0.6d0			! Uniform flow depth
	u0 = dis/(2d0*b0*h0)								! uniform flow velocity
	beta0 = b0/h0										! aspect ratio (half channel width - water depth ratio)
	fr = u0/dsqrt(g*h0)					! Froude number
	cf = g*nm**2.d0/h0**0.333333d0		! friction coefficient
	ts0 = cf*u0**2.d0/(spec*g*diam)		! Shields number
	call usc(diam,tsc,spec,g,snu0)		! Critical Shields number
	phi0 = ts0-tsc
	gam0 = dsqrt(tsc/ts0)/mu_s	        ! coefs. of bed slope effect in n direction (Hasegawa)
	
	phi = tsc/ts0
	phi10 = (3.d0+phi**0.5d0)/(1.d0-phi)        !  some coefs of Ashida - Michiue bedload formula
	
	
	write(*,*) "h0= ", h0
	write(*,*) "u0= ", u0
	write(*,*) "Fr= ", fr
	write(*,*) "beta= ", beta
	write(*,*) "T*0= ", ts0
	write(*,*) "T*c= ", tsc
	write(*,*) "gamma= ", gam0

    !  U-shape instability diagram

	open(12,file='contor.dat',status='unknown')
        
        nb = 400
        nl = 400

        write(12,*) beta0
        
        write(12,*) nl, nb

        do n = 1, nb
        
            beta = dble(n)*0.25d0   ! aspect ratio
                            
            do l = 1, nl
            
                lam = l*0.005d0		! dimensionless wavenumber

                call LinearStabilityAnalysis(gro,cel,lam,beta,m,cf,fr,phi10,gam0)
                        
                !  non-dimensional wavenumber, aspect ratio, growrth rate, celerity
                            
                write(12,'(4e16.6)') lam, beta, gro, cel

            end do
            
        end do

	close(12)

    !  grwoth rate of given aspect ratio

    open(13,file='Omega-lam.dat',status='unknown')

        write(13,*) beta0
        
        write(13,*) nl

        do l = 1, nl
        
            lam = l*0.005d0		! dimensionless wavenumber

            call LinearStabilityAnalysis(gro,cel,lam,beta0,m,cf,fr,phi10,gam0)
                    
            !  non-dimensional wavenumber, aspect ratio, growrth rate, celerity
                        
            write(13,'(3e16.6)') lam, gro, cel

        end do
        
    close(13)

end program bar_instability

