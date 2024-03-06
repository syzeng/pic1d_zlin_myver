    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !***********************************************************************************
    !* Modules for program PIC1D, Unmagnetized plasma use electron unit: n,m,e,T=1     *
    !* Basic unit: Debye length sqrt(T/ne^2)=1, plasma frequency sqrt(ne^2/m)=1        *
    !* compile on dragon: "pgf90 -mp pic1d.f90"                 
    !* compile on linux sequentially: "ifort pic1d.f90"        *
    !* compile on linux in parallel:  "ifort -qopenmp pic1d.f90"        *
    !***********************************************************************************
    Module parameter

    implicit none

    ! numerics
    integer,parameter ::  singlep=selected_real_kind(6),doublep=selected_real_kind(12)
    real,parameter :: pi=3.14159265358979,small=1.0e-10

    ! control parameter
    integer :: stdout=1,imarker=1 !stdout=0: output to display; 1: to file "PIC1D.out"
    !1: load physical PDF; 0: uniform for [-vmax,vmax]
    real :: vmax=5.0           !velocity space range [-vmax,vmax], f_m(v=5)=0.000003
    real :: linear=0.0         !1.0 for linear simulation; 0.0 for nonlinear.
    real :: deltaf=1.0         !1.0 for delta_f simulation, 0.0 for full_f

    ! particle array: # of species, particles, velocity grids. First specie is electron
    integer,parameter :: nspecie=2,nparticle=60000,nvelocity=128
    integer,parameter :: coll=0!coll > 0 turns on collision operator
    real,parameter :: nu=0.0   !collision frequency normalized by gyrofrequency
    real,dimension(nspecie) ::  qspecie,aspecie,tspecie
    data qspecie/-1.0,1./         !charge
    data aspecie/1.0,100./          !mass
    data tspecie/1.0,1./          !temperature
    real,dimension(nparticle,nspecie) :: x,xa,v,va,p0,mu,w,wa
    ! x: position; v: velocity; w=delta_f/g; p0=f/g; f=f_0+delta_f, g is marker PDF,
    ! mu: magnetic moment

    ! field array: # of grid points, must be power of 2 due to FFT
    integer,parameter :: ngrid=64
    real :: deltax=0.245437
    real,dimension(0:ngrid) :: charge,phi,efield
    ! charge density, potential, electric field

    ! time array: # of time steps, interval of steps, modes for diagnostics
    integer,parameter :: ntime=500,ndiag=2,nmode=2
    integer,parameter :: nd=(ntime-1)/ndiag+1
    integer mode(nmode),nt
    data mode/1,2/             ! mode # for diagnostics
    real :: tstep=0.1          ! time step size
    real,dimension(nd) :: field_ene,potential ! field energy, potential RMS
    real,dimension(nd,2,nmode) :: phi_mode ! mode amplitudes
    real,dimension(nd,nspecie) :: density,entropy,flow,kin_ene,flow_mar,kin_mar
    ! volume average density, entropy, flow, kinetic energy, marker flow and energy

    end Module parameter

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !***********************************************************************************
    !* PIC1D main program for 1-D ES particle simulation of Vlasov plasmas             *
    !*                        Zhihong Lin @ 2006                                       *
    !*                 University of California, Irvine                                *
    !* Contributors: Hongpeng Qu, Onnie Luk, Joseph McClenaghan                        *
    !***********************************************************************************
    Program PIC1D

    use parameter
    implicit none
    integer nrk
    real ctime0,ctime1

    ! obtain starting CPU time
    call CPU_TIME(ctime0)
    if(stdout==1)open(stdout,file='PIC1D.out',status='replace')

    ! load particles
    call load

    ! 2nd Runge-Kutta time integration
    do nt=1,ntime
        do nrk=1,2

            ! calculate field quatities, time diagnostics
            call field(nrk)

            ! integrate particle orbits
            call particle(nrk)

        enddo

        if(coll>0)call collision
    enddo

    ! write history data
    call history

    ! obtain ending CPU time
    call CPU_TIME(ctime1)
    write(stdout,*)'Program CPU TIME=', ctime1-ctime0

    if(stdout==1)close(stdout)

    end Program PIC1D
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Subroutine load

    use parameter
    implicit none
    integer i,j,k,isize,init,nseed
    integer,dimension(:),allocatable :: iget,iput
    real :: c0=2.515517,c1=0.802853,c2=0.010328,d1=1.432788,d2=0.189269,d3=0.001308
    real winit,xsize,xfreq,tfreq,vth_inv
    real tmp

    ! initial weight
    winit=0.1
    ! first harmonics k=2pi/L, L is system size
    xsize=real(ngrid)*deltax
    xfreq=2.0*pi/xsize

    ! initialize Fortran random number generator
    call random_seed(size=isize)
    allocate(iget(isize),iput(isize))
    call random_seed(get=iget)
    do i=1,isize
        !     call system_clock(init)   !random initialization
        init=1111                  !same initialization
        iput(i)=1111*init+iget(i)
    enddo
    call random_seed(put=iput)

    ! load particles
    do k=1,nspecie

        call random_number(v(:,k))
        if(imarker==1)then ! Maxwellian distribution

            ! marker weight defined as p == f/g, f is physical PDF, g is marker PDF
            p0(:,k)=1.0

            ! transform uniform random number to Gaussian with variance=1
            v(:,k)=v(:,k)-0.5
            va(:,k)=sign(1.0,v(:,k))
            v(:,k)=sqrt(max(small,log(1.0/max(small,v(:,k)*v(:,k)))))
            v(:,k)=v(:,k)-(c0+c1*v(:,k)+c2*v(:,k)*v(:,k))/&
                (1.0+d1*v(:,k)+d2*v(:,k)*v(:,k)+d3*v(:,k)*v(:,k)*v(:,k))
            v(:,k)=va(:,k)*v(:,k)

        elseif(imarker==0)then ! uniform loading for [-vmax,vmax]
            ! marker weight p=exp(-v^2/2)
            v(:,k)=2.0*vmax*(v(:,k)-0.5)
            p0(:,k)=exp(-0.5*v(:,k)**2)

        endif

        ! normalize velocity to thermal velocity
        v(:,k)=v(:,k)*sqrt(tspecie(k)/aspecie(k))

        ! normalize total marker # to nparticle
        tmp=sum(p0(:,k))
        p0(:,k)=p0(:,k)*real(nparticle)/tmp

        ! equal-spacing loading
        do i=1, nparticle
            x(i,k)=xsize*real(i)/real(nparticle)
            !        x(i,k)=x(i,k)-winit*sin(xfreq*x(i,k))/xfreq
            ! periodic BC: out of bound particle recycled
            !        x(i,k)=x(i,k)/(real(ngrid)*deltax)+10.0
            !        x(i,k)=real(ngrid)*deltax*(x(i,k)-aint(x(i,k)))
        enddo

        if(deltaf>0.5)then !delta_f simulation
            ! weight w defined as w=delta_f/g, g is the marker distribution function
            ! initial random weight, white noise.
            !     call random_number(w(:,k))
            !     w(:,k)=winit*p0(:,k)*(w(:,k)-0.5)
            ! single mode standing wave
            do i=1,nparticle
                w(i,k)=winit*p0(i,k)*sin(xfreq*x(i,k))
            enddo

            tmp=sum(w(:,k))/real(nparticle)
            w(:,k)=w(:,k)-tmp*p0(:,k)

            ! eigenfunction of undamped Case-Van Kampen mode
            !        tfreq=1.282
            !        w(:,k)=w(:,k)+winit*p0(:,k)*sin(xfreq*x(:,k))*&
            !             max(-100.0,min(100.0,xfreq*v(:,k)/(tfreq-xfreq*v(:,k))))

            ! marker weight corrected by initial perturbation, p = f(t=0)/g(t=0)
            p0(:,k)=p0(:,k)+(1.0-linear)*w(:,k)

            tmp=sum(p0(:,k)*v(:,k))/real(nparticle)
            v(:,k)=v(:,k)-tmp

            tmp=sqrt(sum(p0(:,k)*v(:,k)*v(:,k))/real(nparticle))
            v(:,k)=v(:,k)/tmp

        elseif(deltaf<0.5)then ! full_f simulation
            w(:,k)=p0(:,k)
        endif

        call exponential(mu(:,k),nparticle,nseed+2)

    enddo

    end Subroutine load
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Subroutine field(nrk)

    use parameter
    implicit none
    integer,intent(in) :: nrk
    integer,parameter :: nfield=11
    integer i,j,jx,jv,k,ndt
    real dx_inv,dv_inv,delj,xsize,x_inv,xpart,xfreq,vtmp,&
        vth_inv,xvdelt(nvelocity,ngrid),xvf0(nvelocity,ngrid), tempcharge(0:ngrid)
    real entropytemp,densitytemp,flowtemp,kin_enetemp,flow_martemp,kin_martemp
    complex,dimension(ngrid/2+1) :: filter,charge_k,phi_k,field_k

    save xvf0

    ! book keeping
    dv_inv=real(nvelocity)/(2.0*vmax)
    dx_inv=1.0/deltax
    xsize=real(ngrid)*deltax
    x_inv=1.0/xsize
    xfreq=2.0*pi*x_inv

    ! mode filtering: k=0,1,...,ngrid/2
    filter=0.0
    !  filter(2:ngrid/8)=1.0
    filter(2)=1.0		! filter out all modes except k=2*pi/xsize

    charge=0.0
    do k=1,nspecie
        ! loop each thread has its own private array "tempcharge()" which is used to
        ! prevent the thread from writing into the shared array "charge()". We use tempcharge
        ! to avoid contentions between threads.
        !$omp parallel private(tempcharge)
        tempcharge = 0.0
        !$omp do private(i,j,delj,xpart)
        do i=1,nparticle

            ! periodic BC
            xpart=x(i,k)*x_inv+10.0
            xpart=xsize*(xpart-aint(xpart))

            ! particle to grid charge scattering
            j=int(xpart*dx_inv)
            delj=xpart*dx_inv-j
            j=max(1,min(ngrid,j+1))
            tempcharge(j)=tempcharge(j)+delj*w(i,k)*qspecie(k)
            tempcharge(j-1)=tempcharge(j-1)+(1.0-delj)*w(i,k)*qspecie(k)
        enddo
        !$omp end do
        !$omp critical
        ! add all private tempcharges to shared charge variable. The loop is enclosed
        ! in a critical section so that one thread at a time updates charge().
        do i=0, ngrid
            charge(i) = charge(i) + tempcharge(i)
        enddo
        !$omp end critical
        !$omp end parallel
    enddo
    charge(ngrid)=charge(ngrid)+charge(0)     !boundary grid
    charge(0)=charge(ngrid)                     !periodic BC
    charge=charge*real(ngrid)/real(nparticle)

    ! Fourier transform of charge density
    call r2cfft(ngrid,charge(0:ngrid-1),charge_k)

    ! Poisson equation for potential, grad^2_phi=-charge, k*k*phi_k=charge_k
    phi_k(1)=0.0
    do k=1,ngrid/2
        phi_k(k+1)=filter(k+1)*charge_k(k+1)/(xfreq*real(k))**2
    enddo
    ! transform back to configuration space
    call c2rfft(ngrid,phi_k,phi(0:ngrid-1))
    phi(ngrid)=phi(0)                         !periodic BC

    ! filtered electric field, grad E=charge, -ik*efield_k=charge_k*filter
    field_k(1)=0.0
    do k=1,ngrid/2
        field_k(k+1)=(0.0,1.0)*filter(k+1)*charge_k(k+1)/(xfreq*real(k))
    enddo
    ! transform back to configuration space
    call c2rfft(ngrid,field_k,efield(0:ngrid-1))
    efield(ngrid)=efield(0)                   !periodic BC

    if(nrk==1 .and. mod(nt-1,ndiag)==0)then

        ndt=(nt-1)/ndiag+1
        ! field diagnostic: field energy and potential RMS, mode amplitudes
        field_ene(ndt)=0.5*sum(phi(1:ngrid)*charge(1:ngrid))/real(ngrid)
        !     field_ene(ndt)=0.5*sum(efield(1:ngrid)*efield(1:ngrid))/real(ngrid)
        potential(ndt)=sqrt(sum(phi(1:ngrid)**2)/real(ngrid))
        do j=1,nmode
            phi_mode(ndt,1,mode(j))=real(phi_k(mode(j)+1))
            phi_mode(ndt,2,mode(j))=aimag(phi_k(mode(j)+1))
        enddo

        ! particle diagnostic: volume average
        do k=1,nspecie
            ! initializing variables to zero
            density(ndt,k)= 0.0
            entropy(ndt,k) =0.0
            flow(ndt,k) =   0.0
            kin_ene(ndt,k)= 0.0
            flow_mar(ndt,k)=0.0
            kin_mar(ndt,k)= 0.0

            ! loop each thread has its own private array which is used to prevent
            ! the thread from writing into the shared array. We use private variables
            ! to avoid contentions between threads.
            !$omp parallel private(entropytemp,densitytemp,flowtemp,kin_enetemp,&
            !$omp& flow_martemp,kin_martemp)
            ! initializing variables to zero
            densitytemp=0.0
            entropytemp= 0.0
            flowtemp=0.0
            kin_enetemp=0.0
            flow_martemp=0.0
            kin_martemp=0.0

            !$omp do private(i)
            do i=1, nparticle
                densitytemp= densitytemp + w(i,k)
                entropytemp= entropytemp + w(i,k)*w(i,k)
                flowtemp=    flowtemp	 + w(i,k)*v(i,k)
                kin_enetemp= kin_enetemp + w(i,k)*v(i,k)*v(i,k)
                flow_martemp=flow_martemp+ p0(i,k)*v(i,k)
                kin_martemp= kin_martemp + p0(i,k)*v(i,k)*v(i,k)
            enddo
            ! add all private variables to shared variable
            !$omp critical
            density(ndt,k)= density(ndt,k)+densitytemp
            entropy(ndt,k)= entropy(ndt,k)+entropytemp
            flow(ndt,k) =   flow(ndt,k)+flowtemp
            kin_ene(ndt,k)= kin_ene(ndt,k)+kin_enetemp
            flow_mar(ndt,k)=flow_mar(ndt,k)+flow_martemp
            kin_mar(ndt,k)= kin_mar(ndt,k)+kin_martemp
            !$omp end critical

            !$omp end parallel
            density(ndt,k)= density(ndt,k)/real(nparticle)
            entropy(ndt,k)= entropy(ndt,k)/real(nparticle)
            flow(ndt,k) =   flow(ndt,k)/real(nparticle)
            kin_ene(ndt,k)= kin_ene(ndt,k)*0.5*aspecie(k)/real(nparticle)
            flow_mar(ndt,k)=flow_mar(ndt,k)/real(nparticle)
            kin_mar(ndt,k)= kin_mar(ndt,k)*0.5*aspecie(k)/real(nparticle)-0.5
        enddo

        ! total perturbed energy
        field_ene(ndt)=field_ene(ndt)+sum(kin_ene(ndt,:))-(1.0-deltaf)*sum(kin_ene(1,:))
        write(stdout,*)'ntime=',nt,' total energy=',field_ene(ndt)

        if(nt==1)then
            ! Maxwellian f_0, \int(f_0)dv=nparticle
            do jv=1,nvelocity
                vtmp=-vmax+(real(jv)-0.5)*2.0*vmax/real(nvelocity)
                xvf0(jv,:)=exp(-0.5*vtmp*vtmp)
            enddo
            xvf0=xvf0*real(nparticle)/sum(xvf0)

            ! output file of field and particle data for specie=1
            open(nfield,file='movie.out',status='replace')
            write(nfield,101)nvelocity/2,ngrid,(ntime-1)/ndiag+1
            write(nfield,102)2.0*vmax/real(nvelocity-1),deltax,tstep*real(ndiag)
        endif

        k=1
        vth_inv=sqrt(aspecie(k)/tspecie(k))
        xvdelt=0.0
        do i=1,nparticle

            ! velocity space grid
            jv=max(1,min(nvelocity,int((v(i,k)*vth_inv+vmax)*dv_inv+1)))
            ! configuration space grid
            jx=max(1,min(ngrid,int(x(i,k)*dx_inv)+1))

            xvdelt(jv,jx)=xvdelt(jv,jx)+w(i,k)
        enddo

        ! output data is delta_f/f_0, f=f_0+delta_f, \int(f_0)dv=1
        xvdelt=xvdelt/xvf0
        if(deltaf<0.5)xvdelt=xvdelt-1.0
        do j=1,ngrid
            write(nfield,103)phi(j),xvdelt(nvelocity/2:nvelocity-1,j)
        enddo
    endif

101 format(i6)
102 format(e14.5)
103 format(e15.5)

    end Subroutine field
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Subroutine particle(nrk)

    use parameter
    implicit none
    integer,intent(in) :: nrk
    integer i,j,k
    real delt,dx_inv,delj,epart,xsize,x_inv,xpart,temp_inv(nspecie),ratio_qm(nspecie)

    if(nrk==1)then
        ! 1st step of Runge-Kutta method
        delt=0.5*tstep
        do k=1, nspecie
            !$omp parallel do private(i)
            do i=1,nparticle
                xa(i,k)=x(i,k)
                va(i,k)=v(i,k)
                wa(i,k)=w(i,k)
            enddo
        enddo
    else
        ! 2nd step of Runge-Kutta method
        delt=tstep
    endif

    ! book keeping
    dx_inv=1.0/deltax
    ratio_qm=qspecie/aspecie
    temp_inv=1.0/tspecie
    xsize=real(ngrid)*deltax
    x_inv=1.0/xsize

    do k=1,nspecie
        !$omp parallel do private(i,xpart,j,delj,epart)
        do i=1,nparticle

            ! periodic BC
            xpart=x(i,k)*x_inv+10.0
            xpart=xsize*(xpart-aint(xpart))

            ! grid to particle field gathering
            j=int(xpart*dx_inv)
            delj=xpart*dx_inv-j
            j=max(1,min(ngrid,j+1))
            epart=delj*efield(j)+(1.0-delj)*efield(j-1)

            ! update particle information
            x(i,k)=xa(i,k)+delt*v(i,k)
            w(i,k)=wa(i,k)+deltaf*(p0(i,k)-(1.0-linear)*w(i,k))*&
                delt*temp_inv(k)*qspecie(k)*v(i,k)*epart
            v(i,k)=va(i,k)+(1.0-linear)*delt*ratio_qm(k)*epart
        enddo
    enddo

    if(nrk==2)then
        do k=1,nspecie
            !$omp parallel do private(i)
            do i=1,nparticle

                ! periodic BC: out of bound particle recycled
                x(i,k)=x(i,k)*x_inv+10.0
                x(i,k)=xsize*(x(i,k)-aint(x(i,k)))
            enddo
        enddo
    endif

    end Subroutine particle
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Subroutine history

    use parameter
    implicit none
    integer,parameter :: nhistory=22,ndata=9+2*nmode
    integer k
    real time(nd)
    character(len=12) data_name(ndata)
    data data_name/'time*omega_p','<delt_n>/n_0','<(delt_f)^2>','<delt_v/v_T>',&
        '<del_v^2/2T>','<tolt_v/v_T>','<tol_v^2/2T>','<total_enen>',&
        '<pote_RMS/T>','mode1_real/T','mode1_imag/T','mode2_real/T',&
        'mode2_imag/T'/

    do k=1,nd
        time(k)=real((k-1)*ndiag)*tstep
    enddo

    ! open and write history data file
    open(nhistory,file='history.out',status='replace')
    write(nhistory,101)nd,ndata
    write(nhistory,100)data_name
    write(nhistory,102)time,density,entropy,flow,kin_ene,&
        flow_mar,kin_mar,field_ene,potential,phi_mode
    close(nhistory)

100 format(a12)
101 format(i6)
102 format(e14.5)

    end Subroutine history
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Subroutine r2cfft(n,rin,cout)
    ! real to complex transform. n: array size, must be power of 2
    ! rin: real array to be transformed,cout=\int{rin*exp(i*k*x)}dx
    ! cout: complex array of output for Fourier modes=[0,n/2+1], normalized by n

    implicit none
    integer,intent(in) :: n
    real,intent(in) :: rin(n)
    complex,intent(out) :: cout(n/2+1)
    integer i
    real fdata(2*n)

    do i=1,n
        fdata(2*i-1)=rin(i)
        fdata(2*i)=0.0
    enddo

    call fft(1,n,fdata)

    fdata=fdata/n
    do i=1,n/2+1
        cout(i)=cmplx(fdata(2*i-1),fdata(2*i))
    enddo

    end Subroutine r2cfft
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Subroutine c2rfft(n,cin,rout)
    ! complex to real transform. n: array size, must be power of 2
    ! cin: complex array of Fourier modes=[0,n/2+1] to be transformed
    ! rout: real array for output, rout=\int{cout*exp(-ikx)}dk

    implicit none
    integer,intent(in) :: n
    complex,intent(in) :: cin(n/2+1)
    real,intent(out) :: rout(n)
    integer i
    real fdata(2*n)

    do i=1,n/2+1
        fdata(2*i-1)=real(cin(i))
        fdata(2*i)=aimag(cin(i))
    enddo
    do i=1,n/2-1
        fdata(n+1+2*i)=fdata(n+1-2*i)
        fdata(n+2+2*i)=-fdata(n+2-2*i) !complex conjugate
    enddo

    call fft(-1,n,fdata)

    do i=1,n
        rout(i)=fdata(2*i-1)
    enddo

    end Subroutine c2rfft
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Subroutine c2cfft(ifft,n,cin,cout)
    ! complex to complex transform. ifft: direction of transform
    ! n: array size, must be power of 2
    ! cin: complex array to be transformed, cout=\int{cin*exp(ifft*i*k*x)}dx
    ! cout: complex array of output for Fourier modes=[0,n-1], normalized by n

    implicit none
    integer,intent(in) :: ifft,n
    complex,intent(in) :: cin(n)
    complex,intent(out) :: cout(n)
    integer i
    real fdata(2*n)

    do i=1,n
        fdata(2*i-1)=real(cin(i))
        fdata(2*i)=aimag(cin(i))
    enddo

    call fft(ifft,n,fdata)

    if(ifft==1)fdata=fdata/n
    do i=1,n
        cout(i)=cmplx(fdata(2*i-1),fdata(2*i))
    enddo

    end
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Subroutine fft(isign,nn,fdata)
    ! FFT basic program, Numerical Recipes, Page 501

    implicit none
    integer,intent(in) :: isign,nn
    real fdata(2*nn)
    integer i,istep,j,m,mmax,n
    real tempi,tempr
    !  real(kind=8) theta,wi,wr,wpi,wpr,wtemp,two_pi
    double precision theta,wi,wr,wpi,wpr,wtemp,two_pi

    two_pi=6.28318530717959d0
    n=2*nn
    j=1
    do i=1,n,2
        if(j.gt.i)then
            tempr=fdata(j)
            tempi=fdata(j+1)
            fdata(j)=fdata(i)
            fdata(j+1)=fdata(i+1)
            fdata(i)=tempr
            fdata(i+1)=tempi
        endif
        m=n/2
        do while((m.ge.2) .and. (j.gt.m))
            j=j-m
            m=m/2
        enddo
        j=j+m
    enddo
    mmax=2
    do while(n.gt.mmax)
        istep=2*mmax
        theta=two_pi/(isign*mmax)
        wpr=-2.0d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.0d0
        wi=0.0d0
        do m=1,mmax,2
            do i=m,n,istep
                j=i+mmax
                tempr=real(wr)*fdata(j)-real(wi)*fdata(j+1)
                tempi=real(wr)*fdata(j+1)+real(wi)*fdata(j)
                fdata(j)=fdata(i)-tempr
                fdata(j+1)=fdata(i+1)-tempi
                fdata(i)=fdata(i)+tempr
                fdata(i+1)=fdata(i+1)+tempi
            enddo
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
        enddo
        mmax=istep
    enddo

    end Subroutine fft
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !***********************************************************************************
    !*										   *
    !* 			Implimentation of collisional effects			   *
    !* 		  exponential subroutine was written by Hongpeng Qu		   *
    !* 		    collision subroutine was written by Onnie Luk		   *
    !*										   *
    !***********************************************************************************

    subroutine exponential(y,n,a)
    ! Generate an array of exponential distribution

    implicit none

    real, dimension(1:n,1):: x,y	!nparticle=n, nspecie=1
    integer n, i,a,isize,init
    integer,dimension(:),allocatable::iput
    intrinsic random_seed,random_number
    !  real :: total

    !total=0.0

    !! initialize Fortran random number generator
    !  call random_seed(size=isize)
    !  allocate(iget(isize),iput(isize))
    !  call random_seed(get=iget)
    !  do i=1,isize
    !!     call system_clock(init)   !random initialization
    !     init=1111                  !same initialization
    !     iput(i)=1111*init+iget(i)
    !  enddo
    !  call random_seed(put=iput)
    call random_number(x)

    do i=1,n
        y(i,1)=-log(x(i,1))  !  pay attention!!!
    enddo

    return
    end subroutine exponential

    !***********************
    subroutine collision
    use parameter
    implicit none

    real :: c,vpara, vperp, vtotal, pitch
    integer :: i, nseed

    real,dimension(1:nparticle) :: R

    c=nu*tstep
    call random_number(R(:))


    !$omp parallel do private(i,vpara,vperp,vtotal,pitch)
    do i=1,nparticle

        vpara = v(i,1)							!vpara  represents the physical quantity of v_parallel/v_thermal
        vperp = sqrt(max(small,mu(i,1)))					!vperp  represents the physical quantity of v_perpend./v_thermal
        vtotal= sqrt(vpara*vpara + vperp*vperp)				!vtotal represents the physical quantity of |velocity|/v_thermal
        pitch = vpara/vtotal						!old pitch angle
        pitch= max(-1.0,min( 1.0, &					!new pitch angle
            pitch*(1.0-c) + (R(i)-0.5)*sqrt(12.0*c*(1.0-pitch*pitch))))
        vpara = vtotal*pitch						!new vpara
        vperp = sqrt(vtotal*vtotal - vpara*vpara)	   			!new vperp

        v(i,1)     = vpara					!convert vpara back to v
        mu(i,1)    = max(small,vperp*vperp)			!convert vperp back to mu

    enddo

    return
    end subroutine collision

    !*************************
