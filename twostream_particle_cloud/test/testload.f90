!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!选择从排版较好的那版出发，测试full_f下运行
!***********************************************************************************
!* Modules for program PIC1D, Unmagnetized plasma use electron unit: n,m,e,T=1     *
!* Basic unit: Debye length sqrt(T/ne^2)=1, plasma frequency sqrt(ne^2/m)=1        *
!* compile on dragon: "pgf90 -mp pic1d.f90"                                        *
! use ifort
!* compile on linux sequentially: "ifort pic1d.f90"        *
!* compile on linux in parallel:  "ifort -qopenmp pic1d.f90"        *
! use gfortran 
!* compile on gfortran sequentially: "gfortran pic1d.f90"        *
!* compile on gfortran in parallel:  "gfortran -fopenmp pic1d.f90"        *
!* compile on gfortran in parallel:  "gfortran -fopenmp -O3 pic1d_diy_2stream_fullf.f90"        *
!*********************************************************************************** 
    Module parameter

        implicit none
    
        ! numerics
        integer,parameter ::  singlep=selected_real_kind(6),doublep=selected_real_kind(12)
        real,parameter :: pi=3.14159265358979,small=1.0e-10
    
        ! control parameter
        integer :: stdout=1,imarker=1 !stdout=0: output to display; 1: to file "PIC1D.out"
        integer :: parafile=114,field_E_file=514,totenergy_file=1919 !parafile 指代保存参数的文件
        !1: load physical PDF; 0: uniform for [-vmax,vmax]
        real :: vmax=5.0  ,vbeam=1.0          !velocity space range [-vmax,vmax], f_m(v=5)=0.000003
        real :: linear=0.0         !1.0 for linear simulation; 0.0 for nonlinear.
        real :: deltaf=0.0         !1.0 for delta_f simulation, 0.0 for full_f
    
        ! particle array: # of species, particles, velocity grids. First specie is electron
        integer,parameter :: nspecie=2,nparticle=60000,nvelocity=128,nbeam=60000
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
        integer,parameter :: ntime=400,ndiag=5,nmode=2
        integer,parameter :: nd=(ntime-1)/ndiag+1
        integer mode(nmode),nt
        data mode/1,2/             ! mode # for diagnostics
        real :: tstep=0.1          ! time step size
        real,dimension(nd) :: field_ene,potential ! field energy, potential RMS
        real,dimension(nd,2,nmode) :: phi_mode ! mode amplitudes
        real,dimension(nd,nspecie) :: density,entropy,flow,kin_ene,flow_mar,kin_mar
        ! volume average density, entropy, flow, kinetic energy, marker flow and energy
    
        end Module parameter


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
        ! 每个粒子的速度信息保存在二维数组v里，`v(i,k)`表示第k种第i个粒子的速度。子程序`load`中，源代码在`imarker==1`时进行速度的麦克斯韦分布：
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
    
            !add delta function to electrons,k==1
            if(k==2)then
                ! equal-spacing loading with random velocity
                do i=1, nparticle
                    x(i,k)=xsize*real(i)/real(nparticle)    !在空间上均分，分在0到xsize之间
                    !        x(i,k)=x(i,k)-winit*sin(xfreq*x(i,k))/xfreq
                    ! periodic BC: out of bound particle recycled
                    !        x(i,k)=x(i,k)/(real(ngrid)*deltax)+10.0
                    !        x(i,k)=real(ngrid)*deltax*(x(i,k)-aint(x(i,k)))
                enddo
            elseif(k==1)then
                ! equal-spacing loading with some specialized velocity
                do i=nbeam+1, nparticle
                    x(i,k)=xsize*real(i-nbeam)/real(nparticle-nbeam)
                    !        x(i,k)=x(i,k)-winit*sin(xfreq*x(i,k))/xfreq
                    ! periodic BC: out of bound particle recycled
                    !        x(i,k)=x(i,k)/(real(ngrid)*deltax)+10.0
                    !        x(i,k)=real(ngrid)*deltax*(x(i,k)-aint(x(i,k)))
                enddo
                do i=1,nbeam/2
                    v(i,1)=vbeam
                    x(i,1)=xsize*real(i)/real(nbeam/2)
                end do
                
                do i=nbeam/2+1,nbeam
                    v(i,1)=-vbeam
                    x(i,1)=xsize*real(i-nbeam/2)/real(nbeam/2)
                end do
            endif
    
            ! normalize velocity to thermal velocity
            v(:,k)=v(:,k)*sqrt(tspecie(k)/aspecie(k))
    
            ! normalize total marker # to nparticle
            tmp=sum(p0(:,k))
            p0(:,k)=p0(:,k)*real(nparticle)/tmp
    
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