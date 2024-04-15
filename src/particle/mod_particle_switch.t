module mod_particle_switch
  use mod_particle_base

  private

  integer, parameter :: Boris=1, Vay=2, HC=3, LM=4
  !> Variable index for gradient B, with relativistic correction 1/kappa
  !> where kappa = 1/sqrt(1 - E_perp^2/B^2)
  integer, protected, allocatable      :: grad_kappa_B(:)
  !> Variable index for (B . grad)B (curvature B drift)
  integer, protected, allocatable      :: b_dot_grad_b(:)

  ! ExB related drifts (vE = ExB/B^2)
  !> Variable index for curvature drift
  integer, protected, allocatable      :: vE_dot_grad_b(:)
  !> Variable index for polarization drift
  integer, protected, allocatable      :: b_dot_grad_vE(:)
  !> Variable index for polarization drift
  integer, protected, allocatable      :: vE_dot_grad_vE(:)
  !> Variable index for gradient |B|
  integer, protected, allocatable      :: grad_b(:)

  public :: switch_init
  public :: switch_create_particles
  integer, parameter :: RK4=1, ARK4=2

  ! Variables
  public :: bp, ep, grad_kappa_B, b_dot_grad_b, grad_b
  public :: vE_dot_grad_b, b_dot_grad_vE, vE_dot_grad_vE
  public :: vp, jp

contains

  subroutine switch_init()
    use mod_global_parameters
    integer :: idir, nwx

    if (physics_type/='mhd') call mpistop("GCA/Lorentz particles need magnetic field!")
    if (ndir/=3) call mpistop("GCA/Lorentz particles need ndir=3!")

    nwx = 0
    allocate(bp(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      bp(idir) = nwx
    end do
    allocate(ep(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      ep(idir) = nwx
    end do
    allocate(grad_kappa_B(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      grad_kappa_B(idir) = nwx
    end do
    allocate(b_dot_grad_b(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      b_dot_grad_b(idir) = nwx
    end do
    allocate(vE_dot_grad_b(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      vE_dot_grad_b(idir) = nwx
    end do
    allocate(b_dot_grad_vE(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      b_dot_grad_vE(idir) = nwx
    end do
    allocate(vE_dot_grad_vE(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      vE_dot_grad_vE(idir) = nwx
    end do
    allocate(vp(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      vp(idir) = nwx
    end do
    allocate(jp(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      jp(idir) = nwx
    end do
    allocate(grad_b(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      grad_b(idir) = nwx
    end do
    ngridvars=nwx

    particles_fill_gridvars => gca_fill_gridvars

    if (associated(particles_define_additional_gridvars)) then
      call particles_define_additional_gridvars(ngridvars)
    end if

    ! Only the Lorentz integraor can be specified, GCA integrator is default ARK4 for speed and stability
    select case(integrator_type_particles)
    case('Boris','boris')
      integrator = Boris
    case('Vay','vay')
      integrator = Vay
    case('HC','hc','higueracary')
      integrator = HC
    case('LM','lm','lapentamarkidis')
      integrator = LM
    case default
      integrator = Boris
    end select

    particles_integrate => switch_integrate_particles

  end subroutine switch_init

  subroutine switch_create_particles()
    ! initialise the particles
    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles, usr_update_payload, usr_check_particle

    double precision :: b(ndir), u(ndir), magmom
    double precision :: bnorm, lfac, vnorm, vperp, vpar
    integer          :: igrid, ipe_particle
    integer          :: n, idir, nparticles_local
    double precision :: x(3, num_particles)
    double precision :: v(3, num_particles)
    double precision :: q(num_particles)
    double precision :: m(num_particles)
    double precision :: rrd(num_particles,ndir)
    double precision :: defpayload(ndefpayload)
    double precision :: usrpayload(nusrpayload)
    logical          :: follow(num_particles), check
    double precision :: tp(num_particles)

    tp = 0.0d0

    if (mype==0) then
      if (.not. associated(usr_create_particles)) then
        ! Randomly distributed
        do idir=1,ndir
          do n = 1, num_particles
            rrd(n,idir) = rng%unif_01()
          end do
        end do
        do n=1, num_particles
          {^D&x(^D,n) = xprobmin^D + rrd(n,^D) * (xprobmax^D - xprobmin^D)\}
        end do
      else
        call usr_create_particles(num_particles, x, v, q, m, follow, tp)
      end if
    end if

    call MPI_BCAST(x,3*num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(v,3*num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(q,num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(m,num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(follow,num_particles,MPI_LOGICAL,0,icomm,ierrmpi)
    call MPI_BCAST(tp,num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)

    nparticles_local = 0

    ! first find ipe and igrid responsible for particle
    do n = 1, num_particles
      call find_particle_ipe(x(:, n),igrid,ipe_particle)

      particle(n)%igrid = igrid
      particle(n)%ipe   = ipe_particle

      if(ipe_particle == mype) then
        check = .true.

        ! Check for user-defined modifications or rejection conditions
        if (associated(usr_check_particle)) call usr_check_particle(igrid, x(:,n), v(:,n), q(n), m(n), follow(n), check)
        if (check) then
          call push_particle_into_particles_on_mype(n)
        else
          cycle
        end if

        nparticles_local = nparticles_local + 1

        allocate(particle(n)%self)
        particle(n)%self%x      = x(:, n)
        particle(n)%self%q      = q(n)
        particle(n)%self%m      = m(n)
        particle(n)%self%follow = follow(n)
        particle(n)%self%index  = n
        particle(n)%self%time   = global_time + tp(n)
        particle(n)%self%dt     = 0.0d0

        ! The momentum vector u(1:3) is filled with the following components
        ! gca_creation is used for inject tparticle>tmhd GCA particles for the unreachable future electromagnetic field in the initial phase, the velocity injected is in the GCA frame
        if (m(n)<0) then
          particle(n)%self%u(:) = v(:,n) 
        else
          call get_lfac_from_velocity(v(:, n), lfac)
          call get_vec(bp, igrid, x(:, n), particle(n)%self%time, b)

          bnorm = norm2(b(:))
          vnorm = norm2(v(:, n))
          vpar  = sum(v(:, n) * b/bnorm)
          vperp = sqrt(vnorm**2 - vpar**2)
          ! parallel momentum component (gamma v||)
          particle(n)%self%u(1) = lfac * vpar

          ! Mr: the conserved magnetic moment
          magmom = m(n) * (vperp * lfac)**2 / (2.0d0 * bnorm)
          particle(n)%self%u(2) = magmom

          ! Lorentz factor
          particle(n)%self%u(3) = lfac
        end if

        ! initialise payloads for GCA module
        allocate(particle(n)%payload(npayload))
        call switch_update_payload(igrid,x(:,n),particle(n)%self%u(:),q(n),m(n),defpayload,ndefpayload,0.d0)
        particle(n)%payload(1:ndefpayload) = defpayload
        if (associated(usr_update_payload)) then
          call usr_update_payload(igrid,x(:,n),particle(n)%self%u(:),q(n),m(n),usrpayload,nusrpayload,0.d0)
          particle(n)%payload(ndefpayload+1:npayload) = usrpayload
        end if
      end if
    end do

    call MPI_ALLREDUCE(nparticles_local,nparticles,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)

  end subroutine switch_create_particles

  subroutine gca_fill_gridvars
    use mod_global_parameters
    use mod_usr_methods, only: usr_particle_fields
    use mod_geometry

    integer                                   :: igrid, iigrid, idir, idim
    double precision, dimension(ixG^T,1:ndir) :: beta
    double precision, dimension(ixG^T,1:nw)   :: w,wold
    double precision                          :: current(ixG^T,7-2*ndir:3)
    integer                                   :: idirmin
    double precision, dimension(ixG^T,1:ndir) :: vE, bhat
    double precision, dimension(ixG^T)        :: kappa, kappa_B, absB, tmp

    call fill_gridvars_default()

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      w(ixG^T,1:nw) = ps(igrid)%w(ixG^T,1:nw)
      call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
      ! fill with velocity:
      gridvars(igrid)%w(ixG^T,vp(:)) = w(ixG^T,iw_mom(:))

      ! grad(kappa B)
      absB(ixG^T) = sqrt(sum(gridvars(igrid)%w(ixG^T,bp(:))**2,dim=ndim+1))
      vE(ixG^T,1) = gridvars(igrid)%w(ixG^T,ep(2)) * gridvars(igrid)%w(ixG^T,bp(3)) &
            - gridvars(igrid)%w(ixG^T,ep(3)) * gridvars(igrid)%w(ixG^T,bp(2))
      vE(ixG^T,2) = gridvars(igrid)%w(ixG^T,ep(3)) * gridvars(igrid)%w(ixG^T,bp(1)) &
            - gridvars(igrid)%w(ixG^T,ep(1)) * gridvars(igrid)%w(ixG^T,bp(3))
      vE(ixG^T,3) = gridvars(igrid)%w(ixG^T,ep(1)) * gridvars(igrid)%w(ixG^T,bp(2)) &
            - gridvars(igrid)%w(ixG^T,ep(2)) * gridvars(igrid)%w(ixG^T,bp(1))
      do idir=1,ndir
        where (absB(ixG^T) .gt. 0.d0)
          vE(ixG^T,idir) = vE(ixG^T,idir) / absB(ixG^T)**2
        elsewhere
          vE(ixG^T,idir) = 0.d0
        end where
      end do
      if (any(sum(vE(ixG^T,:)**2,dim=ndim+1) .ge. c_norm**2) .and. relativistic) then
        call mpistop("GCA FILL GRIDVARS: vE>c! ABORTING...")
      end if
      if (any(vE .ne. vE)) then
        call mpistop("GCA FILL GRIDVARS: NaNs IN vE! ABORTING...")
      end if

      if (relativistic) then
        kappa(ixG^T) = 1.d0/sqrt(1.0d0 - sum(vE(ixG^T,:)**2,dim=ndim+1)/c_norm**2)
      else
        kappa(ixG^T) = 1.d0
      end if
      kappa_B(ixG^T) = absB(ixG^T) / kappa(ixG^T)

      if (any(kappa_B .ne. kappa_B)) then
        call mpistop("GCA FILL GRIDVARS: NaNs IN kappa_B! ABORTING...")
      end if

      tmp=0.d0
      do idim=1,ndim
        call gradient(kappa_B,ixG^LL,ixG^LL^LSUB1,idim,tmp)
        gridvars(igrid)%w(ixG^T,grad_kappa_B(idim)) = tmp(ixG^T)
      end do

      ! gradient absB
      tmp = 0.d0
      do idim=1,ndim
        call gradient(absB,ixG^LL,ixG^LL^LSUB1,idim,tmp)
        gridvars(igrid)%w(ixG^T,grad_B(idim)) = tmp(ixG^T)
      end do

      ! bhat
      do idir=1,ndir
        where (absB(ixG^T) .gt. 0.d0)
          bhat(ixG^T,idir) = gridvars(igrid)%w(ixG^T,bp(idir)) / absB(ixG^T)
        elsewhere
          bhat(ixG^T,idir) = 0.d0
        end where
      end do
      if (any(bhat .ne. bhat)) then
        call mpistop("GCA FILL GRIDVARS: NaNs IN bhat! ABORTING...")
      end if

      do idir=1,ndir
        ! (b dot grad) b and the other directional derivatives
        do idim=1,ndim
          call gradient(bhat(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
          gridvars(igrid)%w(ixG^T,b_dot_grad_b(idir)) = gridvars(igrid)%w(ixG^T,b_dot_grad_b(idir)) &
               + bhat(ixG^T,idim) * tmp(ixG^T)
          gridvars(igrid)%w(ixG^T,vE_dot_grad_b(idir)) = gridvars(igrid)%w(ixG^T,vE_dot_grad_b(idir)) &
               + vE(ixG^T,idim) * tmp(ixG^T)
          call gradient(vE(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
          gridvars(igrid)%w(ixG^T,b_dot_grad_vE(idir)) = gridvars(igrid)%w(ixG^T,b_dot_grad_vE(idir)) &
               + bhat(ixG^T,idim) * tmp(ixG^T)
          gridvars(igrid)%w(ixG^T,vE_dot_grad_vE(idir)) = gridvars(igrid)%w(ixG^T,vE_dot_grad_vE(idir)) &
               + vE(ixG^T,idim) * tmp(ixG^T)
        end do
      end do

      if (time_advance) then
        ! Fluid velocity
        w(ixG^T,1:nw) = pso(igrid)%w(ixG^T,1:nw)
        call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
        gridvars(igrid)%wold(ixG^T,vp(:)) = w(ixG^T,iw_mom(:))

        ! grad(kappa B)
        absB(ixG^T) = sqrt(sum(gridvars(igrid)%wold(ixG^T,bp(:))**2,dim=ndim+1))
        vE(ixG^T,1) = gridvars(igrid)%wold(ixG^T,ep(2)) * gridvars(igrid)%wold(ixG^T,bp(3)) &
             - gridvars(igrid)%wold(ixG^T,ep(3)) * gridvars(igrid)%wold(ixG^T,bp(2))
        vE(ixG^T,2) = gridvars(igrid)%wold(ixG^T,ep(3)) * gridvars(igrid)%wold(ixG^T,bp(1)) &
             - gridvars(igrid)%wold(ixG^T,ep(1)) * gridvars(igrid)%wold(ixG^T,bp(3))
        vE(ixG^T,3) = gridvars(igrid)%wold(ixG^T,ep(1)) * gridvars(igrid)%wold(ixG^T,bp(2)) &
             - gridvars(igrid)%wold(ixG^T,ep(2)) * gridvars(igrid)%wold(ixG^T,bp(1))
        do idir=1,ndir
          where (absB(ixG^T) .gt. 0.d0) 
            vE(ixG^T,idir) = vE(ixG^T,idir) / absB(ixG^T)**2
          elsewhere
            vE(ixG^T,idir) = 0.d0
          end where
        end do
        if (any(sum(vE(ixG^T,:)**2,dim=ndim+1) .ge. c_norm**2) .and. relativistic) then
          call mpistop("GCA FILL GRIDVARS: vE>c! ABORTING...")
        end if
        if (any(vE .ne. vE)) then
          call mpistop("GCA FILL GRIDVARS: NaNs IN vE! ABORTING...")
        end if

        if (relativistic) then
          kappa(ixG^T) = 1.d0/sqrt(1.0d0 - sum(vE(ixG^T,:)**2,dim=ndim+1)/c_norm**2)
        else
          kappa(ixG^T) = 1.d0
        end if
        kappa_B(ixG^T) = absB(ixG^T) / kappa(ixG^T)
        if (any(kappa_B .ne. kappa_B)) then
          call mpistop("GCA FILL GRIDVARS: NaNs IN kappa_B! ABORTING...")
        end if

        do idim=1,ndim
          call gradient(kappa_B,ixG^LL,ixG^LL^LSUB1,idim,tmp)
          gridvars(igrid)%wold(ixG^T,grad_kappa_B(idim)) = tmp(ixG^T)
        end do

        do idir=1,ndir
          where (absB(ixG^T) .gt. 0.d0)
            bhat(ixG^T,idir) = gridvars(igrid)%wold(ixG^T,bp(idir)) / absB(ixG^T)
          elsewhere
            bhat(ixG^T,idir) = 0.d0
          end where
        end do
        if (any(bhat .ne. bhat)) then
          call mpistop("GCA FILL GRIDVARS: NaNs IN bhat! ABORTING...")
        end if
  
        do idir=1,ndir
          ! (b dot grad) b and the other directional derivatives
          do idim=1,ndim
            call gradient(bhat(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
            gridvars(igrid)%wold(ixG^T,b_dot_grad_b(idir)) = gridvars(igrid)%wold(ixG^T,b_dot_grad_b(idir)) &
                 + bhat(ixG^T,idim) * tmp(ixG^T)
            gridvars(igrid)%wold(ixG^T,vE_dot_grad_b(idir)) = gridvars(igrid)%wold(ixG^T,vE_dot_grad_b(idir)) &
                 + vE(ixG^T,idim) * tmp(ixG^T)
            call gradient(vE(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
            gridvars(igrid)%wold(ixG^T,b_dot_grad_vE(idir)) = gridvars(igrid)%wold(ixG^T,b_dot_grad_vE(idir)) &
                 + bhat(ixG^T,idim) * tmp(ixG^T)
            gridvars(igrid)%wold(ixG^T,vE_dot_grad_vE(idir)) = gridvars(igrid)%wold(ixG^T,vE_dot_grad_vE(idir)) &
                 + vE(ixG^T,idim) * tmp(ixG^T)
          end do
        end do
      end if
    end do

  end subroutine gca_fill_gridvars

  function check_switch_condition(ipart_working, igrid_working) result(DoSwitch)
    use mod_global_parameters
    integer, intent(in)  :: igrid_working, ipart_working
    logical              :: DoSwitch 

    double precision, dimension(ndir)  :: x, u, b, e
    double precision                   :: m,q,t,gamma
    double precision, dimension(ndir)  :: bhat, gradb, rlvector, vE, vp
    double precision                   :: rl, absb, ugyro2
    double precision, dimension(ndir)  :: rg, tmp1, e1, e2
    double precision                   :: theta

    double precision, dimension(ndir)  :: bdotgradb, vEdotgradb, gradkappaB, bdotgradvE, vEdotgradvE, ud, utmp1, utmp2, utmp3
    double precision                   :: Mr, upar, epar, kappa, rcr

    DoSwitch = .false.

    x = particle(ipart_working)%self%x
    u = particle(ipart_working)%self%u
    m = particle(ipart_working)%self%m
    q = particle(ipart_working)%self%q
    t = particle(ipart_working)%self%time

    call get_vec(bp, igrid_working, x, t, b)
    call get_vec(ep, igrid_working, x, t, e)
    call get_vec(grad_b, igrid_working, x, t, gradb)
    absb = norm2(b)
    bhat = b/absb
    if (absb .gt. 0.0d0) then
      bhat = b / absb
    else
      bhat = 0.d0
    end if
    call cross(e,bhat,vE)
    if (absb .gt. 0.0d0) then
      vE = vE/absb
    else 
      vE = 0.d0
    end if
    ! GCA -> Full Lorentz
    if (m > 0) then
      rl = sqrt(u(2)*m*2/absb)/abs(q)
      rcr = dsqrt(rl**2+(u(1)*m/absb/q)**2)/absb*norm2(gradb)
      if (rcr > Kcr) then
        DoSwitch = .true.
        print*, "Particle change from gca to florentz: ", ipart_working
        Mr = u(2); upar = u(1) ; epar=sum(e*bhat)
        if (relativistic) then
          kappa = 1.d0/sqrt(1.0d0 - sum(vE(:)**2)/c_norm**2)
          gamma = sqrt(1.0d0+upar**2/c_norm**2+2.0d0*Mr*absb/m/c_norm**2)*kappa
        else
          kappa = 1.d0
          gamma = 1.d0
        end if
        ! call get_vec(b_dot_grad_b, igrid_working,x,t,bdotgradb)
        ! call get_vec(vE_dot_grad_b, igrid_working,x,t,vEdotgradb)
        ! call get_vec(grad_kappa_B, igrid_working,x,t,gradkappaB)
        ! call get_vec(b_dot_grad_vE, igrid_working,x,t,bdotgradvE)
        ! call get_vec(vE_dot_grad_vE, igrid_working,x,t,vEdotgradvE)
        ! if (absb .gt. 0.d0) then
        !   utmp1(1:ndir) = bhat(1:ndir)/(absb/kappa**2)
        ! else
        !   utmp1 = 0.d0
        ! end if
        ! utmp2(1:ndir) = Mr/(gamma*q)*gradkappaB(1:ndir) &
        !     + m/q* (upar**2/gamma*bdotgradb(1:ndir) + upar*vEdotgradb(1:ndir) &
        !             + upar*bdotgradvE(1:ndir) + gamma*vEdotgradvE(1:ndir))
        ! if (relativistic) then
        !   utmp2(1:ndir) = utmp2(1:ndir) + upar*epar/(gamma)*vE(1:ndir)
        ! end if
        ! call cross(utmp1,utmp2,utmp3)
        ud = vE(1:ndir) ! + utmp3(1:ndir)

        if (norm2(x) .eq. 0) then
          rg = [0,0,1]
        else 
          rg = x/norm2(x)
        end if
        call cross(rg,bhat,tmp1)
        call cross(tmp1,bhat,e1)
        call cross(bhat,e1,e2)
        theta = rng%unif_01() * dpi*2
        particle(ipart_working)%self%x = x+rl*(e2*cos(theta)+e1*(sin(theta)))
        particle(ipart_working)%self%u = u(1)*bhat+sqrt(u(2)*2*absb/abs(m))*(e1*cos(theta)+e2*sin(theta))+gamma*ud
        particle(ipart_working)%self%m = -m
      end if
    ! Full Lorentz -> GCA
    else
      if (relativistic) then
        call get_lfac(u,gamma)
      else
        gamma = 1.d0
      end if
      vp = u/gamma - vE
      call cross(vp,bhat,rlvector)
      if (absb .gt. 0) then
        rlvector = gamma*rlvector*abs(m)/absb/abs(q)
      else
        call mpistop("switch condition founds that absb lt 0")
      end if
      rcr = norm2(vp)*gamma*abs(m)/absb/abs(q)*norm2(gradb)/absb
      if (rcr < 0.9*Kcr) then
        DoSwitch = .true.
        print *, "rcr: ", rcr
        print*, "particle change from florentz to gca: ", ipart_working 
        particle(ipart_working)%self%x = x-rlvector
        particle(ipart_working)%self%u(1) = sum(u*bhat)
        ugyro2 = sum(u**2)-sum((u*bhat)**2)-gamma**2*sum(vE**2)
        particle(ipart_working)%self%u(2) = gamma**2*abs(m)*ugyro2/2/absb
        particle(ipart_working)%self%u(3) = gamma
        particle(ipart_working)%self%m = -m
      end if
    end if

  end function check_switch_condition

  subroutine switch_integrate_particles(end_time)
    use mod_odeint
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_create_particles, usr_update_payload
    double precision, intent(in)        :: end_time

    double precision                    :: lfac, absS
    double precision                    :: defpayload(ndefpayload)
    double precision                    :: usrpayload(nusrpayload)
    double precision                    :: dt_p, tloc, y(ndir+2),dydt(ndir+2),ytmp(ndir+2), euler_cfl, int_factor
    double precision                    :: tk, k1(ndir+2),k2(ndir+2),k3(ndir+2),k4(ndir+2)
    double precision                    :: xp(ndir), xpm(ndir), xpc(ndir), xpcm(ndir)
    double precision                    :: up(ndir), upc(ndir), tp, rho, rhoold, td
    double precision, dimension(1:ndir) :: vE, e, b, bc, ec, bhat, vfluid, current
    double precision, dimension(1:ndir) :: drift1, drift2
    double precision, dimension(1:ndir) :: drift3, drift4, drift5, drift6, drift7
    double precision, dimension(1:ndir) :: bdotgradb, vEdotgradb, gradkappaB
    double precision, dimension(1:ndir) :: bdotgradvE, vEdotgradvE
    double precision, dimension(1:ndir) :: gradBdrift, reldrift, bdotgradbdrift
    double precision, dimension(1:ndir) :: vEdotgradbdrift, bdotgradvEdrift
    double precision, dimension(1:ndir) :: vEdotgradvEdrift
    double precision                    :: kappa, Mr, upar, m, absb, gamma, q, mompar, vpar, vEabs
    double precision                    :: gradBdrift_abs, reldrift_abs, epar
    double precision                    :: bdotgradbdrift_abs, vEdotgradbdrift_abs
    double precision                    :: bdotgradvEdrift_abs, vEdotgradvEdrift_abs
    double precision                    :: momentumpar1, momentumpar2, momentumpar3, momentumpar4
    ! Precision of time-integration:
    double precision,parameter          :: eps=1.0d-6
    ! for odeint:
    double precision                    :: h1, hmin, h_old
    integer                             :: nok, nbad, ic1^D, ic2^D, ierror, nvar
    integer                             :: ipart, iipart, seed, ic^D,igrid_particle, ipe_particle 

    logical                             :: DoSwitch

    nvar=ndir+2

    do iipart=1,nparticles_active_on_mype

      ipart                   = particles_active_on_mype(iipart)
      igrid_working           = particle(ipart)%igrid
      ipart_working           = particle(ipart)%self%index

      ! check whether to switch between GCA and Lorentz
      DoSwitch = check_switch_condition(ipart_working, igrid_working)

      tloc  = particle(ipart)%self%time
      q     = particle(ipart)%self%q
      m     = particle(ipart)%self%m
      xp    = particle(ipart)%self%x
      up    = particle(ipart)%self%u
      tp    = particle(ipart)%self%time  

      if (m>0) then
        dt_p = gca_get_particle_dt(particle(ipart), end_time)
        particle(ipart)%self%dt = dt_p
        ! Adaptive stepwidth RK4:
        ! initial solution vector:
        y(1:ndir) = xp(1:ndir) ! position of guiding center
        y(ndir+1) = particle(ipart)%self%u(1) ! parallel momentum component (gamma v||)
        y(ndir+2) = particle(ipart)%self%u(2) ! conserved magnetic moment Mr
      ! y(ndir+3) = particle(ipart)%self%u(3) ! Lorentz factor of particle
  
        ! we temporarily save the solution vector, to replace the one from the euler
        ! timestep after euler integration
        ytmp=y

        ! euler step check is not neccesary 
        y(1:ndir+2) = ytmp(1:ndir+2)
  
        particle(ipart)%self%x(1:ndir) = ytmp(1:ndir) ! position of guiding center
        particle(ipart)%self%u(1)      = ytmp(ndir+1) ! parallel momentum component (gamma v||)
        particle(ipart)%self%u(2)      = ytmp(ndir+2) ! conserved magnetic moment
  
        ! specify a minimum step hmin. If the timestep reaches this minimum, multiply by
        ! a factor 100 to make sure the RK integration doesn't crash
        h1 = dt_p/2.0d0; hmin=1.0d-9; h_old=dt_p/2.0d0
  
        if(h1 .lt. hmin)then
          h1=hmin
          dt_p=2.0d0*h1
        endif
  
        if (any(y .ne. y)) then
          call mpistop("NaNs DETECTED IN GCA_INTEGRATE BEFORE ODEINT CALL! ABORTING...")
        end if
  
        ! RK4 integration with adaptive stepwidth
        call odeint(y,nvar,tloc,tloc+dt_p,eps,h1,hmin,nok,nbad,derivs_gca_rk,rkqs,ierror)
  
        if (ierror /= 0) then
           print *, "odeint returned error code", ierror
           print *, "1 means hmin too small, 2 means MAXSTP exceeded"
           print *, "Having a problem with particle", iipart
        end if
  
        if (any(y .ne. y)) then
          call mpistop("NaNs DETECTED IN GCA_INTEGRATE AFTER ODEINT CALL! ABORTING...")
        end if
  
        ! original RK integration without interpolation in ghost cells
  !       call odeint(y,nvar,tloc,tloc+dt_p,eps,h1,hmin,nok,nbad,derivs_gca,rkqs)
  
        ! final solution vector after rk integration
        particle(ipart)%self%x(1:ndir) = y(1:ndir)
        particle(ipart)%self%u(1)      = y(ndir+1)
        particle(ipart)%self%u(2)      = y(ndir+2)
        !particle(ipart)%self%u(3)      = y(ndir+3)

        ! now calculate other quantities, mean Lorentz factor, drifts, perpendicular velocity:
        call get_vec(bp, igrid_working,y(1:ndir),tloc+dt_p,b)
        call get_vec(vp, igrid_working,y(1:ndir),tloc+dt_p,vfluid)
        call get_vec(jp, igrid_working,y(1:ndir),tloc+dt_p,current)
        call get_vec(ep, igrid_working,y(1:ndir),tloc+dt_p,e)

        absb         = sqrt(sum(b(:)**2))
        bhat(1:ndir) = b(1:ndir) / absb

        epar         = sum(e(:)*bhat(:))

        call cross(e,bhat,vE)

        vE(1:ndir) = vE(1:ndir) / absb
        vEabs = sqrt(sum(vE(:)**2))
        if (relativistic) then
          kappa = 1.d0/sqrt(1.0d0 - sum(vE(:)**2)/c_norm**2)
        else
          kappa = 1.d0
        end if
        Mr = y(ndir+2); upar = y(ndir+1)
        if (relativistic) then
          gamma = sqrt(1.0d0+upar**2/c_norm**2+2.0d0*Mr*absb/m/c_norm**2)*kappa
        else
          gamma = 1.d0
        end if
        particle(ipart)%self%u(3) = gamma
      else
        dt_p = Lorentz_get_particle_dt(particle(ipart), end_time)
        particle(ipart)%self%dt = dt_p

        ! Push particle over first half time step
        ! all pushes done in Cartesian coordinates
        call partcoord_to_cartesian(xp,xpc)
        call partvec_to_cartesian(xp,up,upc)
        call get_lfac(upc,lfac)
        xpcm = xpc + 0.5d0 * dt_p * upc/lfac
        call partcoord_from_cartesian(xpm,xpcm)
        ! Fix xp if the 0,2*pi boundary was crossed in cylindrical/spherical coords
        call fix_phi_crossing(xpm,igrid_working)
        
        ! Get E, B at n+1/2 position
        call get_vec(bp, igrid_working, xpm, tp+dt_p/2.d0, b)
        call get_vec(vp, igrid_working, xpm, tp+dt_p/2.d0, vfluid)
        call get_vec(jp, igrid_working, xpm, tp+dt_p/2.d0, current)
        call get_vec(ep, igrid_working, xpm ,tp+dt_p/2.d0, e)
        if (particles_etah > zero) then
          call interpolate_var(igrid_working,ixG^LL,ixM^LL,ps(igrid_working)%w(ixG^T,1),ps(igrid_working)%x,xpm,rho)
          if (time_advance) then
            td = (tp+dt_p/2.d0 - global_time) / dt
            call interpolate_var(igrid_working,ixG^LL,ixM^LL,pso(igrid_working)%w(ixG^T,1),ps(igrid_working)%x,xpm,rhoold)
            rho = rhoold * (1.d0-td) + rho * td
          end if
          select case (coordinate)
          case (Cartesian,Cartesian_stretched,spherical)
            e(1) = e(1) + particles_etah/rho * (current(2)*b(3) - current(3)*b(2))
            e(2) = e(2) + particles_etah/rho * (-current(1)*b(3) + current(3)*b(1))
            e(3) = e(3) + particles_etah/rho * (current(1)*b(2) - current(2)*b(1))
          case (cylindrical)
            e(r_) = e(r_) + particles_etah/rho * (current(phi_)*b(z_) - current(z_)*b(phi_))
            e(phi_) = e(phi_) + particles_etah/rho * (-current(r_)*b(z_) + current(z_)*b(r_))
            e(z_) = e(z_) + particles_etah/rho * (current(r_)*b(phi_) - current(phi_)*b(r_))
          end select
        end if
        ! Convert fields to Cartesian frame
        call partvec_to_cartesian(xpm,b,bc)
        call partvec_to_cartesian(xpm,e,ec)

        ! 'Kick' particle (update velocity) based on the chosen integrator
        call Lorentz_kick(upc,ec,bc,q,abs(m),dt_p)

        ! Push particle over second half time step
        ! all pushes done in Cartesian coordinates
        call get_lfac(upc,lfac)
        xpc = xpcm + 0.5d0 * dt_p * upc/lfac
        call partcoord_from_cartesian(xp,xpc)
        ! Fix xp if the 0,2*pi boundary was crossed in cylindrical/spherical coords
        call fix_phi_crossing(xp,igrid_working)
        call partvec_from_cartesian(xp,up,upc)

        ! Store updated x,u
        particle(ipart)%self%x = xp
        particle(ipart)%self%u = up
      end if

      ! Time update
      particle(ipart)%self%time = particle(ipart)%self%time + dt_p

      ! Update payload
      call switch_update_payload(igrid_working,&
             particle(ipart)%self%x,particle(ipart)%self%u,q,m,defpayload,ndefpayload,particle(ipart)%self%time)
      particle(ipart)%payload(1:ndefpayload) = defpayload
      if (associated(usr_update_payload)) then
        call usr_update_payload(igrid_working,&
             particle(ipart)%self%x,particle(ipart)%self%u,q,m,usrpayload,nusrpayload,particle(ipart)%self%time)
        particle(ipart)%payload(ndefpayload+1:npayload) = usrpayload
      end if

    end do

  end subroutine switch_integrate_particles

  subroutine derivs_gca_rk(t_s,y,dydt)
    use mod_global_parameters

    double precision                :: t_s, y(ndir+2)
    double precision                :: dydt(ndir+2)

    double precision,dimension(ndir):: vE, b, e, x, bhat, bdotgradb, vEdotgradb, gradkappaB, vfluid, current
    double precision,dimension(ndir):: bdotgradvE, vEdotgradvE, u, utmp1, utmp2, utmp3
    double precision                :: upar, Mr, gamma, absb, q, m, epar, kappa
    integer                         :: ic^D

    ! Here the terms in the guiding centre equations of motion are interpolated for
    ! the RK integration. The interpolation is also done in the ghost cells such
    ! that the RK integration does not give an error

    q = particle(ipart_working)%self%q
    m = particle(ipart_working)%self%m

    x(1:ndir) = y(1:ndir)
    upar      = y(ndir+1) ! gamma v||
    Mr        = y(ndir+2)
    !gamma     = y(ndir+3)

    if (any(x .ne. x)) then
      write(*,*) "ERROR IN DERIVS_GCA_RK: NaNs IN X OR Y!"
      write(*,*) "x",x
      write(*,*) "y",y(ndir+1:ndir+2)
      call mpistop("ABORTING...")
    end if

    call get_vec(bp, igrid_working,x,t_s,b)
    call get_vec(vp, igrid_working,x,t_s,vfluid)
    call get_vec(jp, igrid_working,x,t_s,current)
    call get_vec(ep, igrid_working,x,t_s,e)
    call get_vec(b_dot_grad_b, igrid_working,x,t_s,bdotgradb)
    call get_vec(vE_dot_grad_b, igrid_working,x,t_s,vEdotgradb)
    call get_vec(grad_kappa_B, igrid_working,x,t_s,gradkappaB)
    call get_vec(b_dot_grad_vE, igrid_working,x,t_s,bdotgradvE)
    call get_vec(vE_dot_grad_vE, igrid_working,x,t_s,vEdotgradvE)

    if (any(b .ne. b) .or. any(e .ne. e) &
        .or. any(bdotgradb .ne. bdotgradb) .or. any(vEdotgradb .ne. vEdotgradb) &
        .or. any(gradkappaB .ne. gradkappaB) .or. any(bdotgradvE .ne. bdotgradvE) &
        .or. any(vEdotgradvE .ne. vEdotgradvE)) then
      write(*,*) "ERROR IN DERIVS_GCA_RK: NaNs IN FIELD QUANTITIES!"
      write(*,*) "b",b
      write(*,*) "e",e
      write(*,*) "bdotgradb",bdotgradb
      write(*,*) "vEdotgradb",vEdotgradb
      write(*,*) "gradkappaB",gradkappaB
      write(*,*) "bdotgradvE",bdotgradvE
      write(*,*) "vEdotgradvE",vEdotgradvE
      call mpistop("ABORTING...")
    end if

    absb  = sqrt(sum(b(:)**2))
    if (absb .gt. 0.d0) then
      bhat(1:ndir) = b(1:ndir) / absb
    else
      bhat = 0.d0
    end if
    epar = sum(e(:)*bhat(:))

    call cross(e,bhat,vE)
    if (absb .gt. 0.d0) then
      vE(1:ndir) = vE(1:ndir) / absb
    else
      vE(1:ndir) = 0.d0
    end if

    if (relativistic) then
      kappa = 1.d0/sqrt(1.0d0 - sum(vE(:)**2)/c_norm**2)
      gamma = sqrt(1.0d0+upar**2/c_norm**2+2.0d0*Mr*absb/m/c_norm**2)*kappa
    else
      kappa = 1.d0
      gamma = 1.d0
    end if

    if (absb .gt. 0.d0) then
      utmp1(1:ndir) = bhat(1:ndir)/(absb/kappa**2)
    else
      utmp1 = 0.d0
    end if
    utmp2(1:ndir) = Mr/(gamma*q)*gradkappaB(1:ndir) &
         + m/q* (upar**2/gamma*bdotgradb(1:ndir) + upar*vEdotgradb(1:ndir) &
                 + upar*bdotgradvE(1:ndir) + gamma*vEdotgradvE(1:ndir))
    if (relativistic) then
      utmp2(1:ndir) = utmp2(1:ndir) + upar*epar/(gamma)*vE(1:ndir)
    end if

    call cross(utmp1,utmp2,utmp3)
    u(1:ndir) = vE(1:ndir) + utmp3(1:ndir)

    ! done assembling the terms, now write rhs:
    dydt(1:ndir) = ( u(1:ndir) + upar/gamma * bhat(1:ndir) )
    dydt(ndir+1) = q/m*epar - Mr/(m*gamma) * sum(bhat(:)*gradkappaB(:)) &
                   + sum(vE(:)*(upar*bdotgradb(:)+gamma*vEdotgradb(:)))
    dydt(ndir+2) = 0.0d0 ! magnetic moment is conserved

  end subroutine derivs_gca_rk

  subroutine derivs_gca(t_s,y,dydt)
    use mod_global_parameters

    double precision                :: t_s, y(ndir+2)
    double precision                :: dydt(ndir+2)

    double precision,dimension(ndir):: vE, b, e, x, bhat, bdotgradb, vEdotgradb, gradkappaB, vfluid, current
    double precision,dimension(ndir):: bdotgradvE, vEdotgradvE, u, utmp1, utmp2, utmp3
    double precision                :: upar, Mr, gamma, absb, q, m, epar, kappa

    ! Here the normal interpolation is done for the terms in the GCA equations of motion

    q = particle(ipart_working)%self%q
    m = particle(ipart_working)%self%m

    x(1:ndir) = y(1:ndir)
    upar      = y(ndir+1) ! gamma v||
    Mr        = y(ndir+2)
    !gamma     = y(ndir+3)

    if (any(x .ne. x)) then
      call mpistop("ERROR IN DERIVS_GCA: NaNs IN X! ABORTING...")
    end if

    call get_vec(bp, igrid_working,x,t_s,b)
    call get_vec(vp, igrid_working,x,t_s,vfluid)
    call get_vec(jp, igrid_working,x,t_s,current)
    call get_vec(ep, igrid_working,x,t_s,e)
    call get_vec(b_dot_grad_b, igrid_working,x,t_s,bdotgradb)
    call get_vec(vE_dot_grad_b, igrid_working,x,t_s,vEdotgradb)
    call get_vec(grad_kappa_B, igrid_working,x,t_s,gradkappaB)
    call get_vec(b_dot_grad_vE, igrid_working,x,t_s,bdotgradvE)
    call get_vec(vE_dot_grad_vE, igrid_working,x,t_s,vEdotgradvE)

    absb         = sqrt(sum(b(:)**2))
    if (absb .gt. 0.d0) then
      bhat(1:ndir) = b(1:ndir) / absb
    else
      bhat = 0.d0
    end if

    epar         = sum(e(:)*bhat(:))
    call cross(e,bhat,vE)
    if (absb .gt. 0.d0) vE(1:ndir)   = vE(1:ndir) / absb

    if (relativistic) then
      kappa = 1.d0/sqrt(1.0d0 - sum(vE(:)**2)/c_norm**2)
      gamma = sqrt(1.0d0+upar**2/c_norm**2+2.0d0*Mr*absb/m/c_norm**2)*kappa
    else
      kappa = 1.d0
      gamma = 1.d0
    end if
    if (absb .gt. 0.d0) then
      utmp1(1:ndir) = bhat(1:ndir)/(absb/kappa**2)
    else
      utmp1 = 0.d0
    end if
    utmp2(1:ndir) = Mr/(gamma*q)*gradkappaB(1:ndir) &
         + m/q* (upar**2/gamma*bdotgradb(1:ndir) + upar*vEdotgradb(1:ndir) &
                 + upar*bdotgradvE(1:ndir) + gamma*vEdotgradvE(1:ndir))
    if (relativistic) then
      utmp2(1:ndir) = utmp2(1:ndir) + upar*epar/(gamma)*vE(1:ndir)
    end if

    call cross(utmp1,utmp2,utmp3)
    u(1:ndir) = vE(1:ndir) + utmp3(1:ndir)

    ! done assembling the terms, now write rhs:
    dydt(1:ndir) = ( u(1:ndir) + upar/gamma * bhat(1:ndir) )
    dydt(ndir+1) = q/m*epar - Mr/(m*gamma) * sum(bhat(:)*gradkappaB(:)) &
                   + sum(vE(:)*(upar*bdotgradb(:)+gamma*vEdotgradb(:)))
    dydt(ndir+2) = 0.0d0 ! magnetic moment is conserved

  end subroutine derivs_gca

  !> Momentum update subroutine for full Lorentz dynamics
  subroutine Lorentz_kick(upart,e,b,q,m,dtp)
    use mod_global_parameters
    use mod_geometry
    double precision, intent(in)      :: e(ndir), b(ndir), q, m, dtp
    double precision, intent(inout)   :: upart(ndir)
    double precision                  :: lfac, cosphi, sinphi, phi1, phi2, r, re, sigma, td
    double precision, dimension(ndir) :: emom, uprime, tau, s, tmp, uplus, xcart1, xcart2, ucart2, radmom
    double precision, dimension(ndir) :: upartk, vbar, Fk, C1, C2, dupartk
    double precision                  :: abserr, tol, lfack, J11, J12, J13, J21, J22, J23, J31, J32, J33
    double precision                  :: iJ11, iJ12, iJ13, iJ21, iJ22, iJ23, iJ31, iJ32, iJ33, Det
    integer                           :: nk, nkmax

    ! Perform momentum update based on the chosen integrator
    select case(integrator)

    ! Boris integrator (works in Cartesian and cylindrical)
    case(Boris)
      ! Momentum update
      emom = q * e * dtp /(2.0d0 * m)
      uprime = upart + emom
      !!!!!! TODO: Adjust and document losses
!      if (losses) then
!        call get_lfac(particle(ipart)%self%u,lfac)
!        re = abs(q)**2 / (m * const_c**2)
!        call cross(upart,b,tmp)
!        radmom = - third * re**2 * lfac &
!             * ( sum((e(:)+tmp(:)/lfac)**2)  &
!             -  (sum(e(:)*upart(:))/lfac)**2 ) &
!             * particle(ipart)%self%u / m / const_c * dt_p
!        uprime = uprime + radmom
!      end if

      call get_lfac(uprime,lfac)
      tau = q * b * dtp / (2.0d0 * lfac * m)
      s = 2.0d0 * tau / (1.0d0+sum(tau(:)**2))

      call cross(uprime,tau,tmp)
      call cross(uprime+tmp,s,tmp)
      uplus = uprime + tmp

      upart = uplus + emom
      !!!!!! TODO: Adjust and document losses
!      if(losses) then
!        call cross(uplus,b,tmp)
!        radmom = - third * re**2 * lfac &
!             * ( sum((e(:)+tmp(:)/lfac)**2)  &
!             -  (sum(e(:)*uplus(:))/lfac)**2 ) &
!             * uplus / m / const_c * dt_p
!        upart = upart + radmom
!      end if

    ! Vay integrator
    case(Vay)
      call get_lfac(upart,lfac)
      emom = q * e * dtp /(2.0d0 * m)
      tau = q * b * dtp / (2.0d0 * m)

      call cross(upart,tau,tmp)
      uprime = upart + 2.d0*emom + tmp/lfac

      call get_lfac(uprime,lfac)
      sigma = lfac**2 - sum(tau(:)*tau(:))
      lfac = sqrt((sigma + sqrt(sigma**2 &
           + 4.d0 * (sum(tau(:)*tau(:)) + sum(uprime(:)*tau(:)/c_norm)**2))) / 2.d0)
      
      call cross(uprime,tau,tmp)
      upart = (uprime + sum(uprime(:)*tau(:))*tau/lfac**2 + tmp/lfac) / (1.d0+sum(tau(:)*tau(:))/lfac**2)

    ! Higuera-Cary integrator
    case(HC)
      call get_lfac(upart,lfac)
      emom = q * e * dtp /(2.0d0 * m)
      tau = q * b * dtp / (2.0d0 * m)
      uprime = upart + emom

      call get_lfac(uprime,lfac)
      sigma = lfac**2 - sum(tau(:)*tau(:))
      lfac = sqrt((sigma + sqrt(sigma**2 &
           + 4.d0 * (sum(tau(:)*tau(:)) + sum(uprime(:)*tau(:)/c_norm)**2))) / 2.d0)
      
      call cross(uprime,tau,tmp)
      upart = (uprime + sum(uprime(:)*tau(:))*tau/lfac**2 + tmp/lfac) / (1.d0+sum(tau(:)*tau(:))/lfac**2) &
              + emom + tmp/lfac

    ! Lapenta-Markidis integrator
    case(LM)
      ! Initialise iteration quantities
      call get_lfac(upart,lfac)
      upartk = upart

      ! START OF THE NONLINEAR CYCLE
      abserr = 1.d0
      tol=1.d-14
      nkmax=10
      nk=0
      do while(abserr > tol .and. nk < nkmax)

        nk=nk+1

        call get_lfac(upartk,lfack)
        vbar = (upart + upartk) / (lfac + lfack)
        call cross(vbar,b,tmp)

        ! Compute residual vector
        Fk = upartk - upart - q*dtp/m * (e + tmp)

        ! Compute auxiliary coefficients
        C1 = (lfack + lfac - upartk(1:ndim) / lfack / c_norm**2 * (upartk + upart)) / (lfack + lfac)**2
        C2 = - upartk / lfack / c_norm**2 / (lfack + lfac)**2

        ! Compute Jacobian
        J11 = 1. - q*dtp/m * (C2(1) * (upartk(2) + upart(2)) * b(3) - C2(1) * (upartk(3) + upart(3)) * b(2))
        J12 = - q*dtp/m * (C1(2) * b(3) - C2(2) * (upartk(3) + upart(3)) * b(2))
        J13 = - q*dtp/m * (C2(3) * (upartk(2) + upart(2)) * b(3) - C1(3) * b(2))
        J21 = - q*dtp/m * (- C1(1) * b(3) + C2(1) * (upartk(3) + upart(3)) * b(1))
        J22 = 1. - q*dtp/m * (- C2(2) * (upartk(1) + upart(1)) * b(3) + C2(2) * (upartk(3) + upart(3)) * b(1))
        J23 = - q*dtp/m * (- C2(3) * (upartk(1) + upart(1)) * b(3) + C1(3) * b(1))
        J31 = - q*dtp/m * (C1(1) * b(2) - C2(1) * (upartk(2) + upart(2)) * b(1))
        J32 = - q*dtp/m * (C2(2) * (upartk(1) + upart(1)) * b(2) - C1(2) * b(1))
        J33 = 1. - q*dtp/m * (C2(3) * (upartk(1) + upart(1)) * b(2) - C2(3) * (upartk(2) + upart(2)) * b(1))

        ! Compute inverse Jacobian
        Det = J11*J22*J33 + J21*J32*J13 + J31*J12*J23 - J11*J32*J23 - J31*J22*J13 - J21*J12*J33
        iJ11 = (J22*J33 - J23*J32) / Det
        iJ12 = (J13*J32 - J12*J33) / Det
        iJ13 = (J12*J23 - J13*J22) / Det
        iJ21 = (J23*J31 - J21*J33) / Det
        iJ22 = (J11*J33 - J13*J31) / Det
        iJ23 = (J13*J21 - J11*J23) / Det
        iJ31 = (J21*J32 - J22*J31) / Det
        iJ32 = (J12*J31 - J11*J32) / Det
        iJ33 = (J11*J22 - J12*J21) / Det

        ! Compute new upartk = upartk - J^(-1) * F(upartk)
        dupartk(1) = - (iJ11 * Fk(1) + iJ12 * Fk(2) + iJ13 * Fk(3))
        dupartk(2) = - (iJ21 * Fk(1) + iJ22 * Fk(2) + iJ23 * Fk(3))
        dupartk(3) = - (iJ31 * Fk(1) + iJ32 * Fk(2) + iJ33 * Fk(3))

        ! Check convergence
        upartk=upartk+dupartk
        abserr=sqrt(sum(dupartk(:)*dupartk(:)))

      end do
      ! END OF THE NONLINEAR CYCLE

      ! Update velocity
      upart = upartk

    end select

  end subroutine Lorentz_kick

  !> Update payload subroutine
  subroutine switch_update_payload(igrid,xpart,upart,qpart,mpart,mypayload,mynpayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,mynpayload
    double precision, intent(in)  :: xpart(1:ndir),upart(1:ndir),qpart,mpart,particle_time
    double precision, intent(out) :: mypayload(mynpayload)
    double precision, dimension(1:ndir) :: vE, e, b, bhat, vfluid, current
    double precision, dimension(1:ndir) :: drift1, drift2
    double precision, dimension(1:ndir) :: drift3, drift4, drift5, drift6, drift7
    double precision, dimension(1:ndir) :: bdotgradb, vEdotgradb, gradkappaB
    double precision, dimension(1:ndir) :: bdotgradvE, vEdotgradvE
    double precision, dimension(1:ndir) :: gradBdrift, reldrift, bdotgradbdrift
    double precision, dimension(1:ndir) :: vEdotgradbdrift, bdotgradvEdrift
    double precision, dimension(1:ndir) :: vEdotgradvEdrift
    double precision                    :: kappa, upar, absb, gamma, vpar, vEabs
    double precision                    :: gradBdrift_abs, reldrift_abs, epar
    double precision                    :: bdotgradbdrift_abs, vEdotgradbdrift_abs
    double precision                    :: bdotgradvEdrift_abs, vEdotgradvEdrift_abs
    double precision                    :: momentumpar1, momentumpar2, momentumpar3, momentumpar4
    ! used for upart transformation to GCA frame
    double precision                    :: lfac, magmom, vperp
    double precision, dimension(1:ndir) :: upart1

    call get_vec(bp, igrid,xpart(1:ndir),particle_time,b)
    call get_vec(vp, igrid,xpart(1:ndir),particle_time,vfluid)
    call get_vec(jp, igrid,xpart(1:ndir),particle_time,current)
    call get_vec(ep, igrid,xpart(1:ndir),particle_time,e)

    absb         = sqrt(sum(b(:)**2))
    bhat(1:ndir) = b(1:ndir) / absb
    epar         = sum(e(:)*bhat(:))
    call cross(e,bhat,vE)
    vE(1:ndir)   = vE(1:ndir) / absb
    vEabs = sqrt(sum(vE(:)**2))
    if (relativistic) then
      kappa = 1.d0/sqrt(1.0d0 - sum(vE(:)**2)/c_norm**2)
    else
      kappa = 1.d0
    end if

    if (mpart<0) then
      call get_lfac(upart, lfac)
      vpar = sum(upart/lfac * bhat)
      vperp = sqrt(norm2(upart/lfac)**2-vpar**2)
      upart1(1) = lfac * vpar
      upart1(2) = abs(mpart) * (vperp * lfac)**2 / (2.0d0 * absb)
      upart1(3) = lfac
    else
      upart1 = upart
    end if

    vpar = upart1(1)/upart1(3)
    upar = upart1(1)

    call get_vec(b_dot_grad_b, igrid,xpart(1:ndir),particle_time,bdotgradb)
    call get_vec(vE_dot_grad_b, igrid,xpart(1:ndir),particle_time,vEdotgradb)
    call get_vec(grad_kappa_B, igrid,xpart(1:ndir),particle_time,gradkappaB)
    call get_vec(b_dot_grad_vE, igrid,xpart(1:ndir),particle_time,bdotgradvE)
    call get_vec(vE_dot_grad_vE, igrid,xpart(1:ndir),particle_time,vEdotgradvE)

    drift1(1:ndir) = bhat(1:ndir)/(absb/kappa**2)
    drift2(1:ndir) = upart1(2)/(upart1(3)*qpart)*gradkappaB(1:ndir)

    call cross(drift1,drift2,gradBdrift)
    gradBdrift_abs = sqrt(sum(gradBdrift(:)**2))

    drift3(1:ndir) = upar*epar/upart1(3)*vE(1:ndir)
    call cross(drift1,drift3,reldrift)
    reldrift_abs = sqrt(sum(reldrift(:)**2))

    drift4(1:ndir) = abs(mpart)/qpart* ( upar**2/upart1(3)*bdotgradb(1:ndir))
    call cross(drift1,drift4,bdotgradbdrift)
    bdotgradbdrift_abs = sqrt(sum(bdotgradbdrift(:)**2))

    drift5(1:ndir) = abs(mpart)/qpart* ( upar*vEdotgradb(1:ndir))
    call cross(drift1,drift5,vEdotgradbdrift)
    vEdotgradbdrift_abs = sqrt(sum(vEdotgradbdrift(:)**2))

    drift6(1:ndir) = abs(mpart)/qpart* ( upar*bdotgradvE(1:ndir))
    call cross(drift1,drift6,bdotgradvEdrift)
    bdotgradvEdrift_abs = sqrt(sum(bdotgradvEdrift(:)**2))

    drift7(1:ndir) = abs(mpart)/qpart* (upart1(3)*vEdotgradvE(1:ndir))
    call cross(drift1,drift7,vEdotgradvEdrift)
    vEdotgradvEdrift_abs = sqrt(sum(vEdotgradvEdrift(:)**2))

    momentumpar1 = qpart/abs(mpart)*epar
    momentumpar2 = -(upart1(2)/abs(mpart)/upart1(3))*sum(bhat(:)*gradkappaB(:))
    momentumpar3 = upar*sum(vE(:)*bdotgradb(:))
    momentumpar4 = upart1(3)*sum(vE(:)*vEdotgradb(:))

    ! Payload update
    if (mynpayload > 0) then
      ! current gyroradius
      mypayload(1) = sqrt(2.0d0*abs(mpart)*upart1(2)*absb)/abs(qpart*absb)
    end if
    if (mynpayload > 1) then
      ! pitch angle
      mypayload(2) = atan2(sqrt((2.0d0*upart1(2)*absb)/(abs(mpart)*upart1(3)**2)),vpar)
    end if
    if (mynpayload > 2) then
      ! particle v_perp
      mypayload(3) = sqrt((2.0d0*upart1(2)*absb)/(abs(mpart)*upart1(3)**2))
    end if
    if (mynpayload > 3) then
      ! particle parallel momentum term 1
      mypayload(4) = momentumpar1
    end if
    if (mynpayload > 4) then
      ! particle parallel momentum term 2
      mypayload(5) = momentumpar2
    end if
    if (mynpayload > 5) then
      ! particle parallel momentum term 3
      mypayload(6) = momentumpar3
    end if
    if (mynpayload > 6) then
      ! particle parallel momentum term 4
      mypayload(7) = momentumpar4
    end if
    if (mynpayload > 7) then
      ! particle ExB drift
      mypayload(8) = vEabs
    end if
    if (mynpayload > 8) then
      ! relativistic drift
      mypayload(9) = gradBdrift_abs
    end if
    if (mynpayload > 9) then
      ! gradB drift
      mypayload(10) = reldrift_abs
    end if
    if (mynpayload > 10) then
      ! bdotgradb drift
      mypayload(11) = bdotgradbdrift_abs
    end if
    if (mynpayload > 11) then
      ! vEdotgradb drift
      mypayload(12) = vEdotgradbdrift_abs
    end if
    if (mynpayload > 12) then
      ! bdotgradvE drift
      mypayload(13) = bdotgradvEdrift_abs
    end if
    if (mynpayload > 13) then
      ! vEdotgradvE drift
      mypayload(14) = vEdotgradvEdrift_abs
    end if

  end subroutine switch_update_payload

  function gca_get_particle_dt(partp, end_time) result(dt_p)
    use mod_odeint
    use mod_global_parameters
    type(particle_ptr), intent(in) :: partp
    double precision, intent(in)   :: end_time
    double precision               :: dt_p

    double precision            :: tout, dt_particles_mype, dt_cfl0, dt_cfl1, dt_a
    double precision            :: dxmin, vp, a, gammap
    double precision            :: v(ndir), y(ndir+2),ytmp(ndir+2), dydt(ndir+2), v0(ndir), v1(ndir), dydt1(ndir+2)
    double precision            :: ap0, ap1, dt_cfl_ap0, dt_cfl_ap1, dt_cfl_ap
    double precision            :: dt_euler, dt_tmp
    ! make these particle cfl conditions more restrictive if you are interpolating out of the grid
    double precision            :: cfl, uparcfl
    double precision, parameter :: uparmin=1.0d-8*const_c
    integer                     :: ipart, iipart, nout, ic^D, igrid_particle, ipe_particle, ipe

    if (const_dt_particles > 0) then
      dt_p = const_dt_particles
      return
    end if

    cfl = particles_cfl
    uparcfl = particles_cfl

    igrid_working = partp%igrid
    ipart_working = partp%self%index
    dt_tmp = (end_time - partp%self%time)
    if(dt_tmp .le. 0.0d0) then
      dt_p = smalldouble
      return
    end if
    ! make sure we step only one cell at a time, first check CFL at current location
    ! then we make an Euler step to the new location and check the new CFL
    ! we simply take the minimum of the two timesteps.
    ! added safety factor cfl:
    dxmin  = min({rnode(rpdx^D_,igrid_working)},bigdouble)*cfl
    ! initial solution vector:
    y(1:ndir) = partp%self%x(1:ndir) ! position of guiding center
    y(ndir+1) = partp%self%u(1) ! parallel momentum component (gamma v||)
    y(ndir+2) = partp%self%u(2) ! conserved magnetic moment
    ytmp=y
    !y(ndir+3) = partp%self%u(3) ! Lorentz factor of guiding centre

    call derivs_gca(partp%self%time,y,dydt)
    v0(1:ndir) = dydt(1:ndir)
    ap0        = dydt(ndir+1)

    ! guiding center velocity:
    v(1:ndir) = abs(dydt(1:ndir))
    vp = sqrt(sum(v(:)**2))

    dt_cfl0    = 1/2 * dxmin / max(vp, smalldouble)
    dt_cfl_ap0 = uparcfl * abs(max(y(ndir+1), uparmin)) / max(abs(ap0), smalldouble)
    ! loose the constraint over the low velocity and low acc, do not let velocity change too much at one step
    !dt_cfl_ap0 = min(dt_cfl_ap0, uparcfl * sqrt(abs(unit_length*dxmin/(ap0+smalldouble))) )

    ! make an Euler step with the proposed timestep:
    ! new solution vector:
    dt_euler = min(dt_tmp,dt_cfl0,dt_cfl_ap0)
    y(1:ndir+2) = y(1:ndir+2) + dt_cfl0 * dydt(1:ndir+2)

    partp%self%x(1:ndir) = y(1:ndir) ! position of guiding center
    partp%self%u(1)      = y(ndir+1) ! parallel momentum component (gamma v||)
    partp%self%u(2)      = y(ndir+2) ! conserved magnetic moment

    ! first check if the particle is outside the physical domain or in the ghost cells
    if(.not. particle_in_igrid(ipart_working,igrid_working)) then
      y(1:ndir+2) = ytmp(1:ndir+2)
    end if

    call derivs_gca_rk(partp%self%time+dt_euler,y,dydt)
    !call derivs_gca(partp%self%time+dt_euler,y,dydt)

    v1(1:ndir) = dydt(1:ndir)
    ap1        = dydt(ndir+1)

    ! guiding center velocity:
    v(1:ndir) = abs(dydt(1:ndir))
    vp = sqrt(sum(v(:)**2))

    dt_cfl1    = dxmin / max(vp, smalldouble)
    dt_cfl_ap1 = uparcfl * abs(max(y(ndir+1), uparmin)) / max(abs(ap1), smalldouble)
    !dt_cfl_ap1 = min(dt_cfl_ap1, uparcfl * sqrt(abs(unit_length*dxmin/(ap1+smalldouble))) )

    dt_tmp = min(dt_euler, dt_cfl1, dt_cfl_ap1)

    partp%self%x(1:ndir) = ytmp(1:ndir) ! position of guiding center
    partp%self%u(1)      = ytmp(ndir+1) ! parallel momentum component (gamma v||)
    partp%self%u(2)      = ytmp(ndir+2) ! conserved magnetic moment
    !dt_tmp = min(dt_cfl1, dt_cfl_ap1)

    ! time step due to parallel acceleration:
    ! The standard thing, dt=sqrt(dx/a) where we compute a from d(gamma v||)/dt and d(gamma)/dt
    ! dt_ap = sqrt(abs(dxmin*unit_length*y(ndir+3)/( dydt(ndir+1) - y(ndir+1)/y(ndir+3)*dydt(ndir+3) ) ) )
    ! vp = sqrt(sum(v(1:ndir)**))
    ! gammap = sqrt(1.0d0/(1.0d0-(vp/const_c)**2))
    ! ap = const_c**2/vp*gammap**(-3)*dydt(ndir+3)
    ! dt_ap = sqrt(dxmin*unit_length/ap)

    !dt_a = bigdouble
    !if (dt_euler .gt. smalldouble) then
    !   a = sqrt(sum((v1(1:ndir)-v0(1:ndir))**2))/dt_euler
    !   dt_a = min(sqrt(dxmin/a),bigdouble)
    !end if

    !dt_p = min(dt_tmp , dt_a)
    dt_p = dt_tmp

    ! Make sure we do not advance beyond end_time
    call limit_dt_endtime(end_time - partp%self%time, dt_p)

  end function gca_get_particle_dt

  ! function gca_get_particle_dt(partp, end_time) result(dt_p)
  !   use mod_odeint
  !   use mod_global_parameters
  !   type(particle_ptr), intent(in) :: partp
  !   double precision, intent(in)   :: end_time
  !   double precision               :: dt_p

  !   double precision            :: tout, dt_particles_mype, dt_cfl0, dt_cfl1, dt_a
  !   double precision            :: dxmin, vp, a, gammap
  !   double precision            :: v(ndir), y(ndir+2),ytmp(ndir+2), dydt(ndir+2), v0(ndir), v1(ndir), dydt1(ndir+2)
  !   double precision            :: ap0, ap1, dt_cfl_ap0, dt_cfl_ap1, dt_cfl_ap
  !   double precision            :: dt_euler, dt_tmp
  !   ! make these particle cfl conditions more restrictive if you are interpolating out of the grid
  !   double precision            :: cfl, uparcfl
  !   double precision, parameter :: uparmin=1.0d-8*const_c
  !   integer                     :: ipart, iipart, nout, ic^D, igrid_particle, ipe_particle, ipe

  !   if (const_dt_particles > 0) then
  !     dt_p = const_dt_particles
  !     return
  !   end if

  !   cfl = particles_cfl
  !   uparcfl = particles_cfl

  !   igrid_working = partp%igrid
  !   ipart_working = partp%self%index
  !   dt_tmp = (end_time - partp%self%time)
  !   if(dt_tmp .le. 0.0d0) dt_tmp = smalldouble
  !   ! make sure we step only one cell at a time, first check CFL at current location
  !   ! then we make an Euler step to the new location and check the new CFL
  !   ! we simply take the minimum of the two timesteps.
  !   ! added safety factor cfl:
  !   dxmin  = min({rnode(rpdx^D_,igrid_working)},bigdouble)*cfl
  !   ! initial solution vector:
  !   y(1:ndir) = partp%self%x(1:ndir) ! position of guiding center
  !   y(ndir+1) = partp%self%u(1) ! parallel momentum component (gamma v||)
  !   y(ndir+2) = partp%self%u(2) ! conserved magnetic moment
  !   ytmp=y
  !   !y(ndir+3) = partp%self%u(3) ! Lorentz factor of guiding centre

  !   call derivs_gca(partp%self%time,y,dydt)
  !   v0(1:ndir) = dydt(1:ndir)
  !   ap0        = dydt(ndir+1)

  !   ! guiding center velocity:
  !   v(1:ndir) = abs(dydt(1:ndir))
  !   vp = sqrt(sum(v(:)**2))

  !   dt_cfl0    = dxmin / max(vp, smalldouble)
  !   dt_cfl_ap0 = uparcfl * abs(max(abs(y(ndir+1)),uparmin) / max(abs(ap0), smalldouble))
  !   !dt_cfl_ap0 = min(dt_cfl_ap0, uparcfl * sqrt(abs(unit_length*dxmin/(ap0+smalldouble))) )

  !   ! make an Euler step with the proposed timestep:
  !   ! new solution vector:
  !   dt_euler = min(dt_tmp,dt_cfl0,dt_cfl_ap0)
  !   y(1:ndir+2) = y(1:ndir+2) + dt_euler * dydt(1:ndir+2)

  !   partp%self%x(1:ndir) = y(1:ndir) ! position of guiding center
  !   partp%self%u(1)      = y(ndir+1) ! parallel momentum component (gamma v||)
  !   partp%self%u(2)      = y(ndir+2) ! conserved magnetic moment

  !   ! first check if the particle is outside the physical domain or in the ghost cells
  !   if(.not. particle_in_igrid(ipart_working,igrid_working)) then
  !     y(1:ndir+2) = ytmp(1:ndir+2)
  !   end if

  !   call derivs_gca_rk(partp%self%time+dt_euler,y,dydt)
  !   !call derivs_gca(partp%self%time+dt_euler,y,dydt)

  !   v1(1:ndir) = dydt(1:ndir)
  !   ap1        = dydt(ndir+1)

  !   ! guiding center velocity:
  !   v(1:ndir) = abs(dydt(1:ndir))
  !   vp = sqrt(sum(v(:)**2))

  !   dt_cfl1    = dxmin / max(vp, smalldouble)
  !   dt_cfl_ap1 = uparcfl * abs(max(abs(y(ndir+1)),uparmin) / max(abs(ap1), smalldouble))
  !   !dt_cfl_ap1 = min(dt_cfl_ap1, uparcfl * sqrt(abs(unit_length*dxmin/(ap1+smalldouble))) )

  !   dt_tmp = min(dt_euler, dt_cfl1, dt_cfl_ap1)

  !   partp%self%x(1:ndir) = ytmp(1:ndir) ! position of guiding center
  !   partp%self%u(1)      = ytmp(ndir+1) ! parallel momentum component (gamma v||)
  !   partp%self%u(2)      = ytmp(ndir+2) ! conserved magnetic moment
  !   !dt_tmp = min(dt_cfl1, dt_cfl_ap1)

  !   ! time step due to parallel acceleration:
  !   ! The standard thing, dt=sqrt(dx/a) where we compute a from d(gamma v||)/dt and d(gamma)/dt
  !   ! dt_ap = sqrt(abs(dxmin*unit_length*y(ndir+3)/( dydt(ndir+1) - y(ndir+1)/y(ndir+3)*dydt(ndir+3) ) ) )
  !   ! vp = sqrt(sum(v(1:ndir)**))
  !   ! gammap = sqrt(1.0d0/(1.0d0-(vp/const_c)**2))
  !   ! ap = const_c**2/vp*gammap**(-3)*dydt(ndir+3)
  !   ! dt_ap = sqrt(dxmin*unit_length/ap)

  !   !dt_a = bigdouble
  !   !if (dt_euler .gt. smalldouble) then
  !   !   a = sqrt(sum((v1(1:ndir)-v0(1:ndir))**2))/dt_euler
  !   !   dt_a = min(sqrt(dxmin/a),bigdouble)
  !   !end if

  !   !dt_p = min(dt_tmp , dt_a)
  !   dt_p = dt_tmp

  !   ! Make sure we do not advance beyond end_time
  !   call limit_dt_endtime(end_time - partp%self%time, dt_p)

  ! end function gca_get_particle_dt
  
  function Lorentz_get_particle_dt(partp, end_time) result(dt_p)
    use mod_global_parameters
    use mod_geometry
    type(particle_ptr), intent(in)   :: partp
    double precision, intent(in)     :: end_time
    double precision                 :: dt_p
    integer                          :: ipart, iipart, nout
    double precision,dimension(ndir) :: b,v
    double precision                 :: lfac,absb,dt_cfl
    double precision                 :: tout
    double precision, parameter      :: cfl = 10.d0

    if (const_dt_particles > 0) then
      dt_p = const_dt_particles
      return
    end if

    call get_vec(bp, partp%igrid,partp%self%x,partp%self%time,b)
    absb = sqrt(sum(b(:)**2))
    call get_lfac(partp%self%u,lfac)

    ! CFL timestep
    ! make sure we step only one cell at a time:
    v(:) = abs(partp%self%u(:) / lfac)

!    ! convert to angular velocity:
!    if(coordinate ==cylindrical.and.phi_>0) then
!      v(phi_) = abs(v(phi_)/partp%self%x(r_))
!    end if

    dt_cfl = min(bigdouble, &
         {rnode(rpdx^D_,partp%igrid)/max(v(^D), smalldouble)})

    if(coordinate ==cylindrical.and.phi_>0) then
      ! phi-momentum leads to radial velocity:
      if(phi_ .gt. ndim) dt_cfl = min(dt_cfl, &
           sqrt(rnode(rpdx1_,partp%igrid)/partp%self%x(r_)) &
           / v(phi_))
      ! limit the delta phi of the orbit (just for aesthetic reasons):
      dt_cfl = min(dt_cfl,0.1d0/max(v(phi_), smalldouble))
      ! take some care at the axis:
      dt_cfl = min(dt_cfl,(partp%self%x(r_)+smalldouble)/max(v(r_), smalldouble))
    end if

    dt_cfl = dt_cfl * cfl

    ! bound by gyro-rotation:
    dt_p = abs( dtheta * partp%self%m * lfac &
         / (partp%self%q * absb) ) * 20

    dt_p = min(dt_p, dt_cfl)
    ! Make sure we don't advance beyond end_time
    call limit_dt_endtime(end_time - partp%self%time, dt_p)

  end function Lorentz_get_particle_dt

end module mod_particle_switch
