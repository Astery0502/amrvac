!> MODULE with Functions which are used in mod_usr.t mainly about TDM and PARTICLE
!> 
!> Each function can initialize a procedure pointer in mod_usr.t
module mod_usr_funcs
    use mod_mhd
    use mod_particles

    implicit none
    public 

    logical :: filter_on
    logical :: uniform_pitch_angle
    logical :: nonuniform_e
    double precision :: vthermal
    double precision :: xparmin^D, xparmax^D
    double precision :: mass_e
    double precision :: mass_p
    double precision :: charge
    double precision :: field_cr
    double precision :: check_cr
    integer :: rhop_, pp_

contains

    subroutine usr_funcs_init
        use mod_global_parameters
        ! default value of parameters
        filter_on = .False.
        uniform_pitch_angle = .False.
        nonuniform_e = .False.
        vthermal = 4.0e8
        mass_e = 9.10956d-28/unit_mass
        mass_p = 1.67262d-24/unit_mass
        charge = 4.8032d-10/unit_charge
        field_cr = 0.0d0
        check_cr = 0.0d0
        call usr_funcs_read(par_files)
    end subroutine usr_funcs_init

    subroutine usr_funcs_read(files)
        use mod_global_parameters, only: unitpar
        character(len=*), intent(in) :: files(:)
        integer                      :: n

        namelist /ufuncs_list/ filter_on, uniform_pitch_angle, nonuniform_e, &
                                xparmin^D, xparmax^D, vthermal, check_cr, field_cr
        do n = 1, size(files)
        open(unitpar, file=trim(files(n)), status="old")
        read(unitpar, ufuncs_list, end=111)
    111   close(unitpar)
        end do

    end subroutine usr_funcs_read

    subroutine usrp_define_gridvars(ngridvars)
        use mod_global_parameters
        integer, intent(inout) :: ngridvars
        integer                :: idir
        ngridvars = ngridvars + 1
        rhop_ = ngridvars
        ngridvars = ngridvars + 1
        pp_ = ngridvars
    end subroutine usrp_define_gridvars

    subroutine usrp_fill_gridvars
        use mod_global_parameters
        use mod_usr_methods, only: usr_particle_fields
        use mod_geometry

        integer                                   :: igrid, iigrid, idir, idim
        double precision, dimension(ixG^T,1:nw)   :: w
        double precision, dimension(ixG^T)        :: pth

        do iigrid=1,igridstail; igrid=igrids(iigrid)
            ! fill with density and temperature
            w(ixG^T,1:nw) = ps(igrid)%w(ixG^T,1:nw)
            !call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
            call mhd_get_pthermal(w,ps(igrid)%x,ixG^LL,ixG^LL,pth)
            gridvars(igrid)%w(ixG^T,rhop_) = w(ixG^T,rho_)
            gridvars(igrid)%w(ixG^T,pp_) = pth(ixG^T)
            if (time_advance) then
                ! fill with old ones 
                w(ixG^T,1:nw) = pso(igrid)%w(ixG^T,1:nw)
                !call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
                call mhd_get_pthermal(w,ps(igrid)%x,ixG^LL,ixG^LL,pth)
                gridvars(igrid)%wold(ixG^T,rhop_) = w(ixG^T,rho_)
                gridvars(igrid)%wold(ixG^T,pp_) = pth(ixG^T)
            end if
        end do
    end subroutine usrp_fill_gridvars

    subroutine usrp_modification(igrid, x, v, q, m, follow, check)  
        use mod_global_parameters
        integer, intent(in)             :: igrid
        double precision, intent(inout) :: x(1:ndir)
        double precision, intent(inout) :: v(1:ndir), q, m
        logical, intent(inout)          :: follow
        logical, intent(out)            :: check

        double precision                :: current(1:ndir), mag(1:ndir), bhat(1:ndir), xaxis(1:ndir), yaxis(1:ndir), umhdi(1:ndir)
        double precision                :: phi, theta, valcr, rhoi, pi
        double precision                :: ue(ndir), up(ndir)
        
        call get_vec(jp, igrid, x(:), global_time, current)
        call get_vec(bp, igrid, x(:), global_time, mag)
        call get_scalar(rhop_, igrid, x(:), global_time, rhoi) 

        if (filter_on) then 
            check = .false.
            if ((norm2(current)/norm2(mag))>(check_cr)) then
            check = .true.
            end if
        end if

        ! reshape the pitch angle distribution to uniform
        if (check .and. uniform_pitch_angle) then 
            phi = rng%unif_01()
            theta = rng%unif_01()
            phi = phi * dpi
            theta = theta*2*dpi
            ! here we assume the zero-phase randomly (not uniform)
            bhat = mag / norm2(mag)
            call random_perpendicular_vector(bhat, xaxis)
            call cross_3d(bhat, xaxis, yaxis)
            yaxis = yaxis/norm2(yaxis)
            v = norm2(v)*(sin(phi)*cos(theta)*xaxis+sin(phi)*sin(theta)*yaxis+cos(phi)*bhat)
        end if 

    end subroutine usrp_modification

    subroutine usrp_custom_field(w,x,E,B,J)
        use mod_global_parameters
        double precision, intent(in)  :: w(ixG^T,1:nw) 
        double precision, intent(in)  :: x(ixG^T,1:ndim) 
        double precision, intent(out) :: E(ixG^T,ndir) 
        double precision, intent(out) :: B(ixG^T,ndir) 
        double precision, intent(out) :: J(ixG^T,ndir) 

        integer                       :: idirmin, idir, ix^D 
        double precision              :: peta
        double precision              :: Btotal(ixG^T), jtotal(ixG^T)
        double precision              :: current(ixG^T,7-2*ndir:3)
        double precision              :: vcr(ixG^T), w_mhd(ixG^T,nw)

        w_mhd(ixG^T,1:nw) = w(ixG^T,1:nw)
        call phys_to_primitive(ixG^LL,ixG^LL,w_mhd,x)

        ! fill with magnetic field:
        if (B0field) then
            call mpistop("user field does not consider b0field")
        else
            B(ixG^T,:) = w_mhd(ixG^T,iw_mag(:))
        end if
        ! fill with current
        current = zero
        call particle_get_current(w_mhd,ixG^LL,ixG^LLIM^D^LSUB1,idirmin,current)
        J(ixG^T,:) = current(ixG^T,:)
        
        E(ixG^T,1) = B(ixG^T,2) * w_mhd(ixG^T,iw_mom(3)) - B(ixG^T,3) * w_mhd(ixG^T,iw_mom(2))
        E(ixG^T,2) = B(ixG^T,3) * w_mhd(ixG^T,iw_mom(1)) - B(ixG^T,1) * w_mhd(ixG^T,iw_mom(3))
        E(ixG^T,3) = B(ixG^T,1) * w_mhd(ixG^T,iw_mom(2)) - B(ixG^T,2) * w_mhd(ixG^T,iw_mom(1))

        Btotal(ixG^T) = sqrt(w_mhd(ixG^T,mag(1))**2+w_mhd(ixG^T,mag(2))**2+w_mhd(ixG^T,mag(3))**2)
        jtotal(ixG^T) = sqrt(current(ixG^T,1)**2+current(ixG^T,2)**2+current(ixG^T,3)**2)
        if (nonuniform_e) then
            where(jtotal(ixG^T)/Btotal(ixG^T)>field_cr) 
            E(ixG^T,1) = E(ixG^T,1)+particles_eta*current(ixG^T,1)
            E(ixG^T,2) = E(ixG^T,2)+particles_eta*current(ixG^T,2)
            E(ixG^T,3) = E(ixG^T,3)+particles_eta*current(ixG^T,3)
            end where
        else
            E(ixG^T,:) = E(ixG^T,:) + particles_eta * current(ixG^T,:)
        end if
    end subroutine usrp_custom_field

    subroutine usrp_custom_e(igrid, x, particle_time, e)
        integer, intent(in)            :: igrid
        double precision, intent(in)   :: x(1:3)
        double precision, intent(in)   :: particle_time
        double precision, intent(out)  :: e(1:3)

        double precision, dimension(ndir) :: b, vfluid, current

        call get_vec(bp, igrid, x, particle_time, b)
        call get_vec(vp, igrid, x, particle_time, vfluid)
        call get_vec(jp, igrid, x, particle_time, current)
        if ((norm2(current)/norm2(b))>field_cr) then
            e(1) = -vfluid(2)*b(3)+vfluid(3)*b(2) + particles_eta*current(1)
            e(2) = vfluid(1)*b(3)-vfluid(3)*b(1) + particles_eta*current(2)
            e(3) = -vfluid(1)*b(2)+vfluid(2)*b(1) + particles_eta*current(3)
        else
            e(1) = -vfluid(2)*b(3)+vfluid(3)*b(2)
            e(2) = vfluid(1)*b(3)-vfluid(3)*b(1)
            e(3) = -vfluid(1)*b(2)+vfluid(2)*b(1)
        end if
    end subroutine usrp_custom_e

    subroutine random_perpendicular_vector(vec, vec_perp)
        double precision, intent(in)  :: vec(1:3)
        double precision, intent(out) :: vec_perp(1:3)
        
        double precision              :: phi, theta, costheta, vechat(1:3), xtmp(1:3)
        ! normalise 
        vechat = vec / norm2(vec)
        phi = rng%unif_01()
        costheta = rng%unif_01()
        phi = phi*2*dpi
        theta = acos(costheta)
        xtmp(1) = sin(theta)*cos(phi)
        xtmp(2) = sin(theta)*sin(phi)
        xtmp(3) = cos(theta)
        call cross_3d(vechat, xtmp, vec_perp)
        vec_perp = vec_perp/norm2(vec_perp)
    end subroutine random_perpendicular_vector

    subroutine get_particle_velocity_base(n_particles, v)
        integer, intent(in) :: n_particles
        double precision, intent(inout) :: v(3, n_particles)

        integer :: itr
        do itr=1, n_particles
            call get_maxwellian_velocity(v(:,itr),vthermal)
        end do
    end subroutine get_particle_velocity_base

    subroutine get_random_position_base(n_particles, x)
        integer, intent(in) :: n_particles
        double precision, intent(inout) :: x(3, n_particles)

        double precision :: xtmp(3)
        integer :: itr
        do itr=1, n_particles
            call get_uniform_pos1(xtmp)
            x(:,itr) = xtmp
        end do
    end subroutine get_random_position_base

    subroutine get_uniform_pos1(x)
        use mod_global_parameters
        use mod_random
        double precision, intent(out) :: x(3)

        call rng%unif_01_vec(x)
        {^D&x(^D) = xparmin^D + x(^D) * (xparmax^D - xparmin^D)\}
        x(ndim+1:) = 0.0d0
    end subroutine get_uniform_pos1

    subroutine quasi_random_particles(n_particles, x)
        integer, intent(in) :: n_particles
        double precision, intent(inout) :: x(3, n_particles)

        ! the number density of particles per unit volume
        double precision :: lx, ly, lz, ix, iy, iz
        integer :: ip, nx, ny, nz, nstack
        double precision :: np 
        double precision, dimension(:,:), allocatable :: xstack
        
        lx = xparmax1-xparmin1
        ly = xparmax2-xparmin2
        lz = xparmax3-xparmin3
        np = (n_particles / (lx * ly * lz))**(1./3.)
        nx = ceiling(lx * np)
        ny = ceiling(ly * np)
        nz = ceiling(lz * np)
        nstack = nx * ny * nz
        allocate(xstack(3, nstack))

        do ip=1, nstack 
            ! ip-1 to confirm the last point in the sequence locate in the same layer
            iz = ceiling(real(ip)/(nx*ny))-1.0 ! to maintain the consistency of all points in xy plane (the first and the last)
            iy = ceiling((ip-iz*nx*ny)/nx)-1.0 ! to avoid the mod(N,N)==0 situation
            ix = ceiling(ip-iz*nx*ny-iy*nx)-1.0
            xstack(3, ip) = (iz+1./2)*lz/nz+xparmin3
            xstack(2, ip) = (iy+1./2)*ly/ny+xparmin2
            xstack(1, ip) = (ix+1./2)*lx/nx+xparmin1
        end do
        x(1:3,1:n_particles) = xstack(1:3,1:n_particles)
    end subroutine quasi_random_particles

    ! three-dimensional vector cross
    subroutine cross_3d(vec1, vec2, vec_cross)
        double precision, intent(in) :: vec1(3), vec2(3)
        double precision, intent(inout) :: vec_cross(3)

        vec_cross(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
        vec_cross(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
        vec_cross(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)
    end subroutine cross_3d

    subroutine get_gcau_from_velocity(b, vm, ugca)
        use mod_global_parameters
        double precision, intent(in) :: b(ndir), vm(ndir)
        double precision, intent(out) :: ugca(ndir)
        double precision :: bnorm, vnorm, vpar, vperp, lfac, magmom

        if (relativistic) then
            lfac = 1.0d0 / sqrt(1.0d0 - sum(vm(:)**2)/c_norm**2)
        else
            lfac = 1.0d0
        end if
        bnorm = norm2(b(:))
        vnorm = norm2(vm(:))
        vpar  = sum(vm * b/bnorm)
        vperp = sqrt(vnorm**2 - vpar**2)
        ! parallel momentum component (gamma v||)
        ugca(1) = lfac * vpar
        ! Mr: the conserved magnetic moment
        magmom = mass_e * (vperp * lfac)**2 / (2.0d0 * bnorm)
        ugca(2) = magmom
        ! Lorentz factor
        ugca(3) = lfac
    end subroutine get_gcau_from_velocity

    subroutine ladder_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
        use mod_global_parameters

        integer, intent(in) :: igrid, level, ixI^L, ixO^L
        double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
        integer, intent(inout) :: refine, coarsen
        double precision :: rho(ixI^S), current(ixI^S,1:ndir), btotal(ixI^S,1:ndir)
        double precision :: gradb2(ixI^S,1:ndir),curlb(ixI^S,1:ndir)
        double precision :: b2(ixI^S), tmp(ixI^S), absj(ixO^S), absb(ixO^S)
        integer :: idirmin, idim, idir
        curlb = 0.d0

        ! fix the bottom layer to the basic highest level
        if (level<refine_max_level-2 .and. any(x(ixO^S,3)<=xprobmin3+0.2d0)) then
            refine=1
            coarsen=-1
        else 
            btotal = w(ixI^S,mag(:))
            rho = w(ixI^S,rho_)
            call curlvector(btotal,ixI^L,ixO^L,current,idirmin,1,ndir)
            b2 = sum(btotal(ixI^S,:)**2,dim=ndim+1)
            absj(ixO^S) = dsqrt(sum(current(ixO^S,:)**2,dim=ndim+1))
            absb(ixO^S) = dsqrt(b2(ixO^S))
            do idir=1,ndir
                call gradient(b2(ixI^S)/2,ixI^L,ixO^L,idir,gradb2(ixI^S,idir))
                gradb2(ixI^S,idir) = gradb2(ixI^S,idir) / rho(ixI^S)
                do idim=1,ndim
                    call gradient(btotal(ixI^S,idir),ixI^L,ixO^L,idim,tmp)
                    curlb(ixI^S,idir) = curlb(ixI^S,idir)+tmp(ixI^S)*btotal(ixI^S,idim)/rho(ixI^S)
                end do
            end do
            if (any((absj/absb)>field_cr) .or. any(gradb2>field_cr*25) .or. any(curlb>field_cr*25)) then
                refine = 1
                coarsen = -1
            else if (level<refine_max_level-2) then
                refine = 0
                coarsen = 0
            else 
                refine = -1
                coarsen = 0
            end if
        end if
    end subroutine ladder_refine_grid
end module mod_usr_funcs
