subroutine simulation(root_dirname, dirname)
  use mt19937
  use fposition_mod
  use fmoment_mod
  use fnoise_mod
  use distance_mod
  use repulsion_mod
  use reflect_mod
  use linkingforce_mod
  use collision_mod
  use electrophoresis_mod
  implicit none

  !1 variable declaration
  !1-1 name
  character(100),intent(in) :: root_dirname, dirname
  character(200) :: filename
  character(60) :: dammy
  character(1) :: slash = '\'
  character(60) :: bondpath = 'C:\Users\Awano\simulation\12A\ishi\movecode\bond.txt'

  !1-2 int
  integer :: i, j, n, nn

  integer :: seed
  integer :: counter = 0 !for debug

  !1-3 parameter.csv
  integer :: f_noise, f_lj, f_wlj, f_lnk, f_clsn, f_elp, f_moment0 , f_mag_bet_dipole, f_mag_bet_ext , bondy
  integer :: f_initial_placement, f_cdview3

  !1-4 parameter
  double precision , parameter :: dt = 1.0d-3 !s delta_t
  integer , parameter :: output_interval = 500 !times output_interval
  integer , parameter :: output_times = 1000 !times output_times

  !1-5 physical constant
  double precision , parameter :: g = 9.80665d0 !m*s^-2 gravitational_acceleration
  double precision , parameter :: pi = dacos(-1.0d0)
  double precision , parameter :: k = 0.01d0!1.38d-23 !J*K-1 (kg*m^2*s^-2*K^-1) Boltzmann_constant

  !1-6 beas chain condition
  integer , parameter :: total_particle_num = 400 !total_particle_num
  double precision , parameter :: m = 1.0d0!1.1d-24 * 17000 !kg mass
  double precision , parameter :: D = 1.0d0!1.0d-6 !m diameter
  double precision , parameter :: sigma = 1.0d0/(2.0d0**(1.0d0/6.0d0))!1.0d-6/(2.0d0**(1.0d0/6.0d0)) !m Collision_diameter shototsuchokkei(=0.891)
  double precision , parameter :: epsilon = 1.0d0 !J (kg*m^2*s^-2) parameter_of_Lennard-Jones_potential
  double precision , parameter :: Cd = 1.0d0!0.5d0 !resistance coefficient depend on particle shape
  double precision , parameter :: phi_max = dacos(-1.0d0)/2.0d0 !rad saidaimagarikaku(=0.523)
  double precision , parameter :: rod = 6.0d0!4.9d-6 !m distance_between_particles ryuusilamkyori

  !1-7 extarnal condition
  double precision , parameter :: wall_x = 11.0d0 !m wall_x
  double precision , parameter :: wall_y = 11.0d0 !m wall_y
  double precision , parameter :: wall_z = 3.05d0 !m wall_z
  double precision , parameter :: temp = 3.0d2 !K temperature
  double precision , parameter :: eta = 0.0d0!8.9d-4 !Pa*s (kg*m^-1*s^-1) viscosity nendo
  double precision , parameter :: rho = 1.0d1!1.0d-1!9.97d2 !kg*m^-3 density mitudo
  double precision , parameter :: qE = 5.0d2!0.5*1.6d-19*17000 * 1.0d7 !C * V*m^-1
  double precision , parameter :: poreR = 0.75d0!50d-9 !m pore_radius
  double precision , parameter :: memT = 20.0d-3!20.0d-9 !m membrane_thickness
  double precision , parameter :: gamma = 40.0d0!3.0d0*pi*D*eta/m !s^-1 coefficient_for_noise nenseiteikoukeisuu

  !1-8 variable
  double precision :: v !m*s^-1 verocity
  double precision , dimension(3 , total_particle_num) :: pos !m position iti
  double precision , dimension(3 , total_particle_num) :: moment !m*s^-1 verocity kakuzikunosokudo
  double precision , dimension(3 , total_particle_num) :: pre_moment !m*s^-1
  double precision , dimension(3 , total_particle_num) :: force !N (kg*m*s^-2)
  double precision , dimension(3 , total_particle_num) :: noise !N (kg*m*s^-2) ramdom_force
  double precision , dimension(3 , total_particle_num) :: lj !N (kg*m*s^-2) Lennard-Jones_potential
  double precision , dimension(3 , total_particle_num) :: wlj !N (kg*m*s^-2) Lennard-Jones_potential_withwall
  double precision , dimension(3 , total_particle_num) :: lnk !N (kg*m*s^-2) linkingforce
  double precision , dimension(3 , total_particle_num) :: clsn !N (kg*m*s^-2) collision with pore
  double precision , dimension(3 , total_particle_num) :: elp !N (kg*m*s^-2) electrophoresis
  double precision , dimension(3 , total_particle_num) :: mag_force
  double precision , dimension(total_particle_num,total_particle_num,4) :: distances !m 4:rx, ry, rz, r_2

  !1-9 magnet
  double precision , parameter :: mag_permeability = 0.10d0  !magnetic permeability
  double precision , parameter :: H = 100.0d0 !external magnetic fields
  double precision , parameter :: r_H = 20.0d0 !imaginally length from particle to external magnet
  double precision , parameter :: magnetization_coef = 50.0d0 !magnetization coefficient of particle
  double precision , parameter :: theta_H = 0.0d0 !angle of dipole with z
  double precision , parameter :: phi_H = 0.0d0 !angle of dipole with x
  double precision , parameter :: theta_mu = 0.0d0 !angle of dipole with z
  double precision , parameter :: phi_mu = 0.0d0 !angle of dipole with x
  double precision , parameter :: mu1 = 1.5d0!strength of dipole vector
  double precision , parameter :: mu2 = 1.0d0 !strength of dipole vector
  double precision , parameter :: r_th = 10.0d0 ! threshold distance of interaction
  double precision , parameter :: xg = 0.0d0! x center of mass of system
  double precision , parameter :: yg = 0.0d0! y center of mass of system
  double precision , parameter :: frequency = 30.0d0
  double precision , parameter :: m2pm1 = 0.0d0 !

  double precision :: mu1_v(3) 
  double precision :: mu2_v(3) 
  double precision :: r_H_v(3) 

  !initialization
  v  = 0.0d0
  pos = 0.0d0
  moment = 0.0d0
  pre_moment = 0.0d0
  force = 0.0d0
  noise = 0.0d0
  lj = 0.0d0
  wlj = 0.0d0
  lnk = 0.0d0
  clsn = 0.0d0
  elp = 0.0d0
  distances = 0.0d0

  mu1_v = 0
  mu1_v(1) = mu1 * sin(theta_mu) * cos(phi_mu)
  mu1_v(2) = mu1 * sin(theta_mu) * sin(phi_mu)
  mu1_v(3) = mu1 * cos(theta_mu)


  r_H_v(1) = sin(theta_H) * cos(phi_H)
  r_H_v(2) = sin(theta_H) * sin(phi_H)
  r_H_v(3) = cos(theta_H)
  ! mu2_v = 0
  ! mu2_v(1) = mu2 * sin(theta_mu) * cos(phi_mu)
  ! mu2_v(2) = mu2 * sin(theta_mu) * sin(phi_mu)
  ! mu2_v(3) = mu2 * cos(theta_mu)

  call system_clock(count=seed)
  call sgrnd(seed)

  !read parameters(1-3)
  open(10, file = 'parameter.csv')
    read(10, *) dammy
    read(10, *) f_noise
    read(10, *) f_lj
    read(10, *) f_wlj
    read(10, *) f_lnk
    read(10, *) f_clsn
    read(10, *) f_elp
    read(10, *) f_moment0
    read(10, *) f_initial_placement
    read(10, *) f_cdview3
    read(10, *) f_mag_bet_dipole
    read(10, *) f_mag_bet_ext
    read(10, *) bondy
  close(10)

  !2 outputs
  !2-1 condition
  filename = ''
  write(filename, '(a, a, "condition.csv")') trim(root_dirname), slash
  open(21,file=filename)
    write(21, *)'dt = ', dt, 'sec'
    write(21, *)'output_interval = ', output_interval, 'times'
    write(21, *)'output_times = ', output_times, 'times'
    write(21, *)'total_particle_num = ', total_particle_num, 'beads'
    write(21, *)'mass = ', m, 'kg'
    write(21, *)'diameter = ', D, 'm'
    write(21, *)'sigma = ', sigma, 'm'
    write(21, *)'epsilon = ', epsilon, 'J'
    write(21, *)'Cd = ', Cd, '-'
    write(21, *)'phi_max = ', phi_max, 'rad'
    write(21, *)'rod = ', rod, 'm'
    write(21, *)'wall_x = ', wall_x, 'm'
    write(21, *)'wall_y = ', wall_y, 'm'
    write(21, *)'wall_z = ', wall_z, 'm'
    write(21, *)'temperature = ', temp, 'K'
    write(21, *)'eta = ', eta, 'Pa*s'
    write(21, *)'rho = ', rho, 'kg*m^-3'
    write(21, *)'gamma = ', gamma, 's^-1'
    write(21, *)'qE = ', qE, 'C * V*m^-1'
    write(21, *)'pore_radius = ', poreR, 'm'
    write(21, *)'membrane_thickness = ', memT, 'm'
    write(21, *)'f_noise = ', f_noise, ''
    write(21, *)'f_lj = ', f_lj, ''
    write(21, *)'f_wlj = ', f_wlj, ''
    write(21, *)'f_lnk = ', f_lnk, ''
    write(21, *)'f_clsn = ', f_clsn, ''
    write(21, *)'f_elp = ', f_elp, ''
    write(21, *)'f_moment0 = ', f_moment0, ''
    write(21, *)'f_initial_placement = ', f_initial_placement, ''
  close(21)

  !2-2 bond
  open(12, file = bondpath)
  close(12,status = 'delete')
  open(12, file = bondpath)
    !beas chain
  if(bondy == 1) then 
    do i = 1 , total_particle_num - 1
      write(12,'(i3," ",i3," 0")')i, i + 1
    end do
  end if
    !square
    write(12,'(i3," ",i3," 0")')total_particle_num + 1, total_particle_num +2
    write(12,'(i3," ",i3," 0")')total_particle_num + 2, total_particle_num +3
    write(12,'(i3," ",i3," 0")')total_particle_num + 3, total_particle_num +4
    write(12,'(i3," ",i3," 0")')total_particle_num + 4, total_particle_num +1
  close(12)

  !3 initial setup
  !3-1 placement
  call fposition(wall_x,wall_y,wall_z&
  ,D,total_particle_num,phi_max,rod,pos,f_initial_placement)

  !3-2 force
  !distance0
  if (f_lj == 1 .or. f_lnk == 1)then
    distances = 0.0d0
    do i=1,total_particle_num-1
      do j=i+1,total_particle_num
        call distance(total_particle_num,i,j,pos,distances)
      end do
    end do
  endif
  !noise0
  if (f_noise == 1) then
    noise = 0.0d0
    call system_clock(count=seed)
    call sgrnd(seed)
    do i=1,total_particle_num
      call fnoise(total_particle_num,i,gamma,k,temp,m,dt,noise,seed)
    end do
  end if
  !repulsion0
  if (f_lj == 1)then
    lj = 0.0d0
    do i=1,total_particle_num-1
      do j=i+1,total_particle_num
        call repulsion(total_particle_num,i,j,sigma,&
        epsilon,lj,distances)
      end do
    end do
  end if
  !relrect0
  if (f_wlj == 1)then
    wlj = 0.0d0
    do i = 1, total_particle_num
      call reflect(total_particle_num,i,wall_x,wall_y,&
      wall_z,sigma,epsilon,pos,wlj,counter)
    end do
  end if
  !linkingforce0
  if (f_lnk == 1)then
    lnk = 0.0d0
    do i = 1, total_particle_num-1
      call linkingforce(total_particle_num,i,rod,lnk,distances)
    end do
  end if

  !collision0
  if (f_clsn == 1)then
    clsn = 0.0d0
    do i = 1, total_particle_num
      call collision(total_particle_num,i,sigma,epsilon,poreR,memT,pos,clsn,counter)
    end do
  end if

  !electrophoresis0
  if (f_elp == 1)then
    elp = 0.0d0
    do i = 1, total_particle_num
      call electrophoresis(total_particle_num,i,qE,poreR,memT,pos,elp)
    end do
  end if

  !moment0 v(0)
  if(f_moment0 == 1)then
    moment = 0.0d0
    do i=1,total_particle_num
      call fmoment(total_particle_num,i,k,temp,m,moment)
    end do
  end if

  !F(0)
  do i=1,total_particle_num
    v = sqrt(moment(1,i)**2.0d0 + moment(2,i)**2.0d0 + moment(3,i)**2.0d0)
    force(1,i) = - rho*Cd*pi*((D/2.0d0)**2.0d0)*v*moment(1,i)/2.0d0 - 3.0d0*pi*D*eta*moment(1,i) &
    + noise(1,i) + lj(1,i) + wlj(1,i) + lnk(1,i) + clsn(1,i) + elp(1,i)
    force(2,i) = - rho*Cd*pi*((D/2.0d0)**2.0d0)*v*moment(2,i)/2.0d0 - 3.0d0*pi*D*eta*moment(2,i) &
    + noise(2,i) + lj(2,i) + wlj(2,i) + lnk(2,i) + clsn(2,i) + elp(2,i)
    force(3,i) = - rho*Cd*pi*((D/2.0d0)**2.0d0)*v*moment(3,i)/2.0d0 - 3.0d0*pi*D*eta*moment(3,i) &
    + noise(3,i) + lj(3,i) + wlj(3,i) + lnk(3,i) + clsn(3,i) + elp(3,i)
  end do

  !3-3 outputs
  filename = ''
  write(filename, '(a, a, "particle_colloid_00000.cdv")') trim(dirname), slash
  filename = trim(filename)

  open(21,file=filename)
    !command
    write(21,'("#bond_file=", a)') bondpath
    write(21,*)'#r8=0.1' !radius of marker baes (pore and square)
    !beas chain
    do i = 1 , total_particle_num
      write(21,'(i4,i3,3e14.4)')i , 0 , pos(1, i) , pos(2, i) , pos(3, i)
    end do
    !marker baes (square)
    write(21,'(i4,i3,3e14.4)')total_particle_num + 1 , 8 , wall_x , -wall_y , 0.0d0
    write(21,'(i4,i3,3e14.4)')total_particle_num + 2 , 8 , wall_x , wall_y , 0.0d0
    write(21,'(i4,i3,3e14.4)')total_particle_num + 3 , 8 , -wall_x , wall_y , 0.0d0
    write(21,'(i4,i3,3e14.4)')total_particle_num + 4 , 8 , -wall_x , -wall_y , 0.0d0
  close(21)

  !4 calculation
  print *, 'start calculation'
  do nn = 1 , output_times
    do n = 1 , output_interval
      !4-1 velocity(n-1/2)
      do i=1,total_particle_num!v(n-1) + F(n-1)*dt/2m
        moment(1,i) = moment(1,i) + force(1,i) * dt / (2.0d0 * m)
        moment(2,i) = moment(2,i) + force(2,i) * dt / (2.0d0 * m)
        moment(3,i) = moment(3,i) + force(3,i) * dt / (2.0d0 * m)
      end do

      !4-2 position
      do i=1,total_particle_num!r(n) = [r(n-1)] + [v(n-1)+F(n-1)*dt/2m] * dt
        pos(1,i) = pos(1,i) + moment(1,i) * dt
        pos(2,i) = pos(2,i) + moment(2,i) * dt
        pos(3,i) = pos(3,i) + moment(3,i) * dt
        ! if (pos(1,i)>wall_x .or. pos(1,i)<(-wall_x)) then
        !   pos(1,i)=-pos(1,i)/abs(pos(1,i))*wall_x + (pos(1,i)-wall_x))
        ! end if
        ! if (pos(2,i)>wall_y .or. pos(2,i)<(-wall_y)) then
        !   pos(2,i)=-pos(2,i)/abs(pos(2,i))*wall_y + (pos(2,i)-wall_y))
        ! end if
        ! if (pos(3,i)>wall_z .or. pos(3,i)<(-wall_z)) then
        !   pos(3,i)=-pos(3,i)/abs(pos(3,i))*wall_z + (pos(3,i)-wall_z))
        ! end if
      end do

      !4-3 force
      !distance
      if (f_lj == 1 .or. f_lnk == 1 .or. f_mag_bet_dipole == 1)then
        distances = 0.0d0
        do i=1,total_particle_num-1
          do j=i+1,total_particle_num
            call distance(total_particle_num,i,j,pos,distances)
          end do
        end do
      endif
      !noise
      if (f_noise == 1)then
        noise = 0.0d0
        do i=1,total_particle_num
          call fnoise(total_particle_num,i,gamma,k,temp,m,dt,noise,seed)
        end do
      end if
      !repulsion
      if (f_lj == 1)then
        lj = 0.0d0
        do i=1,total_particle_num-1
          do j=i+1,total_particle_num
            call repulsion(total_particle_num,i,j,sigma,&
            epsilon,lj,distances)
          end do
        end do
      end if
      !reflect
      if (f_wlj == 1)then
        wlj = 0.0d0
        do i = 1, total_particle_num
          call reflect(total_particle_num,i,wall_x,wall_y,&
          wall_z,sigma,epsilon,pos,wlj,counter)
        end do
      end if
      !linkingforce
      if (f_lnk == 1)then
        lnk = 0.0d0
        do i = 1, total_particle_num-1
          call linkingforce(total_particle_num,i,rod,lnk,distances)
        end do
      end if
      !collision
      if (f_clsn == 1)then
        clsn = 0.0d0
        do i = 1, total_particle_num
          call collision(total_particle_num,i,sigma,epsilon,poreR,memT,pos,clsn,counter)
        end do
      end if

      !mag_bet
      if (f_mag_bet_dipole == 1)then
        mag_force=0.0d0 
        do i = 1, total_particle_num-1
          do j = i+1, total_particle_num
            if (distances(i,j,4)<r_th**2.0d0)then
                call interaction(total_particle_num,i,j,distances,&
                mag_permeability,mag_force,mu1_v,mu1_v)
            end if
          end do
        end do
      end if

      !ext\mag\field
      if (f_mag_bet_ext == 1)then
        call  ext_mag_field(total_particle_num,mu1,r_H,r_H_v,H,&
        theta_mu,phi_mu,mag_permeability,mag_force)
      end if

      !electrophoresis
      if (f_elp == 1)then
        elp = 0.0d0
        do i = 1, total_particle_num
          call electrophoresis(total_particle_num,i,qE,poreR,memT,pos,elp)
        end do
      end if
      !F(n) use r(n),v(n-1)
      do i=1,total_particle_num
        v = sqrt(pre_moment(1,i)**2.0d0 + pre_moment(2,i)**2.0d0 + pre_moment(3,i)**2.0d0)
        force(1,i) = - rho*Cd*pi*((D/2.0d0)**2.0d0)*v*pre_moment(1,i)/2.0d0 - 3.0d0*pi*D*eta*pre_moment(1,i) &
        + noise(1,i) + lj(1,i) + wlj(1,i) + lnk(1,i) + clsn(1,i) + elp(1,i)
        force(2,i) = - rho*Cd*pi*((D/2.0d0)**2.0d0)*v*pre_moment(2,i)/2.0d0 - 3.0d0*pi*D*eta*pre_moment(2,i) &
        + noise(2,i) + lj(2,i) + wlj(2,i) + lnk(2,i) + clsn(2,i) + elp(2,i)
        force(3,i) = - rho*Cd*pi*((D/2.0d0)**2.0d0)*v*pre_moment(3,i)/2.0d0 - 3.0d0*pi*D*eta*pre_moment(3,i) &
        + noise(3,i) + lj(3,i) + wlj(3,i) + lnk(3,i) + clsn(3,i) + elp(3,i)
      end do

      !4-4 velocity
      !v(n) = [v(n-1)+F(n-1)*dt/2m] + F(n)*dt/2m (= velocity(n-1/2) + F(n)*dt/2m)
      do i=1,total_particle_num
        moment(1,i) = moment(1,i) + force(1,i) * dt / (2.0d0 * m)
        moment(2,i) = moment(2,i) + force(2,i) * dt / (2.0d0 * m)
        moment(3,i) = moment(3,i) + force(3,i) * dt / (2.0d0 * m)
        pre_moment(1,i) = moment(1,i)
        pre_moment(2,i) = moment(2,i)
        pre_moment(3,i) = moment(3,i)
      end do
    end do

    !4-5 output
    write(*,fmt='(I5)',advance='no')nn
    filename = ''
    write(filename , '(a, a, "particle_colloid_", i5.5, ".cdv")')trim(dirname), slash, nn
    filename = trim(filename)
    open(22,file=filename)
      !command
      write(22,'("#bond_file=", a)') bondpath
      write(22,*)'#r8=0.1' !radius of marker baes (pore and square)
      !beas chain
      do i=1,total_particle_num
        if( pos(3,i) >  0.0d0 )then
          write(22,'(i4,i3,3e14.4)')i , 3 , pos(1,i) , pos(2,i) , pos(3,i)
        else if(pos(3,i) == 0.0d0 .and. pos(3,i) == -0.0d0)then
          write(22,'(i4,i3,3e14.4)')i , 1 , pos(1,i) , pos(2,i) , pos(3,i)
        else
          write(22,'(i4,i3,3e14.4)')i , 2 , pos(1,i) , pos(2,i) , pos(3,i)
        end if
      end do
      !marker baes (square)
      write(22,'(i4,i3,3e14.4)')total_particle_num + 1 , 8 , wall_x , -wall_y , 0.0d0
      write(22,'(i4,i3,3e14.4)')total_particle_num + 2 , 8 , wall_x , wall_y , 0.0d0
      write(22,'(i4,i3,3e14.4)')total_particle_num + 3 , 8 , -wall_x , wall_y , 0.0d0
      write(22,'(i4,i3,3e14.4)')total_particle_num + 4 , 8 , -wall_x , -wall_y , 0.0d0
    close(22)

  end do

end subroutine simulation
