subroutine fposition(wall_x,wall_y,wall_z,D,&
  total_particle_num,phi_max,rod,pos,f_initial_placement)

  use mt19937
  implicit none

  integer,intent(in)::total_particle_num,f_initial_placement
  double precision,intent(in)::wall_x,wall_y,wall_z,D
  double precision,intent(in)::phi_max, rod
  double precision,intent(out) :: pos(3,total_particle_num)

  integer::i,j,k,l
  integer :: seed
  double precision , parameter :: pi = dacos(-1.0d0)
  double precision :: phi, theta
  double precision :: dr_xy
  double precision :: dr(3)
  double precision :: altx(3),alty(3),altz(3)

  call system_clock(count=seed)
  call sgrnd(seed)

  if(f_initial_placement == 1)then
    do i = 1, int(wall_x*2.0d0/D)-1
      do j = 1, int(wall_y*2.0d0/D)-1
        l = int(wall_y*2-1)*int(i-1) + j
        pos(1,l) = -wall_x + dble(i) * D
        pos(2,l) = -wall_y + dble(j) * D
        pos(3,l) = 0.0d0
        if (l == total_particle_num)exit
      end do
      if (l == total_particle_num)exit
    end do

    if (l < total_particle_num)then
      stop 'too many particles'
    end if
  else
    phi = 0.0d0
    theta = 0.0d0
    pos(1,1) = 0.0d0!5d-7*sin(phi)*cos(theta)
    pos(2,1) = 0.0d0!5d-7*sin(phi)*sin(theta)
    pos(3,1) = 8.0d0!5d-7*cos(phi)

    pos(1,2) = 0.0d0!5d-7*sin(phi)*cos(theta)
    pos(2,2) = 0.0d0!5d-7*sin(phi)*sin(theta)
    pos(3,2) = 2.0d0!5d-7*cos(phi)
  endif
end subroutine fposition


subroutine fmoment(total_particle_num,i,k,temp,m,moment)
  use mt19937
  implicit none

  integer,intent(in)::total_particle_num,i
  double precision,intent(in)::k,temp,m
  double precision,intent(out) :: moment(3 , total_particle_num)

  double precision , parameter :: pi = dacos(-1.0d0)
  integer :: seed

  call system_clock(count=seed)
  call sgrnd(seed)

  moment(1,i)=sqrt(k*temp/m) * sqrt(-2.0d0*log(grnd())) * cos(2.0d0*pi*grnd())
  moment(2,i)=sqrt(k*temp/m) * sqrt(-2.0d0*log(grnd())) * cos(2.0d0*pi*grnd())
  moment(3,i)=sqrt(k*temp/m) * sqrt(-2.0d0*log(grnd())) * cos(2.0d0*pi*grnd())
end subroutine fmoment


subroutine fnoise(total_particle_num,i,gamma,k,temp,m,dt,noise,seed)
  use mt19937
  implicit none

  integer,intent(in)::seed
  integer,intent(in)::total_particle_num,i
  double precision,intent(in)::k,temp,m,gamma,dt
  double precision,intent(out) :: noise(3 , total_particle_num)

  double precision , parameter :: pi = dacos(-1.0d0)
  double precision :: variance !bunsan
  variance = 2.0d0*gamma*k*temp*m/dt !bunsan

  !Box–Muller's method
  noise(1,i) = sqrt(variance) * sqrt(-2.0d0*log(grnd())) * cos(2.0d0*pi*grnd())
  noise(2,i) = sqrt(variance) * sqrt(-2.0d0*log(grnd())) * cos(2.0d0*pi*grnd())
  noise(3,i) = sqrt(variance) * sqrt(-2.0d0*log(grnd())) * cos(2.0d0*pi*grnd())
end subroutine fnoise


subroutine distance(total_particle_num,i,j,pos,distances)
  implicit none
  integer,intent(in)::total_particle_num,i,j
  double precision,intent(in)::pos(3,total_particle_num)
  double precision,intent(out)::distances(total_particle_num,total_particle_num,4)
  distances(i,j,1) = pos(1,i) - pos(1,j)
  distances(i,j,2) = pos(2,i) - pos(2,j)
  distances(i,j,3) = pos(3,i) - pos(3,j)
  distances(i,j,4) = distances(i,j,1)**2.0d0 + distances(i,j,2)**2.0d0 + distances(i,j,3)**2.0d0
end subroutine distance


subroutine repulsion(total_particle_num,i,j,sigma,epsilon,lj,distances)
  implicit none
  integer,intent(in)::total_particle_num,i,j
  double precision,intent(in)::sigma,epsilon
  double precision,intent(in)::distances(total_particle_num,total_particle_num,4)
  double precision,intent(out)::lj(3, total_particle_num) !N (kg*m*s^-1)
  double precision::flj !N*m^-1 (kg*s^-2) F(r)/r = -dU(r)/dr /r
  double precision::sigma_2 !m^2 equilibrium_interatomic_distance^2 heikougensikankyori
  sigma_2 = ((2.0d0**(1.0d0/6.0d0)) * sigma)**2.0d0

  if(distances(i,j,4) < sigma_2)then
    !Lennard-Jones potential
    flj = 24.0d0 * epsilon* ( 2.0d0 * sigma**12.0d0 / distances(i,j,4)**7.0d0 - sigma**6.0d0 / distances(i,j,4)**4.0d0 )

    lj(1,i) = lj(1,i) + flj * distances(i,j,1)
    lj(2,i) = lj(2,i) + flj * distances(i,j,2)
    lj(3,i) = lj(3,i) + flj * distances(i,j,3)

    lj(1,j) = lj(1,j) - flj * distances(i,j,1)
    lj(2,j) = lj(2,j) - flj * distances(i,j,2)
    lj(3,j) = lj(3,j) - flj * distances(i,j,3)
  end if
end subroutine repulsion


subroutine reflect(total_particle_num,i,wall_x,wall_y,&
  wall_z,sigma,epsilon,pos,wlj,counter)
  implicit none
  integer,intent(in)::total_particle_num,i
  double precision,intent(in)::wall_x,wall_y,wall_z
  double precision,intent(in)::sigma,epsilon
  double precision,intent(in)::pos(3,total_particle_num)
  double precision,intent(out)::wlj(3,total_particle_num)
  integer,intent(out)::counter
  double precision::sigma2 !m equilibrium_interatomic_distance heikougensikankyori
  sigma2 = (2.0d0**(1.0d0/6.0d0)) * sigma

  if((pos(1,i))**2.0d0 > (wall_x+sigma/2.0d0)**2.0d0)then
    counter = counter +1
    if (counter ==4)then
      stop "particle out of wall x"
    endif
  elseif ((pos(2,i))**2.0d0 > (wall_y+sigma/2.0d0)**2.0d0)then
    counter = counter +1
    if (counter ==4)then
      stop "particle out of wall y"
    endif
  elseif ((pos(3,i))**2.0d0 > (wall_z+sigma/2.0d0)**2.0d0)then
    counter = counter +1
    if (counter ==4)then
      stop "particle out of wall z"
    endif
  end if

  if(pos(1,i) > (wall_x-sigma2))then
    wlj(1,i) = 24.0d0*epsilon*(2.0d0*sigma**12.0d0/((pos(1,i)-wall_x)**13.0d0)&
    -sigma**6.0d0/((pos(1,i)-wall_x)**7.0d0))
  elseif(pos(1,i) < -(wall_x-sigma2))then
    wlj(1,i) = 24.0d0*epsilon*(2.0d0*sigma**12.0d0/((pos(1,i)+wall_x)**13.0d0)&
    -sigma**6.0d0/((pos(1,i)+wall_x)**7.0d0))
  else
    wlj(1,i) = 0.0d0
  end if

  if(pos(2,i) > (wall_y-sigma2))then
    wlj(2,i) = 24.0d0*epsilon*(2.0d0*sigma**12.0d0/((pos(2,i)-wall_y)**13.0d0)&
    -sigma**6.0d0/((pos(2,i)-wall_y)**7.0d0))
  elseif(pos(2,i) < -(wall_y-sigma2))then
    wlj(2,i) = 24.0d0*epsilon*(2.0d0*sigma**12.0d0/((pos(2,i)+wall_y)**13.0d0)&
    -sigma**6.0d0/((pos(2,i)+wall_y)**7.0d0))
  else
    wlj(2,i) = 0.0d0
  end if

  if(pos(3,i) > (wall_z-sigma2))then
    wlj(3,i) = 24.0d0*epsilon*(2.0d0*sigma**12.0d0/((pos(3,i)-wall_z)**13.0d0)&
    -sigma**6.0d0/((pos(3,i)-wall_z)**7.0d0))
  elseif(pos(3,i) < -(wall_z-sigma2))then
    wlj(3,i) = 24.0d0*epsilon*(2.0d0*sigma**12.0d0/((pos(3,i)+wall_z)**13.0d0)&
    -sigma**6.0d0/((pos(3,i)+wall_z)**7.0d0))
  else
    wlj(3,i) = 0.0d0
  end if
end subroutine reflect


subroutine linkingforce(total_particle_num,i,rod,lnk,distances)
  implicit none
  integer,intent(in)::total_particle_num,i
  double precision,intent(in)::rod
  double precision,intent(in)::distances(total_particle_num,total_particle_num,4)
  double precision,intent(out)::lnk(3, total_particle_num) !N (kg*m*s^-1)
  double precision::lnkf, k
  ! double precision::R0, R1, R2, R3
  k = 250d-1
  lnkf = k*(sqrt(distances(i,i+1,4)) - rod)/sqrt(distances(i,i+1,4))

  lnk(1,i) = lnk(1,i) - lnkf*distances(i,i+1,1)
  lnk(2,i) = lnk(2,i) - lnkf*distances(i,i+1,2)
  lnk(3,i) = lnk(3,i) - lnkf*distances(i,i+1,3)

  lnk(1,i+1) = lnk(1,i+1) + lnkf*distances(i,i+1,1)
  lnk(2,i+1) = lnk(2,i+1) + lnkf*distances(i,i+1,2)
  lnk(3,i+1) = lnk(3,i+1) + lnkf*distances(i,i+1,3)

  ! if(distances(i,i+1,4) < R0**2)then
  !   lnk(1,i)
  !   lnk(2,i)
  !   lnk(3,i)
  ! elseif(distances(i,i+1,4) < R1**2)then
  !   lnk(1,i)
  !   lnk(2,i)
  !   lnk(3,i)
  ! elseif(distances(i,i+1,4) < R2**2)then
  !   lnk(1,i) = 0.0d0
  !   lnk(2,i) = 0.0d0
  !   lnk(3,i) = 0.0d0
  ! elseif(distances(i,i+1,4) < R3**2)then
  !   lnk(1,i)
  !   lnk(2,i)
  !   lnk(3,i)
  ! else
  !   lnk(1,i)
  !   lnk(2,i)
  !   lnk(3,i)
  ! end if
end subroutine linkingforce


subroutine collision(total_particle_num,i,sigma,epsilon,poreR,memT,pos,clsn,counter)
  implicit none
  integer,intent(in)::total_particle_num,i
  double precision,intent(in)::sigma,epsilon,poreR,memT
  double precision,intent(in)::pos(3,total_particle_num)
  double precision,intent(out)::clsn(3,total_particle_num)
  integer,intent(out)::counter
  double precision::P(3)
  double precision::r
  double precision::flj !N*m^-1 (kg*s^-2) F(r)/r = -dU(r)/dr /r
  double precision::sigma2 !m equilibrium_interatomic_distance heikougensikankyori
  sigma2 = (2.0d0**(1.0d0/6.0d0)) * sigma

  if(pos(3,i)**2.0d0 <= (memT/2.0d0)**2.0d0)then
    if(pos(1,i)**2.0d0 + pos(2,i)**2.0d0 >= poreR**2.0d0)then
      counter = counter +1
      if (counter ==4)then
        stop "particle in the pore"
      endif
    endif
  end if

  P = 0.0d0
  r = 0.0d0
  if (pos(1,i)**2.0d0 + pos(2,i)**2.0d0 + pos(3,i)**2.0d0 == 0.0d0)then
    clsn(1,i) = 0.0d0
    clsn(2,i) = 0.0d0
    clsn(3,i) = 0.0d0
  elseif (pos(3,i)**2.0d0 < ((memT/2.0d0)+sigma2)**2.0d0)then
    if (pos(1,i)**2.0d0 + pos(2,i)**2.0d0 > poreR**2.0d0)then
      P(1) = pos(1,i)
      P(2) = pos(2,i)
      P(3) = memT/2.0d0 * pos(3,i)/abs(pos(3,i))
    elseif (pos(1,i)**2.0d0 + pos(2,i)**2.0d0 <= poreR**2.0d0)then
      P(1) = pos(1,i) * poreR/sqrt(pos(1,i)**2.0d0+pos(2,i)**2.0d0)
      P(2) = pos(2,i) * poreR/sqrt(pos(1,i)**2.0d0+pos(2,i)**2.0d0)
      P(3) = 0.0d0
    endif
    r = sqrt((pos(1,i)-P(1))**2.0d0 + (pos(2,i)-P(2))**2.0d0 + (pos(3,i)-P(3))**2.0d0)
    !Lennard-Jones potential
    flj = 24.0d0 * epsilon* ( 2.0d0 * sigma**12.0d0/r**13.0d0 - sigma**6.0d0/r**7.0d0 )
    clsn(1,i) = flj*(pos(1,i)-P(1))/r
    clsn(2,i) = flj*(pos(2,i)-P(2))/r
    clsn(3,i) = flj*(pos(3,i)-P(3))/r
  endif
end subroutine collision


subroutine electrophoresis(total_particle_num,i,qE,poreR,memT,pos,elp)
  implicit none
  integer,intent(in)::total_particle_num,i
  double precision,intent(in)::qE,poreR,memT
  double precision,intent(in)::pos(3,total_particle_num)
  double precision,intent(out)::elp(3, total_particle_num) !N (kg*m*s^-1)
  double precision::R0
  double precision::r !m kyori
  R0 = sqrt(poreR**2.0d0 + (memT/2.0d0)**2.0d0)
  r = sqrt(pos(1,i)**2.0d0 + pos(2,i)**2.0d0 + pos(3,i)**2.0d0)

  if(pos(3,i) == 0.0d0)then
    elp(1,i) = 0.0d0
    elp(2,i) = 0.0d0
    elp(3,i) = - qE
  elseif(r > R0)then
    elp(1,i) = qE/((r/R0)**2.0d0) *(pos(1,i)/r)*(-pos(3,i))/abs(pos(3,i))
    elp(2,i) = qE/((r/R0)**2.0d0) *(pos(2,i)/r)*(-pos(3,i))/abs(pos(3,i))
    elp(3,i) = qE/((r/R0)**2.0d0) *(pos(3,i)/r)*(-pos(3,i))/abs(pos(3,i))
  elseif(r <= R0)then
    elp(1,i) = qE *(pos(1,i)/r)*(-pos(3,i))/abs(pos(3,i))
    elp(2,i) = qE *(pos(2,i)/r)*(-pos(3,i))/abs(pos(3,i))
    elp(3,i) = qE *(pos(3,i)/r)*(-pos(3,i))/abs(pos(3,i))
  else
    elp(1,i) = 0.0d0
    elp(2,i) = 0.0d0
    elp(3,i) = 0.0d0
  end if
end subroutine electrophoresis
