program main

  use simulation_mod
  implicit none
  character(10) :: dammy
  character(100) :: start_date
  character(100) :: root_dirname
  character(100) :: dirname
  character(150) :: mkdir
  character(120) :: cdview3
  character(1) :: slash = '\'
  integer :: initial_date(8)
  integer :: f_noise, f_lj, f_wlj, f_lnk, f_clsn, f_elp, f_moment0 , f_mag_bet_dipole , f_mag_bet_ext , bondy
  integer :: f_initial_placement,f_cdview3

  call date_and_time(dammy, dammy, dammy, initial_date)
  write (start_date,'(i4.4,"_",2i2.2,"_",2i2.2)')initial_date(1),&
  initial_date(2),initial_date(3),initial_date(5),initial_date(6)!1:year 2:month 3:day 5:hour 6:minutes ex. 1997_1226_1030

  write (root_dirname,'("C:\Users\Awano\simulation\12A\ishi\result",a,a)')&
  slash,trim(start_date)!ex. C:\Users\saiki60\Documents\ishii\Simulation\results\1997_1226_1030
  write (mkdir,'("mkdir ",a)')trim(root_dirname)
  call system(mkdir)

  write(dirname,'(a,a,"cdvfiles")')trim(root_dirname),slash!ex. C:\Users\saiki60\Documents\ishii\Simulation\results\1997_1226_1030\cdvfiles
  write(mkdir,'("mkdir ",a)')trim(dirname)
  call system(mkdir)

  write(cdview3, '("cdview3 ",a ,a,"particle_colloid_00000.cdv")')trim(dirname),slash

  !read csv
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
  close(10,status = 'delete')

    f_noise = 1
    f_lj = 1
    f_wlj = 1
    f_lnk = 0
    f_clsn = 0
    f_elp = 0
    f_moment0 = 0
    f_initial_placement = 1
    f_cdview3 = 1
    f_mag_bet_dipole = 1
    f_mag_bet_ext = 1
    bondy=0

  !rewrite csv
  open(10, file = 'parameter.csv')
    write(10, *) 'debug character'
    write(10, *) f_noise
    write(10, *) f_lj
    write(10, *) f_wlj
    write(10, *) f_lnk
    write(10, *) f_clsn
    write(10, *) f_elp
    write(10, *) f_moment0
    write(10, *) f_initial_placement
    write(10, *) f_cdview3
    write(10, *) f_mag_bet_dipole
    write(10, *) f_mag_bet_ext
    write(10, *) bondy
    write(10, *) 'debug character'
  close(10)

  call simulation(root_dirname, dirname)
  if (f_cdview3 == 1)then
    call system(cdview3)
  endif

end program main
