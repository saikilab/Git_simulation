module simulation_mod
  interface
    subroutine simulation(root_dirname, dirname)
      character(100),intent(in) :: root_dirname, dirname
    end subroutine simulation
  end interface
end module simulation_mod

module fposition_mod
  interface
    subroutine fposition(wall_x,wall_y,wall_z,&
      D,total_particle_num,phi_max,rod,pos,f_initial_placement)
      integer,intent(in)::total_particle_num,f_initial_placement
      double precision,intent(in)::wall_x,wall_y,wall_z,D
      double precision,intent(in)::phi_max,rod
      double precision,intent(out) :: pos(3,total_particle_num)
    end subroutine fposition
  end interface
end module fposition_mod

module fmoment_mod
  interface
    subroutine fmoment(total_particle_num,i,k,temp,m,moment)
      integer,intent(in)::total_particle_num,i
      double precision,intent(in)::k,temp,m
      double precision,intent(out) :: moment(3 , total_particle_num)
    end subroutine fmoment
  end interface
end module fmoment_mod

module fnoise_mod
  interface
    subroutine fnoise(total_particle_num,i,gamma,k,temp,m,dt,noise,seed)
      integer,intent(in) :: seed
      integer,intent(in)::total_particle_num
      double precision,intent(in)::k,temp,m,gamma,dt
      double precision,intent(out) :: noise(3 , total_particle_num)
    end subroutine fnoise
  end interface
end module fnoise_mod

module distance_mod
  interface
    subroutine distance(total_particle_num,i,j,pos,distances)
      integer,intent(in)::total_particle_num,i,j
      double precision,intent(in)::pos(3,total_particle_num)
      double precision,intent(out)::distances(total_particle_num,total_particle_num,4)
    end subroutine distance
  end interface
end module distance_mod

module repulsion_mod
  interface
    subroutine repulsion(total_particle_num,i,j,sigma,epsilon,lj,distances)
      integer,intent(in)::total_particle_num,i,j
      double precision,intent(in)::sigma,epsilon
      double precision,intent(in)::distances(total_particle_num,total_particle_num,4)
      double precision,intent(out)::lj(3, total_particle_num)
    end subroutine repulsion
  end interface
end module repulsion_mod

module reflect_mod
  interface
    subroutine reflect(total_particle_num,i,wall_x,wall_y,&
      wall_z,sigma,epsilon,pos,wlj,counter)
      integer,intent(in)::total_particle_num,i
      double precision,intent(in)::wall_x,wall_y,wall_z
      double precision,intent(in)::sigma,epsilon
      double precision,intent(in)::pos(3,total_particle_num)
      double precision,intent(out)::wlj(3,total_particle_num)
      integer,intent(out)::counter
    end subroutine reflect
  end interface
end module reflect_mod

module linkingforce_mod
  interface
    subroutine linkingforce(total_particle_num,i,rod,lnk,distances)
      integer,intent(in)::total_particle_num,i
      double precision,intent(in)::rod
      double precision,intent(in)::distances(total_particle_num,total_particle_num,4)
      double precision,intent(out)::lnk(3, total_particle_num)
    end subroutine linkingforce
  end interface
end module linkingforce_mod

module collision_mod
  interface
    subroutine collision(total_particle_num,i,sigma,epsilon,poreR,memT,pos,clsn,counter)
      integer,intent(in)::total_particle_num,i
      double precision,intent(in)::sigma,epsilon,poreR,memT
      double precision,intent(in)::pos(3,total_particle_num)
      double precision,intent(out)::clsn(3,total_particle_num)
      integer,intent(out)::counter
    end subroutine collision
  end interface
end module collision_mod

module electrophoresis_mod
  interface
    subroutine electrophoresis(total_particle_num,i,qE,poreR,memT,pos,elp)
      integer,intent(in)::total_particle_num,i
      double precision,intent(in)::qE,poreR,memT
      double precision,intent(in)::pos(3,total_particle_num)
      double precision,intent(out)::elp(3, total_particle_num) !N (kg*m*s^-1)
    end subroutine electrophoresis
  end interface
end module electrophoresis_mod
