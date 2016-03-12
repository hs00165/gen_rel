!General Relativity

PROGRAM gen_rel_three_c
  IMPLICIT NONE

DOUBLE PRECISION::M_sun, M_earth, M_merc, rp_earth, ra_earth, rp_merc, ra_merc, a_earth, e_earth, a_merc, e_merc, &
step, M_total_earth, M_total_merc, reduced_mass_earth, reduced_mass_merc, E_i_earth, E_f_earth, E_i_merc, &
E_f_merc, earth_period, merc_period, G, dE_earth, dE_merc, time, ax, ay, az, vx_temp, vy_temp, vz_temp, r
DOUBLE PRECISION, DIMENSION(1:3) :: pos_earth, pos_merc, v_earth, v_merc


!  INTERFACE
!    SUBROUTINE leap_frog_int(x_pos, y_pos, z_pos, vx, vy, vz, M, t)
!      DOUBLE PRECISION :: x_pos, y_pos, z_pos, vx, vy, vz, M, r, G, t, ax, ay, az
!    END SUBROUTINE leap_frog_int
!  END INTERFACE




!=============================
G=1.0
M_sun = 1.0
M_earth = 3.0024584e-6
M_merc = 1.65956463e-7
M_total_earth = M_sun + M_earth
M_total_merc = M_sun + M_merc
reduced_mass_earth = (M_sun*M_earth)/M_total_earth
reduced_mass_merc = (M_sun*M_merc)/M_total_merc



a_earth = 1
e_earth = 0.0167
a_merc = 0.387098
e_merc = 0.205635


rp_earth = 0.9833
ra_earth = 1.0167
rp_merc = 0.307497
ra_merc = 0.466699
!=============================



!WRITE(6,*) "========================================================================="
!WRITE(6,*) "============================INITIAL CONDITIONS==========================="
!WRITE(6,*) "========================================================================="
!WRITE(6,*)

!NOTE: start position is at the pericentre of the orbit, therefore
!      the only velocity is in the y component. This means that the 
!      initial velocity in the y axis is just the square root of Vp^2.



!Time period for Earth:
earth_period = sqrt((39.4784176*(a_earth**3))/(M_total_earth*G))

!Time period for Mercury:
merc_period = sqrt((39.4784176*(a_merc**3))/(M_total_merc*G))

WRITE(6,*) "Earth and Mercury periods:"
WRITE(6,*) earth_period, merc_period

pos_earth(1) = -rp_earth
pos_earth(2) = 0.0
pos_earth(3) = 0.0

v_earth(1) = 0.0
v_earth(2) = sqrt((1.0/a_earth)*((1+e_earth)/(1-e_earth)))
v_earth(3) = 0.0


pos_merc(1) = -rp_merc
pos_merc(2) = 0.0
pos_merc(3) = 0.0

v_merc(1) = 0.0
v_merc(2) = sqrt((1.0/a_merc)*((1+e_merc)/(1-e_merc)))
v_merc(3) = 0.0

OPEN(13, file='Earth.ini')
WRITE(13,*) pos_earth(1), pos_earth(2), pos_earth(3), v_earth(1), v_earth(2), v_earth(3)
  CLOSE(13)

OPEN(14, file='Mercury.ini')
WRITE(14,*) pos_merc(1), pos_merc(2), pos_merc(3), v_merc(1), v_merc(2), v_merc(3)
  CLOSE(14)

WRITE(6,*) "========================================================================="
WRITE(6,*) "============================LEAP FROG INTEGRATION========================"
WRITE(6,*) "========================================================================="
WRITE(6,*)
WRITE(6,*) "Please enter the step size:"
READ(5,*) step 

  !=======INITIAL ORBIT ENERGIES===========
  !Initial Earth orbit energy
  E_i_earth = 0.5*((v_earth(1)**2)+(v_earth(2)**2)+(v_earth(3)**2)) - ((G*M_total_earth)/(sqrt((pos_earth(1)**2)+&
  (pos_earth(2)**2)+(pos_earth(3)**2))))
  !Initial Mercury orbit energy
  E_i_merc = 0.5*((v_merc(1)**2)+(v_merc(2)**2)+(v_merc(3)**2)) - ((G*M_total_merc)/(sqrt((pos_merc(1)**2)+&
  (pos_merc(2)**2)+(pos_merc(3)**2))))

  time = 0.0

  !Initial setting of r and acceleration 
  r = sqrt((pos_earth(1)**2)+(pos_earth(2)**2)+(pos_earth(3)**2))

  ax = -(G*M_total_earth*pos_earth(1))/(r**3)
  ay = -(G*M_total_earth*pos_earth(2))/(r**3)
  az = -(G*M_total_earth*pos_earth(3))/(r**3)



!=================================================================================
!====================================EARTH ORBIT==================================
!=================================================================================
  OPEN(10, file='earth_orbit.out')
  DO WHILE (time <= earth_period)

    pos_earth(1) = pos_earth(1) + v_earth(1)*step + (ax*(step**2)/2.0)
    pos_earth(2) = pos_earth(2) + v_earth(2)*step + (ay*(step**2)/2.0)
    pos_earth(3) = pos_earth(3) + v_earth(3)*step + (az*(step**2)/2.0)

    vx_temp = v_earth(1) + ((ax*step)/2.0)
    vy_temp = v_earth(2) + ((ay*step)/2.0)
    vz_temp = v_earth(3) + ((az*step)/2.0)

    r = sqrt((pos_earth(1)**2)+(pos_earth(2)**2)+(pos_earth(3)**2))
    ax = -(G*M_total_earth*pos_earth(1))/(r**3)
    ay = -(G*M_total_earth*pos_earth(2))/(r**3)
    az = -(G*M_total_earth*pos_earth(3))/(r**3)

    v_earth(1) = vx_temp + ((ax*step)/2.0)
    v_earth(2) = vy_temp + ((ay*step)/2.0)
    v_earth(3) = vz_temp + ((az*step)/2.0)



  
    WRITE(10,*) pos_earth(1), pos_earth(2), pos_earth(3), v_earth(1), v_earth(2), v_earth(3)

    time = time + step

  END DO
   CLOSE(10)

  !Final Earth orbit energy
  E_f_earth = 0.5*((v_earth(1)**2)+(v_earth(2)**2)+(v_earth(3)**2)) - ((G*M_total_earth)/(sqrt((pos_earth(1)**2)+&
  (pos_earth(2)**2)+(pos_earth(3)**2))))

!=================================================================================
!===================================MERCURY ORBIT=================================
!=================================================================================
  time = 0.0

  r = sqrt((pos_merc(1)**2)+(pos_merc(2)**2)+(pos_merc(3)**2))

  ax = -(G*M_total_merc*pos_merc(1))/(r**3)
  ay = -(G*M_total_merc*pos_merc(2))/(r**3)
  az = -(G*M_total_merc*pos_merc(3))/(r**3)

  OPEN(10, file='merc_orbit.out')
  DO WHILE (time <= merc_period)

  
    pos_merc(1) = pos_merc(1) + v_merc(1)*step + (ax*(step**2)/2.0)
    pos_merc(2) = pos_merc(2) + v_merc(2)*step + (ay*(step**2)/2.0)
    pos_merc(3) = pos_merc(3) + v_merc(3)*step + (az*(step**2)/2.0)

    vx_temp = v_merc(1) + ((ax*step)/2.0)
    vy_temp = v_merc(2) + ((ay*step)/2.0)
    vz_temp = v_merc(3) + ((az*step)/2.0)

    r = sqrt((pos_merc(1)**2)+(pos_merc(2)**2)+(pos_merc(3)**2))
    ax = -(G*M_total_merc*pos_merc(1))/(r**3)
    ay = -(G*M_total_merc*pos_merc(2))/(r**3)
    az = -(G*M_total_merc*pos_merc(3))/(r**3)

    v_merc(1) = vx_temp + ((ax*step)/2.0)
    v_merc(2) = vy_temp + ((ay*step)/2.0)
    v_merc(3) = vz_temp + ((az*step)/2.0)



    
    WRITE(10,*) pos_merc(1), pos_merc(2), pos_merc(3), v_merc(1), v_merc(2), v_merc(3)

    time = time + step

  END DO
   CLOSE(10)

  !Final Mercury orbit energy
  E_f_merc = 0.5*((v_merc(1)**2)+(v_merc(2)**2)+(v_merc(3)**2)) - ((G*M_total_merc)/(sqrt((pos_merc(1)**2)+&
  (pos_merc(2)**2)+(pos_merc(3)**2))))



  !=============ENERGY LOSS=============
  dE_earth = ABS((E_f_earth - E_i_earth)/E_i_earth)
  dE_merc = ABS((E_f_merc - E_i_merc)/E_i_merc)


  WRITE(6,*)
  WRITE(6,*) "The Energy loss is:"
  WRITE(6,*) "Earth:", dE_earth, "         Mercury:", dE_merc



END PROGRAM gen_rel_three_c



















