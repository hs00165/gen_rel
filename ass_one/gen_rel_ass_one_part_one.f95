!General Relativity

PROGRAM gen_rel_three_c
  IMPLICIT NONE

DOUBLE PRECISION::M_sun, M_earth, M_merc, rp_earth, ra_earth, rp_merc, ra_merc, a_earth, e_earth, a_merc, e_merc, &
step, M_total_earth, M_total_merc, reduced_mass_earth, reduced_mass_merc, E_i_earth, E_f_earth, E_i_merc, &
E_f_merc, earth_period, merc_period, G, dE_earth, dE_merc, time, ax, ay, az, vx_temp, vy_temp, vz_temp, r, &
dE_earth_rel, dE_merc_rel
DOUBLE PRECISION, DIMENSION(1:3) :: pos_earth, pos_merc, v_earth, v_merc
INTEGER :: j

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

WRITE(6,*)
WRITE(6,*) "Part 3:"
WRITE(6,*) "This part of the program will show the energy loss as a function of the step size."
WRITE(6,*) "Written to screen will be the step size, initial energy, final energy, absolute"
WRITE(6,*) "error and relative error after one orbital period" 
WRITE(6,*)
WRITE(6,*)
WRITE(6,*)

step = 0.1

OPEN(11, file='e_loss_step_size.out')
DO j=1,6,1

!Time period for Earth:
earth_period = sqrt((39.4784176*(a_earth**3))/(M_total_earth*G))

!Time period for Mercury:
merc_period = sqrt((39.4784176*(a_merc**3))/(M_total_merc*G))

!WRITE(6,*) "Earth and Mercury periods:"
!WRITE(6,*) earth_period, merc_period

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


!WRITE(6,*) "========================================================================="
!WRITE(6,*) "============================LEAP FROG INTEGRATION========================"
!WRITE(6,*) "========================================================================="
!WRITE(6,*)


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

IF (j==6) THEN
WRITE(6,*) "Writing out the final orbit to file, please wait..."
END IF

!=================================================================================
!====================================EARTH ORBIT==================================
!=================================================================================

IF (j==6) THEN
OPEN(10, file='earth_orbit.out')
END IF
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



IF (j==6) THEN
WRITE(10,*) pos_earth(1), pos_earth(2), pos_earth(3), v_earth(1), v_earth(2), v_earth(3)
END IF

    time = time + step

  END DO

IF (j==6) THEN
   CLOSE(10)
END IF

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

IF (j==6) THEN
  OPEN(10, file='merc_orbit.out')
END IF

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



IF (j==6) THEN    
WRITE(10,*) pos_merc(1), pos_merc(2), pos_merc(3), v_merc(1), v_merc(2), v_merc(3)
END IF

    time = time + step

  END DO
IF (j==6) THEN
   CLOSE(10)
END IF

  !Final Mercury orbit energy
  E_f_merc = 0.5*((v_merc(1)**2)+(v_merc(2)**2)+(v_merc(3)**2)) - ((G*M_total_merc)/(sqrt((pos_merc(1)**2)+&
  (pos_merc(2)**2)+(pos_merc(3)**2))))



  !=============ENERGY LOSS=============
  dE_earth_rel = ABS((E_f_earth - E_i_earth)/E_i_earth)
  dE_merc_rel = ABS((E_f_merc - E_i_merc)/E_i_merc)
  dE_earth = ABS(E_f_earth - E_i_earth)
  dE_merc = ABS(E_f_merc - E_i_merc)


WRITE(6,*) "==================================================================="
WRITE(6,*) "Step size:", step
WRITE(6,*) "                    EARTH:                  ", "MERCURY:"
WRITE(6,*) "Initial energy: ", E_i_earth, E_i_merc
WRITE(6,*) "Final energy:   ", E_f_earth, E_f_merc
WRITE(6,*) "Absolute error:", dE_earth, dE_merc
WRITE(6,*) "Relative error:", dE_earth_rel, dE_merc_rel
WRITE(6,*) "==================================================================="
WRITE(6,*)




!WRITE(6,'(a5, f7.4, f7.4, f7.4, e7.4, e7.4)') "Earth", step, E_i_earth, E_f_earth, dE_earth, dE_earth_rel
!WRITE(6,'(a7, f7.4, f7.4, f7.4, e7.4, e7.4)') "Mercury", step, E_i_merc, E_f_merc, dE_merc, dE_merc_rel
  !WRITE(6,*)
  !WRITE(6,*) "The total energy loss is:"
  !WRITE(6,*) "Earth:", dE_earth, "         Mercury:", dE_merc
  !WRITE(6,*)
  !WRITE(6,*) "The relative energy loss is:"
  !WRITE(6,*) "Earth:", dE_earth_rel, "         Mercury:", dE_merc_rel




  WRITE(11,*) step, dE_earth, dE_merc, dE_earth_rel, dE_merc_rel
  step = step/10

END DO
  CLOSE(11)

WRITE(6,*)
WRITE(6,*) "Three files have now been created."
WRITE(6,*)
WRITE(6,*) "earth_orbit.out         <- Earths orbit"
WRITE(6,*) "merc_orbit.out          <- Murcurys orbit"
WRITE(6,*) "e_loss_step_size.out    <- Energy loss data file"
WRITE(6,*)
WRITE(6,*) "The energy loss data file is written in the following collumns:"
WRITE(6,*) "Step size | abs E loss Earth | abs E loss Mercury | Rel E loss Earth | Rel E loss Mercury"
WRITE(6,*)




END PROGRAM gen_rel_three_c



















