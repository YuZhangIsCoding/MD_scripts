UBROUTINE Surface_area

INCLUDE 'shape_include.f90'

REAL*8 :: ran01,ran02,theta,phi,xcom,ycom,zcom,xsurf,ysurf,zsurf,sigmix1,sigmix2,surf_point,bin_points,area_bin
REAL*8 :: sigma_Au,probe_d,xr,yr,zr,r2,dr,test
INTEGER :: num_neighbors(max_atoms),neighbors(100,max_atoms),dx,dy,dz,dr2,jatom,jjatom,num_surface_atoms,isurf_atom(max_atoms)

 area_tot=0.0d0

 OPEN(UNIT=501,file='surf.xyz')

 DO i=1,nframes_tot

  num_neighbors = 0
  neighbors = 0

  DO j=1,natoms(i)-1
  DO jj=j+1,natoms(i)

    dx = ipos(1,j,i) - ipos(1,jj,i)
    dy = ipos(2,j,i) - ipos(2,jj,i)
    dz = ipos(3,j,i) - ipos(3,jj,i)

    dr2 = dx*dx + dy*dy + dz*dz

    IF(dr2.LE.4) THEN ! NEIGHBOR
      num_neighbors(j) = num_neighbors(j) + 1
      num_neighbors(jj) = num_neighbors(jj) + 1
      neighbors(num_neighbors(j),j) = jj
      neighbors(num_neighbors(jj),jj) = j
    END IF

  END DO
  END DO

  num_surface_atoms = 0
  isurf_atom = 0

  ! Identify the surface atoms
  DO j=1,natoms(i)
   IF(num_neighbors(j).LT.16) THEN
     num_surface_atoms = num_surface_atoms + 1
     isurf_atom(num_surface_atoms) = j
   END IF
  END DO

  probe_d = 0.20d0
  sigma_Au = 2.0d0
  sigmix1 = 0.50d0*(probe_d + sigma_Au)
  !Perform sample tests of the surface atoms
  DO j=1,num_surface_atoms

    area_bin = 0.0d0
    bin_points = 0.0d0

    jatom = isurf_atom(j)

    xcom = DBLE(ipos(1,jatom,i))
    ycom = DBLE(ipos(2,jatom,i))
    zcom = DBLE(ipos(3,jatom,i))

  DO k=1,10000

      ran01 = ran2(iseed)
      ran02 = ran2(iseed)
      theta = 2.0d0*pi*ran01
      phi = ACOS(2.0d0*ran02 - 1.0d0)

      xsurf = xcom + sigmix1*(SIN(phi)*COS(theta))
      ysurf = ycom + sigmix1*(SIN(phi)*SIN(theta))
      zsurf = zcom + sigmix1*(COS(phi))

    surf_point = 1.0d0 ! reset the overlap flag

    test_neighbor : DO jj=1,num_neighbors(jatom)
       jjatom = neighbors(jj,jatom)

       xr = xsurf - DBLE(ipos(1,jjatom,i))
       yr = ysurf - DBLE(ipos(2,jjatom,i))
       zr = zsurf - DBLE(ipos(3,jjatom,i))

       r2 = xr*xr + yr*yr + zr*zr

       dr = SQRT(r2)

       IF(dr.LT.sigmix1) THEN
         surf_point = 0.0d0 ! Reject if an overlap occurs
         EXIT test_neighbor
       END IF

    END DO test_neighbor

    bin_points = bin_points + 1.0d0
    area_bin = area_bin + surf_point*3.14159d0*(sigmix1*2)**2

    !IF(i.EQ.1) THEN
    ! IF(surf_point.GT.0.0d0) THEN ! save the surface points
    !  WRITE(501,*) 'H',xsurf,ysurf,zsurf
    ! END IF
    !END IF

  END DO ! end loop over test vectors

  area_tot(i) = area_tot(i) + area_bin/bin_points

  END DO ! end of surface_atoms

 PRINT *,'Surface area:',area_tot(i)

 END DO ! end frames cycle

CLOSE(501)

END SUBROUTINE Surface_area

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

