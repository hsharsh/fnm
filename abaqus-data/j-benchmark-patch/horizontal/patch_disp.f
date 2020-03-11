      SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION U(3),TIME(2),COORDS(3)
C
      pi=4.D0*DATAN(1.D0)
      cx = 0.33333334300000000
      cy = 0.50000000000000000
      dk = 1
      E = 1
      nu = 0
      G = E/(1+nu)
      mu = 3-4*nu
      x = COORDS(1)-cx
      y = COORDS(2)-cy
      r = SQRT(x**2+y**2)
      theta = atan2(y,x)


      /* write(*,*) COORDS(1), COORDS(2), COORDS(3)
      write(*,*) r, theta */
      If (JDOF == 1) THEN
        U(1) = (dk/(2*G))*(sqrt(r/(2*pi)))*cos(theta/2)
        U(1) = U(1) * (mu - 1 + 2*sin(theta/2)**2)
        write(*,*) "Node", NODE
        write(*,*) COORDS(1), COORDS(2)
        write(*,*) r, theta
        write(*,*) U(1)
      ELSE
        U(1) = (dk/(2*G))*(sqrt(r/(2*pi)))*sin(theta/2)
        U(1) = U(1) * (mu + 1 - 2*cos(theta/2)**2)
        if (NODE == 21) THEN
          U(1) = -U(1)
        END IF
        write(*,*) U(1)
      END IF

      RETURN
      END
