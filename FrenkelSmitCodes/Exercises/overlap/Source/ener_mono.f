      Subroutine Ener_Mono(U,Ipart,Xk,Yk,Zk,Ltake)
      Implicit None

      Include 'maxarray.inc'
      Include 'system.inc'

      Logical Ltake
      Integer I,Ipart
      Double Precision U,Xk,Yk,Zk,Dx,Dy,Dz,R2

Cccccccccccccccccccccccccccccccccccccccccccccccc
C     Calculate Total Energy Of a Monomer      C
C     Exclude Particle Ipart                   C
C                                              C
C     Ltake    = Interactions With Polymer ?   C
C     U        = Energy                        C
C     Xk/Yk/Zk = Coordinates Of Particle       C
Cccccccccccccccccccccccccccccccccccccccccccccccc

      U = 0.0d0

      Do I=1,Npart
         If(I.Ne.Ipart) Then
                   
            Dx = Xx(I) - Xk
            Dy = Yy(I) - Yk
            Dz = Zz(I) - Zk

            Dx = Dx - Box*Dble(Idint(Dx*Ibox + 999.5d0) - 999)
            Dy = Dy - Box*Dble(Idint(Dy*Ibox + 999.5d0) - 999)
            Dz = Dz - Box*Dble(Idint(Dz*Ibox + 999.5d0) - 999)

            R2 = Dx**2 + Dy**2 + Dz**2

            If(R2.Lt.1.0d0) 
     &           U  = U + Alpha*((Dsqrt(R2) - 1.0d0)**2)
         Endif
      Enddo

Ccccccccccccccccccccccccccccccccccc
C     Interactions With Polymer   C
Ccccccccccccccccccccccccccccccccccc

      If(Ltake) Then
         Do I=1,Nlength
            Dx = X(I) - Xk
            Dy = Y(I) - Yk
            Dz = Z(I) - Zk

            Dx = Dx - Box*Dble(Idint(Dx*Ibox + 999.5d0) - 999)
            Dy = Dy - Box*Dble(Idint(Dy*Ibox + 999.5d0) - 999)
            Dz = Dz - Box*Dble(Idint(Dz*Ibox + 999.5d0) - 999)

            R2 = Dx**2 + Dy**2 + Dz**2

            If(R2.Lt.1.0d0) 
     &           U  = U + Alpha*((Dsqrt(R2) - 1.0d0)**2)
         Enddo
      Endif

      Return
      End
