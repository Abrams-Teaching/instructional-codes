      Subroutine Ener_Tot(U)
      Implicit None

      Include 'maxarray.inc'
      Include 'system.inc'

      Integer I,J,Ii
      Double Precision U,Dx,Dy,Dz,R2

Cccccccccccccccccccccccccccccccccccc
C     Calculate Total Energy (U)   C
Cccccccccccccccccccccccccccccccccccc

      U = 0.0d0

      Do I=1,Npart-1
         Do J=(I+1),Npart

            Dx = Xx(I) - Xx(J)
            Dy = Yy(I) - Yy(J)
            Dz = Zz(I) - Zz(J)

            Dx = Dx - Box*Dble(Idint(Dx*Ibox + 999.5d0) - 999)
            Dy = Dy - Box*Dble(Idint(Dy*Ibox + 999.5d0) - 999)
            Dz = Dz - Box*Dble(Idint(Dz*Ibox + 999.5d0) - 999)

            R2 = Dx**2 + Dy**2 + Dz**2

            If(R2.Lt.1.0d0) 
     &           U  = U + Alpha*((Dsqrt(R2) - 1.0d0)**2)
         Enddo
      Enddo

Ccccccccccccccccccccccccccccc
C     Polymer Interactions  C
Ccccccccccccccccccccccccccccc

      If(Lchain) Then

Cccccccccccccccccccccccccccc
C     Intra Chain          C
C     No 1-2 Interactions  C
C     No Pbc's Here        C
Cccccccccccccccccccccccccccc

         If(Nlength.Ge.3) Then
            Do I=1,Nlength-2
               Do Ii = I+2,Nlength

                  Dx = X(Ii) - X(I)
                  Dy = Y(Ii) - Y(I)
                  Dz = Z(Ii) - Z(I)

                  R2 = Dx**2 + Dy**2 + Dz**2

                  If(R2.Lt.1.0d0) 
     &                 U  = U + Alpha*((Dsqrt(R2) - 1.0d0)**2)
               Enddo
            Enddo
         Endif

Ccccccccccccccccccccccccccccccccccccc
C     Particle-Chain Interactions   C
Ccccccccccccccccccccccccccccccccccccc

         Do I=1,Npart
            Do Ii=1,Nlength
               Dx = X(Ii) - Xx(I)
               Dy = Y(Ii) - Yy(I)
               Dz = Z(Ii) - Zz(I)

               Dx = Dx - Box*Dble(Idint(Dx*Ibox + 999.5d0) - 999)
               Dy = Dy - Box*Dble(Idint(Dy*Ibox + 999.5d0) - 999)
               Dz = Dz - Box*Dble(Idint(Dz*Ibox + 999.5d0) - 999)

               R2 = Dx**2 + Dy**2 + Dz**2

               If(R2.Lt.1.0d0) 
     &              U  = U + Alpha*((Dsqrt(R2) - 1.0d0)**2)
            Enddo
         Enddo
      Endif

      Return
      End
