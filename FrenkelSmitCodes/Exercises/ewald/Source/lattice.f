      Subroutine Lattice
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccc
C     Setup Simple Cubic Lattice   C
Cccccccccccccccccccccccccccccccccccc

      Integer Ix,Iy,Iz,I,Ip,In

      I  = 0
      Ip = 0
      In = 0

      Do Ix=1,Ncell
         Do Iy=1,Ncell
            Do Iz=1,Ncell

               I = Iz + (Iy-1)*Ncell + (Ix-1)*Ncell*Ncell

               If(I.Gt.Maxpart) Then
                  Write(6,*) 'Increase Maxpart !!!'
                  Stop
               Endif

C               Rx(I) = Dble(Ix-1) + 0.5d0
C               Ry(I) = Dble(Iy-1) + 0.5d0
C               Rz(I) = Dble(Iz-1) + 0.5d0
               Rx(I) = Dble(Ix-1) + 0.5d0 + 0.6*(RAND()-0.5)
               Ry(I) = Dble(Iy-1) + 0.5d0 + 0.6*(RAND()-0.5)
               Rz(I) = Dble(Iz-1) + 0.5d0 + 0.6*(RAND()-0.5)

               If(Mod((Ix+Iy+Iz),2).Eq.0) Then
                  Z(I) = 1.0d0
                  Ip   = Ip + 1
               Else
                  Z(I) = -1.0d0
                  In   = In + 1
               Endif
            Enddo
         Enddo
      Enddo

Cccccccccccccccccccccccccc
C     Define Variables   C
Cccccccccccccccccccccccccc

      Npart = Ncell*Ncell*Ncell
      Box   = Dble(Ncell)
      Hbox  = 0.5d0*Box
      Ibox  = 1.0d0/Box
      Hbox2 = Hbox**2
      Twopi = 8.0d0*Datan(1.0d0)
      Piesq = 1.0d0/Dsqrt(4.0d0*Datan(1.0d0))

Cccccccccccccccccccc
C     Write Info   C
Cccccccccccccccccccc

      Write(6,*) 'Number Of Cells      : ',Ncell
      Write(6,*) 'Number Of Ions       : ',Npart
      Write(6,*) 'Number Of Pos. Ions  : ',Ip
      Write(6,*) 'Number Of Neg. Ions  : ',In
      Write(6,*) 'Alpha                : ',Alpha
      Write(6,*) 'Kmax                 : ',Kmax
      Write(6,*) 'Boxsize              : ',Box
      Write(6,*) 'Recommended Alpha    : ',3.2d0/Hbox
      Write(6,*)
      
      If(Kmax.Lt.1.Or.Kmax.Gt.Maxkap) Then
         Write(6,*) 'Error Kmax !!!'
         Stop
      Endif

      If(Ip.Ne.In) Then
         Write(6,*) 'Please Choose Ncell Even !!!'
         Stop
      Endif

      Do I=1,Npart
         Write(22,'(4E15.5)') Rx(I),Ry(I),Rz(I),Z(I)
      Enddo

      Return
      End
