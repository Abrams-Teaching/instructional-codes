C     Maxpart = Maximum Number Of Particles
C     Maxkap  = Maximum NumberOf K-vectors
C     Npart   = Number Of Particles (Ncell**3)
C     Ncell   = Number Of Cells
C     Kmax    = Maximum Number Ofk-Vectors
C     Rx/Y/Z  = Coordinates
C     Z       = Charge
C     Alpha   = Ewald Parameter
C     Box     = Boxsize

      Integer Maxpart,Maxkap,Npart,Ncell,Kmax

      Parameter (Maxpart = 10000)
      Parameter (Maxkap  = 30)

      Double precision Rx(Maxpart),Ry(Maxpart),Rz(Maxpart),
     &                 Z(Maxpart),Alpha,Box,Ibox,Hbox,Hbox2,
     &                 Piesq,Twopi

      Common /Ewww/ Rx,Ry,Rz,Z,Alpha,Piesq,Twopi,Box,Ibox,
     &              Hbox,Hbox2,Npart,Ncell,Kmax
