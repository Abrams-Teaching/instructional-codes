Cccccccccccccccccccccccccc
C     Mdvelo.Inc         C
Cccccccccccccccccccccccccc
C
C     Vxx/Vyy = Velocity
C     Rxx/Ryy = Position
C     Rxf/Ryf = Old Position
C     Fxx/Fyy = Force

      Double Precision Vxx(Maxatom),Vyy(Maxatom),
     &                 Rxx(Maxatom),Ryy(Maxatom),
     &                 Fxx(Maxatom),Fyy(Maxatom),
     &                 Rxf(Maxatom),Ryf(Maxatom)

      Common /Mdv1/    Vxx,Vyy,Rxx,Ryy,Fxx,Fyy,Rxf,Ryf
