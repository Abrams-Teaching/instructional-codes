      Subroutine Disp_Mono(Dispa,Dispb)
      Implicit None

      Include 'maxarray.inc'
      Include 'system.inc'

      Logical Ltake
      Integer Ipart
      Double Precision Unew,Uold,Ran_Uniform,
     &     Xk,Yk,Zk,dispa,dispb

Cccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Select Random Particle And Calculate Delta U   C
C     Accept/Reject This Trial Move                  C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Ltake = Lchain
      Ipart = 1     + Idint(Ran_Uniform()*Dble(Npart))
      Dispb = Dispb + 1.0d0

      Xk = Xx(Ipart)
      Yk = Yy(Ipart)
      Zk = Zz(Ipart)

      Call Ener_Mono(Uold,Ipart,Xk,Yk,Zk,Ltake)

      Xk = Xx(Ipart) + (Ran_Uniform()-0.5d0)*Dispmax
      Yk = Yy(Ipart) + (Ran_Uniform()-0.5d0)*Dispmax
      Zk = Zz(Ipart) + (Ran_Uniform()-0.5d0)*Dispmax

Ccccccccccccccccccccccccc
C     Back In The Box   C
Ccccccccccccccccccccccccc

      Xk = Xk - Box*(Dble(Idint(Xk*Ibox + 999.0d0) - 999))
      Yk = Yk - Box*(Dble(Idint(Yk*Ibox + 999.0d0) - 999))
      Zk = Zk - Box*(Dble(Idint(Zk*Ibox + 999.0d0) - 999))

      Call Ener_Mono(Unew,Ipart,Xk,Yk,Zk,Ltake)

Ccccccccccccccccccccccc
C     Accepted ????   C
Ccccccccccccccccccccccc

      If(Ran_Uniform().Lt.Dexp(-Beta*(Unew-Uold))) Then
         Xx(Ipart) = Xk
         Yy(Ipart) = Yk
         Zz(Ipart) = Zk

         Usim  = Usim  + Unew - Uold
         Dispa = Dispa + 1.0d0
      Endif

      Return
      End
