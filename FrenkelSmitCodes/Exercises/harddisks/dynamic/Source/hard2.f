      Program Hard2
      Implicit None

      Include 'system.inc'

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Hard Disks On A Square; Random Placement  C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer Nstep,I,J,Iii,Kkk,Ipart,Pp,
     &     Ninit,Sstmm
      Double Precision Ran_Uniform,R2,Att1,Att2,Xn,Yn,
     &     Dispmax,M1
      
Cccccccccccccccccccccccccccccccccccccccccccc
C     Initialize Random Number Generator   C
Cccccccccccccccccccccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)

Cccccccccccccccccccccccccccccccccccccccccccc
C     Read Info From Disk                  C
Cccccccccccccccccccccccccccccccccccccccccccc

      Write(*,*) 'How Many Cycles       ?        (Example: 1000     )'
      Read(*,*)  Nstep

      Write(*,*) 'How Many Init. Cycles ?        (Example: 100      )'
      Read(*,*)  Ninit
      
      Write(*,*) 'How Many Particles    ?        (Always 2 < I < 80 )'
      Read(*,*)  Npart

      Write(*,*) 'Maximum Displacement  ?        (Disp > 0          )'
      Read(*,*)  Dispmax

      If(Npart.Gt.80.Or.Npart.Lt.2) Stop

      Call Sample(1)

      Att1 = 0.0d0
      Att2 = 0.0d0
      Pp   = 0

Ccccccccccccccccccccccccccccccccccc
C     Put Particles On A Lattice  C
Ccccccccccccccccccccccccccccccccccc

      Kkk = 0

      Do I=1,9
         Do J=1,9

            Kkk    = Kkk + 1
            X(Kkk) = Dble(I)*1.1d0
            Y(Kkk) = Dble(J)*1.1d0
         Enddo
      Enddo

      Do Iii=1,Nstep
         Do Kkk=1,1000

            Att1 = Att1 + 1.0d0 

Ccccccccccccccccccccccccccccccccccc
C     Select Particle At Random   C
Ccccccccccccccccccccccccccccccccccc

            Ipart = 1        + Idint(Ran_Uniform()*Dble(Npart))
            Xn    = X(Ipart) + (Ran_Uniform()-0.5d0)*Dispmax
            Yn    = Y(Ipart) + (Ran_Uniform()-0.5d0)*Dispmax

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Always Reject When Particle Is Out Of The Box   C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            If(Min(Xn,Yn).Lt.0.0d0.Or.Max(Xn,Yn).Gt.10.0d0) Goto 111

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     See If There Is An Overlap With Any Other Particle  C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            Do I=1,Npart
               If(I.Ne.Ipart) Then

                  R2 = (Xn-X(I))**2 + (Yn-Y(I))**2

                  If(R2.Lt.1.0d0) Goto 111
               Endif
            Enddo

Cccccccccccccccccccccccccccccccccc
C     No Overlaps, So Accepted   C
Cccccccccccccccccccccccccccccccccc

            Att2     = Att2 + 1.0d0
            X(Ipart) = Xn
            Y(Ipart) = Yn

 111        Continue

            If(Mod(Kkk,5).Eq.0.And.Iii.Gt.Ninit) Then
               Call Sample(2)
               Pp = Pp + 1
            Endif
         Enddo
         If(Mod(Iii,5).Eq.0) Call Writepdb
      Enddo

      If(Pp.Gt.4) Call Sample(3)
      Write(6,*) 'Fraction Succes : ',Att2/Att1

      Stop
      End
