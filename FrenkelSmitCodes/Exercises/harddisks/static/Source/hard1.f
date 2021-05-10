      Program Hard1
      Implicit None

      Include 'system.inc'

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Hard Disks On A Square; Random Placement  C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer          Nstep,I,J,Iii,Kkk,Sstmm
      Double Precision Ran_Uniform,R2,Att1,Att2,M1
      
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

      Write(*,*) 'How Many Cycles ?             (Example: 1000     )'
      Read(*,*)  Nstep
      
      Write(*,*) 'How Many Particles ?          (Always 2 < I < 80 )'
      Read(*,*)  Npart

      If(Npart.Gt.80.Or.Npart.Lt.2) Stop

      Call Sample(1)

      Att1 = 0.0d0
      Att2 = 0.0d0

      Do Iii=1,Nstep
         Do Kkk=1,1000

            Att1 = Att1 + 1.0d0 

            Do I=1,Npart

Ccccccccccccccccccccc
C     New Particle  C
Ccccccccccccccccccccc

               X(I) = Ran_Uniform()*10.0d0
               Y(I) = Ran_Uniform()*10.0d0

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Check Distances With Previous Placed Particles  C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               Do J=1,(I-1)
                  R2 = (X(I)-X(J))**2 + (Y(I)-Y(J))**2

                  If(R2.Lt.1.0d0) Goto 111
               Enddo
            Enddo

Ccccccccccccccccccccccccccccccccc
C     No Overlaps Occured !!!   C
Ccccccccccccccccccccccccccccccccc

            Att2 = Att2 + 1.0d0

            Call Sample(2)

 111        Continue

         Enddo
      Enddo

      If(Att2.Ge.0.5d0) Call Sample(3)
      Write(6,*) 'Fraction Succes : ',Att2/Att1

      Stop
      End
