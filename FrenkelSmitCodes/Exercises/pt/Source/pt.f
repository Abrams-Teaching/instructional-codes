      Program Pt
      Implicit None

Cccccccccccccccccccccccccccccccccccccccccccccccc
C     Parallel Tempering For A Simple System   C
C     Written By Thijs J.H. Vlugt              C
Cccccccccccccccccccccccccccccccccccccccccccccccc

      Include 'system.inc'

      Integer Sstmm,I,J,Itemp,Jtemp,Nstep,Imol,Kkk,part
      Double Precision M1,Rho,Del,U,Dispa(9),
     &     Dispb(9),Swapa(9),Swapb(9),Ibox,
     &     Ran_Uniform,Ua(9),Ub(9),Xa(9),Xb(9),
     &     Ppt,Temp(9),Beta(9),Rxold,Ryold

Ccccccccccccccccccccccccccccccccccccc
C     Set Random Number Generator   C
Ccccccccccccccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1) 
 
Cccccccccccccccccccccccccccc
C     Read In Data         C
Cccccccccccccccccccccccccccc
 
      Read(21,*)
      Read(21,*) Rho,Ppt,Nstep,Ntemp

      If(Ntemp.Lt.1.Or.Ntemp.Gt.9) Stop

      Read(21,*)
      Read(21,*) (Temp(I),I=1,Ntemp)

      Do I=1,Ntemp
         Beta(I) = 1.0d0/Temp(I)
      Enddo

      Box  = Dsqrt(9.0d0/Rho)
      Del  = Box/4.0d0
      Ibox = 1.0d0/Box
      part = 0

Cccccccccccccccccccccccccccc
C     Generate A Lattice   C
C     Of 9 Particles       C
Cccccccccccccccccccccccccccc

      Do I=1,3
         Do J=1,3
            part = part + 1
            Do Itemp=1,Ntemp
               Rx(I,Itemp) = Dble(I)*Del
               Ry(I,Itemp) = Dble(J)*Del
            Enddo
         Enddo
      Enddo

      Do I=1,9
         Uold(I)  = 0.0d0
         Dispa(I) = 0.0d0
         Dispb(I) = 0.0d0
         Swapa(I) = 0.0d0
         Swapb(I) = 0.0d0
         Ua(I)    = 0.0d0
         Ub(I)    = 0.0d0
         Xa(I)    = 0.0d0
         Xb(I)    = 0.0d0
      Enddo

Cccccccccccccccccccccccccccccccccccccc
C     Calculate The Initial Energy   C
Cccccccccccccccccccccccccccccccccccccc

      Do I=1,Ntemp
         Call Ener(U,I)

         Uold(I) = U
      Enddo

      Call Sample_Energy(1)
      Call Sample_Xcoord(1)

      Write(6,*) 'Density            : ',Rho
      Write(6,*) 'Boxsize            : ',Box
      Write(6,*) 'Fraction Exchanges : ',Ppt
      Write(6,*) 'Number Of Steps    : ',Nstep
      Write(6,*) 'Number Of Temp.    : ',Ntemp

      Do I=1,Ntemp
         Write(6,*) 'Temp. (Energy)     : ',
     &        Temp(I),Uold(I)
      Enddo

Cccccccccccccccccccccccccccccc
C     Loop Over All Cycles   C
Cccccccccccccccccccccccccccccc

      Do I=1,Nstep
         Do Kkk=1,1000*Ntemp

Cccccccccccccccccccccccccccccccccccccccc
C     Sample Energy And X-Coordinate   C
Cccccccccccccccccccccccccccccccccccccccc

            If(I.Gt.10) Then
               Call Sample_Energy(2)
               Call Sample_Xcoord(2)

               Do Itemp=1,Ntemp
                  Ua(Itemp) = Ua(Itemp) + Uold(Itemp)
                  Ub(Itemp) = Ub(Itemp) + 1.0d0

                  Xa(Itemp) = Xa(Itemp) + Rx(1,Itemp)*Ibox
                  Xb(Itemp) = Xb(Itemp) + 1.0d0
               Enddo
            Endif

Ccccccccccccccccccccccccccccccccccccccc
C     Select A Trial Move At Random   C
Ccccccccccccccccccccccccccccccccccccccc

            If(Ran_Uniform().Lt.Ppt.And.Ntemp.Gt.1) Then

Ccccccccccccccccccccccccccccccccccccccc
C     Exchange Trial Move             C
C     Select Temperatures At Random   C
Ccccccccccccccccccccccccccccccccccccccc

C     Start Modification
C     End   Modification

            Else

Cccccccccccccccccccccccccccccccccccccccccc
C     Particle Displacement              C
C     Select System/Particle At Random   C
Cccccccccccccccccccccccccccccccccccccccccc

               Itemp        = 1 + Idint(Dble(Ntemp)*Ran_Uniform())
               Imol         = 1 + Idint(9.0d0*Ran_Uniform())
               Rxold        = Rx(Imol,Itemp)
               Ryold        = Ry(Imol,Itemp)
               Dispb(Itemp) = Dispb(Itemp) + 1.0d0

               Rx(Imol,Itemp) = Rx(Imol,Itemp) + 
     &              (Ran_Uniform()-0.5d0)*0.05d0

               Ry(Imol,Itemp) = Ry(Imol,Itemp) + 
     &              (Ran_Uniform()-0.5d0)*0.05d0
               
Cccccccccccccccccccccccccccccccc
C     Check For Acceptance     C
C     Check Boundaries First   C
Cccccccccccccccccccccccccccccccc

               If(Min(Rx(Imol,Itemp),Ry(Imol,Itemp)).Gt.0.0d0.And.
     &            Max(Rx(Imol,Itemp),Ry(Imol,Itemp)).Lt.Box) Then

                  Call Ener(U,Itemp)

                  If(Ran_Uniform().Lt.
     &                 Dexp(-Beta(Itemp)*(U-Uold(Itemp)))) Then

                     Uold(Itemp)  = U
                     Dispa(Itemp) = Dispa(Itemp) + 1.0d0
                  Else
                     Rx(Imol,Itemp) = Rxold
                     Ry(Imol,Itemp) = Ryold
                  Endif
               Else
                  Rx(Imol,Itemp) = Rxold
                  Ry(Imol,Itemp) = Ryold
               Endif
            Endif
         Enddo
      Enddo

Cccccccccccccccccccccccccccccccccc
C     Write Info To The Screen   C
Cccccccccccccccccccccccccccccccccc
  
      Call Sample_Energy(3)
      Call Sample_Xcoord(3)

      Write(6,*)
      Write(6,*)
      Write(6,*)
      Write(6,*) 'Average Energy'

      Do I=1,Ntemp
         Write(6,*) 'Temp : ',I,' ',
     &        Ua(I)/Ub(I)
      Enddo

      Write(6,*)
      Write(6,*) 'Average X Coordinate'

      Do I=1,Ntemp
         Write(6,*) 'Temp : ',I,' ',
     &        Xa(I)/Xb(I)
      Enddo

      Write(6,*)
      Write(6,*) 'Fraction Accepted Displacements:'
      
      Do I=1,Ntemp
         Write(6,*) 'Temp : ',I,' ',
     &        Dispa(I)/Dispb(I)
      Enddo

      Write(6,*)
      Write(6,*) 'Fraction Accepted Exchanges:'

      Do I=1,Ntemp-1
         If(Swapb(I).Gt.0.5d0) 
     &        Write(6,*) I,I+1,'  ',
     &        Swapa(I)/Swapb(I)
      Enddo

Ccccccccccccccccccccccccccccccccc
C     End Of The Program        C
Ccccccccccccccccccccccccccccccccc
 
      Stop
      End
