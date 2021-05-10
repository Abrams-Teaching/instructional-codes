      Subroutine Md(Ntstart,Nsign,Nintegrate,Av1,Lpath)
      Implicit None

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Molecular Dynamics (Velocity-Verlet Integrator)             C
C                                                                 C
C     Ntstart = Start Of The Trajectory In The Xxtra-Variabele    C
C     Nsign   = Sign Of The Trajectory (-1 Or 1)                  C
C               The Stored Velocities Have To Be Multiplied       C
C               By Nsign                                          C
C     Av1     = Energy Drift During The Simulation                C
C     Lpath   = Do We Use Path-Sampling ?                         C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'mdvelo.inc'
      Include 'traject.inc'
 
      Double Precision Ukin,Utot,Upot,Av1,Av2,Av0,
     &                 Hstep
     
      Integer          I,Iii,Iiii,Nsign,Ntstart,Nintegrate
      Logical          In_A,In_B,Lpath

Ccccccccccccccccccccccccccccc              
C     Loop Over All Steps   C
C     Calc Initial Force    C
Ccccccccccccccccccccccccccccc
 
      Call Force(Utot,Upot,Ukin)
                  
      Av0   = Utot
      Av1   = 0.0d0
      Av2   = 0.0d0
      Hstep = 0.5d0*Tstep

Cccccccccccccccccccccccccccc
C     Initial Point Start  C
Cccccccccccccccccccccccccccc

      If(Lpath) Then
         La(Ntstart)    = In_A()
         Lb(Ntstart)    = In_B()
         Eetra(Ntstart) = Utot

         If(Nsign.Eq.1) Then
            Do I=1,Natom
               Xxtra(I,Ntstart) = Rxx(I)
               Yytra(I,Ntstart) = Ryy(I)
                                          
               Vxtra(I,Ntstart) = Vxx(I)
               Vytra(I,Ntstart) = Vyy(I)
            Enddo
         Else
            Do I=1,Natom
               Xxtra(I,Ntstart) = Rxx(I)
               Yytra(I,Ntstart) = Ryy(I)
                                          
               Vxtra(I,Ntstart) = -Vxx(I)
               Vytra(I,Ntstart) = -Vyy(I)
            Enddo
         Endif
      Else
         Write(6,*)
         Write(6,*)
         Write(6,*) 'Start Of The Md Simulation'
         Write(6,*)
         Write(6,*) 'Total Energy                           : ',Utot
         Write(6,*) 'Potential Energy                       : ',Upot
         Write(6,*) 'Kinetic Energy                         : ',Ukin
         Write(6,*)
         Write(6,*)

         Call Sample_Dist(1)
         Call Sample_Ct(1,Utot)
      Endif

Ccccccccccccccccccccccccccccccccccccccccccc
C     Return If Nothing Has To Be Done    C
Ccccccccccccccccccccccccccccccccccccccccccc

      If(Nintegrate.Eq.0) Return

Cccccccccccccccccccccccccccccc
C     Loop Over All Cycles   C
C     Loop Over Nshort       C
Cccccccccccccccccccccccccccccc

      Do Iii = 1,Nintegrate
         Do Iiii = 1,Nshort
            
            Do I=1,Natom
               Vxx(I) = Vxx(I) + Hstep*Fxx(I)
               Vyy(I) = Vyy(I) + Hstep*Fyy(I)
                                             
               Rxx(I) = Rxx(I) + Tstep*Vxx(I)
               Ryy(I) = Ryy(I) + Tstep*Vyy(I)
            Enddo

            Call Force(Utot,Upot,Ukin)
            
            Ukin = 0.0d0

            Do I=1,Natom
               Vxx(I) = Vxx(I) + Hstep*Fxx(I)
               Vyy(I) = Vyy(I) + Hstep*Fyy(I)
               
               Ukin = Ukin + Vxx(I)**2 + Vyy(I)**2
            Enddo
         Enddo
         
Cccccccccccccccccccccccccccc
C     Temperature..Energy  C
Cccccccccccccccccccccccccccc
 
         Ukin = 0.5d0*Ukin
         Utot = Upot + Ukin
         Av1  = Av1 + Dabs((Av0-Utot)/Av0)
         Av2  = Av2 + 1.0d0

Cccccccccccccccccccccccccccccccccccccc
C     Sample Correlation Functions   C
Cccccccccccccccccccccccccccccccccccccc

         If(.Not.Lpath.And.Iii.Gt.25000) Then
            Call Sample_Dist(2)
            Call Sample_Ct(2,Utot)

            If(Iii.Gt.50000.And.Mod(Iii,100000).Eq.0) Then
               Call Sample_Dist(3)
               Call Sample_Ct(3,Utot) 
            Endif
         Endif
         
Ccccccccccccccccccccccccccccccc
C     Used For Debugging      C
Ccccccccccccccccccccccccccccccc

         If(Lscreen.And.Mod(Iii,250).Eq.0) Then
            Write(6,'(A,3(1x,E20.10))') 'Utot Ukin Upot ',
     &           Utot,Ukin,Upot
         Endif
         
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Store Particle Positions En Velocities;           C
C     Check If The Particle Is In A Or B At Slice Iiii  C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         Iiii = Ntstart + Iii*Nsign

         If((Iiii.Gt.Maxtraject.Or.Iiii.Lt.1).And.Lpath) Then
            Write(6,*) 'Traject Overflow !!! ',Iiii
            Call Exitt(2)
         Endif

         If(Lpath) Then

            La(Iiii)    = In_A()
            Lb(Iiii)    = In_B()
            Eetra(Iiii) = Utot 

            If(Nsign.Eq.1) Then
               Do I=1,Natom
                  Xxtra(I,Iiii) = Rxx(I)
                  Yytra(I,Iiii) = Ryy(I)
                           
                  Vxtra(I,Iiii) = Vxx(I)
                  Vytra(I,Iiii) = Vyy(I)
               Enddo
            Else
              Do I=1,Natom
                  Xxtra(I,Iiii) = Rxx(I)
                  Yytra(I,Iiii) = Ryy(I)
                           
                  Vxtra(I,Iiii) = -Vxx(I)
                  Vytra(I,Iiii) = -Vyy(I)
               Enddo 
            Endif
         Endif
      Enddo

Ccccccccccccccccccc
C     The End     C
Ccccccccccccccccccc

      If(Av2.Gt.0.5d0) Av1 = Av1/Av2
      
      If(.Not.Lpath) Then
         Write(6,*)
         Write(6,*)
         Write(6,*) 'End Of The Md Simulation'
         Write(6,*)
         Write(6,*) 'F Total Energy                         : ',Utot
         Write(6,*) 'F Drift                                : ',Av1
         Write(6,*) 'F Kinetic Energy                       : ',Ukin
         Write(6,*) 'F Potential                            : ',Upot
         Write(6,*)
         Write(6,*) 

         Call Sample_Dist(3)
         Call Sample_Ct(3,Utot)
      Endif

      Return
      End
