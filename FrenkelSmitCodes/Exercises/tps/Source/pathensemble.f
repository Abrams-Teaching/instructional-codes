      Subroutine Pathensemble
      Implicit None

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'pathtopo.inc'
      Include 'traject.inc'

Cccccccccccccccccccccccccccccccccccccccc
C     Generates An Initial Path        C
C     Input: Molecule At An Arbitrary  C
C     Position                         C
Cccccccccccccccccccccccccccccccccccccccc

      Integer          I,J,Ncycle,Nminit

      Logical          Ldisk

      Double Precision Ran_Uniform,Pshift

Cccccccccccccccccccccccccccccccccc
C     Read Path-Sampling Input   C
Cccccccccccccccccccccccccccccccccc

      Open(32,File='input2')

      Read(32,*)
      Read(32,*) Ldisk,Pshift,Deltaphi
      Read(32,*)
      Read(32,*) Ncycle,Nminit,Lumbrella

      Close(32)

      Write(6,*)
      Write(6,*) 'Prob. Shifting                         : ',Pshift
      Write(6,*) 'Number Of Slices                       : ',Nslice
      Write(6,*) 'Number Of Mc Cycles                    : ',Ncycle
      Write(6,*) 'Number Of Init. Cycles                 : ',Nminit
      Write(6,*) 'Maximum Rotation                       : ',Deltaphi
      Write(6,*)
      Write(6,*)
      Write(6,*)

Ccccccccccccccccccccccccccc
C     Set Stuff To Zero   C
Ccccccccccccccccccccccccccc

      Drift1 = 0.0d0
      Drift2 = 0.0d0
      Shift1 = 0.0d0
      Shift2 = 0.0d0

      Disp1 = 0.0d0
      Disp2 = 0.0d0
            
      Do J=1,Maxtraject

         La_Old(J) = .False.
         Lb_Old(J) = .False.
         Eeold(J)  = 0.0d0

         Do I=1,Maxatom
            Xxold(I,J) = 0.0d0
            Yyold(I,J) = 0.0d0
                      
            Vxold(I,J) = 0.0d0
            Vyold(I,J) = 0.0d0
         Enddo
      Enddo
      
      If(Ldisk) Then

Ccccccccccccccccccccccccccccc
C     Read Stuff From Disk  C
Ccccccccccccccccccccccccccccc

         Call Readpath

      Else

Cccccccccccccccccccccccccccc
C     Initial Coordinates  C
Cccccccccccccccccccccccccccc
      
         Call Init_Path

         Write(6,*)
         Write(6,*) 'Succeeded In Finding A Path'
         Write(6,*)
      
         Call Writepath
      Endif

Ccccccccccccccccccccccccccccccccccc
C     Loop Over All Cycles        C
C     Set Calc. Averages To Zero  C
Ccccccccccccccccccccccccccccccccccc

      If(Lumbrella) Then
         Call Sample_Umbrella(1)
      Else
         Call Sample_Avhb(1)
      Endif

      Do I=1,Ncycle
      
         If(Ran_Uniform().Lt.Pshift) Then
            Call Shifting
         Else
            Call Shooting
         Endif
      
Ccccccccccccccccccccccccc
C     Sample Averages   C
Ccccccccccccccccccccccccc            

         If(I.Ge.Nminit) Then
            If(Mod(I,3).Eq.0) Then
               If(Lumbrella) Then
                  Call Sample_Umbrella(2)
               Else
                  Call Sample_Avhb(2)
               Endif
            Endif
         
            If((Mod(I,1000).Eq.0).Or.(I.Eq.Ncycle)) Then

               If(Lumbrella) Then
                  Call Sample_Umbrella(3)
               Else
                  Call Sample_Avhb(3)
               Endif
            Endif
         Endif

Ccccccccccccccccccccccccccccccc
C     Write Info To Screen    C
Ccccccccccccccccccccccccccccccc

         If((Mod(I,1000).Eq.0).Or.(I.Eq.Ncycle).Or.(I.Eq.1)) Then
            Call Writepath

            Write(6,*)
            Write(6,*)
            Write(6,*) 'Cycle                                  : ',
     &           I
            Write(6,*) 'Energy                                 : ',
     &           Eeold(1)

            If(Drift2.Gt.0.5d0) 
     &           Write(6,*) 'Average Energy Drift                   : ',
     &           Drift1/Drift2

            If(Shift2.Gt.0.5d0)
     &           Write(6,*) 'Fraction Acc. Shifting                 : ',
     &           Shift1/Shift2

            If(Disp2.Gt.0.5d0)
     &           Write(6,*) 'Fraction Acc. Shooting                 : ',
     &           Disp1/Disp2
         Endif
      Enddo

      Write(6,*)
      Write(6,*)
      Write(6,*)

      Return
      End
