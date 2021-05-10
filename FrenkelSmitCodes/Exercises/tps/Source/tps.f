      Program Tps
      Implicit None

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'mdvelo.inc'
      Include 'pathtopo.inc'
      Include 'traject.inc'
      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Transition Path Sampling For Lj Dumbell                C
C                                                            C
C     See: Dellago Etal, J.Chem.Phys., 1999, 110(14),6617    C
C                                                            C
C     Written By Thijs J.H. Vlugt, 12-7-2000                 C
C                                                            C
C     Main Program.. Initialize Rng                          C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      Integer          I,J,Sstmm
      Double Precision Du,M1,Ran_Uniform

      J = Sstmm()

      Write(6,*)
      Write(6,*)
      Write(6,*) 'Path Ensemble Sampling; Written By Thijs J.H. Vlugt'
      Write(6,*)
      Write(6,*) 'Initial Seed                           : ',J
      
Cccccccccccccccccccccccccc
C     Initialize Rng     C
Cccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*J),1000))
      
      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0 
      
      Call Genrand(M1)

      J = 10 + Idint(Ran_Uniform()*1000.0d0)

      Do I=1,J
         M1 = Ran_Uniform()
      Enddo

      Write(6,*) 'Random Number                          : ',
     &     Ran_Uniform()
      Write(6,*)
      Write(6,*)
      Write(6,*)

Cccccccccccccccccccc
C     Read Input   C
Cccccccccccccccccccc
      
      Call Readinput

Ccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Do Either Pathensemble Or Normal Nve/Md     C
Ccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      If(Lensemble) Then
         Call Pathensemble
      Else
         Call Md(1,1,Nstep,Du,.False.)
      Endif

      Call Exitt(1)
      End
