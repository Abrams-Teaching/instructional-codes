      Subroutine Mdloop
      Implicit None
 
      Include 'system.inc'

      Integer I,J
      Double Precision Ran_Vel,A,C1,C2

      Call Sample(1,0)

      C1 = 0.0d0
      C2 = 0.0d0

      Do I=1,Ncycle
         Do J=1,1000

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Generate Initial Coordinates              C
C                                               C
C     Perform 1000*Ncycle Md Simulations With   C
C     Different Initial Conditions              C
C                                               C
C     Xpos   = Starting Position                C
C     Vpos   = Starting Velocity                C
C     Thetan = Storage For Initial Velocity     C
C     Qstar  = Place Of The Dividing Surface    C
C                                               C
C     To Program:                               C
C                                               C
C     -Generate Initial Position/Velocity       C
C     -Integrate The Equations Of Motion By A   C
C      Function Call To Subroutine Integrate    C
C     -Beware That In Integrate The Subroutine  C
C      Sample Is Called !!                      C
C                                               C
C     The Averaged Energy Drift (Over All       C
C     Simulations) Equals C1/C2                 C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

C Start Modification

C End Modification

         Enddo
      Enddo

      Call Sample(3,0)

      Write(6,*) 'Av. Energy Drift      :',C1/C2

      Return
      End
