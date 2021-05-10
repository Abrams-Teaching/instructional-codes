      Subroutine Init_Path
      Implicit None

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'mdvelo.inc'
      Include 'traject.inc'

      Double Precision Vxforw(Maxatom),Vyforw(Maxatom),
     &                 Du,Ran_Gauss,Ukin,Factor,Mx,My,
     &                 Udes

      Integer          I,J,Islice
      Logical          L_Ready,Lslice

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Generate An Initial Path Form A Gaussian   C
C     Scale To The Desired Total Energy          C
C     Set Initial Momenta To Zero                C
C     Momentum Is Not Conserved Here !!!!        C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

      Islice = Nslice/2

 100  Ukin = 0.0d0
      Mx   = 0.0d0
      My   = 0.0d0

      Do I=1,Natom
         Rxx(I) = Rxf(I)
         Ryy(I) = Ryf(I)
         
         Vxx(I) = Ran_Gauss()
         Vyy(I) = Ran_Gauss()
            
         Mx = Mx + Vxx(I)
         My = My + Vyy(I)
      Enddo

      Mx = Mx/Dble(Natom)
      My = My/Dble(Natom)

      Do I=1,Natom
         Vxx(I) = Vxx(I) - Mx
         Vyy(I) = Vyy(I) - My
         
         Ukin = Ukin + 0.5d0*(Vxx(I)**2 + Vyy(I)**2)
      Enddo

      Udes   = Upath - H
      Factor = Dsqrt(Udes/Ukin)

      Do I=1,Natom
         Vxx(I) = Vxx(I)*Factor
         Vyy(I) = Vyy(I)*Factor
            
         Vxforw(I) = -Vxx(I)
         Vyforw(I) = -Vyy(I)
      Enddo
            
Ccccccccccccccccccccccccccccccccccccccc
C     Backward Integration First      C
Ccccccccccccccccccccccccccccccccccccccc

      Call Md(Islice,-1,(Islice-1),Du,.True.)

      If(.Not.La(1)) Goto 100

Cccccccccccccccccccccccccccc
C     Forward Integration  C
C     Copy Coordinates     C
Cccccccccccccccccccccccccccc

      Do I=1,Natom
         Rxx(I) = Rxf(I)
         Ryy(I) = Ryf(I)
      
         Vxx(I) = Vxforw(I)
         Vyy(I) = Vyforw(I)
      Enddo

      Call Md(Islice,1,(Nslice-Islice),Du,.True.)

      L_Ready = .False.

      If(Lumbrella) Then
         Do I=1,Natom
            Rxx(I) = Xxtra(I,Nslice)
            Ryy(I) = Yytra(I,Nslice)
         Enddo

         Call In_Slice(Lslice,Mx)
         
         If(Lslice) L_Ready = .True.
      Else
         Do J=1,Nslice
            If(Lb(J)) Then
               L_Ready = .True.
               Goto 10
            Endif
         Enddo
      Endif
 
 10   Continue
           
      If(.Not.L_Ready) Goto 100

Ccccccccccccccccccccccccccccc
C     We Have A Winner !!!  C
Ccccccccccccccccccccccccccccc

      Do J=1,Nslice
            
         La_Old(J) = La(J)
         Lb_Old(J) = Lb(J)
         Eeold(J)  = Eetra(J)
 
         Do I=1,Natom
            Xxold(I,J) = Xxtra(I,J)
            Yyold(I,J) = Yytra(I,J)
                     
            Vxold(I,J) = Vxtra(I,J)
            Vyold(I,J) = Vytra(I,J)
         Enddo
      Enddo
      
      Return
      End
