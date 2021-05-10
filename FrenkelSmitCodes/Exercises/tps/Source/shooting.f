      Subroutine Shooting
      Implicit None

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'mdvelo.inc'
      Include 'pathtopo.inc'
      Include 'traject.inc'

Ccccccccccccccccccccccccccccccccccccccccc
C     Shooting Move For Path Ensemble   C
Ccccccccccccccccccccccccccccccccccccccccc

      Double Precision Vxforw(Maxatom),Vyforw(Maxatom),
     &                 Du,Ran_Uniform,Order,Mag,Phi

      Integer          I,J,Islice,Imol
      Logical          L_Ready,Lslice

      Islice  = 1 + Idint(Ran_Uniform()*Dble(Nslice))
      Imol    = 1 + Idint(Ran_Uniform()*Dble(Natom))
      Disp2   = Disp2 + 1.0d0
      L_Ready = .False.

Cccccccccccccccccccccccccccc
C     New Configuration    C
Cccccccccccccccccccccccccccc

      Do I=1,Natom
         Rxx(I) = Xxold(I,Islice)
         Ryy(I) = Yyold(I,Islice)

         Vxx(I) = -Vxold(I,Islice)
         Vyy(I) = -Vyold(I,Islice)
         
Cccccccccccccccccccccccccccccccccccccccccccccc
C     Rotation Of The Velocity Of One Atom   C
Cccccccccccccccccccccccccccccccccccccccccccccc

         If(I.Eq.Imol) Then
            Mag    = Dsqrt(Vxx(I)**2 + Vyy(I)**2)
            Phi    = Dacos(Vxx(I)/Mag)

            If(Vyy(I).Lt.0.0d0) Phi = - Phi

            Phi    = Phi + Twopi + Deltaphi*(Ran_Uniform()-0.5d0)

            Vxx(I) = Mag*Dcos(Phi)
            Vyy(I) = Mag*Dsin(Phi)
         Endif
         
         Vxforw(I) = -Vxx(I)
         Vyforw(I) = -Vyy(I)
      Enddo

Ccccccccccccccccccccccccccccccccccccccc
C     Backward Integration First      C
C     Rejected If Initially Not In A  C
C                                     C
C     Reject/Accept Path Energy       C
C     Use The First Point Of The      C
C     Path For This                   C
Ccccccccccccccccccccccccccccccccccccccc

      Call Md(Islice,-1,(Islice-1),Du,.True.)

      Drift1 = Drift1 + Du
      Drift2 = Drift2 + 1.0d0
      
      If(.Not.La(1)) Return
           
Cccccccccccccccccccccccccccc
C     Forward Integration  C
C     Copy Coordinates     C
Cccccccccccccccccccccccccccc

      Do I=1,Natom
         Rxx(I) = Xxold(I,Islice)
         Ryy(I) = Yyold(I,Islice)
      
         Vxx(I) = Vxforw(I)
         Vyy(I) = Vyforw(I)
      Enddo

      Call Md(Islice,1,(Nslice-Islice),Du,.True.)

      Drift1 = Drift1 + Du
      Drift2 = Drift2 + 1.0d0
      
      If(Lumbrella) Then

Ccccccccccccccccccccccccccccccccccc
C     Check If Endpoint In Slice  C
Ccccccccccccccccccccccccccccccccccc

         Do I=1,Natom
            Rxx(I) = Xxtra(I,Nslice)
            Ryy(I) = Yytra(I,Nslice)
         Enddo

         Call In_Slice(Lslice,Order)
         
         If(.Not.Lslice) Then
            Return
         Else
            L_Ready = .True.
         Endif
         
      Else

Ccccccccccccccccccccccccccc
C     Check If In B Once  C
Ccccccccccccccccccccccccccc
         
         Do J=1,Nslice
            If(Lb(J)) Then
               L_Ready = .True.
               Goto 10
            Endif
         Enddo
      Endif
      
 10   If (L_Ready) Then

Ccccccccccccccccccccccccccccc
C     Accepted              C
C     Copy Coordinates      C
Ccccccccccccccccccccccccccccc

         Disp1 = Disp1 + 1.0d0
                  
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
      Endif
      
      Return
      End
