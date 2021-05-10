      Subroutine Shifting
      Implicit None

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'mdvelo.inc'
      Include 'pathtopo.inc'
      Include 'traject.inc'

Ccccccccccccccccccccccccccccccccccccccccc
C     Shifting Move For Path Ensemble   C
Ccccccccccccccccccccccccccccccccccccccccc

      Double Precision Du,Ran_Uniform,Order

      Integer I,J,Ishift,Kkk,Nleft,Nright
      Logical L_Ready,Lslice

      Shift2 = Shift2 + 1.0d0

      L_Ready = .False.
      Ishift  = 1 + Idint(Ran_Uniform()*Dble(Nslice)*0.05d0)

      If(Ran_Uniform().Lt.0.5d0) Ishift = -Ishift
      
Cccccccccccccccccccccccccccccccccccc
C     Shift The Path               C
C     Copy Coordinates And Whether C
C     The Path Is In B At Time T   C
Cccccccccccccccccccccccccccccccccccc

      Do J=1,Maxtraject

         Kkk = J + Ishift
            
         If(Kkk.Ge.1.And.Kkk.Le.Maxtraject) Then

            La(Kkk)    = La_Old(J)
            Lb(Kkk)    = Lb_Old(J)
            Eetra(Kkk) = Eeold(J)

            Do I=1,Natom
               Xxtra(I,Kkk) = Xxold(I,J)
               Yytra(I,Kkk) = Yyold(I,J)
                                        
               Vxtra(I,Kkk) = Vxold(I,J)
               Vytra(I,Kkk) = Vyold(I,J)
            Enddo
         Endif
      Enddo

Cccccccccccccccccccccccccccccccccccccc
C     Calculate New Safe Boundaries  C
Cccccccccccccccccccccccccccccccccccccc

      Nleft  = Max(1          , (Ishift + 1))
      Nright = Min(Maxtraject , (Ishift + Maxtraject))

      If(Nleft.Gt.1) Then

Cccccccccccccccccccccccccc
C     Backward Shifting  C
C     Turn Momenta       C
Cccccccccccccccccccccccccc

         Do I=1,Natom
            Rxx(I) =  Xxtra(I,Nleft)
            Ryy(I) =  Yytra(I,Nleft)
                                   
            Vxx(I) = -Vxtra(I,Nleft)
            Vyy(I) = -Vytra(I,Nleft)
         Enddo

         Call Md(Nleft,-1,(Nleft-1),Du,.True.)

         Drift1 = Drift1 + Du
         Drift2 = Drift2 + 1.0d0
      Endif

Ccccccccccccccccccccccccccccccccccccccccccc
C     Check If Initial Point Is In A      C
Ccccccccccccccccccccccccccccccccccccccccccc
      
      If(.Not.La(1)) Return

Cccccccccccccccccccccccccc
C     So Far, So Good    C
Cccccccccccccccccccccccccc

      If(Nright.Lt.Nslice) Then

Ccccccccccccccccccccccccccccccccc
C     Forward Shifting          C
Ccccccccccccccccccccccccccccccccc

         Do I=1,Natom
            Rxx(I) = Xxtra(I,Nright)
            Ryy(I) = Yytra(I,Nright)
                                            
            Vxx(I) = Vxtra(I,Nright)
            Vyy(I) = Vytra(I,Nright)
         Enddo

         Call Md(Nright,1,(Nslice-Nright),Du,.True.)

         Drift1 = Drift1 + Du
         Drift2 = Drift2 + 1.0d0
      Endif

      If(Lumbrella) Then

Ccccccccccccccccccccccccccccccccccc
C     Check If Endpoint In Slice  C
Ccccccccccccccccccccccccccccccccccc

         Do I=1,Natom
            Rxx(I) = Xxtra(I,Nslice)
            Ryy(I) = Yytra(I,Nslice)
         Enddo

         Call In_Slice(Lslice,Order)

         If(.Not.Lslice) Return
         
         L_Ready = .True.
         
      Else

Ccccccccccccccccccccccccccccc
C     Check If Once In B    C
Ccccccccccccccccccccccccccccc
         
         Do J=1,Nslice
            If(Lb(J)) Then
               L_Ready = .True.
               Goto 10
            Endif
         Enddo
      Endif

 10   If(L_Ready) Then

Ccccccccccccccccccccccccccccccccccccccccc
C     Accepted; Energy Is Not Changed   C
C     (Only Due To Drift...)            C
Ccccccccccccccccccccccccccccccccccccccccc

         Shift1 = Shift1 + 1.0d0

         Do J=1,Maxtraject

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
