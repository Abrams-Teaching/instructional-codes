      Subroutine Force(Utot,Upot,Ukin)
      Implicit None

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'mdvelo.inc'
      
Cccccccccccccccccccccccccccccccc
C     Force And Total Energy   C
Cccccccccccccccccccccccccccccccc

      Integer I,J,Kkk

      Double Precision Dx,Dy,R2i,R6i,Upot,Dis,Ri,
     &                 Utot,Ukin,Ff,R,Wdis

Ccccccccccccccccccccccccccc
C     Set Stuff To Zero   C
Ccccccccccccccccccccccccccc

      Upot = 0.0d0
      Ukin = 0.0d0
      
      Do I=1,Natom
         Fxx(I) = 0.0d0
         Fyy(I) = 0.0d0

         Ukin = Ukin + Vxx(I)**2 + Vyy(I)**2
      Enddo

      Ukin = 0.5d0*Ukin
 
Cccccccccccccccccccccccccccccccc
C     Nonbonded Interactions   C
C     Use A List Of Pairs      C
Cccccccccccccccccccccccccccccccc

      Do Kkk=1,Ninter
         I = Iinter(Kkk)
         J = Jinter(Kkk)

         Dx  = Rxx(I) - Rxx(J)
         Dy  = Ryy(I) - Ryy(J)
         R2i = Dx*Dx + Dy*Dy

         If(R2i.Lt.Rwca2) Then
            R2i    = 1.0d0/R2i
            R6i    = R2i*R2i*R2i

            Upot   = Upot + 4.0d0*R6i*(R6i - 1.0d0) + 1.0d0
            Ff     = 48.0d0*R6i*(R6i-0.5d0)*R2i

            Fxx(I) = Fxx(I) + Ff*Dx
            Fyy(I) = Fyy(I) + Ff*Dy
                  
            Fxx(J) = Fxx(J) - Ff*Dx
            Fyy(J) = Fyy(J) - Ff*Dy
                  
         Endif
      Enddo
      
Cccccccccccccccccccc
C     Spring 1-2   C
Cccccccccccccccccccc

      Dx = Rxx(1) - Rxx(2)
      Dy = Ryy(1) - Ryy(2)

      R2i = Dsqrt(Dx*Dx + Dy*Dy)
      R   = (R2i - Rwca - W)*Iw
      R6i = 1.0d0 - R**2
         
      Upot = Upot + H*R6i*R6i
      Ff   = 4.0d0*H*R6i*R*Iw/R2i

      Fxx(1) = Fxx(1) + Ff*Dx
      Fyy(1) = Fyy(1) + Ff*Dy

      Fxx(2) = Fxx(2) - Ff*Dx 
      Fyy(2) = Fyy(2) - Ff*Dy 
      
Cccccccccccccccccccccccccccccccccccc
C     Interactions With Boundary   C
Cccccccccccccccccccccccccccccccccccc

      Do I=1,Natom
         Dx  = Rxx(I)
         Dy  = Ryy(I)
         Dis = Dx*Dx + Dy*Dy

         If(Dis.Gt.Rad2) Then
            Wdis = Dsqrt(Dis)
            R    = Rwca + Rad - Wdis
            Ri   = 1.0d0/R
            R2i  = Ri*Ri
            R6i  = R2i*R2i*R2i

            Upot = Upot + 4.0d0*R6i*(R6i - 1.0d0) + 1.0d0
            Ff   = (48.0d0*R6i*R6i*Ri - 24.0d0*R6i*Ri)/Wdis

            Fxx(I) = Fxx(I) - Ff*Dx
            Fyy(I) = Fyy(I) - Ff*Dy
         Endif
      Enddo

Ccccccccccccccc
C     Utot    C
Ccccccccccccccc

      Utot = Upot + Ukin

      Return
      End
