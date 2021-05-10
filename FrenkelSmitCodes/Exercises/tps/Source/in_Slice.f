      Subroutine In_Slice(Lslice,R1)
      Implicit None

Ccccccccccccccccccccccccccccccccccccccccccc
C     Checks If A Particle Is In A Slice  C
C                                         C
C     Lmin/Lmax = Boundaries              C
C     R1        = Lambda**2               C
Ccccccccccccccccccccccccccccccccccccccccccc

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'mdvelo.inc'

      Logical          Lslice
      Double Precision R1,Dx,Dy

      Dx = Rxx(1) - Rxx(2)
      Dy = Ryy(1) - Ryy(2)

      R1 = Dx*Dx + Dy*Dy

      If(R1.Gt.Lmin2.And.R1.Lt.Lmax2) Then
         Lslice = .True.
      Else
         Lslice = .False.
      Endif

      Return
      End
