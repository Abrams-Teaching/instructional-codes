      Function In_A()
      Implicit None

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'mdvelo.inc'

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Check If A Point Of The Path Is In A   C
Cccccccccccccccccccccccccccccccccccccccccccccc
      
      Logical          In_A
      Double Precision R1,Dx,Dy

      Dx = Rxx(1) - Rxx(2)
      Dy = Ryy(1) - Ryy(2)

      R1 = Dx*Dx + Dy*Dy
      
      If(R1.Lt.Left2) Then
         In_A = .True.
      Else
         In_A = .False.
      Endif

      Return
      End
