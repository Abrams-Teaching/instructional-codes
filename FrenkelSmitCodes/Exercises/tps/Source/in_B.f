      Function In_B()
      Implicit None

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'mdvelo.inc'

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Check If A Point Of The Path Is In B   C
Cccccccccccccccccccccccccccccccccccccccccccccc
      
      Logical          In_B
      Double Precision R1,Dx,Dy

      Dx = Rxx(1) - Rxx(2)
      Dy = Ryy(1) - Ryy(2)

      R1 = Dx*Dx + Dy*Dy
            
      If(R1.Gt.Right2) Then
         In_B = .True.
      Else
         In_B = .False.
      Endif

      Return
      End
