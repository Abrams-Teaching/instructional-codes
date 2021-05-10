      Subroutine Writepath
      Implicit None

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'traject.inc'

Ccccccccccccccccccccccccccc
C     Write Path To Disk  C
C     Binary Format       C
Ccccccccccccccccccccccccccc

      Integer I,J
      
      Open(35,File='pathnew',Form='Unformatted')

      Do J=1,Nslice
         Do I=1,Natom
            Write(35) Xxold(I,J),Yyold(I,J),Vxold(I,J),Vyold(I,J)
         Enddo
      Enddo

      Close(35)

      Return
      End
