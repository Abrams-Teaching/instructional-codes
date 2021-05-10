      Subroutine Readpath
      Implicit None

      Include 'maxarray.inc'
      Include 'mdvelo.inc'
      Include 'mdtop.inc'
      Include 'traject.inc'

Cccccccccccccccccccccccccccc
C     Read Path From Disk  C
Cccccccccccccccccccccccccccc

      Integer          I,J
      Logical          In_A,In_B,Lready
      Double Precision Order,Ukin,M1,Upot,Utot,Av1
      
      Open(32,File='pathold',Form='Unformatted')
        
      Do J=1,Nslice
         Do I=1,Natom

            Read(32) Xxold(I,J),Yyold(I,J),
     &               Vxold(I,J),Vyold(I,J)

            Xxtra(I,J) = Xxold(I,J)
            Yytra(I,J) = Yyold(I,J)
                                  
            Vxtra(I,J) = Vxold(I,J)
            Vytra(I,J) = Vyold(I,J)
                           
            Rxx(I) = Xxold(I,J)
            Ryy(I) = Yyold(I,J)
         Enddo

         La_Old(J) = In_A()
         Lb_Old(J) = In_B()
      Enddo
         
      Close(32)

Cccccccccccccccccccccccccc
C     Check Some Things  C
Cccccccccccccccccccccccccc

      If(.Not.La_Old(1)) Then
         Write(6,*) 'Initial Path Does Not Start In A !!!'
         Write(6,*) 'Error Readpath'
         Call Exitt(2)
      Endif

      If(Lumbrella) Then
          
         Do I=1,Natom
            Rxx(I) = Xxold(I,Nslice)
            Ryy(I) = Yyold(I,Nslice)
         Enddo

         Call In_Slice(Lready,Order)

         If(.Not.Lready) Then
            Write(6,*) 'Path Is Not In Slice !!!'
            Write(6,*) 'Error Readpath For Umbrella Sampling !!!'
            Call Exitt(2)
         Endif

      Else

         Lready = .False.

         Do J=1,Nslice
            If(Lb_Old(J)) Lready = .True.
         Enddo
         
         If(.Not.Lready) Then
            Write(6,*) 'Path Is Never In B !!!'
            Write(6,*) 'Error Readpath For Pathsampling !!!'
            Call Exitt(2)
         Endif
      Endif
      
ccccccccccccccccccccccccccccccccccccccccccc
C     Calculate The Energy Of All Paths   C
Ccccccccccccccccccccccccccccccccccccccccccc

      Do J=1,Nslice
         Do I=1,Natom
            Rxx(I) = Xxold(I,J)
            Ryy(I) = Yyold(I,J)
                              
            Vxx(I) = Vxold(I,J)
            Vyy(I) = Vyold(I,J)
         Enddo
      
         Call Force(Utot,Upot,Ukin)
            
         Eeold(J) = Utot
      Enddo

      M1  = Eeold(1)
      Av1 = 0.0d0
                     
      Do J=1,Nslice
         Av1 = Av1 + Dabs((Eeold(J) - M1)/M1)
      Enddo

      Write(6,*) 'Initial Energy                         : ',
     &     Eeold(1)      
      Write(6,*) 'Path Energy Drift                      : ',
     &     Av1/Dble(Nslice)
     
      Return
      End
