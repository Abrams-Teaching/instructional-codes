      Program Umbrella
      Implicit None

      Integer Sstmm,Ncycle,I,Ii,Max,J

      Parameter(Max=100000)

      Double Precision M1,Lmin,Lmax,Xold,Xnew,Uold,Unew,
     &     Av1,Av2,Ran_Uniform,Dispm,Gg(-Max:Max),
     &     Width,Dddx

      Parameter(Width=0.001d0)

Ccccccccccccccccccccccccccccccccc
C     Written By Thijs Vlugt    C
C     Initialize Rng            C
Ccccccccccccccccccccccccccccccccc
 
      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)

Cccccccccccccccccccccc
C     Read Data      C
Cccccccccccccccccccccc
 
      Write(6,*) 'Number Of Cycles      ?'
      Read(*,*) Ncycle

      Write(6,*) 'Maximum Displacement  ?'
      Read(*,*) Dispm

      Write(6,*) 'Minimum Value Of X    ?'
      Read(*,*) Lmin

      Write(6,*) 'Maximum Value Of X   ? '
      Read(*,*) Lmax

      If(Lmin.Ge.Lmax) Stop

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Generate Initial Position And Energy   C
Cccccccccccccccccccccccccccccccccccccccccccccc

      Xold  = 0.5d0*(Lmax+Lmin)
      Uold  = Xold**2
      Xnew  = 0.0d0
      Unew  = 0.0d0
      Av1   = 0.0d0
      Av2   = 0.0d0
      Dddx  = 1.0d0/Width
      
      Do I=-Max,Max
         Gg(I) = 0.0d0
      Enddo
         
      Write(6,*)
      Write(6,*) 'Calculating......'
      Write(6,*)

Cccccccccccccccccccccccccccccc
C     Loop Over All Cycles   C
Cccccccccccccccccccccccccccccc

      Do I=1,Ncycle
         Do Ii=1,1000000
      
            Xnew = Xold + (Ran_Uniform()-0.5d0)*Dispm
            Unew = Xnew**2
            Av2  = Av2 + 1.0d0
            J    = Idint(Dddx*Xold)

            If(J.Lt.-Max.Or.J.Gt.Max) Then
               Write(6,*) 'Out Of Range : J !!!'
               Stop
            Else
               If(J.Eq.0) Then
                  Gg(J) = Gg(J) + 0.5d0
               Else
                  Gg(J) = Gg(J) + 1.0d0
               Endif
            Endif

Ccccccccccccccccccccccccccccccccccc
C     Reject When Outside Slice   C
Ccccccccccccccccccccccccccccccccccc

            If(Xnew.Gt.Lmin.And.Xnew.Lt.Lmax) Then
               If(Ran_Uniform().Lt.Dexp(Uold-Unew)) Then
                  Xold = Xnew
                  Uold = Unew
                  Av1  = Av1 + 1.0d0
               Endif
            Endif
         Enddo
      Enddo

Ccccccccccccccccccccccc
C     Write Results   C
Ccccccccccccccccccccccc

      Open(32,File='Umbrella.dat')
         
      Do I=-Max,Max
         If(Gg(I).Gt.0.5d0) Then
            Write(32,'(I7,E20.10)') I,Gg(I)/Av2
         Endif
      Enddo

      Close(32)

      Write(6,*)
      Write(6,*) 'Fraction Of Accepted Displacements : ',Av1/Av2
      Write(6,*)
      Write(6,*) 'Histogram Written To Umbrella.dat'
      Write(6,*)
 
      Stop
      End
