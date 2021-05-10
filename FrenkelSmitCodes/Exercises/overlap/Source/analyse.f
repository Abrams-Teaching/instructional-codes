      Program Analyse
      Implicit None

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Analyse The Histograms Of The Overlapping Dist. Method    C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer I,K,Nhisto

      Parameter(Nhisto = 20000)

      Double Precision Binsize,Dummy,
     &     B22(Nhisto),B23(Nhisto)

      Logical L22(Nhisto),L23(Nhisto)

      Do I=1,Nhisto
         B22(I)  = 0.0d0
         B23(I)  = 0.0d0
         L22(I)  = .False.
         L23(I)  = .False.
      Enddo

      Do I=22,23
         Read(I,*) Binsize

 98      Read(I,*,End=99) K,Dummy
         
         If(K.Le.0.Or.K.Gt.Nhisto) Then
            Write(6,*) 'Change Nhisto'
            Stop
         Endif

         If(I.Eq.22) Then
            B22(K) = Dummy
            L22(K) = .True.
         Else
            B23(K) = Dummy
            L23(K) = .True.
         Endif

         Goto 98

 99      Continue
      Enddo

Ccccccccccccccccccccccc
C     Write Results   C
Ccccccccccccccccccccccc

      Do I=1,Nhisto
         If(L22(I).And.L23(I)) Then
            Write(34,*) Binsize*(Dble(I)-0.5d0),
     &           B22(I)-B23(I)
         Endif

         If(L22(I)) 
     &        Write(32,*) Binsize*(Dble(I)-0.5d0),B22(I)
         
         If(L23(I)) 
     &        Write(33,*) Binsize*(Dble(I)-0.5d0),B23(I)
      Enddo

      Stop
      End
