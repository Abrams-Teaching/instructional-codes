      Program Overdist
      Implicit None

      Include 'maxarray.inc'
      Include 'system.inc'

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Overlapping Distribution Of Polymers Using Cbmc   C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer I,Sstmm,Ncycle,Ii,Iii,Nhisto,Ip
      Double Precision M1,U,Ran_Uniform,Dispa,
     &     Dispb,W1,W2,Weight,Dx,Dy,Dz,Deltaw,
     &     Dispc,Dispd,Av1,Av2

      Parameter (Nhisto = 10000)
      Parameter (Deltaw = 0.05d0)

      Double Precision Disw(Nhisto),Ideltaw

Ccccccccccccccccccccccccccccc
C     Initialize Rng        C
Ccccccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)

Ccccccccccccccccccccccccccccc
C     Read Info From Disk   C
C     Initialize Histogram  C
Ccccccccccccccccccccccccccccc

      Do I=1,Nhisto
         Disw(I) = 0.0d0
      Enddo

      Ideltaw = 1.0d0/Deltaw

      Read(21,*) Npart,Ncycle,Box,Alpha,Beta,
     &     Dispmax,Nchoi,Nlength,Lchain,Pchain

      Write(6,*)
      Write(6,*) 'Number Of Particles   : ',Npart
      Write(6,*) 'Number Of Cycles      : ',Ncycle
      Write(6,*) 'Boxsize               : ',Box
      Write(6,*) 'Density               : ',Dble(Npart)/(Box**3)
      Write(6,*) 'Repulsion Parameter   : ',Alpha
      Write(6,*) 'Beta                  : ',Beta
      Write(6,*) 'Maximum Displacement  : ',Dispmax
      Write(6,*) 'Number Of Trial Dirs  : ',Nchoi
      Write(6,*) 'Chain Length          : ',Nlength
      Write(6,*) 'Prob. Cbmc            : ',Pchain
      Write(6,*)

      If(Lchain) Then
         Write(6,*) 'Chain Is Added (N+1) Particles'
      Else
         Write(6,*) 'Only Test Chain (N) Particles'
      Endif

      If(Pchain.Le.0.0d0.Or.Pchain.Ge.1.0d0) Then
         Write(6,*) 'Pchain Out Of Rangle !!'
         Stop
      Endif

      If(Nlength.Lt.1.Or.Nlength.Gt.Maxlength) Then
         Write(6,*) 'Nlength Out Of Rangle !!'
         Stop
      Endif

      If(Nchoi.Lt.1.Or.Nchoi.Gt.Maxtrial) Then
         Write(6,*) 'Nchoi Out Of Rangle !!'
         Stop
      Endif

      If(Beta.Le.0.0d0) Then
         Write(6,*) 'Beta Out Of Rangle !!'
         Stop
      Endif

      If(Npart.Le.1.Or.Npart.Gt.Maxpart) Then
         Write(6,*) 'Npart Out Of Rangle !!'
         Stop
      Endif

      If(Box.Le.2.0d0) Then
         Write(6,*) 'Box Out Of Rangle !!'
         Stop
      Endif 

      Ibox = 1.0d0/Box

      If(Alpha.Le.0.0d0.Or.Alpha.Gt.50.0d0) Then
         Write(6,*) 'Alpha Out Of Rangle !!'
         Stop
      Endif 

      If(Dispmax.Le.0.001d0.Or.Dispmax.Gt.1.0d0) Then
         Write(6,*) 'Dispmax Out Of Rangle !!'
         Stop
      Endif 
      
Ccccccccccccccccccccccccccccccccccccccc
C     Generate Initial Coordinates    C
Ccccccccccccccccccccccccccccccccccccccc

      Do I=1,Npart
         Xx(I) = Box*Ran_Uniform()
         Yy(I) = Box*Ran_Uniform()
         Zz(I) = Box*Ran_Uniform()
      Enddo

      X(1) = Box*Ran_Uniform()
      Y(1) = Box*Ran_Uniform()
      Z(1) = Box*Ran_Uniform()

      Do I=2,Nlength
         Call Ran_Sphere(Dx,Dy,Dz)

         X(I) = X(I-1) + Dx
         Y(I) = Y(I-1) + Dy
         Z(I) = Z(I-1) + Dz
      Enddo
      
Ccccccccccccccccccccccccccccccc
C     Loop Over Cycles        C
Ccccccccccccccccccccccccccccccc

      Dispa = 0.0d0
      Dispb = 0.0d0
      Dispc = 0.0d0
      Dispd = 0.0d0
      W1    = 0.0d0
      W2    = 0.0d0
      Av1   = 0.0d0
      Av2   = 0.0d0

      Call Ener_Tot(U)

      Write(6,*)
      Write(6,*) 'Initial Energy        : ',U

      Usim = U

      Do I=1,Ncycle
         Do Ii=1,Npart
            
Cccccccccccccccccccccccccccccccccc
C     Calculate Average Energy   C
Cccccccccccccccccccccccccccccccccc

            If(I.Gt.(Ncycle/2)) Then
               Av1 = Av1 + Usim
               Av2 = Av2 + 1.0d0
            Endif

Cccccccccccccccccccccccccccc
C     Select A Trialmove   C
Cccccccccccccccccccccccccccc

            If(Lchain) Then
               If(Ran_Uniform().Lt.Pchain) Then
                  Call Regrow(Dispc,Dispd)
               Else
                  Call Disp_Mono(Dispa,Dispb)
               Endif
            Else
               Call Disp_Mono(Dispa,Dispb)
            Endif

            If(Lchain.And.I.Gt.(Ncycle/2)) Then
               Call Widom(Weight)
               W1 = W1 + Weight
               W2 = W2 + 1.0d0

               If(Weight.Gt.1.0d-200) Then
                  Weight = -Dlog(Weight)
                  Ip     = 1 + Idint(Ideltaw*Weight)

                  If(Ip.Lt.Nhisto.And.Ip.Ge.1) 
     &                 Disw(Ip) = Disw(Ip) + 1.0d0
               Endif

            Endif
         Enddo

         If(.Not.Lchain.And.I.Gt.(Ncycle/2)) Then
            Do Iii=1,Npart
               Call Widom(Weight)
               W1 = W1 + Weight
               W2 = W2 + 1.0d0

               If(Weight.Gt.1.0d-200) Then
                  Weight = -Dlog(Weight)
                  Ip     = 1 + Idint(Ideltaw*Weight)

                  If(Ip.Lt.Nhisto.And.Ip.Ge.1) 
     &                 Disw(Ip) = Disw(Ip) + 1.0d0
               Endif
            Enddo
         Endif
      Enddo

Cccccccccccccccccccccccccccccc
C     End Of The Simulation  C
C     Calculate Averages     C
Cccccccccccccccccccccccccccccc

      Call Ener_Tot(U)

      Write(6,*) 'Final Energy          : ',Usim
      Write(6,*) 'Final Energy (Rigor)  : ',U
      Write(6,*) 'Energy Difference     : ',Dabs(U-Usim)
      Write(6,*) 'Average Energy        : ',Av1/Av2
      Write(6,*) 'Frac. Acc. Displ.     : ',Dispa/Dispb

      If(Dispd.Gt.0.5d0) 
     &     Write(6,*) 'Frac. Acc. Regr.      : ',Dispc/Dispd
      
Cccccccccccccccccccccccccccccccccccccccccccccc
C     Write Histograms To Disk               C
C     The Histogram Is Made For -Ln(W) !!!   C
Cccccccccccccccccccccccccccccccccccccccccccccc

C     Start Modification
C     End   Modification

C     See the file analyse.f and ../Run/runhow to format the output

      Stop
      End
