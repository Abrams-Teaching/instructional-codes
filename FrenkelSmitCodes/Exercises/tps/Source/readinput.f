      Subroutine Readinput
      Implicit None

      Include 'maxarray.inc'
      Include 'mdtop.inc'
      Include 'mdvelo.inc'
      
CCccccccccccccccccccccccccccc
C     Read Info From Disk   C
Ccccccccccccccccccccccccccccc
       
      Integer          I,J

      Double Precision Ran_Gauss,Ukin,Factor,Mx,My,Udes
 
      Open(32,File='input1')

      Read(32,*)
      Read(32,*) Nstep,Lscreen
      Read(32,*)
      Read(32,*) Lensemble,Tstep
      Read(32,*)
      Read(32,*) Upath,H,Lmin,Lmax

      Close(32)

Cccccccccccccccccccccccccccc
C     Set Some Constants   C
Cccccccccccccccccccccccccccc

      Nslice = Maxtraject
      Natom  = Maxatom
      Ninter = 0
      Iw     = 1.0d0/W
      Twopi  = 8.0d0*Datan(1.0d0)
      Rwca   = 2.0d0**(1.0d0/6.0d0)
      Rwca2  = Rwca**2
      Rad2   = Rad**2
      Left2  = Left**2
      Right2 = Right**2
      Lmin2  = Lmin**2
      Lmax2  = Lmax**2

Ccccccccccccccccccccccccccccccccccc
C     Make Interaction Table      C
C     Exclude 1-2 Interactions    C
Ccccccccccccccccccccccccccccccccccc

      Do I=1,Natom-1
         Do J=(I+1),Natom

            If(.Not.(I.Eq.1.And.J.Eq.2)) Then
               Ninter         = Ninter + 1
               Iinter(Ninter) = I
               Jinter(Ninter) = J
            Endif
         Enddo
      Enddo

      If(Upath.Le.H) Then
         Write(6,*) 'Upath Is Too Low !!!'
         Call Exitt(2)
      Endif

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Read Coordinates From Disk; Generate New Momenta    C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Open(32,File='initial-config')

      Do I=1,Natom
         Read(32,*) Rxx(i),Ryy(i)

         Rxf(I) = Rxx(I)
         Ryf(I) = Ryy(I)
      Enddo

      Close(32)

      Ukin = 0.0d0
      Mx   = 0.0d0
      My   = 0.0d0

      Do I=1,Natom
         Vxx(I) = Ran_Gauss()
         Vyy(I) = Ran_Gauss()

         Mx = Mx + Vxx(I)
         My = My + Vyy(I)
      Enddo

      Mx = Mx/Dble(Natom)
      My = My/Dble(Natom)
 
      Do I=1,Natom
         Vxx(I) = Vxx(I) - Mx
         Vyy(I) = Vyy(I) - My
         Ukin   = Ukin + 0.5d0*(Vxx(I)**2 + Vyy(I)**2)
      Enddo

      Udes = Upath - H

      Factor = Dsqrt(Udes/Ukin)

      Do I=1,Natom
         Vxx(I) = Vxx(I)*Factor
         Vyy(I) = Vyy(I)*Factor
      Enddo

Cccccccccccccccccccc
C     Print Info   C
Cccccccccccccccccccc

      Udes = Dble(Natom)/(4.0d0*Datan(1.0d0)*Rad2)

      Write(6,*)
      Write(6,*) 'Total Number Of Atoms                  : ',Natom
      Write(6,*) 'Number Of Interactions                 : ',Ninter
      Write(6,*) 'Radius                                 : ',Rad
      Write(6,*) 'Density                                : ',Udes
      Write(6,*) 'W                                      : ',W
      Write(6,*) 'H                                      : ',H
      Write(6,*) 'Total Energy                           : ',Upath
      Write(6,*) 'Total Number Of Cycles                 : ',Nstep
      Write(6,*) 'Steps In A Cycle                       : ',Nshort
      Write(6,*) 'Time Step                              : ',Tstep
      Write(6,*) 'Left  (Region A)                       : ',Left
      Write(6,*) 'Right (Region B)                       : ',Right
      Write(6,*) 'Lmin                                   : ',Lmin
      Write(6,*) 'Lmax                                   : ',Lmax
      Write(6,*)

      Return
      End
