      Subroutine Fourierspace(Ukappa,Uself)
      Implicit None

      Include 'system.inc'

      Double Precision Ialpha,Ckc,Cks,Clm(Maxpart),
     &                 Prff,Slm(Maxpart),Elc(Maxpart,0:1),
     &                 Els(Maxpart,0:1),Cs,
     &                 Emc(Maxpart,0:Maxkap),
     &                 Ems(Maxpart,0:Maxkap),
     &                 Enc(Maxpart,0:Maxkap),
     &                 Ens(Maxpart,0:Maxkap),
     &                 Rkxx,Rkyy,Rkzz,Sx,Sy,Sz,
     &                 Ukappa,Uself,Rksq
 
      Integer I,J,Ll,Mm,Nn,L,M,N,Mmin,Nmin
 
Ccccccccccccccccccccccccccccccccccccccccccccc
C     Fourier Part Of The Ewald Summation   C
Ccccccccccccccccccccccccccccccccccccccccccccc

      Ialpha = -0.25d0/(Alpha*Alpha)
      Ukappa = 0.0d0
      Uself  = 0.0d0
      Prff   = 2.0d0*Twopi*Ibox*Ibox*Ibox
 
      Do I = 1,Npart
         Elc(I,0) = 1.0d0
         Emc(I,0) = 1.0d0
         Enc(I,0) = 1.0d0
         Els(I,0) = 0.0d0
         Ems(I,0) = 0.0d0
         Ens(I,0) = 0.0d0
         Sx       = Twopi*Ibox*Rx(I)
         Sy       = Twopi*Ibox*Ry(I)
         Sz       = Twopi*Ibox*Rz(I)
         Elc(I,1) = Dcos(Sx)
         Emc(I,1) = Dcos(Sy)
         Enc(I,1) = Dcos(Sz)
         Els(I,1) = Dsin(Sx)
         Ems(I,1) = Dsin(Sy)
         Ens(I,1) = Dsin(Sz)
      Enddo
 
      Mmin = 0
      Nmin = 1
 
      Do I = 1,Npart
         Do L = 2,Kmax
            Emc(I,L) = 
     &           Emc(I,L-1)*Emc(I,1) - Ems(I,L-1)*Ems(I,1)
            Ems(I,L) = 
     &           Ems(I,L-1)*Emc(I,1) + Emc(I,L-1)*Ems(I,1)
         Enddo
         Do L = 2,Kmax
            Enc(I,L) = 
     &           Enc(I,L-1)*Enc(I,1) - Ens(I,L-1)*Ens(I,1)
            Ens(I,L) = 
     &           Ens(I,L-1)*Enc(I,1) + Enc(I,L-1)*Ens(I,1)
         Enddo
      Enddo
 
      Do Ll = 0,Kmax
 
         Rkxx = Twopi*Dble(Ll)*Ibox
 
         If (Ll.Ge.1) Then
            Do I = 1,Npart
               Cs       = Elc(I,0)
               Elc(I,0) = Elc(I,1)*Cs - Els(I,1)*Els(I,0)
               Els(I,0) = Elc(I,1)*Els(I,0) + Els(I,1)*Cs
            Enddo
         Endif
 
         Do Mm = Mmin,Kmax
 
            M    = Iabs(Mm)
            Rkyy = Twopi*Dble(Mm)*Ibox
 
            If (Mm.Ge.0) Then
               Do I = 1,Npart
                  Clm(I) = 
     &                 Elc(I,0)*Emc(I,M) - Els(I,0)*Ems(I,M)
                  Slm(I) = 
     &                 Els(I,0)*Emc(I,M) + Ems(I,M)*Elc(I,0)
               Enddo
            Else
               Do I = 1,Npart
                  Clm(I) = 
     &                 Elc(I,0)*Emc(I,M) + Els(I,0)*Ems(I,M)
                  Slm(I) = 
     &                 Els(I,0)*Emc(I,M) - Ems(I,M)*Elc(I,0)
               Enddo
            Endif
 
            Do Nn = Nmin,Kmax
               N     = Iabs(Nn)
               Rkzz  = Twopi*Dble(Nn)*Ibox
               Rksq  = Rkxx*Rkxx + Rkyy*Rkyy + Rkzz*Rkzz
               Ckc   = 0.0d0
               Cks   = 0.0d0
 
               If (Nn.Ge.0) Then
                  Do I = 1,Npart
                     Ckc = Ckc + 
     &                    Z(I)*(Clm(I)*Enc(I,N)-Slm(I)*Ens(I,N))
                     Cks = Cks + 
     &                    Z(I)*(Slm(I)*Enc(I,N)+Clm(I)*Ens(I,N))
                  Enddo
               Else
                  Do I = 1,Npart
                     Ckc = Ckc + 
     &                    Z(I)*(Clm(I)*Enc(I,N)+Slm(I)*Ens(I,N))
                     Cks = Cks + 
     &                    Z(I)*(Slm(I)*Enc(I,N)-Clm(I)*Ens(I,N))
                  Enddo
               Endif
 
               Ukappa = Ukappa + 
     &              (Prff*Dexp(Ialpha*Rksq)/Rksq)*(Ckc*Ckc+Cks*Cks)
  
            Enddo
            Nmin = -Kmax
         Enddo
         Mmin = -Kmax
      Enddo
 
Ccccccccccccccccccccccccccccc
C     Add The Self Energy   C
Ccccccccccccccccccccccccccccc

      Do J = 1,Npart
         Uself = Uself - Alpha*Piesq*Z(J)*Z(J)
      Enddo
 
      return
      end
