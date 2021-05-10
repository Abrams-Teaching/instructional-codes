      Subroutine Grow(Weight,Uchain,Xf,Yf,Zf,Lold)
      Implicit None

      Include 'maxarray.inc'
      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccc
C     Grow A Chain Using Cbmc        C
C                                    C
C     Xf/Yf/Zf = Coordinates         C
C     Lold     = Old Config ?        C
C     Uchain   = Energy              C
C     Weight   = Rosenbluth Weight   C
Cccccccccccccccccccccccccccccccccccccc

      Integer Ip,I,Ii,Ichoi
      Logical Lold

      Double Precision R2,Dx,Dy,Dz,Xf(Maxlength),
     &     Yf(Maxlength),Zf(Maxlength),Xt(Maxtrial),
     &     Yt(Maxtrial),Zt(Maxtrial),Ran_Uniform,
     &     Cumw,Sumw,Cmmm,Bfac(Maxtrial),U,Weight,
     &     Uchain,Ut(Maxtrial)

Ccccccccccccccccccccccccccccccccccccccccc
C     First Atom Of The Chain           C
C     Old Config Or Random In The Box   C
Ccccccccccccccccccccccccccccccccccccccccc

      If(Lold) Then
         Xf(1) = X(1)
         Yf(1) = Y(1)
         Zf(1) = Z(1)
      Else
         Xf(1) = Ran_Uniform()*Box
         Yf(1) = Ran_Uniform()*Box
         Zf(1) = Ran_Uniform()*Box
      Endif

      Call Ener_Mono(U,0,Xf(1),Yf(1),Zf(1),.False.)

      Weight = Dexp(-Beta*U)
      Uchain = U

      Do Ip = 2,Nlength

         Sumw = 0.0d0

Cccccccccccccccccccccccccccccccccccc
C     Loop Over Trial Directions   C
Cccccccccccccccccccccccccccccccccccc

         Do I=1,Nchoi

            If(Lold.And.I.Eq.1) Then
               Xt(I) = X(Ip)
               Yt(I) = Y(Ip)
               Zt(I) = Z(Ip)
            Else
               Call Ran_Sphere(Dx,Dy,Dz)

               Xt(I) = Xf(Ip-1) + Dx
               Yt(I) = Yf(Ip-1) + Dy
               Zt(I) = Zf(Ip-1) + Dz
            Endif

            Call Ener_Mono(U,0,Xt(I),Yt(I),Zt(I),.False.)

Ccccccccccccccccccccccccccccccccccccc
C     Intramolecular Interactions   C
Ccccccccccccccccccccccccccccccccccccc

            If(Ip.Ge.3) Then
               Do Ii=1,Ip-2
                  Dx = Xf(Ii) - Xt(I)
                  Dy = Yf(Ii) - Yt(I)
                  Dz = Zf(Ii) - Zt(I)

                  R2 = Dx**2 + Dy**2 + Dz**2

                  If(R2.Lt.1.0d0) 
     &                 U  = U + 
     &                 Alpha*((Dsqrt(R2) - 1.0d0)**2)
               Enddo
            Endif

            Bfac(I) = Dexp(-Beta*U)
            Sumw    = Sumw + Bfac(I)
            Ut(I)   = U
         Enddo

Ccccccccccccccccccccccccccccccccc
C     Select A Configuration    C
Ccccccccccccccccccccccccccccccccc

         Weight = Weight*Sumw/Dble(Nchoi)
         Ichoi  = 1

         If(.Not.Lold) Then
            Cumw   = Bfac(1)
            Cmmm   = Sumw*Ran_Uniform()

            Do While(Cumw.Lt.Cmmm)
               Ichoi = Ichoi + 1
               Cumw  = Cumw  + Bfac(Ichoi)
            Enddo
         Endif
         
         Xf(Ip) = Xt(Ichoi)
         Yf(Ip) = Yt(Ichoi)
         Zf(Ip) = Zt(Ichoi)
         
         Uchain = Uchain + Ut(Ichoi)
      Enddo

      Return
      End
