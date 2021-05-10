      Subroutine Sample_Ct(Switch,Energy)
      Implicit None

      Include 'maxarray.inc'
      Include 'mdtop.inc'

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Compute C(t) By Direct Molecular Dynamics             C
C                                                           C
C     Tmax  = Maximum Timesteps For The Correlation Time    C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer   Tmax

      Parameter (Tmax = 200000)

      Integer Ttel,Tt0(Tmax),T0,Switch,I,T,Tvacf,Dt
      Logical Ha0(Tmax),In_A,In_B,Lla,Llb

      Double Precision Vacf(Tmax),Nvacf(Tmax),
     &     Av1,Av2,Av3,M1,M3,Energy
       
      Save Av1,Av2,Av3,Vacf,Nvacf,Tvacf,Tt0,T0,Ha0

      If(Switch.Eq.1) Then

         Tvacf = 0
         T0    = 0
         Av1   = 0.0d0
         Av2   = 0.0d0
         Av3   = 0.0d0

         Do I = 1, Tmax
            Nvacf(I) = 0.0d0
            Vacf(I)  = 0.0d0
            Tt0(I)   = 0
            Ha0(I)   = .False.
         Enddo

      Elseif(Switch.Eq.2) Then

         Lla = In_A()
         Llb = In_B()

         If(Lla) Av1 = Av1 + 1.0d0
         If(llb) Av3 = Av3 + 1.0d0
         
         Av2   = Av2   + 1.0d0
         Tvacf = Tvacf + 1

         If (Mod(Tvacf,50).Eq.0) Then
            T0         = T0 + 1
            Ttel       = Mod(T0-1, Tmax) + 1
            Tt0(Ttel)  = Tvacf
            Ha0(Ttel)  = Lla
         Endif

         Do T = 1, Min(T0, Tmax)
            Dt = Tvacf - Tt0(T) + 1

            If (Dt.Lt.Tmax) Then
               Nvacf(Dt) = Nvacf(Dt) + 1.0d0
            
               If(Llb.And.Ha0(T)) Vacf(Dt) = Vacf(Dt) + 1.0d0
            Endif
         Enddo

      Else

         M1 = Av1/Av2
         M3 = Av3/Av2

         Write(6,'(A,3E15.8)') 
     &        ' Fraction Spend In A,B, Energy          : ',
     &        M1,M3,Energy

         Open(32,File='ct.dat')

         Do I = 1, Tmax-1
            If (Nvacf(I).Ge.0.5d0) 
     &           Write(32,*) Dble(I-1)*Tstep*Dble(Nshort),
     &           Vacf(I)/(Nvacf(I)*M1)
         Enddo

         Close(32)

      Endif

      Return
      End
