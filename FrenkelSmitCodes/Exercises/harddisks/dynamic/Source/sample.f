      Subroutine Sample(Ichoise)
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Samples The Radial Distribution Function   C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer I,J,Maxx,Ichoise,A

      Parameter(Maxx = 500)

      Double Precision Ggt,Gg(Maxx),Delta,R2

      Save Ggt,Gg,Delta

      If(Ichoise.Eq.1) Then
         Do I=1,Maxx
            Gg(I) = 0.0d0
         Enddo

         Ggt   = 0.0d0
         Delta = Dble(Maxx-1)/5.0d0
      
      Elseif(Ichoise.Eq.2) Then

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Sample The Radial Distribution Function   C
C     Loop Over All Particle Pairs              C
C     See Frenkel/Smit P. 77                    C
C                                               C
C     Delta  = 1/Binsize                        C
C     Ggt/Gg = Counters                         C
C     Maxx   = Maximum Number Of Bins           C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

C     Start Modification

C     End   Modification

      Else

Ccccccccccccccccccccccccccccccc
C     Write Results To Disk   C
Ccccccccccccccccccccccccccccccc

         Ggt   = 1.0d0/(Ggt*Dble(Npart))
         Delta = 1.0d0/Delta

         Do I=1,Maxx-1
            R2 = 4.0d0*Datan(1.0d0)*
     &           Dble(Npart-1)*0.01d0*Delta*Delta*
     &           ((Dble(I+1))**2 - (Dble(I)**2))
            
            Write(21,*) ((Dble(I)-0.5d0)*Delta),Gg(I)*Ggt/R2
         Enddo
      Endif

      Return
      End
