      Subroutine Readdat
      Implicit None
 
      Include 'mpif.h'
      Include 'mpifcomms.inc'
      Include 'maxarray.inc'
      Include 'system.inc'
 
      If(Iirank.Eq.0) Then

C     Read In System Information
 
         Read (21,*) Box,Npart,Nstep,Temp,Tstep,Ninit,Nprint
 
         If (Npart.Gt.Maxpart) Then
            Write (6,*) 'Maximum No. Of Particles Is : ',Maxpart
            Call Mpi_Finalize(Iierr)
         Endif
 
C     Calculate Some Parameters
 
         Hbox = 0.5d0*Box
 
         Rcutsq = (0.49999d0*Box)**2
         Ecut = 4.0d0*((Rcutsq**(-6.0d0)) - 
     &        (Rcutsq**(-3.0d0)))
 
         Tstep2 = Tstep*Tstep
         I2tstep = 0.5d0/Tstep
 
C     Print Information To The Screen
 
         Write(6,*) 'Number Of Particles   : ',Npart
         Write(6,*) 'Boxlength             : ',Box
         Write(6,*) 'Density               : ',Dble(Npart)/(Box*Box*Box)
         Write(6,*) 'Temperature           : ',Temp
         Write(6,*) 'Cut-Off Radius        : ',Dsqrt(Rcutsq)
         Write(6,*) 'Cut-Off Energy        : ',Ecut
         Write(6,*) 'Number Of Steps       : ',Nstep
         Write(6,*) 'Number Of Init Steps  : ',Ninit
         Write(6,*) 'Timestep              : ',Tstep
      Endif

C     Copy To All Processors

      Call Mpi_Bcast(Box,1,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Hbox,1,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Rcutsq,1,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Ecut,1,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Tstep,1,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Tstep2,1,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(I2tstep,1,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Temp,1,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Nstep,1,
     &     Mpi_Integer,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Ninit,1,
     &     Mpi_Integer,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Npart,1,
     &     Mpi_Integer,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Nprint,1,
     &     Mpi_Integer,0,Mpi_Comm_World,Iierr)
      Call Mpi_Barrier(Mpi_Comm_World,Iierr)

      Return
      End
