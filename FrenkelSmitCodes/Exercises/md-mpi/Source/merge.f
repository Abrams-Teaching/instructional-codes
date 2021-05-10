      Subroutine Merge
      Implicit None

      Include 'mpif.h'
      Include 'mpifcomms.inc'
      Include 'maxarray.inc'
      Include 'system.inc'

Ccccccccccccccccccccccccccccccccccccccccccccccc
C     Merge All Cordinates Between The Nodes  C
Ccccccccccccccccccccccccccccccccccccccccccccccc

      Call Mpi_Bcast(Rxx,Npart,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Ryy,Npart,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Rzz,Npart,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Rxf,Npart,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Ryf,Npart,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Rzf,Npart,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Vxx,Npart,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Vyy,Npart,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Bcast(Vzz,Npart,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
      Call Mpi_Barrier(Mpi_Comm_World,Iierr)

      Return
      End
