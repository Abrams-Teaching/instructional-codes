      Program Md
      Implicit None

      Include 'mpif.h'
      Include 'mpifcomms.inc'
      Include 'maxarray.inc'
      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccccccccccccc
C     Molecular Dynamics Program Of Argon,     C
C     MPI version                              C
Cccccccccccccccccccccccccccccccccccccccccccccccc

      Double Precision A,Dummy(Maxproc),Ran_Uniform,T
      Integer          I,Sstmm

      Call Mpi_Init(Iierr)
      Call Mpi_Comm_Size(Mpi_Comm_World,Iisize,Iierr)
      Call Mpi_Comm_Rank(Mpi_Comm_World,Iirank,Iierr)

Ccccccccccccccccccccccccccc
C     Initialize Rng      C
Ccccccccccccccccccccccccccc

      A = (Dble(Iirank)+0.5d0)/(Dble(Iisize + Mod(Iabs(Sstmm()),3)))

      Call Genrand(A)

      Do I=1,(5 + 500*Iirank + Iabs(Sstmm()))
         A = Ran_Uniform()
      Enddo

      Do I=1,Maxproc
         Dummy(I) = 0.0d0
      Enddo
      
      Call Mpi_Gather(A,1,Mpi_Double_Precision,Dummy,1,
     &     Mpi_Double_Precision,0,Mpi_Comm_World,Iierr)
            
      If (Iirank.Eq.0) Then
         Write(6,*)
         Write(6,*)
         Write(6,*) 'Parallel Molecular Dynamics Program'
         Write(6,*)
         Write(6,*) 'Number Of Processors  : ',Iisize
         Write(6,*)
         
         Do I=1,Iisize
            Write(6,'(A,I5,A,E20.10)') ' Processor   ',I,
     &           ' Random Number          : ',Dummy(I)
         Enddo

         Write(6,*)
         Write(6,*)
      Endif
             
Cccccccccccccccccccccccccccccccccccccccccccccccc
C     Read Data From Disk                      C
Cccccccccccccccccccccccccccccccccccccccccccccccc
 
      Call Readdat
 
Cccccccccccccccccccccccccccccccccccccccccccccccc
C     Generate Initial Coordinates/Velocities  C
Cccccccccccccccccccccccccccccccccccccccccccccccc
 
      Call Init
 
Cccccccccccccccccccccccccccccccccccccccccccccccc
C     Finally, Perform An Md Simulation        C
Cccccccccccccccccccccccccccccccccccccccccccccccc
 
      T = Mpi_Wtime()

      Call Mdloop
 
      Call Mpi_Barrier(Mpi_Comm_World,Iierr)
            
      If(Iirank.Eq.0) Then
         Write(6,*) 'Elapsed time (s)      : ',Mpi_Wtime() - T
         Write(6,*)
      Endif

Cccccccccccccccccccccccccccccccccccccccccccccccc
C     End Of The Program                       C
Cccccccccccccccccccccccccccccccccccccccccccccccc
 
      Call Mpi_Finalize(Iierr)
      End
