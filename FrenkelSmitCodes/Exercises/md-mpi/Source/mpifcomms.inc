Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Mpifcomms.Inc                                      C
C                                                        C
C     Variables For Mpi                                  C
C                                                        C
C     Iierr    = Error Code                              C
C     Iisize   = Number Of Processors                    C
C     Iirank   = Index Of A Processor (0...Iisize-1)     C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer           Iierr,Iisize,Iirank

      Common /Cmns11/   Iierr,Iisize,Iirank
