Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Npart    = Number Of Particles                  C
C     Nchoi    = Number Of Trial Dirs.                C
C     Nlength  = Length Of The Chain                  C
C     Lchain   = Trial Chain Or Real Chain            C
C     Box      = Boxlength                            C
C     Ibox     = 1/Box                                C
C     Xx/Yy/Zz = Coordinates Of Particles             C
C     X/Y/Z    = Coordinates Of Chain                 C
C     Alpha    = Repulsion Parameter                  C
C     Beta     = 1/(K_B T)                            C
C     Usim     = Running Energy                       C
C     Dispmax  = Maximum Displacement                 C
C     Pchain   = Prob. To Perform A Cbmc Trialmove    C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      Integer Npart,Nchoi,Nlength
      Logical Lchain

      Double Precision Box,Ibox,Xx(Maxpart),Yy(Maxpart),
     &     Zz(Maxpart),X(Maxlength),Y(Maxlength),
     &     Z(Maxlength),Alpha,Beta,Usim,Dispmax,Pchain

      Common /Ssyss/ Pchain,Dispmax,Beta,Usim,Alpha,
     &     X,Y,Z,Xx,Yy,Zz,Box,Ibox,Npart,Nchoi,
     &     Nlength,Lchain
