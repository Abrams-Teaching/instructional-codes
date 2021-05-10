Cccccccccccccccccccccccccccccc
C     Mdtop.Inc              C
Cccccccccccccccccccccccccccccc
C
C     Rwca      = 2**(1/6)
C     Rwca2     = Rwca**2
C     H         = Height Of The Barrier
C     Tstep     = Timestep Md
C     Left2     = Left**2
C     Right2    = Right**2
C     Iw        = 1/W
C     Upath     = Total Energy
C     Lmin      = Left Boundary Slice
C     Lmax      = Right Boundary Slice
C     Lmin2     = Lmin**2
C     Lmax2     = Lmax**2
C     Twopi     = 2*PI
C     Deltaphi  = Rotation Angle
C     Nstep     = Number Of Steps
C     Nslice    = Length Of Thepath
C     Natom     = Numberof Atoms
C     Ninter    = Number Of Interactions
C     Iinter    = Identity I
C     Jinter    = Identity J
C     Lensemble = Md Of Tps ?
C     Lumbrella = Umbrella Sampling ?
C     Lscreen   = Write To Screen During Md ?

      Double Precision Rwca,Rwca2,H,Rad2,Tstep,Left2,Right2,
     &                 Lmin,Lmax,Iw,Upath,Lmin2,Lmax2,Twopi,
     &                 Deltaphi

      Integer          Nstep,Nslice,Natom,Ninter,Iinter(Maxinter),
     &                 Jinter(Maxinter)

      Logical          Lensemble,Lumbrella,Lscreen

      Common /Mdtp1/   Rwca,Rwca2,H,Rad2,Tstep,Left2,Right2,Lmin,
     &                 Twopi,Deltaphi,Lmax,Iw,Upath,Lmin2,Lmax2,
     &                 Nstep,Nslice,Natom,Ninter,Iinter,Jinter,
     &                 Lensemble,Lumbrella,Lscreen
