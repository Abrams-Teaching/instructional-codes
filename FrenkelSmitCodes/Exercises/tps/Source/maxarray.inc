Cccccccccccccccccccccccccccccc
C     Maxarray.Inc           C
Cccccccccccccccccccccccccccccc
C
C     Maxatom    = Number Of Particles
C     Maxtraject = Length Of The Path
C     Maxumbrel  = Maximum Number Of Slices In Histogram
C     Maxdslice  = Resolution
C     Nshort     = Number Of Md Steps Between Two Points
C     W          = Potential Parameter
C     Rad        = Radius Of The Circle
q
C     Left/Right = Identity Of A And B

      Integer Maxatom,Maxtraject,Maxumbrel,Nshort,Maxinter

      Double Precision Maxdslice,W,Left,Right,Rad

      Parameter (Maxatom    = 15              )
      Parameter (Maxinter   = Maxatom*Maxatom )
      Parameter (Maxtraject = 1000            )
      Parameter (Maxumbrel  = 2000            )
      Parameter (Maxdslice  = 0.002d0         )
      Parameter (Nshort     = 5               )
      Parameter (Rad        = 3.0d0           )
      Parameter (W          = 0.25d0          )
      Parameter (Left       = 1.30d0          )
      Parameter (Right      = 1.45d0          )
