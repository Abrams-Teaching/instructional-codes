Ccccccccccccccccccccccccc
C     Traject.Inc       C
Ccccccccccccccccccccccccc
c
C     Xxtra/Yytra   = Position Of New Path
C     Vxtra/Vytra   = Velocity Of New Path
C     Xxold/Yyold   = Position Of Old Path
C     Vxold/Vyold   = Velocity Of Old Path
C     Eeold         = Energy Oldpath
C     Eetra         = Energy New Path
C     La/Lb         = Identity New Path
C     La_Old/Lb_Old = Identity Old Path


      Double Precision Xxtra(Maxatom,Maxtraject),
     &                 Yytra(Maxatom,Maxtraject),
     &                 Vxtra(Maxatom,Maxtraject),
     &                 Vytra(Maxatom,Maxtraject),
     &                 Xxold(Maxatom,Maxtraject),
     &                 Yyold(Maxatom,Maxtraject),
     &                 Vxold(Maxatom,Maxtraject),
     &                 Vyold(Maxatom,Maxtraject),
     &                 Eeold(Maxtraject),
     &                 Eetra(Maxtraject)

      Logical          Lb(Maxtraject),Lb_Old(Maxtraject),
     &                 La(Maxtraject),La_Old(Maxtraject)

      Common /Trassw/  Xxtra,Yytra,Vxtra,Vytra,Xxold,Yyold,
     &                 Vxold,Vyold,Eeold,Eetra,La,Lb,La_Old,
     &                 Lb_Old
