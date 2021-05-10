      Function Dderfc(X)
      Implicit None
 
Cccccccccccccccccccccccccccccccccccccc
C     Complementary Error Function   C
Cccccccccccccccccccccccccccccccccccccc

      Double Precision Pa,P0,P1,P2,P3,P4,P5,P6,
     &     P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,
     &     P17,P18,P19,P20,P21,P22,T,U,Y,X,
     &     Dderfc
 
      Parameter (Pa=3.97886080735226000d+00)
      Parameter (P0=2.75374741597376782d-01)
      Parameter (P1=4.90165080585318424d-01)
      Parameter (P2=7.74368199119538609d-01)
      Parameter (P3=1.07925515155856677d+00)
      Parameter (P4=1.31314653831023098d+00)
      Parameter (P5=1.37040217682338167d+00)
      Parameter (P6=1.18902982909273333d+00)
      Parameter (P7=8.05276408752910567d-01)
      Parameter (P8=3.57524274449531043d-01)
      Parameter (P9=1.66207924969367356d-02)
      Parameter (P10=-1.19463959964325415d-01)
      Parameter (P11=-8.38864557023001992d-02)
      Parameter (P12=2.49367200053503304d-03)
      Parameter (P13=3.90976845588484035d-02)
      Parameter (P14=1.61315329733252248d-02)
      Parameter (P15=-1.33823644533460069d-02)
      Parameter (P16=-1.27223813782122755d-02)
      Parameter (P17=3.83335126264887303d-03)
      Parameter (P18=7.73672528313526668d-03)
      Parameter (P19=-8.70779635317295828d-04)
      Parameter (P20=-3.96385097360513500d-03)
      Parameter (P21=1.19314022838340944d-04)
      Parameter (P22=1.27109764952614092d-03)
 
      T = Pa/(Pa+Dabs(X))
 
      U = T - 0.5d0
 
      Y = (((((((((P22*U+P21)*U+P20)*U+P19)*U+
     &     P18)*U+P17)*U+P16)*U+P15)*U+P14)
     &     *U+P13)*U + P12
 
      Y = ((((((((((((Y*U+P11)*U+P10)*U+P9)*U+P8)*U+
     &     P7)*U+P6)*U+P5)*U+P4)*U+P3)*
     &     U+P2)*U+P1)*U+P0)*T*Dexp(-X*X)
 
      If (X.Lt.0.0d0) Y = 2.0d0 - Y
 
      Dderfc = Y
 
      Return
      End
