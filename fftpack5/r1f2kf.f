CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: r1f2kf.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE R1F2KF (IDO,L1,CC,IN1,CH,IN2,WA1)
      REAL       CH(IN2,IDO,2,L1) ,CC(IN1,IDO,L1,2) , WA1(IDO)
C
      DO 101 K=1,L1
         CH(1,1,1,K) = CC(1,1,K,1)+CC(1,1,K,2)
         CH(1,IDO,2,K) = CC(1,1,K,1)-CC(1,1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(1,I,1,K) = CC(1,I,K,1)+(WA1(I-2)*CC(1,I,K,2)-
     1       WA1(I-1)*CC(1,I-1,K,2))
            CH(1,IC,2,K) = (WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*
     1       CC(1,I-1,K,2))-CC(1,I,K,1)
            CH(1,I-1,1,K) = CC(1,I-1,K,1)+(WA1(I-2)*CC(1,I-1,K,2)+
     1       WA1(I-1)*CC(1,I,K,2))
            CH(1,IC-1,2,K) = CC(1,I-1,K,1)-(WA1(I-2)*CC(1,I-1,K,2)+
     1       WA1(I-1)*CC(1,I,K,2))
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(1,1,2,K) = -CC(1,IDO,K,2)
         CH(1,IDO,1,K) = CC(1,IDO,K,1)
  106 CONTINUE
  107 RETURN
      END
