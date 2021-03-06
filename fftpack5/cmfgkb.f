CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: cmfgkb.f,v 1.2 2004/06/15 21:08:32 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CMFGKB (LOT,IDO,IP,L1,LID,NA,CC,CC1,IM1,IN1,
     1                                      CH,CH1,IM2,IN2,WA)
      REAL       CH(2,IN2,L1,IDO,IP) ,CC(2,IN1,L1,IP,IDO),
     1                CC1(2,IN1,LID,IP)    ,CH1(2,IN2,LID,IP)  ,
     2                WA(IDO,IP-1,2)
C
C FFTPACK 5.0 auxiliary routine
C
      M1D = (LOT-1)*IM1+1
      M2S = 1-IM2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      DO 110 KI=1,LID
         M2 = M2S
         DO 110 M1=1,M1D,IM1
         M2 = M2+IM2
         CH1(1,M2,KI,1) = CC1(1,M1,KI,1)
         CH1(2,M2,KI,1) = CC1(2,M1,KI,1)
  110 CONTINUE
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 112 KI=1,LID
         M2 = M2S
         DO 112 M1=1,M1D,IM1
         M2 = M2+IM2
            CH1(1,M2,KI,J) =  CC1(1,M1,KI,J)+CC1(1,M1,KI,JC)
            CH1(1,M2,KI,JC) = CC1(1,M1,KI,J)-CC1(1,M1,KI,JC)
            CH1(2,M2,KI,J) =  CC1(2,M1,KI,J)+CC1(2,M1,KI,JC)
            CH1(2,M2,KI,JC) = CC1(2,M1,KI,J)-CC1(2,M1,KI,JC)
  112    CONTINUE
  111 CONTINUE
      DO 118 J=2,IPPH
         DO 117 KI=1,LID
         M2 = M2S
         DO 117 M1=1,M1D,IM1
         M2 = M2+IM2
            CC1(1,M1,KI,1) = CC1(1,M1,KI,1)+CH1(1,M2,KI,J)
            CC1(2,M1,KI,1) = CC1(2,M1,KI,1)+CH1(2,M2,KI,J)
  117    CONTINUE
  118 CONTINUE
      DO 116 L=2,IPPH
         LC = IPP2-L
         DO 113 KI=1,LID
         M2 = M2S
         DO 113 M1=1,M1D,IM1
         M2 = M2+IM2
            CC1(1,M1,KI,L) = CH1(1,M2,KI,1)+WA(1,L-1,1)*CH1(1,M2,KI,2)
            CC1(1,M1,KI,LC) = WA(1,L-1,2)*CH1(1,M2,KI,IP)
            CC1(2,M1,KI,L) = CH1(2,M2,KI,1)+WA(1,L-1,1)*CH1(2,M2,KI,2)
            CC1(2,M1,KI,LC) = WA(1,L-1,2)*CH1(2,M2,KI,IP)
  113    CONTINUE
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = MOD((L-1)*(J-1),IP)
            WAR = WA(1,IDLJ,1)
            WAI = WA(1,IDLJ,2)
            DO 114 KI=1,LID
               M2 = M2S
               DO 114 M1=1,M1D,IM1
               M2 = M2+IM2
               CC1(1,M1,KI,L) = CC1(1,M1,KI,L)+WAR*CH1(1,M2,KI,J)
               CC1(1,M1,KI,LC) = CC1(1,M1,KI,LC)+WAI*CH1(1,M2,KI,JC)
               CC1(2,M1,KI,L) = CC1(2,M1,KI,L)+WAR*CH1(2,M2,KI,J)
               CC1(2,M1,KI,LC) = CC1(2,M1,KI,LC)+WAI*CH1(2,M2,KI,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      IF(IDO.GT.1 .OR. NA.EQ.1) GO TO 136
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 KI=1,LID
         DO 119 M1=1,M1D,IM1
            CHOLD1 = CC1(1,M1,KI,J)-CC1(2,M1,KI,JC)
            CHOLD2 = CC1(1,M1,KI,J)+CC1(2,M1,KI,JC)
            CC1(1,M1,KI,J) = CHOLD1
            CC1(2,M1,KI,JC) = CC1(2,M1,KI,J)-CC1(1,M1,KI,JC)
            CC1(2,M1,KI,J) = CC1(2,M1,KI,J)+CC1(1,M1,KI,JC)
            CC1(1,M1,KI,JC) = CHOLD2
  119    CONTINUE
  120 CONTINUE
      RETURN
  136 DO 137 KI=1,LID
         M2 = M2S
         DO 137 M1=1,M1D,IM1
         M2 = M2+IM2
         CH1(1,M2,KI,1) = CC1(1,M1,KI,1)
         CH1(2,M2,KI,1) = CC1(2,M1,KI,1)
  137 CONTINUE
      DO 135 J=2,IPPH
         JC = IPP2-J
         DO 134 KI=1,LID
         M2 = M2S
         DO 134 M1=1,M1D,IM1
         M2 = M2+IM2
            CH1(1,M2,KI,J) = CC1(1,M1,KI,J)-CC1(2,M1,KI,JC)
            CH1(1,M2,KI,JC) = CC1(1,M1,KI,J)+CC1(2,M1,KI,JC)
            CH1(2,M2,KI,JC) = CC1(2,M1,KI,J)-CC1(1,M1,KI,JC)
            CH1(2,M2,KI,J) = CC1(2,M1,KI,J)+CC1(1,M1,KI,JC)
  134    CONTINUE
  135 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 131 I=1,IDO
         DO 130 K=1,L1
         M2 = M2S
         DO 130 M1=1,M1D,IM1
         M2 = M2+IM2
            CC(1,M1,K,1,I) = CH(1,M2,K,I,1)
            CC(2,M1,K,1,I) = CH(2,M2,K,I,1)
  130    CONTINUE
  131 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
         M2 = M2S
         DO 122 M1=1,M1D,IM1
         M2 = M2+IM2
            CC(1,M1,K,J,1) = CH(1,M2,K,1,J)
            CC(2,M1,K,J,1) = CH(2,M2,K,1,J)
  122    CONTINUE
  123 CONTINUE
      DO 126 J=2,IP
         DO 125 I=2,IDO
            DO 124 K=1,L1
               M2 = M2S
               DO 124 M1=1,M1D,IM1
               M2 = M2+IM2
               CC(1,M1,K,J,I) = WA(I,J-1,1)*CH(1,M2,K,I,J)
     1                      -WA(I,J-1,2)*CH(2,M2,K,I,J)
               CC(2,M1,K,J,I) = WA(I,J-1,1)*CH(2,M2,K,I,J)
     1                      +WA(I,J-1,2)*CH(1,M2,K,I,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
      END
