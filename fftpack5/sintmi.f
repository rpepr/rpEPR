CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: sintmi.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SINTMI (N, WSAVE, LENSAV, IER)
      INTEGER    N, LENSAV, IER
      REAL       WSAVE(LENSAV)
C
      IER = 0
C
      IF (LENSAV .LT. N/2 + N + INT(LOG(REAL(N))) +4) THEN
        IER = 2
        CALL XERFFT ('SINTMI', 3)
        GO TO 300
      ENDIF
C
      PI = 4.*ATAN(1.)
      IF (N .LE. 1) RETURN
      NS2 = N/2
      NP1 = N+1
      DT = PI/FLOAT(NP1)
      DO 101 K=1,NS2
         WSAVE(K) = 2.*SIN(K*DT)
  101 CONTINUE
      LNSV = NP1 + INT(LOG(REAL(NP1))) +4
      CALL RFFTMI (NP1, WSAVE(NS2+1), LNSV, IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINTMI',-5)
      ENDIF
C
  300 CONTINUE
      RETURN
      END
