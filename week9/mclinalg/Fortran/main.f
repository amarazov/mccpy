      PROGRAM GENSSOB
C
C       IT USES A NONSTANDARD TIMING
C       ROUTINE "DTIME"
C
C       THIS PROGRAM TESTS ACCURACY OF
C       NUMERICAL INTEGRATION USING "GOSNIEDX"
C       AND INTEGRAND (2) OF DAVIS AND
C       RABINOWITZ, PAGE 406
C
C      User Define:
C        DIMEN : dimension 
C        ATMOST : sequence length
C        SAM : Number of replications
C        MAX : Maximum Digits of Scrambling Of Owen type Scrambling
C        IFLAG: User Choice of Sequences
C        IFLAG = 0 : No Scrambling
C        IFLAG = 1 : Owen type Scrambling
C        IFLAG = 2 : Faure-Tezuka type Scrambling
C        IFLAG = 3 : Owen + Faure-Tezuka type Scrambling
C      
      LOGICAL FLAG(2)
      INTEGER DIMEN,ATMOST,I,II,J,TAUS
      INTEGER IFLAG,MAX,SAM
      REAL*8 QUASI(40),SUM,F,EI
      REAL*4 TIM, tarray(2)

      SAM = 1
      MAX = 30
      DIMEN = 2
      IFLAG = 1
      ATMOST = 2**12
      DO 22 II = 1,SAM
         WRITE(*,*) 'I = ITERATION NUMBER'
         WRITE(*,*) 'EI = ESTIMATED INTEGRAL'
         TIM = DTIME(tarray)
         CALL INSSOBL(FLAG,DIMEN,ATMOST,TAUS,QUASI,MAX,IFLAG)
         F= 1.0
         DO 55 J = 1,DIMEN
              F = F*ABS(4.0*QUASI(J)-2.0)
 55         CONTINUE
         SUM = F
         DO 100 I=2,ATMOST
             CALL GOSSOBL(QUASI)
             F = 1.0
             DO 50 J=1,DIMEN
                F = F*ABS(4.0*QUASI(J)-2.0)
 50           CONTINUE
                IF (MOD(I,500) .EQ. 0) THEN
                  WRITE(*,*) 'I = ',I
                  WRITE(*,*) 'EI = ',SUM/I
                ENDIF
             SUM = SUM + F
 100       CONTINUE
         TIM = DTIME(tarray)
         WRITE(*,*) 'Total time elapsed = ', tarray(1)
 22    CONTINUE
      STOP
      END
