      PROGRAM GENSFAUR
C
C       THIS PROGRAM TESTS ACCURACY OF
C       NUMERICAL INTEGRATION USING "GOSFAUR"
C       AND INTEGRAND (2) OF DAVIS AND
C       RABINOWITZ, PAGE 406
C
C       IT USES A NONSTANDARD TIMING
C       ROUTINE "DTIME"
C
C       User Define: 
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
      INTEGER*4 DIMEN,ATMOST,II,I,J,QS,SAM,MAX,MAXX
      INTEGER IFLAG
      REAL*8 QUASI(500),F,SUM,EI
      REAL*4 TIM,tarray(2)
      COMMON /FAURE/ QS
 
      IFLAG = 3
      DIMEN = 2
      SAM = 1      
      MAX= 30
      ATMOST = 2**12
      DO 10 II = 1,SAM
          TIM = DTIME(tarray)
          WRITE(*,*) 'I = ITERATION NUMBER'
          WRITE(*,*) 'EI = ESTIMATED INTEGRAL'       
          CALL INSFAUR(FLAG,DIMEN,ATMOST,QUASI,MAX,IFLAG,MAXX)
            IF(.NOT. FLAG(2)) THEN
              WRITE(*,*) 'ATMOST = ',ATMOST
              WRITE(*,*) 'ATMOST IS NOT OK'
              RETURN
            ENDIF
          WRITE(*,*) 'I = ITERATION NUMBER'
          F= 1.0           
          DO 20 J = 1,DIMEN
                F = F*ABS(4.0*QUASI(J)-2.0)
 20            CONTINUE 
          SUM = F             
          DO 30 I=2,ATMOST
               CALL GOSFAUR(QUASI,MAXX)
               F = 1.0 
               DO 40 J=1,DIMEN
                  F = F*ABS(4.0*QUASI(J)-2.0)
40               CONTINUE
               SUM = SUM+F
               IF (MOD(I,500) .EQ. 0) THEN
                  WRITE(*,*) 'I = ',I
                  WRITE(*,*) 'EI = ',SUM/I
               ENDIF   
30          CONTINUE
           TIM = DTIME(tarray)
           WRITE(*,*) 'Total time elapsed = ', tarray(1) 
10       CONTINUE
      STOP
      END
