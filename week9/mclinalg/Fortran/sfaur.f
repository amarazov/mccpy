C      THIS MODIFIED ALGORITHM 659, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 14, NO. 1, P.88.
      
      SUBROUTINE INSFAUR(FLAG,DIMEN,ATMOST,QUASI,MAX,IFLAG,MAXX)
C
C       THIS SUBROUTINE FIRST CHECKS WHETHER
C       THE USER-SUPPLIED DIMENSION "DIMEN" OF THE
C       QUASIRANDOM VECTORS IS ACCEPTABLE
C       (STRICTLY BETWEEN 1 AND 41) : IF SO,
C       FLAG(1)=.TRUE.
C 
C       THEN IT CALCULATES AN UPPER SUMMATION
C       LIMIT "HISUM" BASED ON "DIMEN" AND THE
C       USER-SUPPLIED NUMBER "ATMOST" OF QUASIRANDOM
C       VECTORS REQUIRED. FLAG(2)=.TRUE. IF
C       ATMOST IS OK.
C
C       IF FLAG(1) AND FLAG(2) ARE TRUE,
C       "INFAUR" NEXT PRODUCES THE OTHER
C       OUTPUTS LISTED BELOW PASSED TO
C       SUBROUTINE GOFAUR VIA LABELLED
C       COMMON "FAURE". THESE OUTPUTS ARE
C       IRRELEVANT TO THE USER.
C
C       FIRST CALL INSFAUR. IF FLAG(1) AND
C       FLAG(2) ARE TRUE, EACH (SUBSEQUENT)
C       CALL TO GOSFAUR GENERATES A NEW
C       QUASIRANDOM VECTOR.
C
C       INPUTS : DIMEN, ATMOST, MAX
C
C       OUTPUTS
C          TO USERS CALLING PROGRAM:
C             QUASI :FIRST SEQUENCE
C             FLAG
C             QSS   : SAME AS QS - SEE BELOW
C
C
C          TO GOSFAUR:
C             S      :DIMENSION
C             QS     :SMALLEST PRIME >=S
C             SCOEF  :SCRAMBLING GENERATING MATRICES
C             NEXTN  :THE NUMBER OF THE
C                     NEXT QUASIRANDOM
C                     VECTOR,INITIALIZED
C                     TO TESTN-1 HERE.
C             TESTN  :INITIALIZED TO 0
C             HISUM  :AFTER BEING USED TO
C                     PRODUCE COEF, INITIALIZED
C                     TO 3 FOR GOFAUR.
C             RQS    :1.0/QS.
C
      LOGICAL FLAG(2)
C
      INTEGER*4 S,ATMOST,QS,COEF(1:500,0:29,0:29),NEXTN,
     +          TESTN,HISUM,I,J,PRIMES(500),DIMEN,L,TEMS,
     +          PGTEMP(0:29),ZTEMP(1:500,0:29),MAX,IN,
     +          LSM(1:500,0:29,0:29),SHIFT(1:500,0:29),
     +          SCOEF(1:500,0:29,0:29),USHIFT(0:29),MAXX,
     +          USM(0:29,0:29),TEMS2,TSCOEF(1:500,0:29,0:29)
      INTEGER IFLAG
C
      REAL*8 RQS,QUASI(500),R
C
      COMMON /FAURE/ S,QS,HISUM,SCOEF,NEXTN,TESTN,
     +               RQS,PGTEMP,ZTEMP,COEF
      SAVE /FAURE/
C
      DATA (PRIMES(I),I=1,500)/1,2,3,5,5,7,7,11,11,11,
     +                         11,13,13,17,17,17,17,19,19,23,
     +                         23,23,23,29,29,29,29,29,29,31,
     +                         31,37,37,37,37,37,37,41,41,41,
     +                         41,43,43,47,47,47,47,53,53,53,
     +                         53,53,53,59,59,59,59,59,59,61,
     +                         61,67,67,67,67,67,67,71,71,71,
     +                         71,73,73,79,79,79,79,79,79,83,
     +                         83,83,83,89,89,89,89,89,89,97,
     +                         97,97,97,97,97,97,97,101,101,101,
     +                         101,103,103,107,107,107,107,109,109,113,
     +                         113,113,113,127,127,127,127,127,127,127,
     +                         127,127,127,127,127,127,127,131,131,131,
     +                         131,137,137,137,137,137,137,139,139,149,
     +                         149,149,149,149,149,149,149,149,149,151,
     +                         151,157,157,157,157,157,157,163,163,163,
     +                         163,163,163,167,167,167,167,173,173,173,
     +                         173,173,173,179,179,179,179,179,179,181,
     +                         181,191,191,191,191,191,191,191,191,191,
     +                         191,193,193,197,197,197,197,199,199,211,
     +                         211,211,211,211,211,211,211,211,211,211,
     +                         211,223,223,223,223,223,223,223,223,223,
     +                         223,223,223,227,227,227,227,229,229,233,
     +                         233,233,233,239,239,239,239,239,239,241,
     +                         241,251,251,251,251,251,251,251,251,251,
     +                         251,257,257,257,257,257,257,263,263,263,
     +                         263,263,263,269,269,269,269,269,269,271,
     +                         271,277,277,277,277,277,277,281,281,281,
     +                         281,283,283,293,293,293,293,293,293,293,
     +                         293,293,293,307,307,307,307,307,307,307,
     +                         307,307,307,307,307,307,307,311,311,311,
     +                         311,313,313,317,317,317,317,331,331,331,
     +                         331,331,331,331,331,331,331,331,331,331,
     +                         331,337,337,337,337,337,337,347,347,347,
     +                         347,347,347,347,347,347,347,349,349,353,
     +                         353,353,353,359,359,359,359,359,359,367,
     +                         367,367,367,367,367,367,367,373,373,373,
     +                         373,373,373,379,379,379,379,379,379,383,
     +                         383,383,383,389,389,389,389,389,389,397,
     +                         397,397,397,397,397,397,397,401,401,401,
     +                         401,409,409,409,409,409,409,409,409,419,
     +                         419,419,419,419,419,419,419,419,419,421,
     +                         421,431,431,431,431,431,431,431,431,431,
     +                         431,433,433,439,439,439,439,439,439,443,
     +                         443,443,443,449,449,449,449,449,449,457,
     +                         457,457,457,457,457,457,457,461,461,461,
     +                         461,463,463,467,467,467,467,479,479,479,
     +                         479,479,479,479,479,479,479,479,479,487,
     +                         487,487,487,487,487,487,487,491,491,491,
     +                         491,499,499,499,499,499,499,499,499,503/

C
C       CHECK S
C
      S=DIMEN
      FLAG(1) = S.GT.1 .AND. S.LT.41
      IF (.NOT.FLAG(1)) RETURN
C
      TESTN = 0
      QS=PRIMES(S)

C
C         COMPUTE LOG(ATMOST+TESTN) IN BASE QS
C         USING A RATIO OF NATURAL LOGS TO GET
C         AN UPPER BOUND ON (THE NUMBER OF
C         DIGITS IN THE BASE QS REPRESENTATION
C         OF ATMOST+TESTN) MINUS ONE.
C
      HISUM=NINT(LOG(REAL(ATMOST+TESTN))/LOG(REAL(QS)))
      
      FLAG(2) = HISUM .LT. 30
      IF(.NOT. FLAG(2)) RETURN
C
C        NOW FIND BINOMIAL COEFFICIENTS MOD QS
C        IN A UPPER-TRIANGULAR MATRIX "COEF"
C        USING RECURSION BINOM(I,J)=BINOM(I,J-1)
C        +BINOM(I-1,J-1) AND A=B+C IMPLIES MOD(A,D)=
C        MOD(MOD(B,D)+MOD(C,D),D)
C
C CONSTRUCTING UPPER-TIRANGULAR GENERATOR MATRICES
C UP TO DIMENSION S
      
      DO 555 I = 0,HISUM
          DO 556 J = 0,HISUM
             IF (I .EQ. J) THEN
                COEF(1,I,J) = 1
             ELSE
                COEF(1,I,J) = 0
             ENDIF
556        CONTINUE
555      CONTINUE
               
      COEF(2,0,0)=1
      DO 50 J=1,HISUM
        COEF(2,0,J)=1
        COEF(2,J,J)=1
   50 CONTINUE
   
      DO 200 I=1,HISUM
        DO 100 J=I+1,HISUM
          COEF(2,I,J)=MOD(COEF(2,I,J-1)+COEF(2,I-1,J-1),QS)
  100   CONTINUE
  200 CONTINUE
       
     
      DO 300 K=3,S
        DO 400 J=0,HISUM
         DO 500 I=0,HISUM
               COEF(K,I,J) = 0
              IF (J .LT. I) THEN
                  goto 500
              ELSE   
                 TEMS=0
              ENDIF
              DO 600 L=0,J
                  TEMS=TEMS+COEF(K-1,I,L)*COEF(2,L,J)          
 600             CONTINUE           
               COEF(K,I,J)=MOD(TEMS,QS)
 500         CONTINUE
 400      CONTINUE  
 300    CONTINUE

      
C
C GENERATE USER CHOICE OF GENERATOR MATRICES
C
      IF (IFLAG .EQ. 0) THEN
        MAXX = HISUM
        DO 699 K = 1,S
           DO 799 J = 0, HISUM
              DO 899 I = 0, HISUM
                    SCOEF(K,I,J) = COEF(K,I,J) 
 899            CONTINUE
                 SHIFT(K,J) = 0
 799        CONTINUE
 699     CONTINUE             
      ELSE
        IF ((IFLAG .EQ. 1) .OR. (IFLAG .EQ. 3)) THEN
          MAXX = MAX
          CALL GENSCRML(LSM,SHIFT,MAX)
          DO 700 K=1,S
             DO 800 J=0,HISUM
              DO 900 I=0,MAX
                IF (K. EQ. 1) THEN
                   SCOEF(K,I,J) = LSM(K,I,J)
                  IF (IFLAG .EQ. 3) THEN 
                   TSCOEF(K,I,J) = SCOEF(K,I,J)
                  ENDIF 
                ELSE  
                TEMS = 0
                 IF (I .GE. HISUM) THEN
                   IN = HISUM
                ELSE
                   IN = J   
                ENDIF   
                DO 901 L=0,IN
                   TEMS=TEMS+LSM(K,I,L)*COEF(K,L,J) 
 901             CONTINUE           
                 SCOEF(K,I,J)=MOD(TEMS,QS)
                 IF (IFLAG .EQ. 3) THEN
                   TSCOEF(K,I,J) = SCOEF(K,I,J)
                 ENDIF 
                ENDIF 
 900          CONTINUE
 800        CONTINUE 
 700      CONTINUE
         ENDIF

         IF ((IFLAG .EQ. 2) .OR. (IFLAG .EQ. 3)) THEN
           CALL GENSCRMU(USM,USHIFT,HISUM)
           IF (IFLAG .EQ. 2) THEN
               MAXX = HISUM
           ENDIF
           DO 701 K=1,S
             DO 801 J=0,HISUM
               DO 902 I=0,MAXX
                TEMS = 0
                TEMS2 = 0
                IN = HISUM   
                DO 903 L=0,IN
                   IF(IFLAG .EQ. 2) THEN
                      TEMS=TEMS+COEF(K,I,L) *USM(L,J) 
                     IF (J .EQ. 0) THEN
                      TEMS2 = TEMS2 +COEF(K,I,L)*USHIFT(L)
                     ENDIF 
                   ENDIF
                   IF(IFLAG .EQ. 3) THEN
                     TEMS=TEMS+TSCOEF(K,I,L)*USM(L,J) 
                     IF (J .EQ. 0) THEN
                      TEMS2 = TEMS2 +TSCOEF(K,I,L)*USHIFT(L)
                     ENDIF 
                   ENDIF  
 903             CONTINUE           
                 SCOEF(K,I,J)=MOD(TEMS,QS)
                 IF (J .EQ. 0) THEN
                   TEMS2 = MOD(TEMS2,QS)
                   IF (IFLAG .EQ. 3) THEN
                     SHIFT(K,I) = MOD((TEMS2+SHIFT(K,I)),QS)
                   ELSE
                     SHIFT(K,I) = TEMS2 
                   ENDIF    
                ENDIF 
 902         CONTINUE
 801        CONTINUE 
 701       CONTINUE
         ENDIF
       ENDIF  
C        CALCULATING THESE COEFFICIENTS
C        MOD QS AVOIDS POSSIBLE OVERFLOW
C        PROBLEMS WITH RAW BINOMIAL COEFFICIENTS
C

      RQS=1.0/REAL(QS)      

C     GENERATE FIRST SCRAMBLED QUASI RANDOM VECTOR

      DO 1000 K = 1,S
          R = 0.0
          DO 1100 I = MAXX,0,-1
             ZTEMP(K,I) = SHIFT(K,I)
             R = ZTEMP(K,I)+RQS*R            
1100        CONTINUE
          QUASI(K) = R*RQS
1000     CONTINUE 
C
C     SET UP FIRST VECTOR AND VALUES FOR GOSFAUR
C      
      DO 1200 I = 0,HISUM
          PGTEMP(I) = 0
1200    CONTINUE         

      TESTN= QS
      NEXTN= 1
      HISUM=0

C     NOW COMPLETE INITIALIZATION
C     AS DESCRIBED IN SECTION 4.
C     NEXTN HAS 1 DIGITS IN BASE
C     QS, SO HISUM EQUALS 0.
C      
      RETURN
      END

      INTEGER FUNCTION EXOR (IIN, JIN)
      INTEGER*4 I,J,K,L,IIN,JIN
C
C     THIS FUNCTION CALCULATES THE EXCLUSIVE-OR OF ITS
C     TWO INPUT PARAMETERS
C
      I = IIN
      J = JIN
      K = 0
      L = 1
C
    1 IF (I .EQ. J) THEN
        EXOR = K
        RETURN
      ENDIF
C
C     CHECK THE CURRENT RIGHT-HAND BITS OF I AND J.
C     IF THEY DIFFER, SET THE APPROPRIATE BIT OF K.
C
      IF (MOD(I,2) .NE. MOD(J,2)) K = K + L
      I = I / 2
      J = J / 2
      L = 2 * L
      GOTO 1
      END
      SUBROUTINE GENSCRML(LSM,SHIFT,MAX)
      
      INTEGER*4 LSM(1:500,0:29,0:29),QS,HISUM,P,I,J,QSM1
      INTEGER*4 SHIFT(1:500,0:29),MAX,S
      REAL*8 UNI
      COMMON /FAURE/ S,QS,HISUM
      SAVE /FAURE/
      
      QSM1 = QS-1
      DO 10 P = 1,S
         DO 20 I = 0,MAX
               SHIFT(P,I) = MOD((int(UNI()*1000.0)),QS)
            DO 30 J = 0,HISUM
               IF (J .EQ. I) THEN
                LSM(P,I,J) = MOD((int(UNI()*1000.0)),
     *                          QSM1) +1
               ELSE IF (J .LT. I)  THEN 
                LSM(P,I,J) = MOD((int(UNI()*1000.0)),QS)
               ELSE
                LSM(P,I,J) = 0
               ENDIF
 30           CONTINUE
 20        CONTINUE
 10      CONTINUE
 
      RETURN
      END   
      
      SUBROUTINE GENSCRMU(USM,USHIFT,MAX)

C     GENERATING LOWER TRIANGULAR SCRMABLING MATRICES AND 
C     SHIFT VECTORS.
C     INPUTS :
C       FROM INSSOBL : MAX
C       FROM BLOCK DATA "SOBOL" : S, MAXCOL,
C
C     OUTPUTS :
C       TO INSSOBL : USM, USHIFT

      INTEGER*4 USM(0:29,0:29),QS,HISUM,P,I,J,QSM1
      INTEGER*4 USHIFT(0:29),MAX,S,TEMP,STEMP,L,LL
      REAL*8 UNI
      COMMON /FAURE/ S,QS,HISUM
      SAVE /FAURE/
     
      QSM1 = QS-1 
          DO 20 I = 0,HISUM
               STEMP = MOD((int(UNI()*1000.0)),QS)
C               STEMP = 0
               USHIFT(I) = STEMP               
            DO 30 J = 0,HISUM 
               IF (J .EQ. I) THEN
                 TEMP = MOD((int(UNI()*1000.0)),
     *                          QSM1) +1
C                  TEMP = 1
               ELSE IF (J .GT. I)  THEN 
C                 TEMP = MOD((int(UNI()*1000.0)),QS)
                 TEMP = 0
               ELSE
                TEMP = 0
               ENDIF
               USM(I,J) = TEMP        
 30           CONTINUE
 20        CONTINUE 
      RETURN
      END   
                          
      SUBROUTINE GOSFAUR(QUASI,MAX)
C
C       THIS SUBROUTINE GENERATES A NEW
C       QUASIRANDOM VECTOR WITH EACH CALL BY USING THE ALPHA-ARY
C       GRAY CODE METHOD OF LICHTNER "SIAM J. DISCRETE MATH(1998),
C       381-386". (SEE ESPECIALLY PAGE 381-382).
C   
C       THE USER MUST CALL "INSFAUR" BEFORE
C       CALLING "GOSFAUR".
C       AFTER CALLING "INFAUR", TEST FLAG(1)
C       AND FLAG(2); IF EITHER IS FALSE, DO
C       NOT CALL GOFAUR. READ THE COMMENTS AT
C       THE BEGINNING OF INSFAUR AND THEN
C       THOSE BELOW.
C
C       ALL INPUTS COME FROM "INFAUR" VIA
C       LABELLED COMMON "FAURE"; FOR THEIR
C       DEFINITIONS, SEE "INSFAUR".
C
C       INPUTS:
C         S,QS,COEF,NEXTN,TESTN,HISUM,RQS,MAX
C
C       OUTPUTS:
C         TO USER'S CALLING PROGRAM:
C         QUASI - A NEW SCRAMBLED QUASIRANDOM VECTOR
C
      
      INTEGER*4 S,QS,COEF(1:500,0:29,0:29),NEXTN,TESTN
      INTEGER*4 HISUM,I,J,K,YTEMP(0:19),ZTEMP(1:500,0:29)
      INTEGER*4 GTEMP(0:29),KTEMP,LTEMP,MTEMP
      INTEGER*4 TEM,PGTEMP(0:29),DIFF,MAX,POS
      INTEGER*4 SCOEF(1:500,0:29,0:29)
C
      REAL*8 QUASI(500),RQS,R
C
      COMMON /FAURE/ S,QS,HISUM,SCOEF,NEXTN,TESTN,
     +               RQS,PGTEMP,ZTEMP,COEF
      SAVE /FAURE/
C
C
C       NEXTN HAS A REPRESENTATION IN BASE
C       QS OF THE FORM: SUM OVER J FROM ZERO
C       TO HISUM OF YTEMP(J)*(QS**J)
C
C       WE NOW COMPUTE THE YTEMP(J)'S SAME AS FAURE 
C
C
      KTEMP=TESTN
      LTEMP=NEXTN 
      DO 100 I=HISUM,0,-1
             KTEMP=KTEMP/QS  
             MTEMP=MOD(LTEMP,KTEMP)
             YTEMP(I)=(LTEMP-MTEMP)/KTEMP     
             LTEMP=MTEMP
100    CONTINUE
C
C PERFORM THE CONVERT DIGIT FROM B-ARY TO B-ARY GRAY CODE
C 
       GTEMP(HISUM) = YTEMP(HISUM)
       DO 200 I = HISUM-1,0,-1
           TEM = YTEMP(I)-YTEMP(I+1)
           IF (TEM .LT. 0) THEN
              GTEMP(I) = TEM+QS
           ELSE   
             GTEMP(I) = MOD(TEM,QS)
           ENDIF   
200      CONTINUE
C
C FINDING THE POSITION OF COLUMN THAT IT'S VALUE 
C HAS BEEN CHANGED. 
C
       DO 300 I = 0, HISUM
           IF (GTEMP(I) .NE. PGTEMP(I)) THEN
              POS = I
              DIFF= GTEMP(I)-PGTEMP(I)
           ENDIF
           PGTEMP(I) = GTEMP(I)            
300       CONTINUE           
C
C UPDATE THE NEW VECTOR BY ADDING THE VALUE OF THE COLUMN THAT
C HAS BEEN CHANGED.
C       
      DO 400 K=1,S
          R=0.0
          DO 500 I=MAX,0,-1
              ZTEMP(K,I)=ZTEMP(K,I)+SCOEF(K,I,POS)*DIFF
              R=MOD(ZTEMP(K,I),QS)+RQS*R
500         CONTINUE             
          QUASI(K) = R*RQS
400   CONTINUE

C
C       UPDATE NEXTN AND, IF NEEDED, TESTN AND
C       HISUM
C
      NEXTN=NEXTN+1
      IF(NEXTN.EQ.TESTN) THEN
        TESTN=TESTN*QS
        HISUM=HISUM+1
C
C       SINCE FLAG(2) IS TRUE,
C       HISUM STAYS UNDER 20
C
      ENDIF
C
      RETURN
      END
      
      
      DOUBLE PRECISION FUNCTION UNI()
*
*     Random number generator, adapted from F. James
*     "A Review of Random Number Generators"
*      Comp. Phys. Comm. 60(1990), pp. 329-344.
*
      DOUBLE PRECISION SEEDS(24), TWOM24, CARRY
      PARAMETER ( TWOM24 = 1D0/16777216 )
      INTEGER I, J
      SAVE I, J, CARRY, SEEDS
      DATA I, J, CARRY / 24, 10, 0.0 /
      DATA SEEDS / 
     & 0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635,
     & 0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043,
     & 0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127,
     & 0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/
      
      
      UNI = SEEDS(I) - SEEDS(J) - CARRY
      IF ( UNI .LT. 0 ) THEN 
         UNI = UNI + 1
         CARRY = TWOM24
      ELSE 
         CARRY = 0
      ENDIF
      SEEDS(I) = UNI
      I = 24 - MOD( 25-I, 24 )
      J = 24 - MOD( 25-J, 24 )
      RETURN
      END
