      SUBROUTINE INSNIEDX (FLAG, DIMEN, ATMOST, QUASI,
     *                     MAX,IFLAG)
C
C     THIS INITIALIZE NIED-XING SEQUENCE.
C     FIRST CHECK WHETHER THE USER-SUPPLIED
C     DIMENSION "DIMEN" OF THE QUASI-RANDOM
C     VECTORS IS STRICTLY BETWEEN 1 AND 16.
C     IF SO, FLAG(1) = .TRUE.
C
C     NEXT CHECK "ATMOST", AN UPPER BOUND ON THE NUMBER
C     OF CALLS THE USER INTENDS TO MAKE ON "GOSOBL".  IF
C     THIS IS POSITIVE AND LESS THAN 2**30, THEN FLAG(2) = .TRUE.
C     (WE ASSUME WE ARE WORKING ON A COMPUTER WITH
C     WORD LENGTH AT LEAST 31 BITS EXCLUDING SIGN.)
C     THE NUMBER OF COLUMNS OF THE ARRAY V WHICH
C     ARE INITIALIZED IS
C          MAXCOL = NUMBER OF BITS IN ATMOST.
C     IN "GOSNIEDX" WE CHECK THAT THIS IS NOT EXCEEDED.
C  
C     NX STORED NIEDERREiITER-XING SEQUENCE GENERATORS
C     WHICH IS BINARY FRACTIONS WITH RESPECT TO ROWS.
C     THE GENERATORS CAN BE FOUND URL OF WEBSITE
C     http://www.dismat.oeaw.ac.at/pirs/niedxing.html
C
C
C     THE NUMBERS IN V ARE ACTUALLY BINARY FRACTIONS WITH
C     RESPECT TO COLUMNS.
C
C     LSM ARE LOWER TRIANGULAR SCRAMBLING MATRICES.
C     USM ARE UPPER TRIANGULAR SCRAMBLING MATRICES.
C     SV ARE SCRAMBLING GENERATOR MATRICES AND THE NUMBERS
C     ARE BINARY FRACTIONS.
C     "RECIPD" HOLDS 1/(THE COMMON DENOMINATOR OF ALL
C     OF THEM).
C
C
C     "INSNIEDX" IMPLICITLY COMPUTES THE FIRST SHIFTED 
C     VECTOR "LASTQ", AND RETURN IT TO THE CALLING
C     PROGRAM. SUBSEQUENT VECTORS COME FROM "GOSNIEDX".
C     "LASTQ" HOLDS NUMERATORS OF THE LAST VECTOR GENERATED.
C
C
C     INPUTS :
C       FROM USER'S PROGRAM : DIMEN, ATMOST, MAX, 
C                             IFLAG, QUASI
C
C     OUTPUTS :
C       TO USER'S PROGRAM : FLAG, LASTQ
C       TO "GOSNIEDX" VIA /NIEDX/ :
C         SV, S, MAXCOL, COUNT, LASTQ, RECIPD
C
C
      INTEGER  NX(16,30),MAX
      INTEGER  V(40,30),S,MAXCOL,COUNT,LASTQ(40)
      INTEGER  DIMEN,ATMOST,I,J,K,P,L,M,NEWV,EXOR
      INTEGER  SHIFT(40),LSM(40,31),SV(40,31)
      INTEGER  OUTPUT(22),DIGIT(60),TEMP1,TEMP2,TEMP3,TEMP4
      INTEGER  USM(30,30),USHIFT(30),IFLAG,TV(40,31,31)
      REAL*8   RECIPD,QUASI(40),LL
      LOGICAL  FLAG(2),INCLUD(8)
      CHARACTER*20 name
      EXTERNAL EXOR
      COMMON   /NIEDX/  S,MAXCOL,SV,COUNT,LASTQ,RECIPD
      SAVE     /NIEDX/
 
C
C     CHECK PARAMETERS
C
      S = DIMEN
      FLAG(1) = (S .GE. 1 .AND. S .LE. 16)
      FLAG(2) = (ATMOST .GT. 0 .AND. ATMOST .LT. 2**30)
      IF (.NOT. (FLAG(1) .AND. FLAG(2))) RETURN
      I = ATMOST
      MAXCOL = 0
   10 MAXCOL = MAXCOL + 1
      I = I / 2
      IF (I .GT. 0) GOTO 10
C
C   GET GENERATOR MATRICES
C
      IF (S .LE. 4) THEN
         name = 'nxs4m30'  
      ELSE IF (S. EQ. 5) THEN
         name = 'nxs5m30' 
      ELSE IF (S .EQ. 6) THEN
         name = 'nxs6m30' 
      ELSE IF (S .LE. 7) THEN
         name = 'nxs7m30'  
      ELSE IF (S. EQ. 8) THEN
         name = 'nxs8m30' 
      ELSE IF (S .EQ. 9) THEN
         name = 'nxs9m30'          
      ELSE IF (S .EQ. 10) THEN
         name = 'nxs10m30' 
      ELSE IF (S .LE. 11) THEN
         name = 'nxs11m30'  
      ELSE IF (S. EQ. 12) THEN
         name = 'nxs12m30' 
      ELSE IF (S .EQ. 13) THEN
         name = 'nxs13m30'                   
      ELSE IF (S .LE. 14) THEN
         name = 'nxs14m30'  
      ELSE IF (S. EQ. 15) THEN
         name = 'nxs15m30' 
      ELSE IF (S .EQ. 16) THEN
         name = 'nxs16m30'
      ENDIF        
            
      OPEN(unit=1,file=name,status='UNKNOWN')          
      DO 20 I = 1, S
           READ(1,*) (NX(I,J), J = 1,30)
 21       CONTINUE
 20     CONTINUE      
      CLOSE(unit=1,status='KEEP')      
    
      DO 111 I = 1,S
        DO 121 J = 1, MAXCOL
            V(I,J) = 0
 121      CONTINUE
 111    CONTINUE 
        
      L = 1
      DO 120 J = 30,1,-1
       DO 110 I = 1, S
         DO 105 K =1,MAXCOL
            V(I,K) = V(I,K)+IBITS(NX(I,J),30-K,1)* L
  105      CONTINUE                   
  110   CONTINUE
        L = 2*L
  120 CONTINUE

C
C   COMPUTING  GENERATOR MATRICES OF USER CHOICE
C      
     
      IF (IFLAG .EQ. 0) THEN
         DO 30 I = 1,S
            DO 40 J= 1,MAXCOL 
              SV(I,J) = V(I,J)
 40          CONTINUE
             SHIFT(I) = 0
 30        CONTINUE   
          LL=2.0**(30)          
      ELSE
       IF ((IFLAG .EQ. 1) .OR. (IFLAG .EQ. 3)) THEN
         CALL GENSCRML(MAX,LSM,SHIFT)
         DO 130 I = 1,S
           DO 140 J = 1,MAXCOL
              L = 1
              TEMP2 = 0
              DO 150 P = MAX,1,-1
                   TEMP1 = 0
                   DO 160 K = 1,30
                       TEMP1 = TEMP1+(IBITS(LSM(I,P),K-1,1) *
     *                                             IBITS(V(I,J),K-1,1))
 160               CONTINUE
                   TEMP1 = MOD(TEMP1,2)
                   TEMP2 = TEMP2+TEMP1*L   
                   L = 2 * L
 150          CONTINUE
             SV(I,J) = TEMP2
 140      CONTINUE
 130     CONTINUE
            LL = 2.0 **(MAX)
       ENDIF
      
       IF ((IFLAG .EQ. 2) .OR. (IFLAG .EQ. 3)) THEN
         CALL GENSCRMU(USM,USHIFT) 
         IF (IFLAG .EQ. 2) THEN
             MAX = 30
         ENDIF    
         DO 230 I = 1,S
           DO 240 J = 1,MAXCOL
            P = MAX
            DO 250 K = 1,MAX
               IF (IFLAG .EQ. 2) THEN
                 TV(I,P,J) = IBITS(V(I,J),K-1,1)
               ELSE
                 TV(I,P,J) = IBITS(SV(I,J),K-1,1) 
               ENDIF 
               P = P-1
 250        CONTINUE
 240       CONTINUE
    
          DO 241 PP = 1,MAXCOL 
            TEMP2 = 0 
            TEMP4 = 0
            L = 1
            DO 260 J = MAX,1,-1
                 TEMP1 = 0
                 TEMP3 = 0
                 DO 270 P = 1,MAXCOL
                     TEMP1 = TEMP1 + TV(I,J,P)*USM(P,PP)
                     IF (PP .EQ. 1) THEN
                       TEMP3 = TEMP3 + TV(I,J,P)*USHIFT(P)
                     ENDIF 
270            CONTINUE
                TEMP1 = MOD(TEMP1,2)
                TEMP2 = TEMP2 + TEMP1*L
                IF (PP .EQ. 1) THEN 
                  TEMP3  = MOD(TEMP3,2)
                  TEMP4 = TEMP4 + TEMP3*L
                ENDIF  
                L = 2*L
 260        CONTINUE
              SV(I,PP) = TEMP2
              IF (PP .EQ. 1) THEN
                 IF (IFLAG .EQ. 3) THEN
                   SHIFT(I) = EXOR(TEMP4, SHIFT(I))           
                 ELSE
                   SHIFT(I) = TEMP4
                 ENDIF  
              ENDIF
 241        CONTINUE 
 230       CONTINUE
            LL= 2.0**(MAX)
       ENDIF
      ENDIF
C
C     RECIPD IS 1/(COMMON DENOMINATOR OF THE ELEMENTS IN V)
C
      RECIPD = 1.0 /LL

C
C     SET UP FIRST VECTOR AND VALUES FOR "GOSNIEDX"
C     AND RETURN THE FIRST SEQUENCE
C
      COUNT = 0
      DO 170 I = 1, S
        LASTQ(I) = SHIFT(I)
        QUASI(I) = LASTQ(I)*RECIPD
 170   CONTINUE
      RETURN
      END

      SUBROUTINE GENSCRML(MAX,LSM,SHIFT)

C     GENERATING LOWER TRIAGULAR SCRMABLING MATRICES 
C     AND SHIFT VECTORS.
C     INPUTS :
C       FROM INSNIEDX : MAX
C       FROM BLOCK DATA "SOBOL" : S, MAXCOL,
C
C     OUTPUTS :
C       TO INSSOBL : LSM, SHIFT

      INTEGER LSM(40,31),MAXCOL,P,I,J
      INTEGER SHIFT(40),MAX,S,TEMP,STEMP,L,LL
      REAL*8 UNI
      COMMON /NIEDX/ S,MAXCOL
      SAVE /NIEDX/
      
      DO 10 P = 1,S
              SHIFT(P) = 0
               L = 1
         DO 20 I = MAX,1,-1
               LSM(P,I) = 0
               STEMP =  MOD((int(UNI()*1000.0)),2)
               SHIFT(P) = SHIFT(P)+STEMP*L
               L = 2 * L
               LL = 1
            DO 30 J = 30,1,-1
               IF (J .EQ. I) THEN
                TEMP = 1
               ELSE IF (J .LT. I)  THEN 
                TEMP = MOD((int(UNI()*1000.0)),2)
               ELSE
                TEMP = 0
               ENDIF
               LSM(P,I) = LSM(P,I)+TEMP*LL
               LL = 2 * LL            
 30           CONTINUE
 20        CONTINUE 
 10      CONTINUE
 
      RETURN
      END   
      
      SUBROUTINE GENSCRMU(USM,USHIFT)

C     GENERATING UPPER TRIANGULAR SCRMABLING MATRICES AND 
C     SHIFT VECTORS.
C     INPUTS :
C       FROM BLOCK DATA "NIEDX" : S, MAXCOL,
C
C     OUTPUTS :
C       TO INSSOBL : USM, USHIFT

      INTEGER USM(30,30),MAXCOL,P,I,J
      INTEGER USHIFT(30),MAX,S,TEMP,STEMP,L,LL
      REAL*8 UNI
      COMMON /NIEDX/ S,MAXCOL
      SAVE /NIEDX/
      
          DO 20 I = 1,MAXCOL
               STEMP =  MOD((int(UNI()*1000.0)),2)
               USHIFT(I) = STEMP               
            DO 30 J = 1,MAXCOL 
               IF (J .EQ. I) THEN
                 TEMP = 1
               ELSE IF (J .GT. I)  THEN 
                TEMP = MOD((int(UNI()*1000.0)),2)
               ELSE
                TEMP = 0
               ENDIF
               USM(I,J) = TEMP        
 30           CONTINUE
 20        CONTINUE 
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

      SUBROUTINE GOSNIEDX (QUASI)
C
C     THIS SUBROUTINE GENERATES A NEW
C     QUASIRANDOM VECTOR WITH EACH CALL
C
C     IT ADAPTS THE IDEAS OF ANTONOV AND SALEEV,
C     USSR COMPUT. MATHS. MATH. PHYS. 19 (1980),
C     252 - 256 
C
C     THE USER MUST CALL "INSNIEDX" BEFORE CALLING
C     "GOSNIEDX".  AFTER CALLING "INSNIEDX", TEST
C     FLAG(1) AND FLAG(2);  IF EITHER IS FALSE,
C     DO NOT CALL "GOSNIEDX".  "GOSNIEDX" CHECKS
C     THAT THE USER DOES NOT MAKE MORE CALLS
C     THAN HE SAID HE WOULD : SEE THE COMMENTS
C     TO "INSNIEDX".
C
C     INPUTS:
C       FROM USER'S CALLING PROGRAM:
C         NONE
C
C       FROM LABELLED COMMON /NIEDX/:
C         SV       TABLE OF DIRECTION NUMBERS
C         S        DIMENSION
C         MAXCOL   LAST COLUMN OF V TO BE USED
C         COUNT    SEQUENCE NUMBER OF THIS CALL
C         LASTQ    NUMERATORS FOR LAST VECTOR GENERATED
C         RECIPD   (1/DENOMINATOR) FOR THESE NUMERATORS
C
      INTEGER  SV(40,31),S,MAXCOL,COUNT,LASTQ(40)
      INTEGER  I,L,EXOR
      REAL*8     RECIPD,QUASI(40)
      COMMON   /NIEDX/ S,MAXCOL,SV,COUNT,LASTQ,RECIPD
      SAVE     /NIEDX/
C
C     FIND THE POSITION OF THE RIGHT-HAND ZERO IN COUNT
C
      L = 0
      I = COUNT
    1 L = L + 1
      IF (MOD(I,2) .EQ. 1) THEN
        I = I / 2
        GOTO 1
      ENDIF
C
C     CHECK THAT THE USER IS NOT CHEATING !
C

      IF (L .GT. MAXCOL) STOP ' TOO MANY CALLS ON GOSOBL'
C
C     CALCULATE THE NEW COMPONENTS OF QUASI,
C     FIRST THE NUMERATORS, THEN NORMALIZED
C
      DO 2 I = 1, S
        LASTQ(I) = EXOR (LASTQ(I), SV(I,L))
C
C     IF A FULL-WORD EXCLUSIVE-OR, SAY .XOR., IS AVAILABLE
C     THEN REPLACE THE PRECEDING STATEMENT BY
C
C         LASTQ(I) = LASTQ(I) .XOR. SV(I,L)
C
C     TO GET A FASTER, EXTENDED FORTRAN PROGRAM
C
        QUASI(I) = LASTQ(I) * RECIPD
    2 CONTINUE
C
      COUNT = COUNT + 1
C
      RETURN
      END
      
      INTEGER FUNCTION EXOR (IIN, JIN)
      INTEGER I,J,K,L,IIN,JIN
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

      PROGRAM GENSNIEDX
C
C       IT USES A NONSTANDARD TIMING
C       ROUTINE "DTIME"
C
C       THIS PROGRAM TESTS ACCURACY OF
C       NUMERICAL INTEGRATION USING "GOSNIEDX"
C       AND INTEGRAND (2) OF DAVIS AND
C       RABINOWITZ, PAGE 406
C
C     
C      User Define: 
C        DIMEN : dimension 
C        ATMOST : sequence length
C        SAM : Number of replications
C        MAX : Maximum Digits of Scrambling
C        IFLAG: User Choice of Sequences
C        IFLAG = 0 : No Scrambling
C        IFLAG = 1 : Owen type Scrambling
C        IFLAG = 2 : Faure-Tezuka type Scrambling
C        IFLAG = 3 : Owen + Faure-Tezuka type Scrambling
C      
      LOGICAL FLAG(2)
      INTEGER DIMEN,ATMOST,II,I,J,SEQL,MAX,SAM
      INTEGER IFLAG
      REAL*8 QUASI(40),SUM,F,EI
      REAL*4 TIM,tarray(2) 
      
      SAM = 1
      MAX = 30
      DIMEN = 2
      IFLAG = 1
      ATMOST = 2**12
      DO 22 II = 1,SAM
         WRITE(*,*) 'I = ITERATION NUMBER'
         WRITE(*,*) 'EI = ESTIMATED INTEGRAL' 
         TIM = DTIME_(tarray)
         CALL INSNIEDX(FLAG,DIMEN,ATMOST,QUASI,MAX,IFLAG)
         F= 1.0
         DO 55 J = 1,DIMEN
             F = F*ABS(4.0*QUASI(J)-2.0) 
 55         CONTINUE 
          SUM = F 
          DO 100 I=2,ATMOST    
             CALL GOSNIEDX(QUASI)
             F = 1.0
             DO 50 J=1,DIMEN
                F = F*ABS(4.0*QUASI(J)-2.0)                 
 50           CONTINUE
             SUM = SUM + F
             IF (MOD(I,500) .EQ. 0) THEN
                  WRITE(*,*) 'I = ',I
                  WRITE(*,*) 'EI = ',SUM/I
             ENDIF                
 100       CONTINUE  
         TIM = DTIME_(tarray)   
         WRITE(*,*) 'Total time elapsed = ', tarray(1)      
 22    CONTINUE
      STOP
      END
