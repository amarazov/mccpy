      BLOCK DATA BDSOBL
C
C     INITIALIZES LABELLED COMMON /SOBDAT/
C     FOR "INSOBL".
C
C     THE ARRAY POLY GIVES SUCCESSIVE PRIMITIVE
C     POLYNOMIALS CODED IN BINARY, E.G.
C          45 = 100101
C     HAS BITS 5, 2, AND 0 SET (COUNTING FROM THE
C     RIGHT) AND THEREFORE REPRESENTS
C          X**5 + X**2 + X**0
C
C     THESE  POLYNOMIALS ARE IN THE ORDER USED BY
C     SOBOL IN USSR COMPUT. MATHS. MATH. PHYS. 16 (1977),
C     236-242. A MORE COMPLETE TABLE IS GIVEN IN SOBOL AND
C     LEVITAN, THE PRODUCTION OF POINTS UNIFORMLY
C     DISTRIBUTED IN A MULTIDIMENSIONAL CUBE (IN RUSSIAN),
C     PREPRINT IPM AKAD. NAUK SSSR, NO. 40, MOSCOW 1976.
C
C         THE INITIALIZATION OF THE ARRAY VINIT IS FROM THE
C     LATTER PAPER. FOR A POLYNOMIAL OF DEGREE M, M INITIAL
C     VALUES ARE NEEDED :  THESE ARE THE VALUES GIVEN HERE.
C     SUBSEQUENT VALUES ARE CALCULATED IN "INSOBL".
C
      INTEGER POLY(2:40),VINIT(2:40,8)
      COMMON  /SOBDAT/ POLY,VINIT
      SAVE    /SOBDAT/
C
      DATA POLY    /3,7,11,13,19,25,37,59,47,
     +              61,55,41,67,97,91,109,103,115,131,
     +              193,137,145,143,241,157,185,167,229,171,
     +              213,191,253,203,211,239,247,285,369,299/
C
      DATA (VINIT(I,1),I=2,40)  /39*1/
      DATA (VINIT(I,2),I=3,40)  /1,3,1,3,1,3,3,1,
     +                           3,1,3,1,3,1,1,3,1,3,
     +                           1,3,1,3,3,1,3,1,3,1,
     +                           3,1,1,3,1,3,1,3,1,3/
      DATA (VINIT(I,3),I=4,40)  /7,5,1,3,3,7,5,
     +                           5,7,7,1,3,3,7,5,1,1,
     +                           5,3,3,1,7,5,1,3,3,7,
     +                           5,1,1,5,7,7,5,1,3,3/
      DATA (VINIT(I,4),I=6,40)  /1,7,9,13,11,
     +                           1,3,7,9,5,13,13,11,3,15,
     +                           5,3,15,7,9,13,9,1,11,7,
     +                           5,15,1,15,11,5,3,1,7,9/
      DATA (VINIT(I,5),I=8,40)  /9,3,27,
     +                           15,29,21,23,19,11,25,7,13,17,
     +                           1,25,29,3,31,11,5,23,27,19,
     +                           21,5,1,17,13,7,15,9,31,9/
      DATA (VINIT(I,6),I=14,40) /37,33,7,5,11,39,63,
     +                           27,17,15,23,29,3,21,13,31,25,
     +                           9,49,33,19,29,11,19,27,15,25/
      DATA (VINIT(I,7),I=20,40) /13,
     +                           33,115,41,79,17,29,119,75,73,105,
     +                           7,59,65,21,3,113,61,89,45,107/
      DATA (VINIT(I,8),I=38,40) /7,23,39/
C
      END
      
      SUBROUTINE INSSOBL (FLAG, DIMEN, ATMOST, TAUS, QUASI,
     *                    MAX,IFLAG)
C
C     THIS IS MODIFIED ROUTINE OF "INSOBL".
C     FIRST CHECK WHETHER THE USER-SUPPLIED
C     DIMENSION "DIMEN" OF THE QUASI-RANDOM
C     VECTORS IS STRICTLY BETWEEN 1 AND 41.
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
C     IN "GOSOBL" WE CHECK THAT THIS IS NOT EXCEEDED.
C
C     THE LEADING ELEMENTS OF EACH ROW OF V ARE
C     INITIALIZED USING "VINIT" FROM "BDSOBL".
C     EACH ROW CORRESPONDS TO A PRIMITIVE POLYNOMIAL
C     (AGAIN, SEE "BDSOBL").  IF THE POLYNOMIAL HAS
C     DEGREE M, ELEMENTS AFTER THE FIRST M ARE CALCULATED.
C
C     THE NUMBERS IN V ARE ACTUALLY BINARY FRACTIONS.
C     LSM ARE LOWER TRIAUGULAR SCRAMBLING MATRICES.
C     USM ARE UPPER TRIAUGULAR SCRMABLING MATRIX.
C     SV ARE SCAMBLING GENERATING MATRICES AND THE NUMBERS
C     ARE BINARY FRACTIONS.
C     "RECIPD" HOLDS 1/(THE COMMON DENOMINATOR OF ALL
C     OF THEM).
C
C
C     "INSSOBL" IMPLICITLY COMPUTES THE FIRST SHIFTED 
C     VECTOR "LASTQ", AND RETURN IT TO THE CALLING
C     PROGRAM. SUBSEQUENT VECTORS COME FROM "GOSSOBL".
C     "LASTQ" HOLDS NUMERATORS OF THE LAST VECTOR GENERATED.
C
C     "TAUS" IS FOR DETERMINING "FAVORABLE" VALUES. AS
C     DISCUSSED IN BRATLEY/FOX, THESE HAVE THE FORM
C     N = 2**K WHERE K .GE. (TAUS+S-1) FOR INTEGRATION
C     AND K .GT. TAUS FOR GLOBAL OPTIMIZATION.
C
C     INPUTS :
C       FROM USER'S PROGRAM : DIMEN, ATMOST, MAX, IFLAG
C       FROM BLOCK DATA "BDSOBL" : POLY, VINIT
C
C     OUTPUTS :
C       TO USER'S PROGRAM : FLAG, TAUS, LASTQ,QUASI
C       TO "GOSOBL" VIA /SOBOL/ :
C         SV, S, MAXCOL, COUNT, LASTQ, RECIPD
C
C
      INTEGER  POLY(2:40),VINIT(2:40,8),MAX
      INTEGER  V(40,31),S,MAXCOL,COUNT,LASTQ(40)
      INTEGER  DIMEN,ATMOST,I,J,K,P,M,NEWV,EXOR,L
      INTEGER  TAU(13),TAUS,SHIFT(40),LSM(40,31),SV(40,31)
      INTEGER  OUTPUT(22),DIGIT(60),TEMP1,TEMP2,TEMP3,TEMP4
      INTEGER  USM(31,31),USHIFT(31),IFLAG,TV(40,31,31),MAXX
      REAL*8   RECIPD,QUASI(40),LL
      LOGICAL  FLAG(2),INCLUD(8)
      EXTERNAL EXOR
      COMMON   /SOBDAT/ POLY,VINIT
      COMMON   /SOBOL/  S,MAXCOL,SV,COUNT,LASTQ,RECIPD
      SAVE     /SOBDAT/,/SOBOL/
      DATA TAU /0,0,1,3,5,8,11,15,19,23,27,31,35/
C
C     CHECK PARAMETERS
C
      S = DIMEN
      FLAG(1) = (S .GE. 1 .AND. S .LE. 40)
      FLAG(2) = (ATMOST .GT.0 .AND. ATMOST .LT. 2**30)
      IF (.NOT. (FLAG(1) .AND. FLAG(2))) RETURN
      IF (S .LE. 13) THEN
        TAUS = TAU(S)
      ELSE
        TAUS = -1
C     RETURN A DUMMY VALUE TO THE CALLING PROGRAM
      ENDIF
*
C
C     FIND NUMBER OF BITS IN ATMOST
C
      I = ATMOST
      MAXCOL = 0
   10 MAXCOL = MAXCOL + 1
      I = I / 2
      IF (I .GT. 0) GOTO 10
      
C
C     INITIALIZE ROW 1 OF V
C
      DO 20 I = 1, MAXCOL
   20 V(1,I) = 1
C
C     INITIALIZE REMAINING ROWS OF V
C
      DO 100 I = 2, S
C
C     THE BIT PATTERN OF POLYNOMIAL I GIVES ITS FORM
C     (SEE COMMENTS TO "BDSOBL")
C     FIND DEGREE OF POLYNOMIAL I FROM BINARY ENCODING
C
        J = POLY(I)
        M = 0
   30   J = J / 2
        IF (J .GT. 0) THEN
          M = M + 1
          GOTO 30
        ENDIF
C
C     WE EXPAND THIS BIT PATTERN TO SEPARATE COMPONENTS
C     OF THE LOGICAL ARRAY INCLUD.
C
        J = POLY(I)
        DO 40 K = M, 1, -1
          INCLUD(K) = (MOD(J,2) .EQ. 1)
          J = J / 2
   40   CONTINUE
C
C     THE LEADING ELEMENTS OF ROW I COME FROM VINIT
C
        DO 50 J = 1, M
          V(I,J) = VINIT(I, J)
   50   CONTINUE
C
C     CALCULATE REMAINING ELEMENTS OF ROW I AS EXPLAINED
C     IN BRATLEY AND FOX, SECTION 2
C
        DO 70 J = M+1, MAXCOL
          NEWV = V(I, J-M)
          L = 1
          DO 60 K = 1, M
            L = 2 * L
            IF (INCLUD(K)) NEWV = EXOR(NEWV, L * V(I, J-K))
C
C     IF A FULL-WORD EXCLUSIVE-OR, SAY .XOR., IS AVAILABLE,
C     THEN REPLACE THE PRECEDING STATEMENT BY
C
C         IF (INCLUD(K)) NEWV = NEWV .XOR. (L * V(I, J-K))
C
C     TO GET A FASTER, EXTENDED FORTRAN PROGRAM
C
   60     CONTINUE
          V(I,J) = NEWV
          
   70   CONTINUE
C
      
  100 CONTINUE
C
C     MULTIPLY COLUMNS OF V BY APPROPRIATE POWER OF 2
C
      L = 1
      DO 120 J = MAXCOL-1, 1, -1
        L = 2 * L
       DO 110 I = 1, S
           V(I,J) = V(I,J) * L
  110   CONTINUE
  120 CONTINUE
C 
C COMPUTING GENERATOR MATRICES OF USER CHOICE
C

      IF (IFLAG .EQ. 0) THEN
         DO 130 I = 1,S
            DO 140 J= 1,MAXCOL 
                  SV(I,J) = V(I,J)
 140          CONTINUE
              SHIFT(I) = 0
 130      CONTINUE
          LL= 2.0**(MAXCOL)
      ELSE             
        IF ((IFLAG .EQ. 1) .OR. (IFLAG .EQ. 3)) THEN
         CALL GENSCRML(MAX,LSM,SHIFT)
         DO  230 I = 1,S
           DO 240 J = 1,MAXCOL
             L = 1
             TEMP2 = 0
             DO 250 P = MAX,1,-1
                TEMP1 = 0
                DO 260 K = 1,MAXCOL
                   TEMP1 = TEMP1+ 
     *                   (IBITS(LSM(I,P),K-1,1)*IBITS(V(I,J),K-1,1))                                         
 260            CONTINUE
                TEMP1 = MOD(TEMP1,2)
                TEMP2 = TEMP2+TEMP1*L   
                L = 2 * L
 250          CONTINUE
              SV(I,J) = TEMP2
 240      CONTINUE
 230     CONTINUE
            LL= 2.0**(MAX)
       ENDIF
       IF ((IFLAG .EQ. 2) .OR. (IFLAG .EQ. 3)) THEN
        CALL GENSCRMU(USM,USHIFT) 
         IF (IFLAG .EQ. 2) THEN
           MAXX = MAXCOL
         ELSE
           MAXX = MAX
         ENDIF    
        DO 330 I = 1,S
         DO 340 J = 1,MAXCOL
          P = MAXX
          DO 350 K = 1,MAXX
             IF (IFLAG .EQ. 2) THEN
              TV(I,P,J) = IBITS(V(I,J),K-1,1)
             ELSE
               TV(I,P,J) = IBITS(SV(I,J),K-1,1) 
             ENDIF 
              P = P-1
 350        CONTINUE
 340       CONTINUE
        
          DO 341 PP = 1,MAXCOL 
            TEMP2 = 0 
            TEMP4 = 0
            L = 1
            DO 360 J = MAXX,1,-1
                TEMP1 = 0
                TEMP3 = 0
             DO 370 P = 1,MAXCOL
                TEMP1 = TEMP1 + TV(I,J,P)*USM(P,PP)
                IF (PP .EQ. 1) THEN
                  TEMP3 = TEMP3 + TV(I,J,P)*USHIFT(P)
                ENDIF 
370           CONTINUE
              TEMP1 = MOD(TEMP1,2)
              TEMP2 = TEMP2 + TEMP1*L
              IF (PP .EQ. 1) THEN 
                TEMP3  = MOD(TEMP3,2)
                TEMP4 = TEMP4 + TEMP3*L
              ENDIF  
              L = 2*L
 360         CONTINUE
              SV(I,PP) = TEMP2
              IF (PP .EQ. 1) THEN
                 IF (IFLAG .EQ. 3) THEN
                   SHIFT(I) = EXOR(TEMP4, SHIFT(I))           
                 ELSE
                   SHIFT(I) = TEMP4
                 ENDIF  
              ENDIF
 341        CONTINUE 
 330       CONTINUE
             LL = 2.0**(MAXX)
       ENDIF
      ENDIF 
C
C     RECIPD IS 1/(COMMON DENOMINATOR OF THE ELEMENTS IN V)
C
      RECIPD = 1.0 /LL

C
C     SET UP FIRST VECTOR AND VALUES FOR "GOSOBL"
C
      COUNT = 0
      DO 400 I = 1, S
        LASTQ(I) = SHIFT(I)
        QUASI(I) = LASTQ(I)*RECIPD
 400   CONTINUE
      RETURN
      END

      SUBROUTINE GENSCRML(MAX,LSM,SHIFT)

C     GENERATING LOWER TRIANULAR SCRMABLING MATRICES AND SHIFT VECTORS.
C     INPUTS :
C       FROM INSSOBL : MAX
C       FROM BLOCK DATA "SOBOL" : S, MAXCOL,
C
C     OUTPUTS :
C       TO INSSOBL : LSM, SHIFT

      INTEGER LSM(40,31),MAXCOL,P,I,J
      INTEGER SHIFT(40),MAX,S,TEMP,STEMP,L,LL
      REAL*8 UNI
      COMMON /SOBOL/ S,MAXCOL
      SAVE /SOBOL/
      
      DO 10 P = 1,S
               SHIFT(P) = 0
               L = 1
         DO 20 I = MAX,1,-1
               LSM(P,I) = 0
               STEMP =  MOD((int(UNI()*1000.0)),2)
               SHIFT(P) = SHIFT(P)+STEMP*L
               L = 2 * L
               LL = 1
            DO 30 J = MAXCOL,1,-1
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
C       FROM BLOCK DATA "SOBOL" : S, MAXCOL,
C
C     OUTPUTS :
C       TO INSSOBL : USM, USHIFT

      INTEGER USM(31,31),MAXCOL,P,I,J
      INTEGER USHIFT(31),S,TEMP,STEMP,L,LL
      REAL*8 UNI
      COMMON /SOBOL/ S,MAXCOL
      SAVE /SOBOL/
      
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

      SUBROUTINE GOSSOBL (QUASI)
C
C     THIS SUBROUTINE GENERATES A NEW
C     QUASIRANDOM VECTOR WITH EACH CALL
C
C     IT ADAPTS THE IDEAS OF ANTONOV AND SALEEV,
C     USSR COMPUT. MATHS. MATH. PHYS. 19 (1980),
C     252 - 256
C
C     THE USER MUST CALL "INSOBL" BEFORE CALLING
C     "GOSOBL".  AFTER CALLING "INSOBL", TEST
C     FLAG(1) AND FLAG(2);  IF EITHER IS FALSE,
C     DO NOT CALL "GOSOBL".  "GOSOBL" CHECKS
C     THAT THE USER DOES NOT MAKE MORE CALLS
C     THAN HE SAID HE WOULD : SEE THE COMMENTS
C     TO "INSOBL".
C
C     INPUTS:
C       FROM USER'S CALLING PROGRAM:
C         NONE
C
C       FROM LABELLED COMMON /SOBOL/:
C         SV        TABLE OF DIRECTION NUMBERS
C         S        DIMENSION
C         MAXCOL   LAST COLUMN OF V TO BE USED
C         COUNT    SEQUENCE NUMBER OF THIS CALL
C         LASTQ    NUMERATORS FOR LAST VECTOR GENERATED
C         RECIPD   (1/DENOMINATOR) FOR THESE NUMERATORS
C
      INTEGER  SV(40,31),S,MAXCOL,COUNT,LASTQ(40)
      INTEGER  I,L,EXOR
      REAL*8     RECIPD,QUASI(40)
      COMMON   /SOBOL/ S,MAXCOL,SV,COUNT,LASTQ,RECIPD
      SAVE     /SOBOL/
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

      SUBROUTINE NSOBOL(N, QUASI)

      INTEGER  I,J,N
      REAL*8    QUASI(40,N)
      REAL*8    Q(40)
      COMMON   /SOBOL/ S
      SAVE     /SOBOL/

      DO 40 I = 1,N
        CALL GOSSOBL(Q)
        DO 50 J = 1,S
            QUASI(J,I) = Q(J)
  50    CONTINUE
  40  CONTINUE

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
      
