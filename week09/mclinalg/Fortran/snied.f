      PROGRAM SNIED
C      
C   This is Modified Program of Niederreiter Sequences Generator  
C   It generates Various Scrambled Niederreiter Sequences.
C    SNIED and its associated subroutines are listed below.
C
C       SNIED
C          SINLO2
C          SGOLO2
C          SCALCC2
C          CALCV
C          CHARAC
C          SETFLD
C          PLYMUL
C          TESTF
C
C    The suffix 2 indicates  routines for use only by
C    the set of programs tailored for base 2.
C
C       IT USES A NONSTANDARD TIMING
C       ROUTINE "DTIME"
C
C       THIS PROGRAM TESTS ACCURACY OF
C       NUMERICAL INTEGRATION USING "SGOLO2"
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
C
      INTEGER MAXDIM
      PARAMETER (MAXDIM=318)
C
      INTEGER II,I,J,DIMEN,ATMOST,IFLAG,SAM,SKIP
      REAL*8 QUASI(MAXDIM),SUM,F,EI
      REAL*4 TIM,tarray(2)   

      IFLAG = 1
      DIMEN = 2
      ATMOST = 2**12
      SAM = 1
      SKIP = 0
      DO 100 II = 1,SAM
        WRITE(*,*) 'I = ITERATION NUMBER'
        WRITE(*,*) 'EI = ESTIMATED INTEGRAL'         
        TIM = DTIME_(tarray)
        CALL SINLO2 (DIMEN, SKIP, IFLAG)
        SUM = 0.0
        DO 110 I = 1, ATMOST
           CALL SGOLO2(QUASI)
           F = 1.0
           DO 120 J=1,DIMEN
               F = F*ABS(4.0*QUASI(J)-2.0)              
  120        CONTINUE
           SUM = SUM+F
           IF (MOD(I,500) .EQ. 0) THEN
             WRITE(*,*) 'I = ',I
             WRITE(*,*) 'EI = ',SUM/I
           ENDIF     
  110     CONTINUE
          TIM = DTIME_(tarray)
          WRITE(*,*) 'Total time elapsed = ', tarray(1)
  100    CONTINUE 


      STOP
      END
C
       SUBROUTINE SCALCC2(IFLAG)
C
C   This is modified routine of SCALCC2.
C
C   This program calculates the values of the constants C(I,J,R).
C   As far as possible, we use Niederreiter's notation.
C   For each value of I, we first calculate all the corresponding
C   values of C :  these are held in the array CI.  All these
C   values are either 0 or 1.  Next we pack the values into the
C   array CJ, in such a way that CJ(I,R) holds the values of C
C   for the indicated values of I and R and for every value of
C   J from 1 to NBITS.  The most significant bit of CJ(I,R)
C   (not counting the sign bit) is C(I,1,R) and the least
C   significant bit is C(I,NBITS,R).
C     When all the values of CJ have been calculated, we return
C   this array to the calling program.
C
C  --------------------------------------------------------------
C
C   We define the common block /COMM2/ and some
C   associated parameters.  These are for use in the base 2
C   version of the generator.
C
       INTEGER MAXDIM, NBITS
       PARAMETER (MAXDIM=318, NBITS=31)
C
C   The parameter MAXDIM is the maximum dimension that will be used.
C   NBITS is the number of bits (not counting the sign) in a
C   fixed-point integer.
C
       INTEGER CJ(MAXDIM, 0:NBITS-1),SCJ(MAXDIM, 0:NBITS-1)
       INTEGER DIMEN, COUNT, NEXTQ(MAXDIM)
       COMMON /COMM2/ SCJ, DIMEN, COUNT, NEXTQ
       SAVE   /COMM2/
C
C   The common block /COMM2/ :
C     CJ    - Contains the packed values of Niederreiter's C(I,J,R)
C     SCJ    - Contains the packed values of User choice of scrambled
C                  Niederreiter's C(I,J,R)
C     DIMEN   - The dimension of the sequence to be generated
C     COUNT - The index of the current item in the sequence,
C             expressed as an array of bits.  COUNT(R) is the same
C             as Niederreiter's AR(N) (page 54) except that N is
C             implicit.
C     NEXTQ - The numerators of the next item in the series.  These
C             are like Niederreiter's XI(N) (page 54) except that
C             N is implicit, and the NEXTQ are integers.  To obtain
C             the values of XI(N), multiply by RECIP (see GOLO2).
C
C   Array SCJ of the common block is set up by subroutine SCALCC2.
C   The other items in the common block are set up by INLO2.
C
C   --------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
       INTEGER MAXQ
       PARAMETER (MAXQ=50)

       INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
       INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
       COMMON /FIELD/ P, Q, ADD, MUL, SUB
       SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
       INTEGER MAXDEG, DEG
       PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   ---------------------------------------------------------------
C
C
C
C   MAXE   - We need MAXDIM irreducible polynomials over Z2.
C            MAXE is the highest degree among these.
C   MAXV   - The maximum possible index used in V.
C
       INTEGER MAXE, MAXV
       PARAMETER (MAXE=11, MAXV=NBITS+MAXE)
C
C INPUT :
C   The array SCJ to be initialised, and DIMEN the number of
C   dimensions we are using, are transmitted through /COMM2/.
C
C OUTPUT :
C   The array SCJ is returned to the calling program.
C
C USES :
C   Subroutine SETFLD is used to set up field arithmetic tables.
C   (Although this is a little heavy-handed for the field of
C   order 2, it makes for uniformity with the general program.)
C   Subroutine CALCV is used to for the auxiliary calculation
C   of the values V needed to get the Cs.
C
       INTEGER PX(-1:MAXDEG), B(-1:MAXDEG)
       INTEGER V(0:MAXV), CI(NBITS, 0:NBITS-1)
       INTEGER E, I, J, R, U, TERM, EXOR
       INTEGER MAXS,TEMP1,TEMP2,TEMP3,L
       INTEGER LSM(MAXDIM,0:NBITS-1),SHIFT(MAXDIM)
       INTEGER PP,TV(MAXDIM,0:NBITS-1,0:NBITS-1),IFLAG
       INTEGER USHIFT(0:NBITS-1),USM(0:NBITS-1,0:NBITS-1)
C
       INTEGER IRRED(MAXDIM, -1:MAXE)
       SAVE IRRED
C
C   This DATA statement supplies the coefficients and the
C   degrees of the first 12 irreducible polynomials over Z2.
C   They are taken from Lidl and Niederreiter, FINITE FIELDS,
C   Cambridge University Press (1984), page 553.
C   The order of these polynomials is the same as the order in
C   file 'irrtabs.dat' used by the general program.
C
C   In this block PX(I, -1) is the degree of the Ith polynomial,
C   and PX(I, N) is the coefficient of x**n in the Ith polynomial.
C
        DATA (IRRED(1,I), I=-1,1)  / 1,0,1 /
        DATA (IRRED(2,I), I=-1,1)  / 1,1,1 /
        DATA (IRRED(3,I), I=-1,2)  / 2,1,1,1 /
        DATA (IRRED(4,I), I=-1,3)  / 3,1,1,0,1 /
        DATA (IRRED(5,I), I=-1,3)  / 3,1,0,1,1 /
        DATA (IRRED(6,I), I=-1,4)  / 4,1,1,0,0,1 /
        DATA (IRRED(7,I), I=-1,4)  / 4,1,0,0,1,1 /
        DATA (IRRED(8,I), I=-1,4)  / 4,1,1,1,1,1 /
        DATA (IRRED(9,I), I=-1,5)  / 5,1,0,1,0,0,1 /
        DATA (IRRED(10,I), I=-1,5) / 5,1,0,0,1,0,1 /
        DATA (IRRED(11,I), I=-1,5) / 5,1,1,1,1,0,1 /
        DATA (IRRED(12,I), I=-1,5) / 5,1,1,1,0,1,1 /
        DATA (IRRED(13,I), I=-1,5) / 5,1,1,0,1,1,1 /
        DATA (IRRED(14,I), I=-1,5) / 5,1,0,1,1,1,1 /
        DATA (IRRED(15,I), I=-1,6) / 6,1,1,0,0,0,0,1 /
        DATA (IRRED(16,I), I=-1,6) / 6,1,0,0,1,0,0,1 /
        DATA (IRRED(17,I), I=-1,6) / 6,1,1,1,0,1,0,1 /
        DATA (IRRED(18,I), I=-1,6) / 6,1,1,0,1,1,0,1 /
        DATA (IRRED(19,I), I=-1,6) / 6,1,0,0,0,0,1,1 /
        DATA (IRRED(20,I), I=-1,6) / 6,1,1,1,0,0,1,1 /
        DATA (IRRED(21,I), I=-1,6) / 6,1,0,1,1,0,1,1 /
        DATA (IRRED(22,I), I=-1,6) / 6,1,1,0,0,1,1,1 /
        DATA (IRRED(23,I), I=-1,6) / 6,1,0,1,0,1,1,1 /
        DATA (IRRED(24,I), I=-1,7) / 7,1,1,0,0,0,0,0,1 /
        DATA (IRRED(25,I), I=-1,7) / 7,1,0,0,1,0,0,0,1 /
        DATA (IRRED(26,I), I=-1,7) / 7,1,1,1,1,0,0,0,1 /
        DATA (IRRED(27,I), I=-1,7) / 7,1,0,0,0,1,0,0,1 /
        DATA (IRRED(28,I), I=-1,7) / 7,1,0,1,1,1,0,0,1 /
        DATA (IRRED(29,I), I=-1,7) / 7,1,1,1,0,0,1,0,1 /
        DATA (IRRED(30,I), I=-1,7) / 7,1,1,0,1,0,1,0,1 /
        DATA (IRRED(31,I), I=-1,7) / 7,1,0,0,1,1,1,0,1 /
        DATA (IRRED(32,I), I=-1,7) / 7,1,1,1,1,1,1,0,1 /
        DATA (IRRED(33,I), I=-1,7) / 7,1,0,0,0,0,0,1,1 /
        DATA (IRRED(34,I), I=-1,7) / 7,1,1,0,1,0,0,1,1 /
        DATA (IRRED(35,I), I=-1,7) / 7,1,1,0,0,1,0,1,1 /
        DATA (IRRED(36,I), I=-1,7) / 7,1,0,1,0,1,0,1,1 /
        DATA (IRRED(37,I), I=-1,7) / 7,1,0,1,0,0,1,1,1 /
        DATA (IRRED(38,I), I=-1,7) / 7,1,1,1,1,0,1,1,1 /
        DATA (IRRED(39,I), I=-1,7) / 7,1,0,0,0,1,1,1,1 /
        DATA (IRRED(40,I), I=-1,7) / 7,1,1,1,0,1,1,1,1 /
        DATA (IRRED(41,I), I=-1,7) / 7,1,0,1,1,1,1,1,1 /
        DATA (IRRED(42,I), I=-1,8) / 8,1,1,0,1,1,0,0,0,1 /
        DATA (IRRED(43,I), I=-1,8) / 8,1,0,1,1,1,0,0,0,1 /
        DATA (IRRED(44,I), I=-1,8) / 8,1,1,0,1,0,1,0,0,1 /
        DATA (IRRED(45,I), I=-1,8) / 8,1,0,1,1,0,1,0,0,1 /
        DATA (IRRED(46,I), I=-1,8) / 8,1,0,0,1,1,1,0,0,1 /
        DATA (IRRED(47,I), I=-1,8) / 8,1,1,1,1,1,1,0,0,1 /
        DATA (IRRED(48,I), I=-1,8) / 8,1,0,1,1,0,0,1,0,1 /
        DATA (IRRED(49,I), I=-1,8) / 8,1,1,1,1,1,0,1,0,1 /
        DATA (IRRED(50,I), I=-1,8) / 8,1,1,0,0,0,1,1,0,1 /
        DATA (IRRED(51,I), I=-1,8) / 8,1,0,1,0,0,1,1,0,1 /
        DATA (IRRED(52,I), I=-1,8) / 8,1,0,0,1,0,1,1,0,1 /
        DATA (IRRED(53,I), I=-1,8) / 8,1,0,0,0,1,1,1,0,1 /
        DATA (IRRED(54,I), I=-1,8) / 8,1,1,1,0,1,1,1,0,1 /
        DATA (IRRED(55,I), I=-1,8) / 8,1,1,0,1,1,1,1,0,1 /
        DATA (IRRED(56,I), I=-1,8) / 8,1,1,1,0,0,0,0,1,1 /
        DATA (IRRED(57,I), I=-1,8) / 8,1,1,0,1,0,0,0,1,1 /
        DATA (IRRED(58,I), I=-1,8) / 8,1,0,1,1,0,0,0,1,1 /
        DATA (IRRED(59,I), I=-1,8) / 8,1,1,1,1,1,0,0,1,1 /
        DATA (IRRED(60,I), I=-1,8) / 8,1,1,0,0,0,1,0,1,1 /
        DATA (IRRED(61,I), I=-1,8) / 8,1,0,0,1,0,1,0,1,1 /
        DATA (IRRED(62,I), I=-1,8) / 8,1,0,0,0,1,1,0,1,1 /
        DATA (IRRED(63,I), I=-1,8) / 8,1,0,1,1,1,1,0,1,1 /
        DATA (IRRED(64,I), I=-1,8) / 8,1,1,0,0,0,0,1,1,1 /
        DATA (IRRED(65,I), I=-1,8) / 8,1,1,1,1,0,0,1,1,1 /
        DATA (IRRED(66,I), I=-1,8) / 8,1,1,1,0,1,0,1,1,1 /
        DATA (IRRED(67,I), I=-1,8) / 8,1,0,1,1,1,0,1,1,1 /
        DATA (IRRED(68,I), I=-1,8) / 8,1,1,1,0,0,1,1,1,1 /
        DATA (IRRED(69,I), I=-1,8) / 8,1,1,0,0,1,1,1,1,1 /
        DATA (IRRED(70,I), I=-1,8) / 8,1,0,1,0,1,1,1,1,1 /
        DATA (IRRED(71,I), I=-1,8) / 8,1,0,0,1,1,1,1,1,1 /
        DATA (IRRED(72,I), I=-1,9) / 9,1,1,0,0,0,0,0,0,0,1 /
        DATA (IRRED(73,I), I=-1,9) / 9,1,0,0,0,1,0,0,0,0,1 /
        DATA (IRRED(74,I), I=-1,9) / 9,1,1,1,0,1,0,0,0,0,1 /
        DATA (IRRED(75,I), I=-1,9) / 9,1,1,0,1,1,0,0,0,0,1 /
        DATA (IRRED(76,I), I=-1,9) / 9,1,0,0,0,0,1,0,0,0,1 /
        DATA (IRRED(77,I), I=-1,9) / 9,1,0,1,1,0,1,0,0,0,1 /
        DATA (IRRED(78,I), I=-1,9) / 9,1,1,0,0,1,1,0,0,0,1 /
        DATA (IRRED(79,I), I=-1,9) / 9,1,1,0,1,0,0,1,0,0,1 /
        DATA (IRRED(80,I), I=-1,9) / 9,1,0,0,1,1,0,1,0,0,1 /
        DATA (IRRED(81,I), I=-1,9) / 9,1,1,1,1,1,0,1,0,0,1 /
        DATA (IRRED(82,I), I=-1,9) / 9,1,0,1,0,0,1,1,0,0,1 /
        DATA (IRRED(83,I), I=-1,9) / 9,1,0,0,1,0,1,1,0,0,1 /
        DATA (IRRED(84,I), I=-1,9) / 9,1,1,1,1,0,1,1,0,0,1 /
        DATA (IRRED(85,I), I=-1,9) / 9,1,1,1,0,1,1,1,0,0,1 /
        DATA (IRRED(86,I), I=-1,9) / 9,1,0,1,1,1,1,1,0,0,1 /
        DATA (IRRED(87,I), I=-1,9) / 9,1,1,1,0,0,0,0,1,0,1 /
        DATA (IRRED(88,I), I=-1,9) / 9,1,0,1,0,1,0,0,1,0,1 /
        DATA (IRRED(89,I), I=-1,9) / 9,1,0,0,1,1,0,0,1,0,1 /
        DATA (IRRED(90,I), I=-1,9) / 9,1,1,0,0,0,1,0,1,0,1 /
        DATA (IRRED(91,I), I=-1,9) / 9,1,0,1,0,0,1,0,1,0,1 /
        DATA (IRRED(92,I), I=-1,9) / 9,1,1,1,1,0,1,0,1,0,1 /
        DATA (IRRED(93,I), I=-1,9) / 9,1,1,1,0,1,1,0,1,0,1 /
        DATA (IRRED(94,I), I=-1,9) / 9,1,0,1,1,1,1,0,1,0,1 /
        DATA (IRRED(95,I), I=-1,9) / 9,1,1,1,1,0,0,1,1,0,1 /
        DATA (IRRED(96,I), I=-1,9) / 9,1,0,0,0,1,0,1,1,0,1 /
        DATA (IRRED(97,I), I=-1,9) / 9,1,1,0,1,1,0,1,1,0,1 /
        DATA (IRRED(98,I), I=-1,9) / 9,1,0,1,0,1,1,1,1,0,1 /
        DATA (IRRED(99,I), I=-1,9) / 9,1,0,0,1,1,1,1,1,0,1 /
        DATA (IRRED(100,I), I=-1,9) / 9,1,0,0,0,0,0,0,0,1,1 /
        DATA (IRRED(101,I), I=-1,9) / 9,1,1,0,0,1,0,0,0,1,1 /
        DATA (IRRED(102,I), I=-1,9) / 9,1,0,1,0,1,0,0,0,1,1 /
        DATA (IRRED(103,I), I=-1,9) / 9,1,1,1,1,1,0,0,0,1,1 /
        DATA (IRRED(104,I), I=-1,9) / 9,1,1,0,0,0,1,0,0,1,1 /
        DATA (IRRED(105,I), I=-1,9) / 9,1,0,0,0,1,1,0,0,1,1 /
        DATA (IRRED(106,I), I=-1,9) / 9,1,1,0,1,1,1,0,0,1,1 /
        DATA (IRRED(107,I), I=-1,9) / 9,1,0,0,1,0,0,1,0,1,1 /
        DATA (IRRED(108,I), I=-1,9) / 9,1,1,1,1,0,0,1,0,1,1 /
        DATA (IRRED(109,I), I=-1,9) / 9,1,1,0,1,1,0,1,0,1,1 /
        DATA (IRRED(110,I), I=-1,9) / 9,1,0,0,0,0,1,1,0,1,1 /
        DATA (IRRED(111,I), I=-1,9) / 9,1,1,0,1,0,1,1,0,1,1 /
        DATA (IRRED(112,I), I=-1,9) / 9,1,0,1,1,0,1,1,0,1,1 /
        DATA (IRRED(113,I), I=-1,9) / 9,1,1,0,0,1,1,1,0,1,1 /
        DATA (IRRED(114,I), I=-1,9) / 9,1,1,1,1,1,1,1,0,1,1 /
        DATA (IRRED(115,I), I=-1,9) / 9,1,0,1,0,0,0,0,1,1,1 /
        DATA (IRRED(116,I), I=-1,9) / 9,1,1,1,1,0,0,0,1,1,1 /
        DATA (IRRED(117,I), I=-1,9) / 9,1,0,0,0,0,1,0,1,1,1 /
        DATA (IRRED(118,I), I=-1,9) / 9,1,0,1,0,1,1,0,1,1,1 /
        DATA (IRRED(119,I), I=-1,9) / 9,1,0,0,1,1,1,0,1,1,1 /
        DATA (IRRED(120,I), I=-1,9) / 9,1,1,1,0,0,0,1,1,1,1 /
        DATA (IRRED(121,I), I=-1,9) / 9,1,1,0,1,0,0,1,1,1,1 /
        DATA (IRRED(122,I), I=-1,9) / 9,1,0,1,1,0,0,1,1,1,1 /
        DATA (IRRED(123,I), I=-1,9) / 9,1,0,1,0,1,0,1,1,1,1 /
        DATA (IRRED(124,I), I=-1,9) / 9,1,0,0,1,1,0,1,1,1,1 /
        DATA (IRRED(125,I), I=-1,9) / 9,1,1,0,0,0,1,1,1,1,1 /
        DATA (IRRED(126,I), I=-1,9) / 9,1,0,0,1,0,1,1,1,1,1 /
        DATA (IRRED(127,I), I=-1,9) / 9,1,1,0,1,1,1,1,1,1,1 /
        DATA (IRRED(128,I), I=-1,10) / 10,1,0,0,1,0,0,0,0,0,0,1 /
        DATA (IRRED(129,I), I=-1,10) / 10,1,1,1,1,0,0,0,0,0,0,1 /
        DATA (IRRED(130,I), I=-1,10) / 10,1,1,0,1,1,0,0,0,0,0,1 /
        DATA (IRRED(131,I), I=-1,10) / 10,1,0,1,1,1,0,0,0,0,0,1 /
        DATA (IRRED(132,I), I=-1,10) / 10,1,1,1,0,0,1,0,0,0,0,1 /
        DATA (IRRED(133,I), I=-1,10) / 10,1,0,1,1,0,1,0,0,0,0,1 /
        DATA (IRRED(134,I), I=-1,10) / 10,1,0,1,0,1,1,0,0,0,0,1 /
        DATA (IRRED(135,I), I=-1,10) / 10,1,1,1,0,0,0,1,0,0,0,1 /
        DATA (IRRED(136,I), I=-1,10) / 10,1,1,0,0,1,0,1,0,0,0,1 /
        DATA (IRRED(137,I), I=-1,10) / 10,1,1,0,0,0,1,1,0,0,0,1 /
        DATA (IRRED(138,I), I=-1,10) / 10,1,0,1,0,0,1,1,0,0,0,1 /
        DATA (IRRED(139,I), I=-1,10) / 10,1,1,1,1,0,1,1,0,0,0,1 /
        DATA (IRRED(140,I), I=-1,10) / 10,1,0,0,0,0,0,0,1,0,0,1 /
        DATA (IRRED(141,I), I=-1,10) / 10,1,1,0,1,0,0,0,1,0,0,1 /
        DATA (IRRED(142,I), I=-1,10) / 10,1,0,0,1,1,0,0,1,0,0,1 /
        DATA (IRRED(143,I), I=-1,10) / 10,1,0,0,1,0,1,0,1,0,0,1 /
        DATA (IRRED(144,I), I=-1,10) / 10,1,1,1,1,0,1,0,1,0,0,1 /
        DATA (IRRED(145,I), I=-1,10) / 10,1,0,1,0,0,0,1,1,0,0,1 /
        DATA (IRRED(146,I), I=-1,10) / 10,1,0,0,1,0,0,1,1,0,0,1 /
        DATA (IRRED(147,I), I=-1,10) / 10,1,1,1,0,1,0,1,1,0,0,1 /
        DATA (IRRED(148,I), I=-1,10) / 10,1,1,1,0,0,1,1,1,0,0,1 /
        DATA (IRRED(149,I), I=-1,10) / 10,1,0,1,1,0,1,1,1,0,0,1 /
        DATA (IRRED(150,I), I=-1,10) / 10,1,1,0,0,1,1,1,1,0,0,1 /
        DATA (IRRED(151,I), I=-1,10) / 10,1,1,1,1,1,1,1,1,0,0,1 /
        DATA (IRRED(152,I), I=-1,10) / 10,1,1,0,1,0,0,0,0,1,0,1 /
        DATA (IRRED(153,I), I=-1,10) / 10,1,0,1,1,0,0,0,0,1,0,1 /
        DATA (IRRED(154,I), I=-1,10) / 10,1,0,0,1,1,0,0,0,1,0,1 /
        DATA (IRRED(155,I), I=-1,10) / 10,1,1,1,1,1,0,0,0,1,0,1 /
        DATA (IRRED(156,I), I=-1,10) / 10,1,1,0,0,0,1,0,0,1,0,1 /
        DATA (IRRED(157,I), I=-1,10) / 10,1,0,0,0,1,1,0,0,1,0,1 /
        DATA (IRRED(158,I), I=-1,10) / 10,1,0,1,1,1,1,0,0,1,0,1 /
        DATA (IRRED(159,I), I=-1,10) / 10,1,1,0,0,0,0,1,0,1,0,1 /
        DATA (IRRED(160,I), I=-1,10) / 10,1,1,1,0,1,0,1,0,1,0,1 /
        DATA (IRRED(161,I), I=-1,10) / 10,1,0,0,0,0,1,1,0,1,0,1 /
        DATA (IRRED(162,I), I=-1,10) / 10,1,1,1,0,0,1,1,0,1,0,1 /
        DATA (IRRED(163,I), I=-1,10) / 10,1,1,0,1,0,1,1,0,1,0,1 /
        DATA (IRRED(164,I), I=-1,10) / 10,1,0,1,0,0,0,0,1,1,0,1 /
        DATA (IRRED(165,I), I=-1,10) / 10,1,1,1,1,0,0,0,1,1,0,1 /
        DATA (IRRED(166,I), I=-1,10) / 10,1,1,1,0,1,0,0,1,1,0,1 /
        DATA (IRRED(167,I), I=-1,10) / 10,1,1,0,1,1,0,0,1,1,0,1 /
        DATA (IRRED(168,I), I=-1,10) / 10,1,0,0,0,0,1,0,1,1,0,1 /
        DATA (IRRED(169,I), I=-1,10) / 10,1,1,0,1,0,1,0,1,1,0,1 /
        DATA (IRRED(170,I), I=-1,10) / 10,1,0,0,1,1,1,0,1,1,0,1 /
        DATA (IRRED(171,I), I=-1,10) / 10,1,0,0,0,0,0,1,1,1,0,1 /
        DATA (IRRED(172,I), I=-1,10) / 10,1,1,1,0,0,0,1,1,1,0,1 /
        DATA (IRRED(173,I), I=-1,10) / 10,1,0,1,0,0,1,1,1,1,0,1 /
        DATA (IRRED(174,I), I=-1,10) / 10,1,1,1,0,1,1,1,1,1,0,1 /
        DATA (IRRED(175,I), I=-1,10) / 10,1,1,0,1,1,1,1,1,1,0,1 /
        DATA (IRRED(176,I), I=-1,10) / 10,1,1,0,0,1,0,0,0,0,1,1 /
        DATA (IRRED(177,I), I=-1,10) / 10,1,0,1,0,1,0,0,0,0,1,1 /
        DATA (IRRED(178,I), I=-1,10) / 10,1,1,0,0,0,1,0,0,0,1,1 /
        DATA (IRRED(179,I), I=-1,10) / 10,1,0,1,0,0,1,0,0,0,1,1 /
        DATA (IRRED(180,I), I=-1,10) / 10,1,0,0,0,1,1,0,0,0,1,1 /
        DATA (IRRED(181,I), I=-1,10) / 10,1,1,1,0,1,1,0,0,0,1,1 /
        DATA (IRRED(182,I), I=-1,10) / 10,1,1,0,0,0,0,1,0,0,1,1 /
        DATA (IRRED(183,I), I=-1,10) / 10,1,1,1,1,0,0,1,0,0,1,1 /
        DATA (IRRED(184,I), I=-1,10) / 10,1,0,0,0,1,0,1,0,0,1,1 /
        DATA (IRRED(185,I), I=-1,10) / 10,1,1,0,1,1,0,1,0,0,1,1 /
        DATA (IRRED(186,I), I=-1,10) / 10,1,0,0,1,1,1,1,0,0,1,1 /
        DATA (IRRED(187,I), I=-1,10) / 10,1,1,1,1,1,1,1,0,0,1,1 /
        DATA (IRRED(188,I), I=-1,10) / 10,1,0,1,0,0,0,0,1,0,1,1 /
        DATA (IRRED(189,I), I=-1,10) / 10,1,0,0,1,0,0,0,1,0,1,1 /
        DATA (IRRED(190,I), I=-1,10) / 10,1,1,1,0,0,1,0,1,0,1,1 /
        DATA (IRRED(191,I), I=-1,10) / 10,1,0,1,1,0,1,0,1,0,1,1 /
        DATA (IRRED(192,I), I=-1,10) / 10,1,0,1,0,1,1,0,1,0,1,1 /
        DATA (IRRED(193,I), I=-1,10) / 10,1,1,1,1,1,1,0,1,0,1,1 /
        DATA (IRRED(194,I), I=-1,10) / 10,1,0,0,0,0,0,1,1,0,1,1 /
        DATA (IRRED(195,I), I=-1,10) / 10,1,0,1,1,0,0,1,1,0,1,1 /
        DATA (IRRED(196,I), I=-1,10) / 10,1,1,0,0,1,0,1,1,0,1,1 /
        DATA (IRRED(197,I), I=-1,10) / 10,1,1,1,1,1,0,1,1,0,1,1 /
        DATA (IRRED(198,I), I=-1,10) / 10,1,1,1,0,1,1,1,1,0,1,1 /
        DATA (IRRED(199,I), I=-1,10) / 10,1,0,1,1,1,1,1,1,0,1,1 /
        DATA (IRRED(200,I), I=-1,10) / 10,1,1,1,1,0,0,0,0,1,1,1 /
        DATA (IRRED(201,I), I=-1,10) / 10,1,0,0,0,1,0,0,0,1,1,1 /
        DATA (IRRED(202,I), I=-1,10) / 10,1,1,1,0,1,0,0,0,1,1,1 /
        DATA (IRRED(203,I), I=-1,10) / 10,1,0,1,1,1,0,0,0,1,1,1 /
        DATA (IRRED(204,I), I=-1,10) / 10,1,0,0,0,0,1,0,0,1,1,1 /
        DATA (IRRED(205,I), I=-1,10) / 10,1,1,0,1,0,1,0,0,1,1,1 /
        DATA (IRRED(206,I), I=-1,10) / 10,1,0,1,0,1,1,0,0,1,1,1 /
        DATA (IRRED(207,I), I=-1,10) / 10,1,0,0,1,1,1,0,0,1,1,1 /
        DATA (IRRED(208,I), I=-1,10) / 10,1,1,1,0,0,0,1,0,1,1,1 /
        DATA (IRRED(209,I), I=-1,10) / 10,1,0,1,1,0,0,1,0,1,1,1 /
        DATA (IRRED(210,I), I=-1,10) / 10,1,0,1,0,1,0,1,0,1,1,1 /
        DATA (IRRED(211,I), I=-1,10) / 10,1,0,0,1,1,0,1,0,1,1,1 /
        DATA (IRRED(212,I), I=-1,10) / 10,1,1,0,0,0,1,1,0,1,1,1 /
        DATA (IRRED(213,I), I=-1,10) / 10,1,1,0,1,1,1,1,0,1,1,1 /
        DATA (IRRED(214,I), I=-1,10) / 10,1,0,1,1,1,1,1,0,1,1,1 /
        DATA (IRRED(215,I), I=-1,10) / 10,1,0,0,0,0,0,0,1,1,1,1 /
        DATA (IRRED(216,I), I=-1,10) / 10,1,1,1,0,0,0,0,1,1,1,1 /
        DATA (IRRED(217,I), I=-1,10) / 10,1,0,1,1,0,0,0,1,1,1,1 /
        DATA (IRRED(218,I), I=-1,10) / 10,1,1,0,0,1,0,0,1,1,1,1 /
        DATA (IRRED(219,I), I=-1,10) / 10,1,0,0,1,0,1,0,1,1,1,1 /
        DATA (IRRED(220,I), I=-1,10) / 10,1,0,0,0,1,1,0,1,1,1,1 /
        DATA (IRRED(221,I), I=-1,10) / 10,1,0,1,0,0,0,1,1,1,1,1 /
        DATA (IRRED(222,I), I=-1,10) / 10,1,1,0,1,1,0,1,1,1,1,1 /
        DATA (IRRED(223,I), I=-1,10) / 10,1,1,0,1,0,1,1,1,1,1,1 /
        DATA (IRRED(224,I), I=-1,10) / 10,1,1,0,0,1,1,1,1,1,1,1 /
        DATA (IRRED(225,I), I=-1,10) / 10,1,0,0,1,1,1,1,1,1,1,1 /
        DATA (IRRED(226,I), I=-1,10) / 10,1,1,1,1,1,1,1,1,1,1,1 /
        DATA (IRRED(227,I), I=-1,11) / 11,1,0,1,0,0,0,0,0,0,0,0,1 /
        DATA (IRRED(228,I), I=-1,11) / 11,1,1,1,0,1,0,0,0,0,0,0,1 /
        DATA (IRRED(229,I), I=-1,11) / 11,1,1,0,1,0,1,0,0,0,0,0,1 /
        DATA (IRRED(230,I), I=-1,11) / 11,1,0,1,1,0,1,0,0,0,0,0,1 /
        DATA (IRRED(231,I), I=-1,11) / 11,1,1,1,0,0,0,1,0,0,0,0,1 /
        DATA (IRRED(232,I), I=-1,11) / 11,1,1,0,0,0,1,1,0,0,0,0,1 /
        DATA (IRRED(233,I), I=-1,11) / 11,1,0,1,0,0,1,1,0,0,0,0,1 /
        DATA (IRRED(234,I), I=-1,11) / 11,1,0,0,0,1,1,1,0,0,0,0,1 /
        DATA (IRRED(235,I), I=-1,11) / 11,1,1,0,1,1,1,1,0,0,0,0,1 /
        DATA (IRRED(236,I), I=-1,11) / 11,1,0,1,1,0,0,0,1,0,0,0,1 /
        DATA (IRRED(237,I), I=-1,11) / 11,1,0,1,0,1,0,0,1,0,0,0,1 /
        DATA (IRRED(238,I), I=-1,11) / 11,1,1,1,1,1,0,0,1,0,0,0,1 /
        DATA (IRRED(239,I), I=-1,11) / 11,1,0,0,1,0,1,0,1,0,0,0,1 /
        DATA (IRRED(240,I), I=-1,11) / 11,1,0,0,0,1,1,0,1,0,0,0,1 /
        DATA (IRRED(241,I), I=-1,11) / 11,1,1,0,0,0,0,1,1,0,0,0,1 /
        DATA (IRRED(242,I), I=-1,11) / 11,1,1,1,1,0,0,1,1,0,0,0,1 /
        DATA (IRRED(243,I), I=-1,11) / 11,1,0,0,0,1,0,1,1,0,0,0,1 /
        DATA (IRRED(244,I), I=-1,11) / 11,1,0,0,0,0,1,1,1,0,0,0,1 /
        DATA (IRRED(245,I), I=-1,11) / 11,1,1,1,0,0,1,1,1,0,0,0,1 /
        DATA (IRRED(246,I), I=-1,11) / 11,1,1,0,1,0,1,1,1,0,0,0,1 /
        DATA (IRRED(247,I), I=-1,11) / 11,1,0,1,0,1,1,1,1,0,0,0,1 /
        DATA (IRRED(248,I), I=-1,11) / 11,1,0,1,1,0,0,0,0,1,0,0,1 /
        DATA (IRRED(249,I), I=-1,11) / 11,1,1,0,0,1,0,0,0,1,0,0,1 /
        DATA (IRRED(250,I), I=-1,11) / 11,1,0,1,0,0,1,0,0,1,0,0,1 /
        DATA (IRRED(251,I), I=-1,11) / 11,1,0,0,1,0,1,0,0,1,0,0,1 /
        DATA (IRRED(252,I), I=-1,11) / 11,1,1,1,0,1,1,0,0,1,0,0,1 /
        DATA (IRRED(253,I), I=-1,11) / 11,1,1,0,1,1,1,0,0,1,0,0,1 /
        DATA (IRRED(254,I), I=-1,11) / 11,1,0,1,1,1,1,0,0,1,0,0,1 /
        DATA (IRRED(255,I), I=-1,11) / 11,1,0,1,0,0,0,1,0,1,0,0,1 /
        DATA (IRRED(256,I), I=-1,11) / 11,1,0,0,1,0,0,1,0,1,0,0,1 /
        DATA (IRRED(257,I), I=-1,11) / 11,1,0,0,0,1,0,1,0,1,0,0,1 /
        DATA (IRRED(258,I), I=-1,11) / 11,1,1,0,1,1,0,1,0,1,0,0,1 /
        DATA (IRRED(259,I), I=-1,11) / 11,1,1,0,0,1,1,1,0,1,0,0,1 /
        DATA (IRRED(260,I), I=-1,11) / 11,1,0,1,0,1,1,1,0,1,0,0,1 /
        DATA (IRRED(261,I), I=-1,11) / 11,1,1,1,1,1,1,1,0,1,0,0,1 /
        DATA (IRRED(262,I), I=-1,11) / 11,1,1,0,0,0,0,0,1,1,0,0,1 /
        DATA (IRRED(263,I), I=-1,11) / 11,1,1,1,1,0,0,0,1,1,0,0,1 /
        DATA (IRRED(264,I), I=-1,11) / 11,1,1,0,1,0,1,0,1,1,0,0,1 /
        DATA (IRRED(265,I), I=-1,11) / 11,1,0,1,1,0,1,0,1,1,0,0,1 /
        DATA (IRRED(266,I), I=-1,11) / 11,1,0,0,1,1,1,0,1,1,0,0,1 /
        DATA (IRRED(267,I), I=-1,11) / 11,1,1,1,0,0,0,1,1,1,0,0,1 /
        DATA (IRRED(268,I), I=-1,11) / 11,1,0,0,1,1,0,1,1,1,0,0,1 /
        DATA (IRRED(269,I), I=-1,11) / 11,1,0,1,0,0,1,1,1,1,0,0,1 /
        DATA (IRRED(270,I), I=-1,11) / 11,1,1,1,1,0,1,1,1,1,0,0,1 /
        DATA (IRRED(271,I), I=-1,11) / 11,1,1,1,0,1,1,1,1,1,0,0,1 /
        DATA (IRRED(272,I), I=-1,11) / 11,1,0,0,0,0,0,0,0,0,1,0,1 /
        DATA (IRRED(273,I), I=-1,11) / 11,1,1,1,0,0,0,0,0,0,1,0,1 /
        DATA (IRRED(274,I), I=-1,11) / 11,1,1,0,0,1,0,0,0,0,1,0,1 /
        DATA (IRRED(275,I), I=-1,11) / 11,1,0,1,0,1,0,0,0,0,1,0,1 /
        DATA (IRRED(276,I), I=-1,11) / 11,1,0,0,1,0,1,0,0,0,1,0,1 /
        DATA (IRRED(277,I), I=-1,11) / 11,1,0,0,1,0,0,1,0,0,1,0,1 /
        DATA (IRRED(278,I), I=-1,11) / 11,1,0,0,0,0,1,1,0,0,1,0,1 /
        DATA (IRRED(279,I), I=-1,11) / 11,1,0,1,1,0,1,1,0,0,1,0,1 /
        DATA (IRRED(280,I), I=-1,11) / 11,1,0,0,1,1,1,1,0,0,1,0,1 /
        DATA (IRRED(281,I), I=-1,11) / 11,1,1,1,1,1,1,1,0,0,1,0,1 /
        DATA (IRRED(282,I), I=-1,11) / 11,1,0,1,0,0,0,0,1,0,1,0,1 /
        DATA (IRRED(283,I), I=-1,11) / 11,1,0,0,0,1,0,0,1,0,1,0,1 /
        DATA (IRRED(284,I), I=-1,11) / 11,1,0,1,1,1,0,0,1,0,1,0,1 /
        DATA (IRRED(285,I), I=-1,11) / 11,1,1,1,0,0,1,0,1,0,1,0,1 /
        DATA (IRRED(286,I), I=-1,11) / 11,1,1,0,1,0,1,0,1,0,1,0,1 /
        DATA (IRRED(287,I), I=-1,11) / 11,1,1,0,0,1,1,0,1,0,1,0,1 /
        DATA (IRRED(288,I), I=-1,11) / 11,1,0,1,0,1,1,0,1,0,1,0,1 /
        DATA (IRRED(289,I), I=-1,11) / 11,1,0,1,0,1,0,1,1,0,1,0,1 /
        DATA (IRRED(290,I), I=-1,11) / 11,1,1,1,1,1,0,1,1,0,1,0,1 /
        DATA (IRRED(291,I), I=-1,11) / 11,1,1,0,0,0,1,1,1,0,1,0,1 /
        DATA (IRRED(292,I), I=-1,11) / 11,1,0,0,1,0,1,1,1,0,1,0,1 /
        DATA (IRRED(293,I), I=-1,11) / 11,1,1,1,1,0,1,1,1,0,1,0,1 /
        DATA (IRRED(294,I), I=-1,11) / 11,1,0,0,0,1,1,1,1,0,1,0,1 /
        DATA (IRRED(295,I), I=-1,11) / 11,1,1,0,1,1,1,1,1,0,1,0,1 /
        DATA (IRRED(296,I), I=-1,11) / 11,1,1,0,0,0,0,0,0,1,1,0,1 /
        DATA (IRRED(297,I), I=-1,11) / 11,1,0,0,1,0,0,0,0,1,1,0,1 /
        DATA (IRRED(298,I), I=-1,11) / 11,1,0,0,0,1,0,0,0,1,1,0,1 /
        DATA (IRRED(299,I), I=-1,11) / 11,1,1,0,0,1,1,0,0,1,1,0,1 /
        DATA (IRRED(300,I), I=-1,11) / 11,1,1,1,1,1,1,0,0,1,1,0,1 /
        DATA (IRRED(301,I), I=-1,11) / 11,1,0,0,0,0,0,1,0,1,1,0,1 /
        DATA (IRRED(302,I), I=-1,11) / 11,1,1,0,1,0,0,1,0,1,1,0,1 /
        DATA (IRRED(303,I), I=-1,11) / 11,1,0,0,1,1,0,1,0,1,1,0,1 /
        DATA (IRRED(304,I), I=-1,11) / 11,1,1,1,1,1,0,1,0,1,1,0,1 /
        DATA (IRRED(305,I), I=-1,11) / 11,1,0,1,0,0,1,1,0,1,1,0,1 /
        DATA (IRRED(306,I), I=-1,11) / 11,1,1,1,1,0,1,1,0,1,1,0,1 /
        DATA (IRRED(307,I), I=-1,11) / 11,1,0,1,1,1,1,1,0,1,1,0,1 /
        DATA (IRRED(308,I), I=-1,11) / 11,1,1,1,0,0,0,0,1,1,1,0,1 /
        DATA (IRRED(309,I), I=-1,11) / 11,1,1,0,1,0,0,0,1,1,1,0,1 /
        DATA (IRRED(310,I), I=-1,11) / 11,1,1,0,0,1,0,0,1,1,1,0,1 /
        DATA (IRRED(311,I), I=-1,11) / 11,1,0,1,0,1,0,0,1,1,1,0,1 /
        DATA (IRRED(312,I), I=-1,11) / 11,1,1,1,1,0,1,0,1,1,1,0,1 /
        DATA (IRRED(313,I), I=-1,11) / 11,1,1,1,0,1,1,0,1,1,1,0,1 /
        DATA (IRRED(314,I), I=-1,11) / 11,1,0,1,1,1,1,0,1,1,1,0,1 /
        DATA (IRRED(315,I), I=-1,11) / 11,1,0,0,1,0,0,1,1,1,1,0,1 /
        DATA (IRRED(316,I), I=-1,11) / 11,1,1,0,1,1,0,1,1,1,1,0,1 /
        DATA (IRRED(317,I), I=-1,11) / 11,1,0,1,1,1,0,1,1,1,1,0,1 /
        DATA (IRRED(318,I), I=-1,11) / 11,1,1,1,0,0,1,1,1,1,1,0,1 /
C
C   all the 10 th degree polynomials for mod 2 are used; n=10
C   the first 92, 11 th degree polyn. are also included; n=11.
C   niederreiter-lidl book, pg 385 - the first 3 columns are done
C
C   Prepare to work in Z2
C
       CALL SETFLD (2)
C
       DO 1000 I = 1, DIMEN
C
C   For each dimension, we need to calculate powers of an
C   appropriate irreducible polynomial :  see Niederreiter
C   page 65, just below equation (19).
C     Copy the appropriate irreducible polynomial into PX,
C   and its degree into E.  Set polynomial B = PX ** 0 = 1.
C   M is the degree of B.  Subsequently B will hold higher
C   powers of PX.
C
         E = IRRED(I, DEG)
         DO 10 J = -1, E
           PX(J) = IRRED(I,J)
10        CONTINUE
         B(DEG) = 0
         B(0) = 1
C
C   Niederreiter (page 56, after equation (7), defines two
C   variables Q and U.  We do not need Q explicitly, but we
C   do need U.
C
         U = 0
C
         DO 90 J = 1, NBITS
C
C   If U = 0, we need to set B to the next power of PX
C   and recalculate V.  This is done by subroutine CALCV.
C
           IF (U .EQ. 0) CALL CALCV (PX, B, V, MAXV)
C
C Now C is obtained from V.  Niederreiter
C obtains A from V (page 65, near the bottom), and then gets
C C from A (page 56, equation (7)).  However this can be done
C in one step.  Here CI(J,R) corresponds to
C Niederreiter's C(I,J,R).
C
           DO 50 R = 0, NBITS-1
             CI(J,R) = V(R+U)
 50         CONTINUE
C
C Increment U.  If U = E, then U = 0 and in Niederreiter's
C paper Q = Q + 1.  Here, however, Q is not used explicitly.
C
           U = U + 1
           IF (U .EQ. E) U = 0
 90      CONTINUE
C
C  The array CI now holds the values of C(I,J,R) for this value
C  of I.  We pack them into array CJ so that CJ(I,R) holds all
C  the values of C(I,J,R) for J from 1 to NBITS.
C
         DO 120 R = 0, NBITS-1
           TERM = 0
           DO 110 J = 1, NBITS
             TERM = 2 * TERM + CI(J,R)
 110       CONTINUE
           CJ(I,R) = TERM
 120      CONTINUE
C
 1000   CONTINUE

C 
C COMPUTING GENERATOR MATRICES OF USER CHOICE
C
      IF (IFLAG .EQ. 0) THEN
        DO 200 I= 1, DIMEN
            DO 210 J= 0,NBITS-1 
                SCJ(I,J) = CJ(I,J)
 210         CONTINUE
             SHIFT(I) = 0
 200     CONTINUE 
      ELSE
        IF ((IFLAG .EQ. 1) .OR. (IFLAG .EQ. 3)) THEN
         CALL GENSCRML(MAXDIM,NBITS,DIMEN,LSM,SHIFT)     
         DO 122 I = 1,DIMEN
           DO 123 J = 0,NBITS-1
             L = 1
             TEMP2 = 0
             DO 124 P = NBITS-1,0,-1
                TEMP1 = 0
                DO 125 K = 0,NBITS-1
                   TEMP1 = TEMP1+
     *                     (IBITS(LSM(I,P),K,1)*IBITS(CJ(I,J),K,1))
 125              CONTINUE
                TEMP1 = MOD(TEMP1,2)
                TEMP2 = TEMP2+TEMP1*L   
                L = 2 * L
 124          CONTINUE
              SCJ(I,J) = TEMP2
 123      CONTINUE
 122    CONTINUE
      ENDIF       

      IF ((IFLAG .EQ. 2) .OR. (IFLAG .EQ. 3)) THEN
        CALL GENSCRMU(NBITS,DIMEN,USM,USHIFT)
        DO 230 I = 1,DIMEN
         DO 240 J = 0,NBITS-1
            P = NBITS-1
            DO 250 K = 0, NBITS-1
                IF (IFLAG .EQ. 2) THEN
                   TV(I,P,J) = IBITS(CJ(I,J),K,1)
                ELSE
                    TV(I,P,J) = IBITS(SCJ(I,J),K,1) 
                ENDIF
                P = P-1
 250          CONTINUE           
 240       CONTINUE    
          DO 241 PP = 0, NBITS-1
             TEMP2 = 0
             TEMP4 = 0
             L = 1
             DO 260 J= NBITS-1,0,-1
                TEMP1 = 0
                TEMP3 = 0
                DO 270 P = 0, NBITS-1
                   TEMP1=TEMP1+TV(I,J,P)*USM(P,PP)
                   IF (PP. EQ. 0) THEN
                    TEMP3 = TEMP3+TV(I,J,P)*USHIFT(P)
                   ENDIF
 270              CONTINUE
                TEMP1 = MOD(TEMP1,2)
                TEMP2 = TEMP2 + TEMP1*L
                IF (PP .EQ. 1) THEN
                  TEMP3= MOD(TEMP1,2)
                  TEMP4 = TEMP4+TEMP3*L
                ENDIF
                   L= 2*L                               
 260           CONTINUE
               SCJ(I,PP) = TEMP2
               IF (PP .EQ. 0) THEN
                 IF (IFLAG .EQ. 3) THEN
                   SHIFT(I) = EXOR(TEMP4,SHIFT(I))
                 ELSE
                   SHIFT(I) = TEMP4
                 ENDIF
               ENDIF
  241       CONTINUE                     
 230      CONTINUE          
        ENDIF      
       ENDIF      

      DO 400 I = 1,DIMEN
          NEXTQ(I) = SHIFT(I)
 400     CONTINUE         
      RETURN
      END
 
C      *****  end of SCALCC2


      SUBROUTINE GENSCRML(MAXDIM,NBITS,DIMEN,LSM,SHIFT)
C     GENERATING LOWER TRIAGULAR SCRMABLING MATRICES 
C     AND SHIFT VECTORS.
C     INPUTS :
C       FROM SCALCC2 : MAXDIM,NBITS,DIMEN
C
C     OUTPUTS :
C       TO SCALCC2 : LSM, SHIFT
      
      INTEGER MAXDIM,NBITS,DIMEN
      INTEGER LSM(MAXDIM,0:NBITS),MAXCOL,P,I,J
      INTEGER SHIFT(MAXDIM),MAXS,S,TEMP,STEMP,L,LL
      REAL*8 UNI

      DO 10 P = 1,DIMEN
               SHIFT(P) = 0
               L = 1
         DO 20 I = NBITS-1,0,-1
                LSM(P,I) = 0
                STEMP =  MOD((int(UNI()*1000.0)),2)
                SHIFT(P) = SHIFT(P)+STEMP*L  
                L = 2 * L
                LL = 1
            DO 30 J = NBITS-1,0,-1
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
      
      SUBROUTINE GENSCRMU(NBITS,DIMEN,USM,USHIFT)

C     GENERATING UPPER TRIAGULAR SCRMABLING MATRIX
C     AND SHIFT VECTOR.
C     INPUTS :
C       FROM SCALCC2: NBITS,DIMEN
C
C     OUTPUTS :
C       TO SCALCC2 : USM, USHIFT
      
      INTEGER NBITS,DIMEN
      INTEGER USM(0:NBITS-1,0:NBITS-1),MAXCOL,P,I,J
      INTEGER USHIFT(0:NBITS-1),MAXS,S,TEMP,STEMP,L,LL
      REAL*8 UNI

         DO 20 I = 0,NBITS-1
               STEMP =  MOD((int(UNI()*1000.0)),2)
               USHIFT(I) = STEMP               
            DO 30 J = 0,NBITS-1 
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
C
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

      SUBROUTINE CALCV (PX, B, V, MAXV)
C
C   This version :  12 February 1991
C
C   See the general comments on implementing Niederreiter's
C   low-discrepancy sequences.
C
C   This program calculates the values of the constants V(J,R) as
C   described in BFN section 3.3.  It is called from either CALCC or
C   SCALCC2.  The values transmitted through common /FIELD/ determine
C   which field we are working in.
C
C INPUT :
C   PX is the appropriate irreducible polynomial for the dimension
C     currently being considered.  The degree of PX will be called E.
C   B is the polynomial defined in section 2.3 of BFN.  On entry to
C     the subroutine, the degree of B implicitly defines the parameter
C     J of section 3.3, by degree(B) = E*(J-1).
C   MAXV gives the dimension of the array V.
C   On entry, we suppose that the common block /FIELD/ has been set
C     up correctly (using SETFLD).
C
C OUTPUT :
C   On output B has been multiplied by PX, so its degree is now E*J.
C   V contains the values required.
C
C USES :
C   The subroutine PLYMUL is used to multiply polynomials.
C
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication, and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
 
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C
C
C
      INTEGER MAXV, E, I, J, KJ, M, BIGM, R, TERM
      INTEGER PX(-1:MAXDEG), B(-1:MAXDEG), V(0:MAXV)
      INTEGER H(-1:MAXDEG)
C
      INTEGER ARBIT, NONZER
      ARBIT() = 1
C
C   We use ARBIT() to indicate where the user can place
C   an arbitrary element of the field of order Q, while NONZER
C   shows where he must put an arbitrary non-zero element
C   of the same field.  For the code,
C   this means 0 <= ARBIT < Q and 0 < NONZER < Q.  Within these
C   limits, the user can do what he likes.  ARBIT is declared as
C   a function as a reminder that a different arbitrary value may
C   be returned each time ARBIT is referenced.
C
C    BIGM is the M used in section 3.3.
C    It differs from the [little] m used in section 2.3,
C    denoted here by M.
C
      NONZER = 1
C
      E = PX(DEG)
C
C   The poly H is PX**(J-1), which is the value of B on arrival.
C   In section 3.3, the values of Hi are defined with a minus sign :
C   don't forget this if you use them later !
C
      DO 10 I = -1, B(DEG)
   10   H(I) = B(I)
      BIGM = H(DEG)
C
C   Now multiply B by PX so B becomes PX**J.
C   In section 2.3, the values of Bi are defined with a minus sign :
C   don't forget this if you use them later !
C
      CALL PLYMUL (PX, B, B)
      M = B(DEG)
C
C   We don't use J explicitly anywhere, but here it is just in case.
C
      J = M / E
C
C   Now choose a value of Kj as defined in section 3.3.
C   We must have 0 <= Kj < E*J = M.
C   The limit condition on Kj does not seem very relevant
C   in this program.
C
      KJ = BIGM
C
C   Now choose values of V in accordance with the conditions in
C   section 3.3
C
      DO 20 R = 0, KJ-1
   20   V(R) = 0
      V(KJ) = 1
C
      IF (KJ .LT. BIGM) THEN
C
        TERM = SUB (0, H(KJ))
C
        DO 30 R = KJ+1, BIGM-1
          V(R) = ARBIT()
C
C         Check the condition of section 3.3,
C         remembering that the H's have the opposite sign.
C
          TERM = SUB (TERM, MUL (H(R), V(R)))
   30   CONTINUE
C
C         Now V(BIGM) is anything but TERM
C
          V(BIGM) = ADD (NONZER, TERM)
C
        DO 40 R = BIGM+1, M-1
   40     V(R) = ARBIT()
C
      ELSE
C       This is the case KJ .GE. BIGM
C
        DO 50 R = KJ+1, M-1
   50     V(R) = ARBIT()
C
      ENDIF
C
C   Calculate the remaining V's using the recursion of section 2.3,
C   remembering that the B's have the opposite sign.
C
      DO 70 R = 0, MAXV-M
        TERM = 0
        DO 60 I = 0, M-1
          TERM = SUB (TERM, MUL (B(I), V(R+I)))
   60   CONTINUE
        V(R+M) = TERM
   70 CONTINUE
C
      RETURN
      END
C
C     ***** end of SUBROUTINE CALCV
      INTEGER FUNCTION CHARAC (QIN)
C
C   This version :  12 December 1991
C
C   This function gives the characteristic for a field of
C   order QIN.  If no such field exists, or if QIN is out of
C   the range we can handle, returns 0.
C
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
 
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C
C
C
      INTEGER QIN, CH(MAXQ)
      SAVE CH
C
      DATA CH  / 0,  2,  3,  2,  5,  0,  7,  2,  3,  0,
     1          11,  0, 13,  0,  0,  2, 17,  0, 19,  0,
     2           0,  0, 23,  0,  5,  0,  3,  0, 29,  0,
     3          31,  2,  0,  0,  0,  0, 37,  0,  0,  0,
     4          41,  0, 43,  0,  0,  0, 47,  0,  7,  0/
C
      IF (QIN .LE. 1 .OR. QIN .GT. MAXQ) THEN
        CHARAC = 0
      ELSE
        CHARAC = CH(QIN)
      ENDIF
C
      END
C
C     ***** end of INTEGER FUNCTION CHARAC
      SUBROUTINE SINLO2 (DIM,SKIP,IFLAG)
C
C   This Modified Routine of INLO2
C
C
C   This subroutine calculates the values of Niederreiter's
C   C(I,J,R) and Various Scrambled Niederreiter's C(I,J,R)
C   performs other initialisation necessary
C   before calling GOLO2.
C
C INPUT :
C   DIMEN - The dimension of the sequence to be generated.
C        {DIMEN is called DIM in the argument of INLO2,
C        because DIMEN is subsequently passed via COMMON
C        and is called DIMEN there.}
C
C   SKIP  - The number of values to throw away at the beginning
C           of the sequence.
C
C   IFLAG -User choice of Genertor matrices.
C
C
C OUTPUT :
C   To SGOLO2, labelled common /COMM2/.
C
C USES :
C   Calls SCALCC2 to calculate the values of SCJ.
C
C
C   ------------------------------------------------------------
C
C
C   This file defines the common block /COMM2/ and some
C   associated parameters.  These are for use in the base 2
C   version of the generator.
C
      INTEGER MAXDIM, NBITS, EXOR
      PARAMETER (MAXDIM=318, NBITS=31)
C
C   The parameter MAXDIM is the maximum dimension that will be used.
C   NBITS is the number of bits (not counting the sign) in a
C   fixed-point integer.
C
      INTEGER SCJ(MAXDIM, 0:NBITS - 1), DIMEN, COUNT
      INTEGER NEXTQ(MAXDIM),SHIFT(MAXDIM)
      COMMON /COMM2/ SCJ, DIMEN, COUNT, NEXTQ
      SAVE   /COMM2/
      
C
C   The common block /COMM2/ :
C     SCJ    - Contains the packed values of User choice of 
C                  Niederreiter's C(I,J,R)
C     DIMEN   - The dimension of the sequence to be generated
C     COUNT - The index of the current item in the sequence,
C             expressed as an array of bits.  COUNT(R) is the same
C             as Niederreiter's AR(N) (page 54) except that N is
C             implicit.
C     NEXTQ - The numerators of the next item in the series.  These
C             are like Niederreiter's XI(N) (page 54) except that
C             N is implicit, and the NEXTQ are integers.  To obtain
C             the values of XI(N), multiply by RECIP (see GOLO2).
C
C   Array SCJ of the common block is set up by subroutine SCALCC2.
C   The other items in the common block are set up by INLO2.
C
C   ------------------------------------------------------------
C
C
C
      INTEGER I, R, DIM, SKIP, GRAY
C
      DIMEN = DIM
C
C       This assignment just relabels the variable for
C       subsequent use.
C
C        IF ((DIMEN .LE. 0) .OR. 
C    *       (DIMEN .GT. MAXDIM)) THEN
C          WRITE (*,*) ' SINLO2 :  Bad dimension'
C          STOP
C       ENDIF
C
       CALL SCALCC2(IFLAG)
C
C   Translate SKIP into Gray code
C
      GRAY = EXOR (SKIP, SKIP/2)
C
C   Now set up NEXTQ appropriately for this value of the Gray code
C
C
      R = 0
   10 IF (GRAY .NE. 0) THEN
        IF (MOD(GRAY,2) .NE. 0) THEN
          DO 20 I = 1, DIMEN
            NEXTQ(I) = EXOR(NEXTQ(I), SCJ(I,R))
   20     CONTINUE
        ENDIF
        GRAY = GRAY / 2
        R = R + 1
        GOTO 10
      ENDIF
C
      COUNT = SKIP
      RETURN
      END
C
C     *****  end of SUBROUTINE SINLO2
      SUBROUTINE PLYMUL (PA, PB, PC)
C
C   This version :  12 December 1991
C
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
 
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C
C
C
      INTEGER I, J, DEGA, DEGB, DEGC, TERM
      INTEGER PA(-1:MAXDEG), PB(-1:MAXDEG), PC(-1:MAXDEG)
      INTEGER PT(-1:MAXDEG)
C
C   Multiplies polynomial PA by PB putting the result in PC.
C   Coefficients are elements of the field of order Q.
C
      DEGA = PA(DEG)
      DEGB = PB(DEG)
      IF (DEGA .EQ. -1 .OR. DEGB .EQ. -1) THEN
        DEGC = -1
      ELSE
        DEGC = DEGA + DEGB
      ENDIF
      IF (DEGC .GT. MAXDEG) THEN
        WRITE (*,*) ' PLYMUL :  Degree of product exceeds MAXDEG'
        STOP
      ENDIF
C
      DO 20 I = 0, DEGC
        TERM = 0
        DO 10 J = MAX(0, I-DEGA), MIN(DEGB, I)
   10     TERM = ADD(TERM, MUL(PA(I-J), PB(J)))
   20   PT(I) = TERM
C
      PC(DEG) = DEGC
      DO 30 I = 0, DEGC
   30   PC(I) = PT(I)
      DO 40 I = DEGC+1, MAXDEG
   40   PC(I) = 0
      RETURN
      END
C
C     *****   end of SUBROUTINE PLYMUL

      SUBROUTINE SGOLO2 (QUASI)
C
C   This modified routine of GOLO2
C
C
C This subroutine generates a new quasi-random vector
C on each call.
C
C INPUT
C   From SINLO2, labelled common /COMM2/, properly initialized.
C
C OUTPUT
C   To the caller, the next vector in the sequence in the
C   array QUASI.
C
C   ------------------------------------------------------------
C
C
C   This file defines the common block /COMM2/ and some
C   associated parameters.  These are for use in the base 2
C   version of the generator.
C
      INTEGER MAXDIM, NBITS, EXOR
      PARAMETER (MAXDIM=318, NBITS=31)
C
C   The parameter MAXDIM is the maximum dimension that will be used.
C   NBITS is the number of bits (not counting the sign) in a
C   fixed-point integer.
C
      INTEGER SCJ(MAXDIM, 0:NBITS-1), DIMEN, COUNT, NEXTQ(MAXDIM)
      COMMON /COMM2/ SCJ, DIMEN, COUNT, NEXTQ
      SAVE   /COMM2/
C
C   The common block /COMM2/ :
C     SCJ    - Contains the packed values of Niederreiter's C(I,J,R)
C     DIMEN   - The dimension of the sequence to be generated
C     COUNT - The index of the current item in the sequence,
C             expressed as an array of bits.  COUNT(R) is the same
C             as Niederreiter's AR(N) (page 54) except that N is
C             implicit.
C     NEXTQ - The numerators of the next item in the series.  These
C             are like Niederreiter's XI(N) (page 54) except that
C             N is implicit, and the NEXTQ are integers.  To obtain
C             the values of XI(N), multiply by RECIP (see GOLO2).
C
C   Array SCJ of the common block is set up by subroutine SCALCC2.
C   The other items in the common block are set up by INLO2.
C
C   ------------------------------------------------------------
C
C
C
      REAL*8 RECIP
      PARAMETER (RECIP=2.0**(-NBITS))
C
C   The parameter RECIP is the multiplier which changes the
C   integers in NEXTQ into the required real values in QUASI.
C
      INTEGER I, R, MAX, TEMP(11,50), J
      REAL*8    QUASI(*),NNEWQ
C
C Multiply the numerators in NEXTQ by RECIP to get the next
C quasi-random vector
C
      DO 5 I = 1, DIMEN
        QUASI(I) = NEXTQ(I) * RECIP
    5 CONTINUE
C
C Find the position of the right-hand zero in COUNT.  This
C is the bit that changes in the Gray-code representation as
C we go from COUNT to COUNT+1.
C
      R = 0
      I = COUNT
   10 IF (MOD(I,2).NE.0) THEN
        R = R + 1
        I = I/2
        GOTO 10
      ENDIF
C
C Check that we have not passed 2**NBITS calls on GOLO2
C
      IF (R .GE. NBITS) THEN
        WRITE (*,*) ' SGOLO2 :  Too many calls'
        STOP
      ENDIF
C
C Compute the new numerators in vector NEXTQ
C
      DO 20 I = 1, DIMEN
       NEXTQ(I) = EXOR(NEXTQ(I), SCJ(I,R))
C       QUASI(I) = NEXTQ(I)*RECIP
   20 CONTINUE
C
      COUNT = COUNT + 1
      RETURN
      END
C
C     *****   end of PROCEDURE GOLO2

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

      REAL*8 FUNCTION TESTF (N, DIMEN, QUASI)
      INTEGER I, N, DIMEN
      REAL*8 X, EXACTF, QUASI(*)
C
C This version :  4 Mar 1992
C
C Provides a variety of test integrals for quasi-random
C sequences.  A call on TESTF computes an estimate of the
C integral ;  a call on EXACTF computes the exact value.
C
      GOTO (100, 200, 300, 400) N
C
      ENTRY EXACTF (N, DIMEN)
      GOTO (1100, 1200, 1300, 1400) N
C
C Test integral 1
C
  100 TESTF = 1.0
      DO 110 I = 1, DIMEN
        TESTF = TESTF * ABS(4 * QUASI(I) - 2)
  110 CONTINUE
      RETURN
C
 1100 EXACTF = 1.0
      RETURN
C
C Test integral 2
C
  200 TESTF = 1.0
      DO 210 I = 1, DIMEN
        TESTF = TESTF * I * COS(I * QUASI(I))
  210 CONTINUE
      RETURN
C
 1200 EXACTF = 1.0
      DO 1210 I = 1, DIMEN
        EXACTF = EXACTF * SIN(FLOAT(I))
 1210 CONTINUE
      RETURN
C
C Test integral 3
C
  300 TESTF = 1.0
      DO 350 I = 1, DIMEN
        X = 2 * QUASI(I) - 1
        GOTO (310, 320, 330, 340) MOD(I, 4)
  310   TESTF = TESTF * X
        GOTO 350
  320   TESTF = TESTF * (2*X*X - 1)
        GOTO 350
  330   TESTF = TESTF * (4*X*X - 3) * X
        GOTO 350
  340   X = X * X
        TESTF = TESTF * (8*X*X - 8*X + 1)
  350 CONTINUE
      RETURN
C
 1300 EXACTF = 0.0
      RETURN
C
C Test integral 4
C
  400        TESTF = 0
             X = 1
             DO 410   I = 1, DIMEN
                  X = - X * QUASI(I)
                  TESTF = TESTF + X
  410        CONTINUE
             RETURN
C
C
 1400 X = 1.0 / (2 ** (DIMEN ))
      IF (MOD(DIMEN, 2) .EQ. 0) THEN
        EXACTF = (X - 1) / 3
      ELSE
        EXACTF = (X + 1) / 3
      ENDIF
      RETURN
C
      END
      SUBROUTINE SETFLD (QIN)
      INTEGER QIN
C
C   This version : 12 December 1991
C
C   This subroutine sets up addition, multiplication, and
C   subtraction tables for the finite field of order QIN.
C   If necessary, it reads precalculated tables from the file
C   'gftabs.dat' using unit 1.  These precalculated tables
C   are supposed to have been put there by GFARIT.
C
C      *****  For the base-2 programs, these precalculated
C      *****  tables are not needed and, therefore, neither
C      *****  is GFARIT.
C
C
C   Unit 1 is closed both before and after the call of SETFLD.
C
C USES
C   Integer function CHARAC gets the characteristic of a field.
C
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
 
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C
C
C
      INTEGER I, J, N, CHARAC
C
      IF (QIN .LE. 1 .OR. QIN .GT. MAXQ) THEN
        WRITE (*,*) ' SETFLD :  Bad value of Q'
        STOP
      ENDIF
C
      Q = QIN
      P = CHARAC(Q)
C
      IF (P .EQ. 0) THEN
        WRITE (*,*) ' SETFLD :  There is no field of order', Q
        STOP
      ENDIF
C
C Set up to handle a field of prime order :  calculate ADD and MUL.
C
      IF (P .EQ. Q) THEN
        DO 10 I = 0, Q-1
          DO 10 J = 0, Q-1
            ADD(I,J) = MOD(I+J, P)
            MUL(I,J) = MOD(I*J, P)
   10   CONTINUE
C
C Set up to handle a field of prime-power order :  tables for
C ADD and MUL are in the file 'gftabs.dat'.
C
      ELSE
        OPEN (UNIT=1, FILE='gftabs.dat', STATUS='old')
C
C    *****  OPEN statements are system dependent.
C
   20   READ (1, 900, END=500) N
  900   FORMAT (20I3)
        DO 30 I = 0, N-1
          READ (1, 900) (ADD(I,J), J = 0, N-1)
   30   CONTINUE
        DO 40 I = 0, N-1
          READ (1, 900) (MUL(I,J), J = 0, N-1)
   40   CONTINUE
        IF (N .NE. Q) GOTO 20
        CLOSE (1)
      ENDIF
C
C Now use the addition table to set the subtraction table.
C
      DO 60 I = 0, Q-1
        DO 50 J = 0, Q-1
          SUB(ADD(I,J), I) = J
   50   CONTINUE
   60 CONTINUE
      RETURN
C
  500 WRITE (*,*) ' SETFLD :  Tables for q =', Q, ' not found'
      STOP
C
      END
C
C     ***** end of SUBROUTINE SETFLD









