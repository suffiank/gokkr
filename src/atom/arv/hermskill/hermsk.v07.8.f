C     HARTREE-FOCK-SLATER SELF-CONSISTENT ATOMIC FIELD PROGRAM
C      ORIGINALLY  WRITTEN BY SHERWOOD SKILLMAN
C      RCA LABORATORIES, PRINCETON, NEW JERSEY, SPRING 1961
C      MODIFIED BY FRANK HERMAN, SUMMER 1961
C      FURTHER MODIFIED BY RICHARD KORTUM,  LOCKHEED RESEARCH
C      LABORATORIES, PALO ALTO, CALIFORNIA,  SUMMER 1962
C      Modified by L F MATTHEISS for CCL Monitor System
C      Modified by S.CRAMPIN and D.K.SALDIN for MUFPOT use
C
C      Rewritten by M. Pauli into F77 code with comments, Nov. 2000
C        Cleaned up the program; determined what was original vs. 
C          modified by prior people.  Fixed some of the prior modifies.
C        Commented all modifications.  Added extra comments to the
C          the original program.
C        Made a few corrections to Herman-Skillman original code.
C
C      MODIFICATIONS:
C        1) Desired level of output is controllable
C           a) CONT = continue run; output hsall.out with headers
C               BRF = continue run; omit hsall.out (creates empty file)
C           b)  POT = print out potential in hspot.out
C           c)  RAD = print out radial wf in hswf.out
C           d)  hsinfo.out always output with energy eigenvalues and 
C               convergence information
C        2) Descriptive headers in output files
C        3) Added output to charge densities
C        4) Initialization of variables
C        5) Ability to set the ionic radius of an ionized atom and 
C           give the branching ratio for charge contained within the
C           ionic radius.
C        6) Fix for when KUT = 1 to fill in the theoretical limiting
C           value of the unmodified potential for MESH > 441.
C        7) Theory improvement adding ALPHA factor for exchange
C           potential.
C        8) Numerical improvement changing simple trapezoidal integrals
C           for atomic coulomb potential into Simpson's rule integral
C           exact up to polynomials of degree 3.
C        9) Various fixes for divide by zero errors
C       10) Fixed the output of the potential to distinguish between
C           choice of modified potential output or unmodified potential
C           output.  NOTE:  Book is incorrect in claiming it gives Cu+1
C           in a modified potential output on pp. 7-7, 7-9, 7-10.
C ITERATION NUMBER,MEASURE OF SELF-CONSISTENCY,AND ATOMIC NUMBER Z
C      ARE printed to hsinfo.out
C
C REMOVED the COMMON and EQUIVALENCE statements and modified the
C          function calls to properly pass the necessary variables.
C MODIFIED SNLO, RU, XI, XJ, RUFNL2 to have a max of 541 values due to 
C          a rare "gotcha" with evaluating those variables to dimension
C          ICUT+20 inside the Pratt Improvement Scheme where ICUT can
C           be from 0 to MESH
C ----------------------------------------------------------------------
      DIMENSION R(521), X(521), SNLO(541), SNL(24,521),
     &          RSCORE(521), RSVALE(521), RSATOM(521),
     &          RU2(521), RU3(521), RU(541), RUEXCH(521), XI(541),
     &          XJ(541), RUINL1(521), RUFNL1(521), RUINL2(521),
     &          RUFNL2(541), V(521), RNUM(521), DENM(521), 
     &          NKKK(24), NNLZ(24), WWNL(24), EE(24), A(4,5)

C ----------------------------------------------------------------------
C Open Files and Set Output Level
C
C c DETAIL = T/F to print combined output with headers
C c TITLE = string to read in header title for input file
C c Q001,2,3,4 = control strings choices for comparison with REC1,2,3
C c REC1 = 'CONT' : continue with run and output all to header file (9)
C          or ' BRF' : continue with run and skip output to header file
C          or anything else : quit program
C c REC2 = ' POT' : output potential to file (10)
C          or anything else : skip potential output
C c REC3 = ' RAD' : output all radial wave functions to file (11)
C          or anything else : skip radial w.f. output
C
C ADDED variables to control desired output level
C ADDED RSATOR to help with new atomic coulomb potential integral
C       method
C ADDED DETAIL to control whether to print all results into a
C       single file that includes descriptive headers (8 hsall)
C ADDED modern opening of files.
C ----------------------------------------------------------------------
      DIMENSION RSATOR(521)
      LOGICAL DETAIL
      CHARACTER*8 TITLE(9)
      CHARACTER*4 Q001, Q002, Q003, Q004, REC1, REC2, REC3
      DATA Q001, Q002, Q003, Q004/'CONT', ' BRF', ' POT', ' RAD'/

      OPEN ( 5, FILE='input.dat',          STATUS='OLD')
      OPEN ( 8, FILE='hsinfo.out',         STATUS='NEW')
 7000 READ (5,1000) REC1, REC2, REC3
 1000 FORMAT (3(A4,4X))
      IF ((REC1 .NE. Q001) .AND. (REC1 .NE. Q002)) STOP

      DETAIL = .TRUE.
      IF (Q002 .EQ. REC1) DETAIL = .FALSE.
      IF (DETAIL) THEN
        OPEN ( 9, FILE='hsall.out',          STATUS='NEW')
      END IF
      IF (Q003 .EQ. REC2) THEN
        OPEN (10, FILE='hspot.out',          STATUS='NEW')
      END IF
      IF (Q004 .EQ. REC3) THEN
        OPEN (11, FILE='hswf.out',           STATUS='NEW')
      END IF


C ----------------------------------------------------------------------
C Initialization of Variables
C   Starting point for each new run defined in the input file
C   Potentials are in Rydbergs... rV_coul(r) = -Z*(e^2) = -Z*(2)
C
C d R = R mesh; actual radial mesh
C d X = X mesh; normalized radial mesh grid
C d SNLO = normalized radial wave function for the orbital being solved
C          by SCHEQ
C d SNL = normalized radial wave functions for all of the possible
C         orbitals rR_nl
C d RSCORE = normalized core charge density 4pi*|rR|^2_core
C d RSVALE = normalized valence charge density 4pi*|rR|^2_valence
C d RSATOM = normalized total charge density = RSCORE + RSVALE 
C d RU2 = temp. atomic potential rV(r) used for reading input file
C d RU3 = temp. atomic potential rV(r) used for reading input file
C         and to extrapolate to trial RU
C d RU = trial (pre-SCHEQ) full 441 pt. atomic potential rV(r)
C d RUEXCH = atomic exchange potential rV_exch(r)
C d XI(I) = temp variable.  Used to calculate atomic potential.
C d XJ(I) = temp variable.  Used to calculate atomic potential.  Ends
C           up becoming the post-SCHEQ potential rV(r).
C d RUINL1 = store prior iteration pre-SCHEQ potential rV(r)
C d RUFNL1 = store prior iteration post-SCHEQ potential rV(r)
C d RUINL2 = store current iteration pre-SCHEQ potential rV(r)
C d RUFNL2 = store current iteration post-SCHEQ potential rV(r)
C d V = full 441,481,521 pt. atomic potential V(r) determined from RU
C       Used by SCHEQ.
C d RNUM = temp. variable.  Numerator in Pratt improvement ALPH calc.
C d DENM = temp. variable.  Denominat. in Pratt improvement ALPH calc.
C i NKKK = max radial mesh pt of normalized radial wave function and 
C          charge density distribution for all of the possible 
C          orbitals nl - by which, the wave func./density is zero
C i NNLZ = 100n + 10l + c; c = dummy number used to distinguish
C          different run configurations - default is 0
C d WWNL = occupation number (# of electrons); 2(2l+1) for filled shells
C d EE = trial energy eigenvalues for the orbital shell nl: Rydbergs
C d A = temp variable used to hold perturbation terms for creating
C       new trial energy eigenvalues
C
C ADDED initialization of variables
C ----------------------------------------------------------------------
      DO I = 1,521
        R(I) = 0.0
        X(I) = 0.0
        SNLO(I) = 0.0
        DO J = 1,24
          SNL(J,I) = 0.0
        END DO
        RSCORE(I) = 0.0
        RSVALE(I) = 0.0
        RSATOM(I) = 0.0
        RU2(I) = 0.0
        RU3(I) = 0.0
        RU(I) = 0.0
        RUEXCH(I) = 0.0
        XI(I) = 0.0
        XJ(I) = 0.0
        RUINL1(I) = 0.0
        RUFNL1(I) = 0.0
        RUINL2(I) = 0.0
        RUFNL2(I) = 0.0
        V(I) = 0.0
        RNUM(I) = 0.0
        DENM(I) = 0.0
      END DO
      DO I = 522, 541
        SNLO(I) = 0.0
        RU(I) = 0.0
        XI(I) = 0.0
        XJ(I) = 0.0
        RUFNL2(I) = 0.0
      END DO
      DO I = 1,24
        NKKK(I) = 0
        NNLZ(I) = 0
        WWNL(I) = 0.0
        EE(I) = 0.0
      END DO
      DO I = 1,4
        DO J = 1,5
          A(I,J) = 0.0
        END DO
      END DO




C ----------------------------------------------------------------------
C Start of Original Herman-Skillman program 
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C Read Heading
C   Set up wavefunction file(11) formatted for use with
C   the mufpot Muffin-Tin program.  Give the potential file(10) a short
C   explanatory heading.  Give the info file(8) a title.
C ----------------------------------------------------------------------
      READ (5,1001) TITLE
      WRITE ( 8,1002) TITLE
 1001 FORMAT (10A8)
 1002 FORMAT (1X, 9A8///)

C ----------------------------------------------------------------------
C Read Control Values.
C
C i KEY = 0,1,2 : type of potential being inputted... see below
C d TOL = Self-consistent criterion tolerance; 0.001 is good
C d THRESH = energy eigenvalue accuracy criterion = abs(dE/E);
C            0.00001 or 0.00002 is good
C i MESH = # of mesh points used for Schroed. Eqtn solver
C          441, 481, or 521.  441 adequate for most atoms.  481 or 521
C          may be need for highly excited configurations to allow
C          larger radial distances for the outermost orbitals.
C i IPRATT = # of consecutive iterations to use the Pratt Improvement
C            scheme, to generate a new trial potential, following each
C            application of the arithmetic average scheme in an 
C            iteration.  1 is usually sufficient
C i MAXIT = max # of iterations allowed before fail Self-consistency
C           convergence; program insures 20 iter.; 30 to 50 is good
C i KUT = choice of H-F-S potential to use.  Controls how the potential
C         tails off at large radius.  0 or 1.
C         0 = V(r); H-F-S-Latter to enforce ideal asymptotic behavior,
C           fixing the tendency of H-F-S to go to zero too fast.
C         1 = V0(r); H-F-S with the free-electron exchange potential
C         V(r) = V0(r) for r < r0; where V0(r0) = -2(Z-N+1)/r0
C         V(r) = -2(Z-N+1)/r for r >= r0
C d RADION = ionic radius.  If given, it defines the radius of the
C            charged atomic sphere used in calculations of the coulomb
C            potential.  Default is 0.0 to ignore this option.
C d BRATIO = branching ratio for ionic radius; 0.0 < RATIO <= 1.0.
C            0.0 defaults to 1.0
C d ALPHA = Xa local-statistical-exchange approximation parameter.
C           ALPHA is about 0.70 for all but the lightest atoms where it
C           rises to about 0.78
C ----------------------------------------------------------------------
      READ (5,1004) KEY, TOL, THRESH, MESH, IPRATT, MAXIT, KUT,
     &              RADION, BRATIO, ALPHA

C ----------------------------------------------------------------------
C Write Control Values.
C
C ADDED output of the control values into the info file(9).
C ADDED RADION, BRATIO, and ALPHA for expanded ability with potential.
C ----------------------------------------------------------------------
      WRITE ( 8,1005) KEY, TOL, THRESH, MESH, IPRATT, MAXIT, KUT,
     &                RADION, BRATIO, ALPHA
 1004 FORMAT (I4, 2F8.6, 4I4, 3F14.9)
 1005 FORMAT ('      KEY   TOL     THRESH   MESH  IPRATT  MAXIT  KUT',
     &        '   RADION     RATIO      ALPHA'/
     &        '      ', I2, F9.4, F10.6, I5, 2I7, I6, F12.7, 2F11.7/)
      IF (RADION .GT. 0.0) THEN
C       Given an ionic radius, make sure we have a valid branching
C       ratio to go with it.
        IF (BRATIO .LT. 0.0) THEN
          WRITE ( 8,*)'RATIO INCORRECT'
          STOP
        ELSE IF (0.0 .EQ. BRATIO) THEN
          BRATIO = 1.0
        END IF
      END IF
      IF (MAXIT .LT. 20) MAXIT = 20

C ----------------------------------------------------------------------
C Read in Atomic Potential.
C   Form based on if KEY = 0,1,2
C
C d ZE3,2 = scaling of full potential; should, ideally, equal atomic #
C ----------------------------------------------------------------------
      IF (0 .EQ. KEY) THEN
C       110pt. Normalized Potential, from point 1 to point 437 with
C       every 4th point only.  It will be extrapolated to a full 441pt.
C       U(r) = -rV(r) / 2Z:  U(0) = 1.00, U(infin) = 1 / Z
        READ (5,1006) (RU2(M), M = 1,437,4)
 1006   FORMAT (F8.5, 9F7.5)
      ELSE IF (1 .EQ. KEY) THEN
C       441pt. Full Potential, possibly unnormalized
C       rV(r): rV(0) = RU3(1) = -2Z, rV(infin) = -2
        READ (5,1007) (RU3(M), M = 1,441)
 1007   FORMAT (1PE15.7, 1P4E14.7)
        ZE3 = -RU3(1) / 2.0
      ELSE IF (2 .EQ. KEY) THEN
C       Two 441pt. Full Potentials near the desired potential.
C       They will be used to extrapolate the desired rV(r) expected.
C       2Z_3 - Z_2 = Z is expected; i.e. Z_3 = Z-1 and Z_2 = Z-2 is good
        READ (5,1007) (RU2(M), M = 1,441)
        READ (5,1007) (RU3(M), M = 1,441)
        ZE2 = -RU2(1) / 2.0
        ZE3 = -RU3(1) / 2.0
      ELSE
        WRITE (*,*)'BAD KEY FOR CONTROLLING INPUT POTENTIAL FORMAT'
        STOP
      END IF

C ----------------------------------------------------------------------
C Read the Atomic Information
C
C d Z = atomic number
C i NCORES = # of atomic shells assumed to be cores = first set of info
C            lines read in
C i NVALES = # of atomic shells assumed to be valences = second set of
C            info lines read in
C d XION = ionicity; i.e. 0 for neutrals, +1, +2, etc. (= Z - # of elec)
C          XION must be 0 or positive.  No method of making the
C          potential energy more negative to bind excess electrons
C          is included.  Allows fractional ionicities.
C i IZ = Z; dummy atomic # as integer for output write statements
C i NCSPVS = Number of CoreS Plus ValenceS
C d TWOION = 2(XION)
C d TWOZZZ = 2(XION + 1)
C d WWW = total occupation, total # of electrons
C d XIONSAV = stored value of XION for use with ionic branching ratio
C
C ADDED output for checks.  Output the energy eigenvalues for the 
C       initial trial potential that is inputted.  
C ADDED XIONSAV for value security with the ionic branching ratio
C ----------------------------------------------------------------------
      READ (5,1008) Z, NCORES, NVALES, XION
 1008 FORMAT (F4.0, 2I4, F10.6)
      IZ = Z
      NCSPVS = NCORES + NVALES
      TWOION = 2.0 * XION
      TWOZZZ = 2.0 * (XION + 1.0)
      READ (5,1010) (NNLZ(I), WWNL(I), EE(I), I = 1,NCSPVS)
      XIONSAV = XION
      WRITE ( 8,1011) (NNLZ(I), WWNL(I), EE(I), I = 1,NCSPVS)
 1010 FORMAT (I4, F10.6, F11.4)
 1011 FORMAT (6x, I4, F10.6, F14.3)
      WWW = 0.0
      DO I = 1,NCSPVS
        WWW = WWW + WWNL(I)
      END DO
      IF (ABS(Z-WWW-XION) .GE. 0.001) THEN
        WRITE (*,*)'TOTAL OCCUPANCY PLUS IONICITY NOT = Z.'
        WRITE ( 8,1012) WWW, XION, Z, NCORES, NVALES, NCSPVS
 1012   FORMAT (' WWW= ', F4.0, ' XION= ', F4.0, '   Z= ', F4.0,
     &        '   NCORES= ', I4, '   NVALES= ', I4, '   NCSPVS= ',
     &        I4/,' CONTROL VALUES INCORRECT.')
        STOP
      END IF

C ----------------------------------------------------------------------
C Construct X Mesh and R Mesh
C.  I = 1..441; R = CMU * X
C
C d CMU = 0.5 * (3pi/4)^2/3 * Z^-1/3
C         also appears in energy eigenvalue perturbation and SCHEQ call
C i NBLOCK = # of blocks of potential mesh.  Each block has a new DELTAX
C d DELTAX = step size of X mesh grid
C ----------------------------------------------------------------------
      NBLOCK = (MESH) / 40
      X(1) = 0.0
      R(1) = 0.0
      I = 1
      DELTAX = 0.0025
      CMU = 0.88534138 / Z**(1.0/3.0)
      DO J = 1,NBLOCK
        DO K = 1,40
          I = I + 1
          X(I) = X(I-1) + DELTAX
          R(I) = CMU * X(I)
        END DO
        DELTAX = 2.0 * DELTAX
      END DO

C ----------------------------------------------------------------------
C Construct the Initial Trial Atomic Potential
C   Take the input data and create the trial atomic potential.  If the
C   trial mesh is chosen to go beyond 441 points or the Latter tail 
C   correction is enacted, define the LIMIT point where the potential
C   is first set to its theoretical value.
C
C d TWOZ, ZOZ = dummy variables used for computational savings
C i LIMIT = mesh point where limiting value of the potential, whichever
C           is chosen, is first used to fill in data in place of input
C i ICUT = mesh point where modified Latter tail correction to the H-F-S
C          potential takes over; i.e. r0 -> r0V(r0) = -2(Z-N + 1)/r0
C ----------------------------------------------------------------------
C     First, create rV(r) from the given input
      IF (0 .EQ. KEY) THEN
        TWOZ = Z + Z
        DO I = 1,437,4
          RU(I) = -RU2(I) * TWOZ
        END DO
        RU(441) = RU(437)
        RU(445) = RU(437)
        M = 9
        DO I = 1,437,4
          M = M - 1
          IF (M .GE. 0) THEN 
C           i mod 36 not = 0  :: so, still within DELTAX block
            RU(I+1) = (21.0*RU(I) + 14.0*RU(I+4) - 3.0*RU(I+8)) / 32.0
            RU(I+2) = ( 3.0*RU(I) +  6.0*RU(I+4) -     RU(I+8)) /  8.0
            RU(I+3) = ( 5.0*RU(I) + 30.0*RU(I+4) - 3.0*RU(I+8)) / 32.0
          ELSE
C           i mod 36 = 0  :: so, set up the 38th, 39th, 40th values
C           This is the end of the values for this DELTAX.
            RU(I+1) = (22.0*RU(I) + 11.0*RU(I+4) - RU(I+8)) / 32.0
            RU(I+2) = (10.0*RU(I) + 15.0*RU(I+4) - RU(I+8)) / 24.0
            RU(I+3) = ( 6.0*RU(I) + 27.0*RU(I+4) - RU(I+8)) / 32.0
            M = 9
          END IF
        END DO
      ELSE IF (1 .EQ. KEY) THEN
        IF (ABS(ZE3-Z) .LE. 0.001) THEN
C         ZE3 from potential value = Z declared; pot. has correct form.
          DO I = 1,441
            RU(I) = RU3(I)
          END DO
        ELSE
C         ZE3 not = Z.  Rescale the input potential for correct form.
          ZOZ = Z / ZE3
          DO I = 1,441
            RU(I) = RU3(I) * ZOZ
          END DO
        END IF
      ELSE
C       KEY = 2 must be true.  Prior if KEY structure insures 0,1,2.
        IF (ABS((2.0*ZE3-ZE2) - Z) .LE. 0.001) THEN 
          DO I = 1,441
            RU(I) = RU3(I) + RU3(I) - RU2(I)
          END DO
        ELSE
          WRITE ( 8,*)'ERROR EXTRAPOLATING POTENTIAL rV(r) FROM GIVEN'
          WRITE ( 8,*)'STARTING POTENTIALS TO GET DESIRED Z = ',Z
          WRITE ( 8,*)'INPUTS GIVE Z = 2 ', ZE3, ' - ', ZE2
          STOP
        END IF
      END IF

C     Second, calculate V(r) from rV(r); pad out to desired mesh size
      V(1) = -9.9E35
C     M should always = 441 as MESH = {441,481,521}.  Maybe this is
C     a left over from possible smaller meshes.
C     MIN0 instead of general MIN to enforce int vs. int comparison
      M = MIN0(441,MESH)
      IF (KUT .NE. 0) THEN
C       KUT = 1 :: Use ideal H-F-S potential
C                  V(r) = -2(Z-N)/r for r > r0
        DO I = 2,M
          V(I) = RU(I) / R(I)
        END DO
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C CORRECTED from Herman and Skillman code which appeared to be backwards
        IF (MESH .GT. M) THEN
C         Requested mesh is greater than input 441 pts.
C         Fill out extra mesh using potential theoretical limit values
          DO I = 442,MESH
            V(I) = -TWOION / R(I)
          END DO
        END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LIMIT = M
        ICUT = MESH
      ELSE
C       KUT = 0 :: Use modified H-F-S + Latter correction potential
C                  V(r) = -2(Z-N + 1)/r for r > r0
        ICUT = 0
        DO I = 2,M
          IF (ICUT .LE. 0) THEN
            IF ((TWOZZZ + RU(I)) .GT. 0.0) THEN
C             Found cut off value r0, start Latter tail correction
              ICUT = I
              V(I) = -TWOZZZ / R(I)
            ELSE
              V(I) = RU(I) / R(I)
            END IF
          ELSE
C           From here on out, use Latter tail correction
            V(I) = -TWOZZZ / R(I)
          END IF
        END DO
        IF (ICUT .LE. 0) ICUT = M
        LIMIT = ICUT
        IF (MESH .GT. M) THEN
C         Requested mesh is greater than input 441 pts.
C         Fill out extra mesh using potential theoretical limit values
          DO I = 442,MESH
            V(I) = -TWOZZZ / R(I)
          END DO
        END IF
      END IF





C ----------------------------------------------------------------------
C Main Iteration Scheme.
C   Runs until achieving self-consistency with the atomic potential
C   or possible failure - no convergence before given iteration limit.
C   DELTA <  TOL; end s-c loop and write out results
C   DELTA >= TOL; if NITER < MAXIT generate a new trial potential rV(r),
C                 a new V(r) for SCHEQ, and new energy eigenvalues
C
C d DELTA = calculated error for meeting self-consistency criterion
C d PDELTA = error from prior iteration.  Used to track the behavior
C            of the error and make sure it is improving (decreasing)
C            with each iteration.
C i NITER = number of self-consistent iterations completed
C i NONMON = the number of times delta is allowed to not monotonically
C            decrease (new delta larger than prior delta); tracked only 
C            by iterations trying to use the Pratt improvement scheme.
C            At the specified limit, Pratt scheme is no longer allowed.
C i IPRSW = Pratt improvment scheme switch.  Tracks whether the
C           current trial potential is to be generated by arithmetic
C           averaging or by Pratt improvement scheme.
C ----------------------------------------------------------------------
      DELTA = 1.0E06
      PDELTA = 1.0E06
      NITER = 0
      NONMON = 3
      IPRSW = 0

      DO WHILE (DELTA .GE. TOL)
C       Failed the self-consistency criterion; prepare a new iteration
C       and run it.


C ----------------------------------------------------------------------
C New Trial Values Generation Block
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C Check the Iteration Number
C   Iterations will alternate between methods to generate the new
C   trial potential between arithmetic averaging and the Pratt 
C   improvement scheme based upon how often Pratt is to be called
C   (IPRATT).
C   e.g. for IPRATT = 2, the new trial value iterations will be arith,
C        pratt, pratt, arith, pratt, pratt, arith, pratt, pratt, arith,
C        arith, arith, ... where at the end I am assuming NONMON has
C        occured and Pratt is then disallowed.
C   NOTE: Initial iteration skips this whole block and jumps straight
C         down to the Iteration Block so the first iteration will use
C         the user-supplied values.
C ----------------------------------------------------------------------
        IF (NITER .GT. 0) THEN

          IF (NITER .GT. MAXIT) THEN
            WRITE ( 8,*) 'Maximum iterations reached without achieving'
            WRITE ( 8,*) 'self consistency criterion tolerance.'
            WRITE ( 8,*) '  Results so far:'
            WRITE ( 8,1013)
 1013       FORMAT ('1'//, 8X, 'I', 6X, 'X', 11X, 'RUINL1', 10X,
     &             'RUFNL1', 10X, 'RUINL2', 10X, 'RUFNL2', 12X, 'RU',/)
            DO I = 1,MESH,5
              WRITE ( 8,1014) I, X(I), RU3(I), RUINL1(I), RUFNL1(I),
     &                        RUINL2(I), RUFNL2(I), RU(I)
 1014         FORMAT (1X, I8, F12.7, 1P6E16.7)
            END DO
            STOP
          END IF

          IF ((IPRSW .GT. 0) .AND. (NONMON .GT. 0)) THEN
C           We are in a possible Pratt scheme iteration.  First, check
C           the behaviour of the error.
            IF (PDELTA .LE. DELTA) NONMON = NONMON - 1
            IF (NONMON .GT. 0) THEN
C             Delta is still behaving well, continue with Pratt scheme.
C ----------------------------------------------------------------------
C Pratt Improvement Scheme for New Trial Potential
C   Fancy method of generating a new trial potential from the old
C   pre- and post-SCHEQ potentials.  It is implemented a limited
C   number of times equal to the value of IPRATT (and controlled by
C   IPRSW) in-between each implementation of the arithmetic averaging
C   scheme as long as the error (DELTA) continues to decrease monot'lly.
C   Due to method of alternating schemes, this can only be called after
C   the second iteration.
C
C d ALPH = Alpha Pratt
C d ADEL = temp. variable used for Latter tail correction tweak
C d ASUM = temp. variable
C ----------------------------------------------------------------------
              ALPH = 0.5
              DO I = 2,ICUT
                RNUM(I) = RUINL1(I) * RUFNL2(I) - RUINL2(I) * RUFNL1(I)
                DENM(I) = RUFNL2(I) - RUFNL1(I) - RUINL2(I) + RUINL1(I)
                IF (ABS(DENM(I)/RUINL2(I)) .LE. 0.0001) THEN
                  ALPH = 0.5
                ELSE 
                  ALPH = (RNUM(I)/DENM(I) - RUFNL2(I)) / SNLO(I)
                  IF (ALPH .LT. 0.0) THEN
                    ALPH = 0.0
                  ELSE IF (ALPH .GT. 0.5) THEN
                    ALPH = 0.5
                  END IF
                END IF
                XI(I) = ALPH
              END DO
              IF (0 .EQ. KUT) THEN
C               Tweak for the Latter tail correction mode.
C               Creates a 20 mesh point bridging zone (half a block)
C               with a simple linear decrease to zero.
C               NOTE: this loop can require an extra 20 elements of the
C                     affected arrays in the case when ICUT = MESH.
C                     Generally, ICUT = MESH will only happen when using
C                     Vo and not for Vlatter.
                ADEL = XI(ICUT) / 20.0
                DO I = ICUT+1,ICUT+20
                  XI(I) = XI(I-1) - ADEL
                  XJ(I) = XI(I)
                END DO
                DO I = ICUT+21,MESH
                  XJ(I) = 0.0
                  RU(I) = RUFNL2(I)
                END DO
              END IF
              XJ(1) = 0.5
              XJ(2) = XI(2)
              ASUM = XI(2) + XI(3) + XI(4) + XI(5)
              DO I = 3,ICUT
                XJ(I) = ASUM * 0.2
                ASUM = ASUM - XI(I-2) + XI(I+3)
              END DO
              DO I = 2,ICUT+20
                RU(I)= RUFNL2(I) + XJ(I)*SNLO(I)
              END DO
              IPRSW = IPRSW - 1
            ELSE
C             Delta made last allowed nonmontonic change, switch to the
C             arithmetic average scheme.  This is just a copy of the
C             arith. avg. below to complete the logic of NONMON IF-ELSE.
              DO I = 2,LIMIT
                RU(I) = 0.5 * (RU(I) + XJ(I))
              END DO 
              IF (MESH .GT. LIMIT) THEN
                IF (XJ(MESH) .NE. XJ(LIMIT)) THEN
                  RUSCALE = (XJ(MESH) - RU(LIMIT)) /
     &                      (XJ(MESH) - XJ(LIMIT))
                ELSE
                  RUSCALE = 0.0
                END IF
                DO I = LIMIT,MESH
                  RU(I) = XJ(MESH) - RUSCALE * (XJ(MESH) - XJ(I))
                END DO
                LIMIT = MESH
              END IF
              IPRSW = IPRATT
            END IF
          ELSE
C ----------------------------------------------------------------------
C Arithmetic Average for New Trial Potential
C   Take an average of the pre- and post-SCHEQ potentials as the new
C   trial potential for the main mesh region.  IPRSW is initialized to
C   0 so this will always be called after the first iteration.  It will
C   again be used each time after the desired # of Pratt improvement
C   iterations has been exhausted and used every time after the Pratt
C   scheme is discounted due to nonmonotonic behavior.
C 
C d RUSCALE = temp variable for computational efficiency
C             Used only after the first iteration to rescale the 
C             tail of the potential mesh that was initialized to
C             theory values. 
C
C ADDED fix for divide by zero TO RUSCALE
C ----------------------------------------------------------------------
            DO I = 2,LIMIT
              RU(I) = 0.5 * (RU(I) + XJ(I))
            END DO 
            IF (MESH .GT. LIMIT) THEN
C             For the first iteration, if the initial mesh was expanded
C             or used a Latter tail correction, use a linear adjustment
C             of the potential for mesh points beyond the limit point.
              IF (XJ(MESH) .NE. XJ(LIMIT)) THEN
                RUSCALE = (XJ(MESH) - RU(LIMIT)) /
     &                    (XJ(MESH) - XJ(LIMIT))
              ELSE
C               If XJ(MESH) = XJ(LIMIT), then we have reached the
C               constant tail of the potential
                RUSCALE = 0.0
              END IF
              DO I = LIMIT,MESH
                RU(I) = XJ(MESH) - RUSCALE * (XJ(MESH) - XJ(I))
              END DO
              LIMIT = MESH
            END IF
            IPRSW = IPRATT
          END IF

C ----------------------------------------------------------------------
C Calculate the Next Iteration V(r) for SCHEQ
C   Uses the new trial potential rV(r) just created
C   If using the Latter tail correction; once the theoretical limit
C     point has been reached (TWOZZZ = RU(I)), set the ICUT point and
C     fill in the rest of the potential with the theoretical value.
C
C d VLAST = temp variable; holds last iteration's V(I) for current R(I)
C ----------------------------------------------------------------------
          IF (KUT .NE. 0) THEN
C           Straight calculation
            ICUT = MESH
            LIMIT = MESH
            DO I = 2,MESH
              VLAST = V(I)
              V(I) = RU(I) / R(I)
              XI(I) = V(I) - VLAST
            END DO
          ELSE
C           Determine the point at which to apply the Latter tail.
            ICUT = 0
            DO I = 2,MESH
              VLAST = V(I)
              IF (ICUT .GT. 0) THEN
                V(I) = -TWOZZZ / R(I)
              ELSE
                IF ((TWOZZZ + RU(I)) .GT. 0.0) THEN
                  ICUT = I
                  V(I) = -TWOZZZ / R(I)
                ELSE
                  V(I) = RU(I) / R(I)
                END IF
              END IF
              XI(I) = V(I) - VLAST
            END DO
          END IF
          XI(1) = 0.0

C ----------------------------------------------------------------------
C Calculate the Next Iteration Trial Energy Eigenvalues
C   Uses perturbation theory
C
C d C = mu from generation of X mesh/R mesh
C i K = # of blocks within mesh of normalized radial wave function
C       that contain pertinent (non-zero) data.  Each block has a new
C       mesh grid step size (double previous block).
C d H = step size for trial energy eigenvalue perturbation
C d A1 = temp variable for trial energy eigenvalue perturbation
C d A2 = temp variable for trial energy eigenvalue perturbation
C
C ADDED output to list energy eigenvalues for the upcoming iteration
C ----------------------------------------------------------------------
          DO M = 1,NCSPVS
            K = (NKKK(M)-1) / 40
C           same C from x = mu r
            H = 0.0025 * CMU
            ASUM = 0.0
            A1 = 0.0
            I = 1
            DO J = 1,K
              DO L = 1,40
                I = I + 1
                A2 = XI(I) * SNL(M,I)**2
                A1 = A1 + A2 * H
              END DO
              ASUM = ASUM + A1 - H * (A2/2.0)
              H = H + H
              A1 = H * (A2/2.0)
            END DO
            EE(M) = EE(M) + ASUM
          END DO
C         NOTE: initial e. eigen. are output just after they are read in
          WRITE ( 8,1011) (NNLZ(I), WWNL(I), EE(I), I = 1,NCSPVS)

        END IF
C         Close the IF (NITER>0)
C ----------------------------------------------------------------------
C End of New Trial Values Generation Block
C ----------------------------------------------------------------------


C ----------------------------------------------------------------------
C Iteration Block
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C Solve the Schroedinger Equation
C   Loop through each orbital state, nl, in turn; from the core states
C   out to the valence states.  Calculate the core and valence charge
C   densities for the current solution.  Store the normalized radial
C   wave functions for each orbital.
C
C d E = energy eigenvalue (BE) of orbital state nl being solved
C i NN = principle quantum number n of orbital being solved derived from
C        NNLZ with truncation by division
C i LAM = angular momentum (orbital QM #) l of orbital being solved
C i KKK = max radial mesh pt for normalized radial wave function and
C         charge density distribution - by which wave func./charge 
C         density is definitely zero; value returned by SCHEQ
C
C MODIFIED the function call to properly pass the necessary variables.
C ----------------------------------------------------------------------
        DO I = 1,MESH
          RSCORE(I) = 0.0
          RSVALE(I) = 0.0
        END DO
        DO M = 1,NCSPVS
          E = EE(M)
          NN = NNLZ(M) / 100
          LAM = NNLZ(M) / 10 - 10 * NN
          CALL SCHEQ(Z, E, LAM, NN, KKK, MESH, CMU, THRESH, V, R, SNLO)
          IF (M .LE. NCORES) THEN
            DO I = 1,KKK
              RSCORE(I) = RSCORE(I) + WWNL(M) * SNLO(I)**2
            END DO
          ELSE
            DO I = 1,KKK
              RSVALE(I) = RSVALE(I) + WWNL(M) * SNLO(I)**2
            END DO
          END IF
          DO I = 1,KKK
            SNL(M,I) = SNLO(I)
          END DO
          NKKK(M) = KKK
          EE(M) = E
        END DO

C ----------------------------------------------------------------------
C Calculate the Total Charge Density and the Atomic Exchange Potential.
C MODIFIED RUEXCH by adding in factor ALPHA
C V_exch = -3 e^2 ((-3/8pi) * rho)^1/3
C 315.82734 = 8pi * 4pi  needed since RSATOM = 4pi*|rR|^2
C ----------------------------------------------------------------------
        DO I = 1,MESH
          RSATOM(I) = RSCORE(I) + RSVALE(I)
          RUEXCH(I) = -6.0 * ALPHA * 
     &               ((3.0 * R(I) * RSATOM(I)) / 315.82734)**(1.0/3.0)
        END DO

C ----------------------------------------------------------------------
C Calculate the Atomic Potential rV(r)
C   rV(r) = rV_en(r) + {rV_ee(r) + rV_ee,self(r)} + rV_exch(r).
C   Calculate the coulombic potential terms and add to the nuclear
C   potential for a neutral atom.  Adjust for ionic atom, then add
C   in the exchange potential.
C
C d RSATOR = RSATOM/R = normalized total charge density / r = 4pi*r*|R|^2
C            used to calculate the self-Coulomb potential
C
C REPLACED the simple trapezoidal integrals of RSATOM and RSATOM/R with
C          Simpson's rule integration subroutine
C ADDED ionic radius ability.  Set radius of charged sphere equal to
C       ionic radius.  This is read in as RADION.
C ----------------------------------------------------------------------
C       No electron density at the exact origin (nucleus) R(1) = 0.0
        RSATOR(1) = 0.0
        DO I = 2,MESH
          RSATOR(I) = RSATOM(I) / R(I)
        END DO
C       Call INTEGR so XI(I)=integral(r^2|R(r)|^2) = coulomb potential
C         with other electrons at mesh position I; rV_ee.
C         XI(0) = 0 and XI(MESH) should equal Z; Assuming neutral atom.
C       Call INTEGR so XJ(I)=integral(r|R(r)|^2) = self-coulomb potential
C         of electron outwards.
        CALL INTEGR(RSATOM, XI, R, NBLOCK)
        CALL INTEGR(RSATOR, XJ, R, NBLOCK)

        DO I = 1,MESH
C         Calculate atomic potential of the neutral atom; i.e. number of
C         electrons equals the number of protons(Z) is assumed.
C                  rV_en        rV_ee         rV_ee,self
          XI(I) = -2.0 * Z + 2.0*(XI(I) + R(I)*(XJ(MESH) - XJ(I)))
          IF (RADION .GT. 0.0) THEN
C           If an ionic radius is declared by radion being non-zero,
C           then redefine the potential so the radius of the charged
C           sphere equals the given ionic radius.
            XION = XIONSAV * BRATIO
            IF (R(I) .LT. RADION) THEN
C             Inside the ionic radius, scale the potential rV(r) so it
C             smoothly matches the theoretical result at  ionic radius.
C             NOTE:  I'm not sure this is completely correct.  It seems
C                    to me that they are not properly accounting for 
C                    the nuclear + coulomb potential.  I think it 
C                    should be + xion(xion+1)*r/radion
              XI(I) = XI(I) + 2.0 * XION * R(I)/RADION
            ELSE
C             Beyond ionic radius, potential rV(r) has a constant
C             additive correction for the ionization... potential does
C             not go to zero.
              XI(I) = XI(I) + 2.0 * XION
            END IF
          END IF
          XJ(I) = XI(I) + RUEXCH(I)
        END DO

        DO I = 1,MESH
C         Shuffle potentials rV(r).
C         Old initial and finals = prior iteration initial and finals.
          RUINL1(I) = RUINL2(I)
          RUFNL1(I) = RUFNL2(I)
C         New initial and finals = current initial and calculated final.
C         Initialized value of RU = input trial potential.
          RUINL2(I) = RU(I)
          RUFNL2(I) = XJ(I)
        END DO

C ----------------------------------------------------------------------
C Calculate the Error for Current Iteration.
C   Error is the largest absolute change in potential from pre-SCHEQ
C   initial value to post-SCHEQ calculated final value at any point
C   of the mesh.
C
C            NOTE: initialized DELTA, PDELTA = 1.0E6 from above
C i IDELTA = mesh point that gives the error value for a given 
C            iteration.  Only used for info file output.  Generally
C            the mesh point moves closer to the center with convergence
C            towards self-consistency.
C ----------------------------------------------------------------------
        PDELTA = DELTA
        DELTA = 0.0
        DO I = 1,LIMIT
          SNLO(I) = RU(I) - XJ(I)
          XI(I) = ABS(SNLO(I))
          IF (XI(I) .GT. DELTA) THEN
            DELTA = XI(I)
            IDELTA = I
          END IF
        END DO

C ----------------------------------------------------------------------
C End of Iteration Block
C ----------------------------------------------------------------------


        NITER = NITER+1
        WRITE ( 8,1015)
        WRITE ( 8,1016) NITER, Z, DELTA,
     &                  IDELTA, X(IDELTA), ICUT, X(ICUT)
 1015   FORMAT (2X, 'ITER', 4X, 'Z', 7X, 'DELTA', 8X,
     &          'I(DEL)    X(DEL)  ', '  I(CUT)    X(CUT)')
 1016   FORMAT (I5, 0PF7.1, 1PE16.8, I8, 0PF11.4, I9, 0PF11.4//)


      END DO
C ----------------------------------------------------------------------
C End of Main Iteration Scheme
C ----------------------------------------------------------------------




C ----------------------------------------------------------------------
C Self-Consistent Results Block
C
C ADDED DETAIL to skip printing atomic info in header file
C ----------------------------------------------------------------------
      IF (DETAIL) THEN
        WRITE ( 9,1018)
        WRITE ( 9,1008) Z, NCORES, NVALES, XION
 1018   FORMAT ('  Z  CORE VAL  ION')
      END IF

C ----------------------------------------------------------------------
C Write out the Modified HFS Potential rV(r)
C
C i NC = counts the output line number for the data being written to
C        file
C
C ADDED a KUT check for the type of limit to apply to the potential in 
C       the output and added modified potential result.
C ADDED DETAIL to skip printing potential in header file
C ADDED choice to print potential in individual file
C ----------------------------------------------------------------------
      IF (KUT .NE. 0) THEN
C       Apply simple coulomb correction to H-F-S potential
        ICUT = 0
        DO I = 1,441
          IF (ICUT .LE. 0) THEN
            IF ((TWOION + RUINL2(I)) .GT. 0.0) THEN
C             Found cut off value r0, start Latter tail correction
              ICUT = I
              RUINL2(I) = -TWOION
            END IF
          ELSE
C           From here on out, use Latter tail correction
            RUINL2(I) = -TWOION
          END IF
        END DO
      ELSE
C       Apply Latter correction to H-F-S potential
        ICUT = 0
        DO I = 1,441
          IF (ICUT .LE. 0) THEN
            IF ((TWOZZZ + RUINL2(I)) .GT. 0.0) THEN
C             Found cut off value r0, start Latter tail correction
              ICUT = I
              RUINL2(I) = -TWOZZZ
            END IF
          ELSE
C           From here on out, use Latter tail correction
            RUINL2(I) = -TWOZZZ
          END IF
        END DO
      END IF

      IF (DETAIL) THEN
        WRITE ( 9,*) ' '
        WRITE ( 9,*) ' '
        WRITE ( 9,*) ' '
        WRITE ( 9,*) 'THE FOLLOWING DATA IS THE MODIFIED HFS',
     &               ' POTENTIAL, RV(R), GIVEN AT THE'
        WRITE ( 9,*) '441 POINTS OF THE MESH, IN RYDBERG AUs,',
     &               ' IDEALLY = -2Z U'
        WRITE ( 9,*) ' '
        NC = 0
        DO MIN = 1,440,5
          MAX = MIN + 4
          NC = NC + 1
          WRITE ( 9,1020)(RUINL2(M), M = MIN,MAX), IZ, NC
 1020     FORMAT (1PE15.7, 1P4E14.7, ' Z', I3, I4)
        END DO
        WRITE ( 9,1022) RUINL2(441), IZ, NC+1
 1022   FORMAT (1PE15.7, 57X, 'Z', I3, I4)
      END IF
      IF (Q003 .EQ. REC2) THEN
C       True if desire to print potential (REC2 = ' POT') in an
C       individual file (for mufpot use, etc.)
        WRITE (10,*) 'HERMAN-S MODIFIED POTENTIAL RV(R) ON 441 pt',
     &               ' GRID, IN HARTREE AUs, IDEALLY = -Z U'
        WRITE (10,*) '441 GRID STARTS WITH dx=0.0025 AND DOUBLES',
     &               ' AFTER EVERY 40 POINTS.'
        WRITE (10,1003) TITLE
        WRITE (10,1009) Z, NCSPVS
 1003   FORMAT (9A8)
 1009 FORMAT (F9.4, /, I4)
        DO MIN = 1,440,5
          MAX = MIN + 4
          WRITE (10,1021) (RUINL2(M)/2.0, M = MIN,MAX)
 1021     FORMAT (1P5E14.7)
        END DO
        WRITE (10,1023) RUINL2(441)/2.0
 1023   FORMAT (1PE14.7)
      END IF

C ----------------------------------------------------------------------
C Write out Normalized Radial Wave Function
C   Loop through orbitals.
C
C ADDED choice to print header info for radial wf in individual file
C ADDED DETAIL to skip printing radial wf in header file
C ADDED choice to print radial wf in individual file
C ----------------------------------------------------------------------
      IF (Q004 .EQ. REC3) THEN
        WRITE (11,1017) 'HERMAN-S'
        WRITE (11,1003) TITLE
        WRITE (11,1009) Z, NCSPVS
      END IF
 1017 FORMAT (A8)

      DO M = 1,NCSPVS
        LAM = NNLZ(M) / 10 - 10 * (NNLZ(M)/100)
        KKK = NKKK(M)
C       Compute first term of series. (SNL(R)/R**(LAM+1) AT R=0)
        DO I = 1,4
          A(I,1) = 1.0
          A(I,2) = R(I+1)
          A(I,3) = R(I+1) * R(I+1)
          A(I,4) = R(I+1) * A(I,3)
          A(I,5) = SNL(M,I+1) / R(I+1)**(LAM+1)
        END DO
        CALL CROSYM (A, 4)

        IF (DETAIL) THEN
          WRITE ( 9,1024) NNLZ(M), KKK
 1024     FORMAT (///,' THE FOLLOWING DATA IS THE NORMALIZED RADIAL',
     &            ' WAVE FUNCTION FOR ORBITAL ', I3,
     &            ' , OUT TO MESH POINT ',I3,/)
          WRITE ( 9,1026) NNLZ(M), KKK, LAM, EE(M), WWNL(M), A(1,5)
 1026     FORMAT (I4, I7, 4X, I7, 7X, 1P3E14.7, /)
          NC = 0
          DO MIN = 1,KKK-1,5
            NC = NC + 1
            MAX = MIN + 4
            WRITE ( 9,1020) (SNL(M,I), I = MIN,MAX), IZ, NC
          END DO
          WRITE ( 9,1022) SNL(M,KKK), IZ, NC+1
        END IF
        IF (Q004 .EQ. REC3) THEN
C         True if desire to print radial wave function (REC3 = ' RAD')
C         in an individual file for MUFPOT use.
          WRITE (11,1025) LAM, KKK, WWNL(M) / (4.0*LAM + 2.0)
 1025     FORMAT (I4, /, I4, /, F9.4)
          DO MIN = 1,KKK-1,5
            MAX = MIN + 4
            WRITE (11,1021) (SNL(M,I), I = MIN,MAX)
          END DO
          WRITE (11,1023) SNL(M,KKK)
        END IF
      END DO

C ----------------------------------------------------------------------
C Write out the Charge Densities
C
C ADDED output of the charge densities
C ADDED DETAIL to skip printing charge densities in header file
C ----------------------------------------------------------------------
      IF (DETAIL) THEN
C ----------------------------------------------------------------------
C   Total Charge Densities write out
C ----------------------------------------------------------------------
        WRITE ( 9,1028)
 1028   FORMAT ('1', ///, 1X, 'TOTAL CHARGE DENSITY OUTPUT', //)
        NC = 0
        DO MIN = 1,436,5
          MAX = MIN + 4
          NC = NC + 1
          WRITE ( 9,1030) (RSATOM(I), I = MIN,MAX), IZ, NC
 1030     FORMAT (1PE15.7, 1P4E14.7, 1X, 'T', I3, I4)
        END DO
        NC = NC + 1 
        WRITE ( 9,1031) RSATOM(441), IZ, NC
 1031   FORMAT (1PE15.7, 57X, 'T', I3, I4)
        IF (NCORES .GT. 0) THEN
C ----------------------------------------------------------------------
C   Core Charge Densities write out
C ----------------------------------------------------------------------
          WRITE ( 9,1032)
 1032     FORMAT ('1', ///, 1X, 'CORE CHARGE DENSITY OUTPUT', //)
          NC = 0
          DO MIN = 1,436,5
            MAX = MIN + 4
            NC = NC + 1
            WRITE ( 9,1033) (RSCORE(I), I = MIN,MAX), IZ, NC
 1033       FORMAT (1PE15.7,1P4E14.7,1X,'C',I3,I4)
          END DO
          NC = NC + 1
          WRITE ( 9,1034) RSCORE(441), IZ, NC
 1034     FORMAT (1PE15.7, 57X, 'C', I3, I4)
        END IF
      END IF

C ----------------------------------------------------------------------
C End Self-Consistent Results Block
C ----------------------------------------------------------------------


C     Set up to start the next run on file?????
      IF (KEY .LE. 1) THEN
        KEY = 1
        DO I = 1,MESH
          RU3(I) = RU(I)
        END DO
        ZE3 = Z
      ELSE
C       Must be that KEY = 2
        DO I = 1,MESH
          RU2(I) = RU3(I)
          RU3(I) = RU(I)
        END DO
        ZE2 = ZE3
        ZE3 = Z
      END IF

     
      GO TO 7000
      STOP
      END
C ----------------------------------------------------------------------
C End of Original Herman Skillman Program
C ----------------------------------------------------------------------





C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C SCHEQ  :: Schroedinger equation subroutine
C   Written by Sherwood Skillman, RCA Laboratories, Princeton, NJ, 
C     Spring 1961.  Modified by Frank Herman, Summer 1961.  Modified
C     by Richard Kortum and Paul Kelly, Lockheed Research Laboratories,
C     Palo Alto, CA, Summer 1962.
C  Compute energy eigenvalue and wave function
C  Defines node in wave function at origin to exist and then searches
C    the rest of the wave function for the appropriate number of nodes.
C    SNLO(1) = 0.0 by definition.
C
C d ZZ = Z; atomic number
C d EE = initially; trial energy eigenvalue.  Finally; solved energy.
C i LL = l; angular (orbital) quantum number
C i NN = n; principle quantum number
C i KKK = number of mesh points used for wave function
C d CMU = scaling factor for mesh; convert between R and X
C d THRESH = eigenvalue accuracy criterion; compare vs DE/E
C d V = potential
C d R = radial mesh
C d SNLO = normalized radial wave function
C i I, J, II = temp variables mainly for loop control
C d E = trial energy of w.f. being manipulated by current iteration
C d PE = prior iteration energy of wave function
C d EPS = shift in energy value from prior iteration; E - PE
C i LTMSLPLS1 = l(l+1)
C d QQ = Effective Kinetic Energy; = E - (V(r) + l(l+1)/r^2)
C        It should equal the negative of the effective potential shifted
C        down by the constant abs(E).
C d P = integration variable; 
C   initialized with a power series expansion
C   P(3) = (1 + A1*DX  + A2*DX^2  + A3*DX^3  + A4*DX^4)  * DX^(l+1)
C   P(4) = (1 + A1*TDX + A2*TDX^2 + A3*TDX^3 + A4*TDX^4) * TDX^(l+1)
C   P(5) = SNLO(I)
C d Q = integration variable; 
C   initialized with a power series expansion
C   Q(3) = (l(l+1) + B1*DX  + B2*DX^2  + B3*DX^3)  / DX^2
C   Q(4) = (l(l+1) + B1*TDX + B2*TDX^2 + B3*TDX^3) / TDX^2
C   Q(5) = -QQ(I)
C d T, D = integration variables; D = difference in T
C i MORE = count of the # of times E is too small; too many wf crossings
C i LESS = count of the # of times E is too large; too few wf crossings
C d EMORE = save best E giving too many crossings
C d ELESS = save best E giving too few crossings
C i NPRINT = count of outward iterations integration has been tried
C i MAXITER = maximum # of iterations to try integration before fail
C i IQQNEG = outer radius mesh point where kin. ener. switches to negat.
C d DX = current step size current block of mesh being accessed
C d TDX = 2 DX
C d DXTOLPLS1 = DX^(l+1); used by initialization of P for outward int.
C d TDXTOLPLS1 = (2 DX)^(l+1); used by initialization of P for out int.
C d B1 = -2 Z
C d B2 = 3 Z/R(2) - E + 2 V(2) - V(3)
C d B3 = (V(3) - V(2))/R(2) - Z/R(2)^2
C d A1 = -Z/(l+1)
C d A2 = (A1*B1 +    B2)        /(4l+6)
C d A3 = (A2*B1 + A1*B2 +    B3)/(6l+12)
C d A4 = (A3*B1 + A2*B2 + A1*B3)/(8l+20)
C i NCROSS = # of times the wf crosses the axis
C d SIGN = +/- 1.0; sign of most recent wf crossing the axis
C i NCOUNT = # of the term of w.f. calculated within the current 5 term 
C            interval of the power series integration.
C i NINT = position of outward integration within current step-size
C          block of the mesh.
C c OUTINTFLAG = control flag for outward integration status:
C                'CONTINU', 'FAILED', 'SUCCESS'
C i IMATCH = mesh point found from outward integration that is 
C            immediately beyond matching radius and starts the second,
C            or later, 5 term interval of a step-size block
C d PMATCH = wave func. of IMATCH radius for outward integration
C d S2 = logarithmic derivative at matching radius for out. integ.
C        SNLO'(Rm)/SNLO(Rm)
C d XIF = integration factor 5*DX(block1)/288/2; DX(block1) = 0.0025CMU
C d SUM1 = temp variable used to find S1
C d SUM2 = temp variable used to find SUM1; 5 interval Newton-Cotes
C          closed quadrature integration.
C d VALUE = prior SNLO(I+5)^2 used to calc SUM2
C d PVALUE = SNLO(I+5)^2 used to calc SUM2
C d S1 = partial normalization factor at matching radius for out. integ.
C        from R = 0 to R = matching radius: int[(SNLO(I))^2/(SNLO(Rm)^2]
C d XINW = radius to start inward integration at; matching radius times
C          (8+l) for n=1 and (5+l) for n>1
C i KKK = mesh point at which inward integration is started
C d SUM3 = temp variable used to find S3
C d SUM4 = temp variable used to find SUM3
C d S3 = partial normalization factor at matching radius for inw. integ.
C        from R = match rad to R = infin.: int[(SNLO(I))^2/(SNLO(Rm)^2]
C d S4 = logarithmic derivative at matching radius for inw. integ.
C        SNLO'(Rm)/SNLO(Rm)
C d DE = delta energy eigenvalue = (S2 - S4) / (S1 - S3)
C
C REMOVED the COMMON and EQUIVALENCE statements and modified the
C         function calls to properly pass the necessary variables.
C ADDED control flag for outward integration status
C REMOVED a bunch of the temp variables used for more extreme speed of 
C         execution savings.  Kept only the variables that had a large
C         affect like ones active in a loop.  Tried to improve
C         readability of program
C ----------------------------------------------------------------------
      SUBROUTINE SCHEQ(ZZ, EE, LL, NN, KKK, MESH, CMU, THRESH,
     &                 V, R, SNLO)
        DIMENSION QQ(521), P(5), Q(5), T(5), D(5),
     &            V(521), R(521), SNLO(541)

        CHARACTER*7 OUTINTFLAG

C ----------------------------------------------------------------------
C SCHEQ Main Iteration Scheme
C   Search for a solution is a two part integration process.  Outward
C     integration from the origin to the matching radius at which the 
C     kinetic energy switches to negative.  Then, if outward int. is
C     valid, start inward integration to the same radius from a point
C     somewhat beyond that.  Integrations work off the effective
C     kinetic energy.
C   Initialize the iteration scheme with the user supplied energy and
C     use it to calculate the initial trial effective kinetic energy.
C     NOTE:  E should be negative to have a bound state.
C   Starting first iteration; NPRINT = 1
C
C ADDED initialization of variables
C ----------------------------------------------------------------------
        DO I = 1,MESH
          SNLO(I) = 0.0
        END DO
        E = EE
        PE = E
        EPS = 0.0
        LTMSLPLS1 = LL*(LL + 1)
        QQ(1) = 0.0
        QQ(2) = 0.0
        QQ(3) = 0.0
        DO I = 4,MESH
          QQ(I) = E - (V(I) + LTMSLPLS1/R(I)/R(I))
        END DO
        DO I = 1,5
          P(I) = 0.0
          Q(I) = 0.0
          T(I) = 0.0
          D(I) = 0.0
        END DO
        MORE = 0
        LESS = 0
        EMORE = 0.0
        ELESS = 0.0
        NPRINT = 1
        MAXITER = 200

        DO WHILE (NPRINT .LE. MAXITER)
C         For the 17 elements, 0N and 1P, I have tried, NPRINT rarely
C         goes beyond 2; then, usually only to 3 and once to 5.

C ----------------------------------------------------------------------
C Find Kinetic Energy Switch Point
C   Work backwards along mesh to find the outer radius point where
C     kinetic energy switches to negative = matching radius; must occur
C     beyond QQ(4) or else the state can not exist - kinetic energy is
C     always negative.
C   Kinetic energy must be zero or negative at least for the last two
C     radial mesh positions or else we do not have a bound state within
C     the scope of the mesh.  Lower the trial energy by an amount equal
C     to the kinetic energy of the mesh point one block prior to the
C     end of mesh and try again.
C ----------------------------------------------------------------------
          I = MESH
          DO WHILE ((QQ(I) .LE. 0.0) .AND. (I .GT.3))
            I = I - 1
          END DO
          IQQNEG = I + 1
          IF (IQQNEG .LT. 5) CALL DUMP(991,E)
          DO WHILE (IQQNEG .GE. MESH)
            EPS = -QQ(MESH-40)
            E = E + EPS
            DO I = 4,MESH
              QQ(I) = QQ(I) + EPS
            END DO

            I = MESH
            DO WHILE ((QQ(I) .LE. 0.0) .AND. (I .GT.3))
              I = I - 1
            END DO
            IQQNEG = I + 1
            IF (IQQNEG .LT. 5) CALL DUMP(991,E)
          END DO


C ----------------------------------------------------------------------
C Outward Integration Block
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C Initialize the Outward Integration
C   Set up integration variables with power series expansion.
C   Second two terms of radial wave function found from power series
C     expansion; SNLO(2), SNLO(3).  Radial w.f. is defined at the
C     origin; SNLO(1) = 0.0.
C ----------------------------------------------------------------------
          DX = R(2)
          TDX =  DX + DX
          DXTOLPLS1 = DX
          TDXTOLPLS1 = TDX
C         Raise R(2) and 2R(2) to the (l+1).  l = 0 done automatically
C         by definition of DXTOLPLS1 and TDXTOLPLS1.
          DO I = 1,LL
            DXTOLPLS1 = DXTOLPLS1 * DX
            TDXTOLPLS1 = TDXTOLPLS1 * TDX
          END DO
          B1 = -2.0*ZZ
          B2 = 3.0*ZZ/DX - E + 2.0*V(2) - V(3)
          B3 = (V(3) - V(2))/DX - ZZ/DX/DX
          A1 = -ZZ / (LL + 1.0)
          A2 = (A1*B1 +    B2)         / (4.0*LL + 6.0)
          A3 = (A2*B1 + A1*B2 +    B3) / (6.0*LL + 12.0)
          A4 = (A3*B1 + A2*B2 + A1*B3) / (8.0*LL + 20.0)
          P(3) = (1.0 + DX*(A1 + DX*(A2 + DX*(A3 + DX*A4)))) * DXTOLPLS1
          P(4) = 
     &      (1.0 + TDX*(A1 + TDX*(A2 + TDX*(A3 + TDX*A4)))) * TDXTOLPLS1
          Q(3) = LTMSLPLS1/DX/DX   + B1/DX  + B2 + B3*DX
          Q(4) = LTMSLPLS1/TDX/TDX + B1/TDX + B2 + B3*TDX
          T(3) = P(3) * (1.0 - Q(3)*DX*DX/12.0)
          T(4) = P(4) * (1.0 - Q(4)*DX*DX/12.0)
          D(4) = T(4) - T(3)
          SNLO(2) = P(3)
          SNLO(3) = P(4)

C ----------------------------------------------------------------------
C Outward Integration
C   Calculate the wave function by integrating outward from the origin.
C     Every time the sign of the w.f. changes, record a node.  Continue
C     until end of mesh or matching radius is reached.  If the matching
C     radius has been reached, verify that the current mesh point fits
C     calculation requirements and then test the results.
C   By excepting the origin, no nodes have yet been found in the wave
C     function; NCROSS = 0.  Radial wave function starts out positive;
C     SIGN = positive.
C   Starting with the current wave function term calculated (SNLO(3))
C     within the current set is the NCOUNT = 3 term. It is at position
C     NINT = 3 of the 40 points within the current step-size block of
C     the mesh. The new term being calculated is the mesh point I = 4
C     term.
C ----------------------------------------------------------------------
          NCROSS = 0
          SIGN = 1.0
          NCOUNT = 3
          NINT = 3
          I = 4
          OUTINTFLAG = 'CONTINU'
          DO WHILE ( (I .LT. MESH) .AND. ('CONTINU' .EQ. OUTINTFLAG) )

C ----------------------------------------------------------------------
C Possible Solution to Outward Integration
C   The matching radius has been reached and the integration has
C     completed 2 terms into a defined 5 term interval that is at least
C     the second such interval within the current step-size block.  The
C     requirements on NCOUNT and NINT are for S2 and to determine a well
C     defined matching radius mesh point for the Newton-Cotes
C     integration used by S1 and S3.
C   If the # of nodes matches theory, check to see that wave is in the
C     damped region (absolute value decreasing and signs alike)
C ----------------------------------------------------------------------
            IF ( (I .GE. IQQNEG) .AND. (2 .EQ. NCOUNT)
     &           .AND. (NINT .GT. 5) ) THEN

              IF ((NN - LL - 1) .EQ. NCROSS) THEN
                IF (ABS(SNLO(I-1)) .LT. ABS(SNLO(I-2))) THEN
                  IF ( ((P(5) .LT. 0.0) .AND. (SNLO(I-2) .LT. 0.0)) .OR.
     &                 ((P(5) .GT. 0.0) .AND. (SNLO(I-2) .GT. 0.0)) )
     &            OUTINTFLAG = 'SUCCESS'
                ELSE
                  IF (ABS(P(5)) .GT. 1.0E25) THEN
C                   Large absolute value of P in what should be the
C                   damped region indicates too few peaks.  Set the
C                   control flags appropriately to raise E.
                    OUTINTFLAG = 'FAILED '
                    NCROSS = 0
                  END IF
                END IF
              ELSE 
                OUTINTFLAG = 'FAILED '
              END IF

            END IF

C ----------------------------------------------------------------------
C Single Step of Outward Integration
C   Record SNLO for this step of the outward integration.  Progress
C     of the integration is tracked as completed sections defined by 
C     5 term intervals; increment or reset NCOUNT.
C   Test this step for the w.f. making a new crossing of the axis.  If
C     true, then increment the node count NCROSS and store the new SIGN
C     of the wave function.
C   Test this step for reaching a new step-size block within the mesh. 
C     If true, reset the integration variables for the new DX.
C     Increment or reset NINT.
C   Move the integration variables P, Q, T, D to start the next step.
C ----------------------------------------------------------------------
            IF ( (ABS(QQ(I)*DX*DX) .LT. 12.0) .AND.
     &           ('CONTINU' .EQ. OUTINTFLAG) ) THEN
C             Successful step in integration
              D(5) = D(4) + Q(4)*P(4)*DX*DX
              T(5) = D(5) + T(4)
              Q(5) = -QQ(I)
              P(5) = T(5) / (1.0 - Q(5)*DX*DX/12.0)
              SNLO(I) = P(5)
              NCOUNT = NCOUNT + 1
              IF (6 .EQ. NCOUNT) THEN
C               The prior 5 term interval of the power series
C               expansion integral has been completed and the current
C               "6th" term of the interval is actually the 1st term
C               of the new interval.  There are 8 intervals to each
C               40 point step-size block of the mesh.
                NCOUNT = 1
              END IF

              IF ( ((SIGN .LT. 0) .AND. (SNLO(I) .GT. 0)) .OR.
     &             ((SIGN .GT. 0) .AND. (SNLO(I) .LT. 0)) ) THEN
                NCROSS = NCROSS + 1
                SIGN = -SIGN
              END IF

              NINT = NINT + 1
              IF (41 .EQ. NINT) THEN
C               Reached new step-size block within the mesh.  
C               Reset info to continue integration over next block.
                NINT = 1
                DX = DX + DX
                T(3) = P(3) * (1.0 - Q(3)*DX*DX/12.0)
                T(5) = P(5) * (1.0 - Q(5)*DX*DX/12.0)
                D(5) = T(5) - T(3)
              END IF

              DO J = 2,4
                P(J) = P(J+1)
                T(J) = T(J+1)
                Q(J) = Q(J+1)
              END DO
              P(1) = P(2)
              D(4) = D(5)
              I = I + 1
            ELSE IF (OUTINTFLAG .NE. 'SUCCESS') THEN
C             Bad step in integration
              OUTINTFLAG = 'FAILED '
            END IF
          END DO
C         Rare case when IQQNEG is only a few points away from end of
C         mesh and so cannot fulfill NCOUNT or NINT requirement without
C         I going beyond MESH
          IF (I .GE. MESH) OUTINTFLAG = 'FAILED '
C ----------------------------------------------------------------------
C End of Outward Integration Block
C ----------------------------------------------------------------------


          IF ('FAILED ' .EQ. OUTINTFLAG) THEN
C ----------------------------------------------------------------------
C NCROSS Failed to Match Theory Value
C   Matching radius or end of mesh has been reached going out, but
C     NCROSS does not equal theoretical value for number of nodes
C     in radial wavefunction, n-l-1, excluding the origin (node
C     at origin is defined to exist).
C   Modify trial energy eigenvalue and return to beginning of outward
C     integration.  Store the energy value if it is the best E for the
C     result, so far.  Shift the energy up or down based on current
C     and prior behavior of wavefunction; has a "binary search-like"
C     approach to convergence on the SCF energy eigenvalue
C ----------------------------------------------------------------------
C           From the 17 elements, 0N and 1P, I have tried, MORE went
C           past 0 once and then only to 3.  LESS never went past 0.
            IF (NCROSS .GT. (NN - LL - 1)) THEN
C             Too many crossings: lower E
              MORE = MORE + 1
              IF ((E .LT. EMORE) .OR. (1 .EQ. MORE)) THEN
                EMORE = E
              END IF
              IF (0 .EQ. LESS) THEN
C               Never had case of too few crossings, so big adjust
                E = 1.25 * PE
              ELSE
                E = 0.5 * (EMORE + ELESS)
              END IF
            ELSE
C             NCROSS <= NN-LL-1
C             Too few crossings: raise E
              LESS = LESS + 1
              IF ((E .GT. ELESS) .OR. (1 .EQ. LESS)) THEN
                ELESS = E
              END IF
              IF (0 .EQ. MORE) THEN
C               Never had case of too many crossings, so big adjust
                E = 0.75 * PE
              ELSE
                E = 0.50 * (EMORE + ELESS)
              END IF
            END IF

          ELSE 
C ----------------------------------------------------------------------
C NCROSS Good and Matching Radius Good Block
C   Matching radius is found and lies in a properly damped region plus
C     NCROSS = theoretical value n-l-1 (excluding origin);
C     OUTINTFLAG = 'SUCCESS'.
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C Outward Integration Boundary Conditions
C   Calculate logarithmic derivative, S2, and partial normalization
C     factor, S1, for outward integration at the matching radius.
C     NOTE: Integration was taken two extra steps to find the values
C           used to calculate S2 that relate to the actual chosen
C           matching radius mesh point, IMATCH (I is left equal to
C           the 3rd mesh point into the 5 term interval by above).
C ----------------------------------------------------------------------
            IMATCH = I - 2
            PMATCH = SNLO(IMATCH)
            S2 = ((T(4) - T(2)) - 0.5*(P(4) - P(2)))/DX / PMATCH

C           Integration is by 8 applications of Newton-Cotes closed
C           quadrature for five intervals on each block.
C           XIFC =(5*DX(block1)/288)/2; DX(block1)=0.0025*CMU
C            XIF = CMU * 0.0125/576
            XIF = 2.0 * CMU * 2.1701389E-5
            SUM1 = 0.0
            SUM2 = 0.0
            VALUE = 0.0
            I = 1
            J = 1
            DO WHILE (IMATCH .NE. I)
ccc              IF (IMATCH .LT. I) CALL DUMP(990,E)
              IF (J .GT. 8) THEN
                J = 1
                SUM1 = SUM1 + XIF  * SUM2
                SUM2 = 0.0
                XIF = XIF + XIF
              END IF
              PVALUE = VALUE
              VALUE = SNLO(I+5)**2
              SUM2 = SUM2 + 19.0*(VALUE + PVALUE) 
     &               + 75.0*(SNLO(I+4)**2 + SNLO(I+1)**2) 
     &               + 50.0*(SNLO(I+2)**2 + SNLO(I+3)**2)
              I = I + 5
              J = J + 1
            END DO
            SUM1 = SUM1 + XIF * SUM2
            S1 = SUM1 / PMATCH**2


C ----------------------------------------------------------------------
C Inward Integration Block
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C Initialize Inward Integration
C   Use XINW to find an acceptable I to start integrating backwards from
C     that starts at the first point of a new step-size block just
C     beyond R = XINW
C ----------------------------------------------------------------------
            IF (1 .EQ. NN) THEN
C             Start inward integration at (8+l)*Rmatch or end of mesh
              XINW = (8 + LL) * R(IMATCH)
            ELSE IF (NN .GT. 1) THEN
C             Start at (5+l)*Rmatch or end of mesh
              XINW = (5 + LL) * R(IMATCH)
            END IF
            I = 41
            DO WHILE ((XINW .GT. R(I)) .AND. (I .LT. MESH))
              I = I + 40
            END DO
C           KKK equals the point found or MESH
            KKK = I

            I = I - 1
            DX = R(I) - R(I+1)
            XIF= 0.17361111E-1 * DX
            Q(3) = -QQ(I+1)
            P(3) = EXP(-R(I+1) * SQRT(Q(3)))
            Q(4) = -QQ(I)
            P(4) = EXP(-R(I) * SQRT(Q(4)))
            SUM3 = P(3) / Q(3)
            DO WHILE (ABS(P(4)) .LE. 1.0E-35)
              KKK = KKK - 40
              IF (KKK .LE. IMATCH) THEN
                WRITE ( 8,9001) ZZ, NN, LL, KKK
 9001           FORMAT  ('0AT Z=',F6.0,'  N,L =',I3,I1,'  KKK =',I5,
     &                   ' IS LESS THAN IMATCH =',I5,
     &                   ' INWARD INTEGRATION WILL BE TRIED AT KKK+40')
                KKK = KKK + 40
                P(3) = 1.0E-35
                P(4) = 1.5E-35
C               Leave the DO WHILE loop
                EXIT
              END IF
C             Reset the testing parameters and loop again
              I = KKK - 1
              DX = R(I) - R(I+1)
              XIF= 0.17361111E-1 * DX
              Q(3) = -QQ(I+1)
              P(3) = EXP(-R(I+1) * SQRT(Q(3)))
              Q(4) = -QQ(I)
              P(4) = EXP(-R(I) * SQRT(Q(4)))
              SUM3 = P(3) / Q(3)
            END DO

            IF (PMATCH .LT. 0.0) THEN
              P(3) = -P(3)
              P(4) = -P(4)
            ELSE IF (0.0 .EQ. PMATCH) THEN 
              CALL DUMP(990,E)
            END IF
            SNLO(I+1) = P(3)
            SNLO(I) = P(4)
            T(3) = P(3) * (1.0 - Q(3)*DX*DX/12.0)
            T(4) = P(4) * (1.0 - Q(4)*DX*DX/12.0)
            D(4) = T(4) - T(3)

C ----------------------------------------------------------------------
C Inward Integration
C ----------------------------------------------------------------------
C           Start loop over I down to IMATCH-1
            J = 2
            I = I - 1
            Q(5) = -QQ(I)
            D(5) = D(4) + Q(4)*P(4)*DX*DX
            T(5) = D(5) + T(4)
            P(5) = T(5) / (1.0 - Q(5)*DX*DX/12.0)
            DO WHILE (I .NE. IMATCH-1)
              IF (I .LT. IMATCH-1) CALL DUMP(990,E)
              SNLO(I) = P(5)
              DO II = 1,4
                P(II) = P(II+1)
                T(II) = T(II+1)
                Q(II) = Q(II+1)
              END DO
              D(4) = D(5)
              J = J + 1
              IF (J .GT. 40) THEN
C               Reached the end of the current mesh block, reset values
C               for the new step size of the new block
                Q(5) = -QQ(I-2)
                D(5) = D(4) + Q(4)*P(4)*DX*DX
                T(5) = D(5) + T(4)
                P(5) = T(5) / (1.0 - Q(5)*DX*DX/12.0)
C               ? perturbation adjustment of P ?
                P(5) = 1.09375*P(4) + 0.2734375*P(5) - 0.546875*P(3)
     &                 + 0.21875*P(2) - 0.0390625*P(1)
                I = I - 1
                DX = DX / 2.0
                Q(5) = -QQ(I)
                T(4) = P(4) * (1.0 - Q(4)*DX*DX/12.0)
                T(5) = P(5) * (1.0 - Q(5)*DX*DX/12.0)
                D(5) = T(5) - T(4)
                SNLO(I) = P(5)
                DO II = 1,4
                  P(II) = P(II+1)
                  T(II) = T(II+1)
                  Q(II) = Q(II+1)
                END DO
                D(4) = D(5)
                J = 2
              END IF
              I = I - 1
              Q(5) = -QQ(I)
              D(5) = D(4) + Q(4)*P(4)*DX*DX
              T(5) = D(5) + T(4)
              P(5) = T(5) / (1.0 - Q(5)*DX*DX/12.0)
            END DO

C           Matching radius has been reached coming in.
C ----------------------------------------------------------------------
C Inward Integration Boundary Conditions
C   Calculate logarithmic derivative, S4, and partial normalization
C     factor, S3, for outward integration at the matching radius.
C     NOTE: Integration was taken two extra steps to find the values
C           used to calculate S2 that relate to the actual
C           matching radius mesh point, IMATCH.
C ----------------------------------------------------------------------
            SUM4 = 0.0
            VALUE = SNLO(KKK)**2
            I = KKK
            J = 1
            DO WHILE (I .NE. IMATCH)
ccc              IF (I .LT. IMATCH) CALL DUMP(990,E)
              IF (J .GT. 8) THEN
                J = 1
                SUM3 = SUM3 + XIF * SUM4
                SUM4 = 0.0
                XIF = 0.5 * XIF
              END IF
              PVALUE = VALUE
              VALUE = SNLO(I-5)**2
              SUM4 = SUM4 + 19.0*(VALUE + PVALUE)
     &               + 75.0*(SNLO(I-1)**2 + SNLO(I-4)**2)
     &               + 50.0*(SNLO(I-2)**2 + SNLO(I-3)**2)
              I = I - 5
              J = J + 1
            END DO

            SUM3 = SUM3 + XIF * SUM4
            S3 = SUM3 / P(4)**2
            S4 = ( (T(5) - T(3)) - 0.5*(P(5) - P(3)) ) /DX/P(4)

C ----------------------------------------------------------------------
C End ofInward Integration Block
C ----------------------------------------------------------------------


C           The value (S1 - S3) tends to scale the vale (S2 - S4) 
C           appropriately for the energy E of the wavefunction; thus,
C           all orbitals will have a similar comparison with thresh.
            DE = (S2 - S4) / (S1 - S3)
            IF (ABS(DE/E) .LT. THRESH) THEN
C ----------------------------------------------------------------------
C Successful SCF Iteration 
C   Calculate the normalized wave functions.
C   Return to program.
C ----------------------------------------------------------------------
              SUM1 = PMATCH / P(4)
              DO J = IMATCH,KKK
                SNLO(J) = SNLO(J) * SUM1
              END DO

C             NOTE: would be more efficient as a DO[statements]WHILE
              I = 1
C              XIF = CMU * 0.0125/576
              XIF = 2.0 * CMU * 2.1701389E-5
              SUM1 = 0.0
              SUM2 = 0.0
              VALUE = 0.0
              DO J = 1,8
                PVALUE = VALUE
                VALUE = SNLO(I+5)**2
                SUM2 = SUM2 + 19.0*(VALUE + PVALUE)
     &                 + 75.0*(SNLO(I+4)**2 + SNLO(I+1)**2)
     &                 + 50.0*(SNLO(I+2)**2 + SNLO(I+3)**2)
                I = I + 5
              END DO
              SUM1 = SUM1 + XIF * SUM2
              DO WHILE (KKK .NE. I)
                IF (KKK .LT. I) CALL DUMP(990,E)
                XIF = XIF + XIF
                SUM2 = 0.0
                DO J = 1,8
                  PVALUE = VALUE
                  VALUE = SNLO(I+5)**2
                  SUM2 = SUM2 + 19.0*(VALUE + PVALUE)
     &                   + 75.0*(SNLO(I+4)**2 + SNLO(I+1)**2)
     &                   + 50.0*(SNLO(I+2)**2 + SNLO(I+3)**2)
                  I = I + 5
                END DO
                SUM1 = SUM1 + XIF * SUM2
              END DO

              SUM1 = SQRT(SUM1)
              IF (0.0 .EQ. SNLO(3)) THEN
                CALL DUMP(990,E)
              ELSE IF (SNLO(3) .LT. 0.0) THEN
                SUM1 = -SUM1
              END IF
              DO I = 1,KKK
                SNLO(I) = SNLO(I) / SUM1
              END DO
              EE = E

              RETURN
            END IF

C ----------------------------------------------------------------------
C Failed SCF Criterion
C   Generate a new trial energy eigenvalue for the wave function from
C     the energy determined by the outward integration and the DE
C     from the boundary matching condition (perturbation theory 
C     adjustment).
C   Ensure that the energy is for a bound state problem.
C ----------------------------------------------------------------------
            E = E + DE
            DO WHILE (E .GE. 0.0)
              E = E - DE
              DE = DE / 2.0
              E = E + DE
            END DO
C ----------------------------------------------------------------------
C End of NCROSS Good and Matching Radius Good Block
C ----------------------------------------------------------------------
          END IF


C ----------------------------------------------------------------------
C New Trial Values Generation
C   Reinitialize the normalized radial wave function for the new trial.
C   Store the difference between the prior energy and the new trial
C     energy; then store the new trial energy. 
C     NOTE: new trial energy comes either from just outward 
C           integration or from failed SCF criterion.
C   Generate the new trial effective kinetic energy to be used by the
C     next round of integrations.
C   Increment iteration count.
C ----------------------------------------------------------------------
          DO I = 1,MESH
            SNLO(I) = 0.0
          END DO
          EPS = E - PE
          PE = E
          DO I = 4,MESH
            QQ(I) = QQ(I) + EPS
          END DO
          NPRINT = NPRINT + 1
        END DO
C ----------------------------------------------------------------------
C End of SCHEQ Main Iteration Scheme
C ----------------------------------------------------------------------


C       Number of iterations is beyond defined limit of trials for scf
        WRITE ( 8,9000) NN,LL ,ZZ
 9000   FORMAT ('   NO CONVERGENCE ON NN=',I4,'    LL=',I2,
     &          '     Z=',F6.2)
        STOP

      END
C ----------------------------------------------------------------------
C End of SCHEQ: Schroedinger Equation Solver
C ----------------------------------------------------------------------


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C INTEGR  :: Simpson's rule integration scheme
C   In general, Y = Integral{X(R) dR} takes X, R, N and returns Y.
C   It is exact for polynomials up to degree 3 due to symmetries.
C   Y(I) = value of finite integral from 0(or R(0)) to R(I).
C    
C
C   N = # of blocks in the mesh grid, each block having a new dR
C
C ADDED simpson's subroutine to improve the integration for the atomic
C   coulomb potential.  NOTE: An individual Simpson's rule integration
C   uses an interval of 2*deltaR.  The integrals are set up so the 
C   area used goes backwards from the point selected; i.e. Y(41), which
C   is right at the break point for starting a new deltaR, uses Y(39) to
C   Y(41) - in which both intervals are the old deltaR.  Y(43), then,
C   uses Y(41) to Y(43) - in which both intervals are the new deltaR.
C   Y(42), however, uses Y(40) to Y(42) - in which both the old and new
C   deltaR are involved, and so it is recalculated using a weighted formula
C   that only uses the new deltaR.
C ----------------------------------------------------------------------
      SUBROUTINE INTEGR(X, Y, R, N)

        DIMENSION X(521), Y(521), R(521)

C       Initializing two separate integrations:
C         Y(1) for the odd integration points Y(3),Y(5),Y(7),...
C         Y(2) for the even integration points Y(4),Y(6),Y(8),...
        H = R(2)
        Y(1) = 0.0
        Y(2) = H*(5.0*X(1) + 8.0*X(2) - X(3)) / 12.0

        DO J = 1, N
          DO K = 1, 40
            I = 40*(J - 1) + K
            IF (I .LT. 40*N) THEN
C             Make sure we only go up to MESH-2 so Y(I+2) doesn't
C             go out of range of the declared array space.  NOTE: 40*N = MESH-1
              Y(I+2) = Y(I) + H*(X(I) + 4.0*X(I+1) + X(I+2)) / 3.0
            END IF
          END DO
          H = H + H
          IF (I .NE. 40*N) THEN
C           If not at the end of the mesh, then only at the end of
C           the current block; reinitialize the even integration points for the
C           next block of 40 mesh points.  The odd integration points work
C           naturally and don't need this.
            Y(I+2) = Y(I+1)
     &             + H*(5.0*X(I+1) + 8.0*X(I+2) - X(I+3)) / 12.0
          END IF
        END DO

        RETURN
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C CROSYM  :: Simultaneous equation solver.
C   Written by I. C. Hanson, Scientific Computation Dept., Lockheed
C     Missiles and Space Co., Sunnyvale, CA.
C   Solve M simultaneous equations by the method of Crout.
C
C i I1, I3, J2 = temp. variables used to track indices of A
C i I, J = looping variables
C d SUM = temp. variable used to work with A
C i M = dimension of array passed to function
C
C REMOVED the COMMON and EQUIVALENCE statements and modified the
C         function calls to properly pass the necessary variables.
C ADDED error check for divide by zero when evaluating A
C ----------------------------------------------------------------------
      SUBROUTINE CROSYM(A, M)

        DIMENSION A(4,5)

        DO I1 = 1, M
          IF (I1 .NE. M) THEN
            I3 = I1
            SUM = ABS(A(I1,I1))
            DO I = I1,M
              IF (SUM.LT.ABS(A(I,I1))) THEN
                I3 = I
                SUM = ABS(A(I,I1))
              END IF
            END DO
            IF (I3.NE.I1) THEN
              DO J = 1,M+1
                SUM = -A(I1,J)
                A(I1,J) = A(I3,J)
                A(I3,J) = SUM
              END DO
            END IF
            I3 = I1 + 1
            DO I = I3,M
              IF (0.0 .EQ. ABS(A(I1,I1))) THEN
                A(I,I1) = 0.0
              ELSE
                A(I,I1) = A(I,I1) / A(I1,I1)
              END IF
            END DO
          END IF

          J2 = I1 - 1
          I3 = I1 + 1
          IF (J2 .NE. 0) THEN 
            DO J = I3,M+1
              DO I = 1,J2
                A(I1,J) = A(I1,J) - A(I1,I) * A(I,J)
              END DO
            END DO
            IF (M .EQ. I1) THEN
              DO I = 1,M
                J2 = M - I
                I3 = J2 + 1
                IF (0.0 .EQ. ABS(A(I3,I3))) THEN
                  A(I3,M+1) = 0.0
                ELSE
                  A(I3,M+1) = A(I3,M+1) / A(I3,I3)
                END IF
                IF (J2 .NE. 0) THEN
                  DO J = 1,J2
                    A(J,M+1) = A(J,M+1) - A(I3,M+1) * A(J,I3)
                  END DO
                ELSE
                  EXIT
                END IF
              END DO

              RETURN
            END IF
          END IF

          J2 = I1
          DO I = I1+1,M
            DO J = 1,J2
              A(I,I1+1) = A(I,I1+1) - A(I,J) * A(J,I1+1)
            END DO
          END DO

        END DO
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C DUMP  :: dummy program for Herman's SCF
C   Kills the program execution when an error is found in SCHEQ.
      SUBROUTINE DUMP(NSTP,E)
        WRITE ( 8,*) '!!!ERROR at NSTOP = ', NSTP, 
     &               'in subroutine SCHEQ.  ', E
        WRITE (*,*) 'Humbug, we have a problem in SCHEQ.'
        STOP
      END
