!!SUBROUTINE SI(X,N,TIME,Z,M,DELTA,SM,S,ETA)!Based on the quantities defined in the paper "Observation and Characterization of Chimera States in Coupled Dynamical Systems with Nonlocal Coupling", by R. Gopal et. al. in PRE 89.052914 (2014)
!!
!! Written by Vagner dos Santos
!! Professor at State Univertsity of Ponta Grossa, Paraná State,  Brazil.
!!
   SUBROUTINE SI(X,!(INPUT) -- X = MATRIX CONTAINING THE TIME SERIES OF EACH NODE
   		 N,!(INPUT) -- N = SIZE OF THE NETWORK
   		 TIME,!(INPUT) -- TIME = TIME EXTENT OF THE SERIES
   		 Z,!(INPUT) -- TIME SERIES OF THE DIFFERENCE VARIABLES -- Z = X(i+1)-X(i)
   		 M,!(INPUT) -- M = NUMBER OF BINS,
   		 DELTA,!(INPUT) -- DELTA = THRESHOLD FOR THE DETECTION OF INCOHERENT NODES
   		 SDEV,!(OUTPUT) -- SDEV(M) = STANDARD DEVIATION OF THE NODES IN EACH BIN
   		 SM,!(OUTPUT) -- SM(M) = BINARY VARIABLE THAT IS EQUAL TO 1 IF THE NODES IN THE BIN ARE COHERENT, OR 0 OTHERWISE.
   		 S,!(OUTPUT) -- S = 1 - SUM(SM)/M, IE THE STRENGHT OF INCOHERENCE.
   		 ETA)!(OUTPUT) -- ETA = USED TO DIFFERENTIATE BETWEEN CHIMERA AND MULTICHIMERA STATES *** DOESN'T WORK PROPERLY ***
   IMPLICIT NONE
   INTEGER, PARAMETER :: BIT = 8
   INTEGER,INTENT(IN) :: N,M,TIME
   REAL(KIND=8), DIMENSION(N,TIME), INTENT(INOUT) :: X,Z
   REAL(KIND=8), INTENT(IN) :: DELTA
   REAL(KIND=8), INTENT(OUT) :: S,ETA!
   REAL(KIND=8), DIMENSION(M), INTENT(OUT) :: SDEV,SM!
   REAL(KIND=8), DIMENSION(TIME,M) :: VAR,AVE!
   REAL(KIND=8) :: TRSHLD
   INTEGER :: I,J,K
!___Transforms the standart X variables into the difference variables Z
   DO J = 1, TIME
    DO I = 1, N-1
     Z(I,J) = X(I+1,J) - X(I,J)
    ENDDO
    Z(N,J) =  X(1,J) - X(N,J)
   ENDDO
!
!--------------------------------------------------------------------
   TRSHLD = DELTA*ABS(MAXVAL(X) - MINVAL(X))
   AVE = 0.E0_BIT
   DO K=1, M
    DO J=1, TIME
     DO I = (K - 1)*N/M + 1, K*N/M
      AVE(J,K) = AVE(J,K) + Z(I,J)*M/N!    Calculates the average over the bin for a fixed time
     ENDDO
    ENDDO
   ENDDO
!
   VAR = 0.E0_BIT
   DO K = 1, M
    DO J = 1, TIME
     DO I = (K - 1)*N/M + 1, K*N/M
      VAR(J,K) = VAR(J,K) + (Z(I,J) - AVE(J,K))**2.E0_BIT!    Calculates the variance over the bin for a fixed time
     ENDDO
     VAR(J,K) = SQRT(VAR(J,K)/(N/M))
    ENDDO
   ENDDO
   SDEV = SUM(VAR,1)/TIME
!
!____________Calculates_S_and_SM_____________________________________
!
   WHERE(SDEV > TRSHLD)
    SM = 0
   ELSEWHERE
    SM = 1
   END WHERE
   S = 1 - SUM(SM)/M
!______________Calculates ETA_________________________________________
    ETA = 0.E0_BIT
    DO I = 1, M-1
     ETA = ETA + ABS(SM(I) - SM(I+1))
    ENDDO
    ETA = (ETA + ABS(SM(M) - SM(1)))/2.E0_BIT
!
   RETURN
  END SUBROUTINE SI
