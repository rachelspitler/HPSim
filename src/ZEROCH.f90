    REAL FUNCTION ZEROCH(AX,F,TOL1,TOL2,DX,FB,IERROR)
    !
    LOGICAL FIRST
    REAL AX,FAX,BX,FBX,F,TOL1,TOL2
    REAL A,B,C,D,E,EPS,FA,FB,FC,TOLX,TOLF,XM,P,Q,R,S
    EXTERNAL F
    DATA FIRST / .TRUE. /
    !
    CALL GUESS3(AX,FAX,BX,FBX,F,DX,TOL2,IERROR)
    IF (IERROR .NE. 0) THEN
        RETURN
    END IF
    IF (FIRST) THEN
        EPS = 1.0
        ONE = 1.0   !VL
        DO WHILE (ONE .GT. 1.0)
            EPS = EPS/2.0
            ONE = 1.0 + EPS0
        END DO 

        FIRST = .FALSE.
    END IF
    !

    A = AX
    B = BX
    FA = FAX
    FB = FBX
    !

    C = A
    FC = FA
    D = B-A
    E = D

    DO WHILE (.TRUE.)

        IF (ABS(FC) .LT. ABS(FB)) THEN
            A = B
            B = C
            C = A
            FA = FB
            FB = FC
            FC = FA        
        END IF

        !40	TOLX = 2.0*EPS*ABS(B) + 0.5*TOL1	!ISI - 05/31/05
        !	TOLF = 2.0*EPS*ABS(FB) + 0.5*TOL2	!ISI - 05/31/05
        !	XM = 0.5*(C - B)				    !ISI - 05/31/05

        !VL: Previously: 40      TOLX = TOL1		 !ISI - 05/31/05 ! all GOTO 40 statements eliminated ....
        TOLX = TOL1		 !ISI - 05/31/05
        TOLF = TOL2		 !ISI - 05/31/05
        XM = 0.5*(C - B) !ISI - 05/31/05 

        IF (ABS(XM) .LE. TOLX) THEN
            EXIT
        END IF
        IF (ABS(FB) .LE. TOLF) THEN
            EXIT
        END IF
        
        IF ((ABS(E) .LT. TOLX) .OR. (ABS(FA) .LE. ABS(FB))) THEN
            D = XM
            E = D
        ELSE

            IF (A .NE. C) THEN

                Q = FA/FC
                R = FB/FC
                S = FB/FA
                P = S*(2.0*XM*Q*(Q - R) - (B - A)*(R - 1.0))
                Q = (Q - 1.0)*(R - 1.0)*(S - 1.0)            

            ELSE

                S = FB/FA
                P = 2.0*XM*S
                Q = 1.0 - S

            END IF

            IF (P .GT. 0.0) THEN
                Q = -Q
            END IF
            P = ABS(P)

            IF (((2.0*P) .GE. (3.0*XM*Q - ABS(TOLX*Q))) .OR. (P .GE. ABS(0.5*E*Q))) THEN
                D = XM
                E = D                
            ELSE
                E = D
                D = P/Q          
            END IF

        END IF

        A = B
        FA = FB
        IF (ABS(D) .GT. TOLX) THEN
            B = B + D
        END IF
        IF (ABS(D) .LE. TOLX) THEN
            B = B + SIGN(TOLX,XM)
        END IF
        FB = F(B,IERR)
        IF ((FB*(FC/ABS(FC))) .GT. 0.) THEN
            C = A
            FC = FA
            D = B-A
            E = D
        END IF

    END DO

    !
    !	SET ERROR CODES
    !	  IERROR = 0, NORMAL RETURN
    !	  IERROR = 1, TOLERANCE ON INDEPENDENT VARIABLE EXCEEDED
    !	  IERROR = 2, TOLERANCE ON FUNCTION VALUE EXCEEDED
    !	  IERROR = 3, TOLERANCES ON INDEPENDENT VARIABLE AND FUNCTION VALUE
    !		      EXCEEDED
    !
    ZEROCH = B
    IERROR = 0
    IF (ABS(XM) .GT. TOLX) THEN
        IERROR = 1
    END IF
    IF (ABS(FB).GT.TOLF) THEN
        IERROR = IERROR + 2
    END IF
    RETURN
    !
    !	COMPUTE MACHINE PRECISION "EPS"
    !
    END
