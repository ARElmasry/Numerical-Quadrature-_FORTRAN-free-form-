MODULE numerical_quadrature
	IMPLICIT NONE
	PRIVATE
	PUBLIC :: adaptive_quadrature
	
CONTAINS 
	
	
	SUBROUTINE adaptive_quadrature (f, a, b, epsilon, &
	                                subdivide_limit, answer, error)
		! This subroutine integrates the function f from a to b
		! using an adaptive method based on the trapezoidal rule.
		! epsilon is the user-specified error tolerance.
		! subdivide_limit is a user-specified samallest interval
		! size to use.
		! answer is the calculated answer success.
	
		! Dummy arguments
		REAL, EXTERNAL :: f
		REAL, INTENT (IN) :: a, b, epsilon, subdivide_limit
		REAL, INTENT (OUT) :: answer
		INTEGER, INTENT (OUT) :: error
	
		! Validity check
		IF (epsilon <= 0.0) THEN
			error = -1
			RETURN
		END IF
		IF (a < b) THEN
			CALL adap_quad (f, a, b, f(a), f(b), subdivide_limit, &
			                epsilon/(b-a), answer, error)        
		ELSE IF (a > b) THEN 
			CALL adap_quad (f, b, a, f(b), f(a), subdivide_limit, & 
			                epsilon/(a-b), answer, error)
			IF (error == 0) answer = -answer
		ELSE
			error = 0
			answer = 0.0
		END IF
	END SUBROUTINE adaptive_quadrature
	
	RECURSIVE SUBROUTINE adap_quad (f, xl, xu, fl, fu,        &
	                                lower, delta, answer, error)
		! This subrouitne performs an adaptive numerical
		! quadratue using the trapezoidal rule.
	
		! Dummy arguments
		REAL, EXTERNAL :: f
		REAL, INTENT (IN) :: xl, xu, fl, fu, lower, delta
		REAL, INTENT (OUT) :: answer
		INTEGER, INTENT (OUT) :: error
		
		! Local variables
		REAL :: h, t, c, xm, fm, e, ans1, ans2
		
		h = xu - xl
		IF (ABS(h) < lower) THEN 
			! Interval has become too small
			error = -2
			answer = HUGE(answer)
			RETURN
		END IF
		t = h*(fl + fu)/2.0
		xm = xl + h/2.0
		fm = f(xm)
		c = h*(fl + 2.0*fm + fu)/4.0
		e = 4.0*(c - t)/3.0
		IF (ABS(e) <= delta*h) THEN 
			! Trapezoidal  rule has achieved required accuracy
			! The PRINT statement is only for during development
			! It will be removed when code is certified as 
			! functional
			PRINT '(1X, ''Interval Used ('', E12.4, '','', E12.4, '')'', &
			        3X, ''h ='', E12.4)', xl, xu, xu-xl
			error = 0
			answer = t
		ELSE
			! Subdivide the interval
			CALL adap_quad (f, xl, xm, fl, fm, lower, delta, &
			                ans1, error)
			IF (error /= 0 ) RETURN
			CALL adap_quad (f, xm, xu, fm, fu, lower, delta, &
			                ans2, error)
			IF (error /= 0 ) RETURN
		answer = ans1 + ans2
		END IF 
	END SUBROUTINE adap_quad
END MODULE numerical_quadrature
