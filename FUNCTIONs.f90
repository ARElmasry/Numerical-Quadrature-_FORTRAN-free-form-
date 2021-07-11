! These are two functions used as an example
!f(x)=1.0/x
!f(x)=cos(x)

REAL FUNCTION f(x)
	IMPLICIT NONE
	REAL, INTENT(IN) :: x
	f = 1.0/x
	RETURN
END FUNCTION f

FUNCTION g(x)
	IMPLICIT NONE
	REAL, INTENT(IN) :: x
	REAL :: g
	g = cos(x)
	RETURN
END FUNCTION g
