!This is an example test run program 
!on how to use the Module numerical_quadrature and Functions

INCLUDE "MODULE_numerical_quadrature.f90"
INCLUDE "Example_functions.f90"		

PROGRAM test_quadrature
	USE numerical_quadrature    
	IMPLICIT NONE
	
	! Declarations
	REAL, PARAMETER :: pi = 3.1415926
	REAL, EXTERNAL :: f, g
	
	REAL :: a, b, accuracy_tolerance, value
	REAL :: smallest_subdivision = 1.E-5
	INTEGER :: error, ios
	CHARACTER ans
	
	!Open a file for results
	OPEN(UNIT=2, FILE="Results.txt", STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND', IOSTAT=ios)
	IF (ios /= 0) THEN
		PRINT*, "Error during opening of file"
		PRINT*, "Enter '1' to continue or '2' to exit program"
		READ*, ans
		IF (ans == '2') STOP
	END IF
	! Calculate integral of f on [0.1, 1.0]
	a = 1.0E-1
	b = 1.0
	accuracy_tolerance = 1.0E-2
	CALL adaptive_quadrature (f, a, b, accuracy_tolerance, &
	                         smallest_subdivision, value, error)
	! Print result or error message, as appropriate
	SELECT CASE (error)
	CASE (0)
		PRINT 100,  a, b,          &
		        accuracy_tolerance, value, -log(a)
		WRITE (UNIT=2, FMT=200)   a, b, accuracy_tolerance, value, -log(a)     
	CASE (-2)
		PRINT*, 'Failed to converge to a solution for first &
		         &problem'
		  
	CASE (-1)
		PRINT *, "Epsilon was less than or equal to zero - should &
		         & be impossible"
	END SELECT
	
	! Calculate integral of g on [0, pi/2]
	a = 0.0
	b = pi/2.0
	accuracy_tolerance = 1.0E-2
	CALL adaptive_quadrature (g, a, b, accuracy_tolerance, &
	                         smallest_subdivision, value, error)
	                         
	! Print result or error message, as appropriate
	SELECT CASE (error)
	CASE (0)
		PRINT 200,      &
		          a, b, accuracy_tolerance, value, 1.0
		WRITE (UNIT=2, FMT=200)  a, b, accuracy_tolerance, value, 1.0
	CASE (-2)
		PRINT*, 'Failed to converge to a solution for first &
		         &problem'
		  
	CASE (-1)
		PRINT *, "Epsilon was less than or equal to zero - should &
		         & be impossible"
	END SELECT
100 FORMAT (//1X,'Value of integral of x**(-1) from ', E9.2, &
		        ' to ', E9.1/ 1x, 'with accuracy tolerance ', &
		        F14.6/, 1x, 'is ', F14.6/ 1x,                   &
		        'Correct answer is ', F14.6//)
200 FORMAT (//1x, 'Value of integral cos(x) from ', F5.1 &
		              ' to ', F10.6/                         &
		          1x, 'with accuracy tolerance ', F14.6/    &
		          1x, 'is ', F14.6/                         &
		          1x, 'Correct answer is ', F14.6//)
STOP
END PROGRAM test_quadrature 
