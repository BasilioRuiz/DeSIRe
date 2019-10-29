	FUNCTION FINT (X,Y1,Y2,Y3,X1,X2,X3)
	implicit real*4 (a-h,o-z)
	D1=X-X1
	D2=X-X2
	D3=X-X3
	D12=X1-X2
	D13=X1-X3
	D23=X2-X3
	FINT=Y1*D2*D3/(D12*D13)-Y2*D1*D3/(D12*D23)+Y3*D1*D2/(D13*D23)
	RETURN
	END
	FUNCTION FINT_db (X,Y1,Y2,Y3,X1,X2,X3)
	implicit real*8 (a-h,o-z)
	real*8 FINT_db
	D1=X-X1
	D2=X-X2
	D3=X-X3
	D12=X1-X2
	D13=X1-X3
	D23=X2-X3
	FINT_db=Y1*D2*D3/(D12*D13)-Y2*D1*D3/(D12*D23)+Y3*D1*D2/(D13*D23)
	RETURN
	END
