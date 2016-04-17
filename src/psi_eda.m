function [psi, psi_h] = psi_eda(eMat, w, q=[], p=[] )
	# w, q, p = preference parameters given as row vectors

	# Precision for rounding of results
	eps = 1e-6;

	[n,m] = size( eMat );

	# Possibly adapting to alternative preference parameter "forms"
	[w,q,p] = check_prefs(w,q,p);

	# Threshold difference on all criteria
	deltaThres = kron( ones( n, 1 ), (p-q)' ); #'

	# Altered evaluations on all criteria:
	ymq = @(y) max(0, y - kron(ones(rows(y),1),q'));   #' ymq = value minus q threshold
	ypq = @(y) min(1, y + kron(ones(rows(y),1),q'));   #' ypq = value plus q threshold
	ymp = @(y) max(0, y - kron(ones(rows(y),1),p'));   #' ymp = value minus p threshold
	ypp = @(y) min(1, y + kron(ones(rows(y),1),p'));   #' ypp = value plus p threshold

	deltaP = @(y,zeta) max( 0, min( 1, 1-(zeta-ymp(y))./deltaThres ) ) + max( -1, min( 0, -(zeta-ypq(y))./deltaThres ) );
	dPh = deltaP(eMat,1);

	# Step for sampling set to have minimum 1000 and maximum 100000
	# Otherwise, 1000*n steps is fine.
	# Here we have to balance **precision** and **calculation time**...
	meshStep = (1/min(max(1000,100*n),100000));

	yMesh = (0:meshStep:1)';  #'

	# Compute the evaluation CDF for each criterion
	CDF_eMat = emp_cdf( yMesh, eMat );

	# **Integrate** evaluation CDF's for each criterion    #'
	intCDF_eMat = cumsum( CDF_eMat.*meshStep, 1 );

	iCDF_eMat = zeros(n, m);
	iintCDF_eMat_ymq = zeros(n, m);
	iintCDF_eMat_ypq = zeros(n, m);
	iintCDF_eMat_ymp = zeros(n, m);
	iintCDF_eMat_ypp = zeros(n, m);

	% ToDo: Time gets lost in the loop (but should not be too critical for few dimensions)
	for h = 1:m
		iCDF_eMat(:,h) = interp1( yMesh, CDF_eMat(:,h), eMat(:,h) );
		iintCDF_eMat_ymp(:,h) = interp1( yMesh, intCDF_eMat(:,h), ymp(eMat)(:,h) );
		iintCDF_eMat_ymq(:,h) = interp1( yMesh, intCDF_eMat(:,h), ymq(eMat)(:,h) );
		iintCDF_eMat_ypq(:,h) = interp1( yMesh, intCDF_eMat(:,h), ypq(eMat)(:,h) );
		iintCDF_eMat_ypp(:,h) = interp1( yMesh, intCDF_eMat(:,h), ypp(eMat)(:,h) );
	endfor

	ucEda = dPh + (( iintCDF_eMat_ymq-iintCDF_eMat_ymp + iintCDF_eMat_ypp-iintCDF_eMat_ypq ) ./ deltaThres);
	
	psi_h = round(ucEda./eps).*eps;
	psi = round(( ucEda * w )./eps).*eps;

end