function [psi, psi_h] = psi_pla(eMat, w, q=[], p=[] )
	# w, q, p = preference parameters given as row vectors

	# Precision for rounding of results
	eps = 1e-6;

	[n,m] = size( eMat );

	# Possibly adapting to alternative preference parameter "forms"
	[w,q,p] = check_prefs(w,q,p);

    lambda = ((q + p)./2)'; #'

    ucPlaFn = @(y,lambda) min( max( ( y + kron(lambda,ones(rows(y),1)) - 1 ), ( 2*y - 1 ) ), ( y - kron(lambda,ones(rows(y),1)) ) );
	ucPla = ucPlaFn(eMat,lambda);

	psi_h = round(ucPla./eps).*eps;
	psi = round(( ucPla * w )./eps).*eps;

end