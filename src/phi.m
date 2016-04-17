function [phi, phi_h] = phi( eMat, w, q=[], p=[] )
    
    # w, q, p = preference parameters given as row vectors

	# Precision for rounding of results
	eps = 1e-6;

	[n,m] = size( eMat );

	# Possibly adapting to alternative preference parameter "forms"
	[w,q,p] = check_prefs(w,q,p);
    
    res_h = zeros( n, m );

    # PROMETHEE piecewise linear preference function (the usual one)
    PFn = @(D, q, p) max(0, min(1, (D - q)./(p - q) ));
    
    for h = 1:m
        d = kron(ones(1,n),eMat(:,h));
        D = d - d';   #'
        P = PFn( D, q(h), p(h) );
        res_h(:,h) = sum(P')' - sum(P)';
    end;

    res_h = res_h ./ (n-1);

    phi_h = round(res_h./eps).*eps;
    phi = round(( res_h * w )./eps).*eps;

endfunction
