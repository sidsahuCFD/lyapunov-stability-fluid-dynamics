void gramSchmidt (double ** a, double ** r, int m, int n, bool full) {
    int i, j;
    double anorm, tol = 10e-7;

    for(i = 0; i < n; i++) {
        r[i][i] = norm(a[i], m);                  // r_ii = ||a_i||

        if(r[i][i] > tol) {
            scalar_div(a[i], r[i][i], m, a[i]);   // a_i = a_i/r_ii
        }
        else if(i == 0) { // set a[0] = [1 0 0 ... 0]^T
            a[i][0] = 1;
            for(j = 1; j < m; j++) {
                a[i][j] = 0;
            }
        }
        else{ // need to choose a_i orthogonal to < a_1, ... a_{i-1} >
            for(j = 0; j < m; j++) {
                a[i][j] = -a[0][i] * a[0][j];
            }
            a[i][i] += 1;

            for(j = 1; j < i; j++) {
                scalar_sub(a[j], a[j][i], m, a[i]);
            }

            anorm = norm(a[i], m);
            scalar_div(a[i], anorm, m, a[i]);
        }

        for(j = i+1; j < n; j++) {
            r[j][i] = dot_product(a[i], a[j], m); // r_ij = a_i*a_j
            scalar_sub(a[i], r[j][i], m, a[j]);   // a_j -= r_ij a_i
        }
    }

    /* if full QR factorization requested, we choose remaining 
       columns of Q so that the m columns of Q form an 
       orthonormal set                                          */
    if(full) {
        for(; i < m; i++) {
            for(j = 0; j < m; j++) {
                    a[i][j] = -a[0][i] * a[0][j];
                }
                a[i][i] += 1;
    
                for(j = 1; j < i; j++) {
                    scalar_sub(a[j], a[j][i], m, a[i]);
                }
    
                anorm = norm(a[i], m);
                scalar_div(a[i], anorm, m, a[i]);
        }
    }
}
