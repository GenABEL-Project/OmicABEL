/*
 * Copyright (c) 2010-2012, Diego Fabregat-Traver and Paolo Bientinesi.
 * All rights reserved.
 *
 * This file is part of OmicABEL.
 * 
 * OmicABEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * OmicABEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with OmicABEL. If not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 * Coded by:
 *   Diego Fabregat-Traver (fabregat@aices.rwth-aachen.de)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "options.h"
#include "blas.h"
#include "lapack.h"
#include "ooc_BLAS.h"
#include "wrappers.h"
#include "utils.h"
#include "databel.h"
#include "statistics.h"
#include "optimization.h"
#include "REML.h"

// Auxiliar functions
static void eigen_mr3(int n, double *Phi, double *Z, double *W);
static void premultiplyXL( int n, int wXL, double *Z, double *XL, double *ZtXL );


void estimates_chol( FGLS_config_t *cf )
{
	estimates_chol_t data;
	int j, k;
#if DEBUG
	FILE *f_ests;
#endif
	
	// Dimensions
	data.n   = cf->n;
	data.wXL = cf->wXL;
	// Load data
	data.Phi = cf->Phi;
	data.X = cf->XL;

	set_multi_threaded_BLAS(cf->num_threads);
	data.Y = fgls_malloc( data.n * sizeof(double) );
	data.beta = fgls_malloc( cf->wXL * sizeof(double) );
	for ( j = 0; j < cf->t; j++ )
	{
		// Should Y_fp
		sync_read( data.Y, cf->Y, cf->n, cf->n * j );
		// Sanity
		/*checkNoNans(cf->n, data.Y, "Y");*/
		average( data.Y, cf->n, 1, cf->threshold, "TRAIT", &cf->Y_fvi->fvi_data[(cf->n+j)*NAMELENGTH], NAMELENGTH, 1 );
		/*printf("var(Y): %.15e\n", variance(data.Y, cf->n));*/

		data.sigma = variance(data.Y, data.n);

		// [a,b] = [0,0.98], tol = 0.01
		minimize( 0.0, 0.99, 1e-8, &polygenic_REML_logLik_chol_wrapper, (void*)&data );
#if 0 //DEBUG
		printf("h: %.15e - sigma: %.15e - res_sigma: %.15e\n", data.h, data.sigma, data.res_sigma);
#endif
		cf->h2[j] = data.h;
		cf->sigma2[j] = data.sigma;
		cf->res_sigma2[j] = data.res_sigma;
		for ( k = 0; k < cf->wXL; k++ )
			cf->beta_ests[cf->wXL*j+k] = data.beta[k];
	}
#if DEBUG
	// Store data - For testing purposes
/*	f_ests = fgls_fopen( cf->ests_path, "wb" );
	sync_write( cf->ests, f_ests, 3 * cf->t, 0 );
	fclose( f_ests );*/
#endif
	free( data.Y );
	free( data.beta );
}

double polygenic_REML_logLik_chol (
        int n, int widthXL,
        double *Phi, double *X, double *Y, 
        double sigma2, double h2,
        double *loglik, double *res_sigma, double *beta
)
{
    double //*beta, //
		   *initX, *initY, // Copy of the initial X and Y
           alpha, gamma, // to scale Phi
           *K, // readibility
           *V, det, // temporary
           ZERO = 0.0, // BLAS / LAPACK
           ONE = 1.0,
           MINUS_ONE = -1.0;

    int    wXL = widthXL,
           info,
           nn = n*n,
           iONE = 1,
           i; //, j, k;

	/*beta  = (double *) fgls_malloc ( wXL * sizeof(double) );*/
    initX = (double *) fgls_malloc ( n * wXL * sizeof(double) );
    initY = (double *) fgls_malloc ( n * sizeof(double) );
    V     = (double *) fgls_malloc ( wXL * wXL * sizeof(double) );

    memcpy( initX, X, n * wXL * sizeof(double) );
    memcpy( initY, Y, n * sizeof(double) );

    /* 1) K := h^2 Phi - (1-h^2) I */
    K = Phi;
	alpha = h2; // * sigma2;
	gamma = (1 - h2); // * sigma2;
    dscal_(&nn, &alpha, K, &iONE);
    for ( i = 0; i < n; i++ )
        K[i*n + i] = K[i*n + i] + gamma;

    /* 2) L * L' = K */
	int *ipiv = (int *) fgls_malloc (n * sizeof(int));
	int lwork = n * 192; // nb
	double *work = (double *) fgls_malloc ( lwork * sizeof(double));
    dsysv_(LOWER, &n, &widthXL, K, &n, ipiv, X, &n, work, &lwork, &info);
	/*dpotrf_(LOWER, &n, K, &n, &info);*/
    if (info != 0)
    {
        fprintf(stderr, __FILE__ ": dsysv_ inv(M)*X - The matrix is singular(%d)\n", info);
        exit(-1);
    }

    /* 3) X := inv(L) X */
	/*dtrsm_(LEFT, LOWER, NO_TRANS, NON_UNIT, &n, &wXL, &ONE, K, &n, X, &n);*/

    /* 4) y := inv(L) y */
	/*dtrsv_(LOWER, NO_TRANS, NON_UNIT, &n, K, &n, Y, &iONE);*/
	/*dtrsv_(LOWER, NO_TRANS, NON_UNIT, &n, K, &n, &Y[j * n], &iONE);*/

    /* 5) beta := X' * y */
    dgemv_(TRANS, &n, &widthXL, &ONE, X, &n, Y, &iONE, &ZERO, beta, &iONE);

    /* 6) V := X' * X */
	/*dsyrk_(LOWER, TRANS, &wXL, &n, &ONE, X, &n, &ZERO, V, &wXL);*/
    dgemm_(TRANS, NO_TRANS, &wXL, &wXL, &n, &ONE, initX, &n, X, &n, &ZERO, V, &wXL);

    /* Chol(V) */
	/*dpotrf_(LOWER, &wXL, V, &wXL, &info);*/
    dsysv_(LOWER, &wXL, &iONE, V, &wXL, ipiv, beta, &wXL, work, &lwork, &info);
    if (info != 0)
    {
        fprintf(stderr, __FILE__ ": [ERROR] inv(V) y failed - info: %d\n", info);
        exit(-1);
    }
    /* beta := inv( V ) beta */
	/*dtrsv_(LOWER, NO_TRANS, NON_UNIT, &wXL, V, &wXL, beta, &iONE);*/
	/*dtrsv_(LOWER,    TRANS, NON_UNIT, &wXL, V, &wXL, beta, &iONE);*/

    // y - X beta
    dgemv_(NO_TRANS, 
            &n, &wXL, 
            &MINUS_ONE, initX, &n, beta, &iONE,
            &ONE, initY, &iONE);
    // initY stores now "y - X beta"

    // residual sigma and loglik
	*res_sigma = variance( initY, n );

    memcpy( Y, initY, n * sizeof(double) );
	dsytrs_( LOWER, &n, &iONE, K, &n, ipiv, initY, &n, &info );
    if (info != 0)
    {
        fprintf(stderr, __FILE__ ": [ERROR] inv(K) y - info: %d\n", info);
        exit(-1);
    }
	/*dtrsv_(LOWER, NO_TRANS, NON_UNIT, */
	/*&n, K, &n, initY, &iONE);*/
    *loglik = ddot_(&n, Y, &iONE, initY, &iONE) / (*res_sigma);

	/*
	 * In principle, the det is just product of diagonal entries (D)
	 * D can be 2x2 block diagonal
	 * Have to check this case and compute accordingly
	 */
    det = 0;
    for ( i = 0; i < n; i++ )
	{
		/*det += log(K[i*n + i]);*/
		if ( ipiv[i] > 0 )
        	det += log(K[i*n + i]);
		else if ( ipiv[i] < 0 && i > 0 && ipiv[i-1] == ipiv[i] )
			det += log( K[(i-1)*n + (i-1)] * K[i*n + i] - 
					    K[(i-1)*n + i] * K[(i-1)*n + i] );
	}
	/*det = det * 2;*/
	*loglik = det + (*loglik) + n * log(*res_sigma);

    // Clean up
	/*free(beta);*/
    free(initX);
    free(initY);
    free(V);

	free( ipiv );
	free( work );
	
	return *loglik;
}

double polygenic_REML_logLik_chol_wrapper ( double h, void *data )
{
	estimates_chol_t *chol_t = (estimates_chol_t*)data;
	size_t n   = chol_t->n, 
	       wXL = chol_t->wXL;
	double *Phi, //  = chol_t->Phi,
		   *Y, //    = chol_t->Y,
		   *X; //    = chol_t->X;
	double sigma = chol_t->sigma;
	double loglik;

	Phi = fgls_malloc( n * n * sizeof(double) );
	X = fgls_malloc( n * wXL * sizeof(double) );
	Y = fgls_malloc( n * sizeof(double) );
	memcpy( Phi, chol_t->Phi, n * n * sizeof(double) );
	memcpy( X, chol_t->X, n * wXL * sizeof(double) );
	memcpy( Y, chol_t->Y, n * sizeof(double) );

	polygenic_REML_logLik_chol (
	        n, wXL,
	        Phi, X, Y,
			sigma, h,
	        &loglik, &chol_t->res_sigma, chol_t->beta
	);
	
	chol_t->h = h;

	free( Phi );
	free( X );
	free( Y );

	return loglik;
}


/*
 * Eigendecomposition-based REML
 */
void estimates_eig( FGLS_config_t *cf )
{
	estimates_eig_t data;
	int j, k;
#if DEBUG
	FILE *f_ests;
#endif
	
	// Dimensions
	data.n   = cf->n;
	data.wXL = cf->wXL;
	// Load data
	data.Phi = cf->Phi;
	// Compute eigenvectors and eigenvalues in data.Z and data.W
	cf->Z = data.Z = fgls_malloc( data.n * data.n * sizeof(double) );
	cf->W = data.W = fgls_malloc( data.n * sizeof(double) );
	set_multi_threaded_BLAS(cf->num_threads); // Probably can be refined
	eigen_mr3(data.n, data.Phi, data.Z, data.W);
	// Phi is not used anymore
	free( cf->Phi ); cf->Phi = NULL;
	// Load XL
	data.X = cf->XL;
	cf->ZtXL = data.ZtX = fgls_malloc( data.n * data.wXL * sizeof(double) );
	// Premultiply ZtXL
	premultiplyXL( data.n, data.wXL, data.Z, data.X, data.ZtX );
	// Premultiply ZtY
	ooc_gemm( data.n, cf->t, cf->ooc_b, data.Z, cf->Y_data_path, cf->ZtY_path,
			cf->threshold, "TRAIT", &cf->Y_fvi->fvi_data[cf->n*NAMELENGTH], NAMELENGTH);

	// Iterate over phenotypes (Y)
	data.Y = fgls_malloc( data.n * sizeof(double) );
	data.ZtY = fgls_malloc( data.n * sizeof(double) );
	data.beta = fgls_malloc( cf->wXL * sizeof(double) );
	cf->ZtY = fgls_fopen( cf->ZtY_path, "rb" );
	for ( j = 0; j < cf->t; j++ )
	{
		sync_read( data.Y, cf->Y, cf->n, cf->n * j );
		// Sanity
		average( data.Y, cf->n, 1, cf->threshold, "TRAIT", &cf->Y_fvi->fvi_data[(cf->n+j)*NAMELENGTH], NAMELENGTH, 0 );

		sync_read( data.ZtY, cf->ZtY, cf->n, cf->n * j );
		data.sigma = variance(data.Y, data.n);

		// [a,b] = [0,0.98], tol = 0.01
		minimize( 0.0, 0.99, 1e-8, &polygenic_REML_logLik_eig_wrapper, (void*)&data );
#if 0 //DEBUG
		printf("h: %.15e - sigma: %.15e - res_sigma: %.15e\n", data.h, data.sigma, data.res_sigma);
#endif

		cf->h2[j] = data.h;
		cf->sigma2[j] = data.sigma;
		cf->res_sigma2[j] = data.res_sigma;
		for ( k = 0; k < cf->wXL; k++ )
			cf->beta_ests[cf->wXL*j+k] = data.beta[k];
	}

#if DEBUG
	// Store data - For testing purposes
/*	f_ests = fgls_fopen( cf->ests_path, "wb" );
	sync_write( cf->ests, f_ests, 3 * cf->t, 0 );
	fclose( f_ests );*/
#endif
	// Encapsulate (maybe reuse)
	free( data.Y );
	free( data.beta );
	free( data.ZtY );
	// Free Phi? Not used anymore, plus wrong data (from the eigensolver overwritting it)
}

double polygenic_REML_logLik_eig (
        int n, int widthXL,
		double *X, double *Y,
        double *Z, double *W, double *ZtX, double *ZtY,
		double sigma, double h2,
        double *loglik, double *res_sigma, double *beta
)
{
    double alpha, gamma, // to scale W
//           *beta,
           *ZtY_upd, *ZtX_upd, *YmXB, *V, *D, det, // temporary
           ZERO = 0.0,
           ONE = 1.0,
           MINUS_ONE = -1.0;

    int    wXL = widthXL,
           iONE = 1, info,
           i, j; //, k;

	/*printf("n: %d\n", n);*/
	/*printf("widthXL: %d\n", widthXL);*/
	/*checkNoNans(widthXL*n, X, "X");*/
	/*checkNoNans(n, Y, "Y");*/
	/*checkNoNans(n*n, Z, "Z");*/
	/*checkNoNans(n, W, "W");*/
	/*checkNoNans(n, ZtY, "ZtY");*/
	/*checkNoNans(widthXL * n, ZtX, "ZtX");*/
	/*printf("sigma: %f\n", sigma);*/
	/*printf("h2: %f\n", h2);*/
	/*beta    = (double *) fgls_malloc ( wXL * sizeof(double) );*/
    V       = (double *) fgls_malloc ( wXL * wXL * sizeof(double) );
	ZtY_upd = (double *) fgls_malloc ( n * sizeof(double) );
	ZtX_upd = (double *) fgls_malloc ( n * wXL * sizeof(double) );
	D = (double *) fgls_malloc ( n * sizeof(double) );
	YmXB = (double *) fgls_malloc ( n * sizeof(double) );

	/* W = sqrt( inv ( W ) ) */
    alpha = h2; // * sigma2;
    gamma = (1 - h2); // * sigma2;
	det = 0.0;
	/*checkNoNans(n, W, "\n[ERROR] W\n");*/
	for ( i = 0; i < n; i++ )
	{
		D[i] = alpha * W[i] + gamma;
		det += log(fabs(D[i]));
		D[i] = 1.0 / D[i];
	}

	/* W * ZtX */
	/* W * ZtY */
#if 0
	/*printf("ZtX: %x\n", ZtX);*/
	for (i = 0; i < wXL; i++)
	{
		for (j = 0; j < wXL; j++)
			printf("%.16e ", ZtX[j*n + i]);
		printf("\n");
	}
	printf("\n");
	for (i = 0; i < wXL; i++)
		printf("%.16e ", D[i]);
	printf("\n\n");
#endif
	/*checkNoNans(n, D, "\n[ERROR] D\n");*/
	// inv( X' inv(M) X) X' inv(M) y
	// inv(M) -> inv( h Phi + (1-h) I )
	// inv(M) -> inv( h Z W Z' + (1-h) I )
	// inv(M) -> Z inv( h W + (1-h) I ) Z'
	// inv( X' Z D Z' X) X' Z D Z' y
	// inv( X' Z D Z' X) X' Z D Z' y
	// inv(V) y
	for (i = 0; i < n; i++ )
	{
		for (j = 0; j < wXL; j++ )
			ZtX_upd[j*n + i] = ZtX[j*n + i] * D[i];
		/*ZtY_upd[i] = ZtY[i] * D[i];*/
		ZtY_upd[i] = ZtY[i];
	}

    /* 5) beta := X' * y */
	/*checkNoNans(n, ZtY_upd, "\n[ERROR] ZtY_upd\n");*/
    dgemv_(TRANS, &n, &widthXL, &ONE, ZtX_upd, &n, ZtY_upd, &iONE, &ZERO, beta, &iONE);

    /* 6) V := X' * X */
#if 0
	for (i = 0; i < wXL; i++)
	{
		for (j = 0; j < wXL; j++)
			printf("%.16e ", ZtX_upd[j*n + i]);
		printf("\n");
	}
	printf("\n");
#endif
	/*checkNoNans(wXL*n, ZtX_upd, "\n[ERROR] ZtX_upd\n");*/
    dgemm_(TRANS, NO_TRANS, &wXL, &wXL, &n, &ONE, ZtX, &n, ZtX_upd, &n, &ZERO, V, &wXL);
	/*dsyrk_(LOWER, TRANS, &wXL, &n, &ONE, ZtX_upd, &n, &ZERO, V, &wXL);*/

#if 0
	printf("wXL: %d\n", wXL);
	for (i = 0; i < wXL; i++)
	{
		for (j = 0; j < wXL; j++)
			printf("%.16e ", V[j*wXL + i]);
		printf("\n");
	}
	printf("\n");
#endif
	/*checkNoNans(wXL*wXL, V, "V");*/
	int *ipiv = (int *) fgls_malloc (wXL * sizeof(int));
	int lwork = wXL * 192; // nb
	double *work = (double *) fgls_malloc ( lwork * sizeof(double));
    dsysv_(LOWER, &wXL, &iONE, V, &wXL, ipiv, beta, &wXL, work, &lwork, &info);
	/*dpotrf_(LOWER, &wXL, V, &wXL, &info);*/
    if (info != 0)
    {
        fprintf(stderr, __FILE__ ": [ERROR] inv(V) y failed - info: %d\n", info);
        exit(-1);
    }
    /* beta := inv( V ) beta */
	/*dtrsv_(LOWER, NO_TRANS, NON_UNIT, &wXL, V, &wXL, beta, &iONE);*/
	/*dtrsv_(LOWER,    TRANS, NON_UNIT, &wXL, V, &wXL, beta, &iONE);*/

    // y - X beta
	memcpy( YmXB, Y, n * sizeof(double) );
    dgemv_(NO_TRANS, 
            &n, &wXL, 
            &MINUS_ONE, X, &n, beta, &iONE,
            &ONE, YmXB, &iONE);

/*	memcpy( YmXB, ZtY, n * sizeof(double) );
    dgemv_(NO_TRANS, 
            &n, &wXL, 
            &MINUS_ONE, ZtX, &n, beta, &iONE,
            &ONE, YmXB, &iONE); // YmXB == y_0
*/
    // residual sigma and loglik
	*res_sigma = variance( YmXB, n );

	// loglik = a + b
	//  a -> log(det(M))
	//  b -> qt' inv(M) qt
	//  qt = YmXB

	// Z' * YmXB
	// 
	
	//  qt = YmXB
	//  b -> qt' 1/sigma Z D Z' qt
	//  b -> 1/sigma (ZtY' D ZtY)
    dgemv_(TRANS, 
            &n, &n, 
            &ONE, Z, &n, YmXB, &iONE,
            &ZERO, ZtY_upd, &iONE);
	// YmXB' * inv( M ) * YmXB
	*loglik = 0.0;
	for (i = 0; i < n; i++ )
	//	*loglik += YmXB[i] * D[i] * YmXB[i];
		*loglik += ZtY_upd[i] * D[i] * ZtY_upd[i];
	/**loglik += ZtY_upd[i] * D[i] * D[i] * ZtY_upd[i];*/
	/**res_sigma = (*loglik) / (n - wXL);*/
	*loglik = (*loglik) / (*res_sigma); // <- b
    *loglik = (*loglik) + det + n * log(*res_sigma);

    // Clean up
	/*free( beta );*/
	free( V );
	free( ZtY_upd );
	free( ZtX_upd );
	free( YmXB );
	free( D );

	free( ipiv );
	free( work );
	
	return *loglik;
}


double polygenic_REML_logLik_eig_wrapper ( double h, void *data )
{
	estimates_eig_t *eig_t = (estimates_eig_t*)data;
	size_t n   = eig_t->n, 
	       wXL = eig_t->wXL;
	double *X = eig_t->X,
		   *Y = eig_t->Y,
	       *Z   = eig_t->Z,
		   *W   = eig_t->W,
	       *ZtX = eig_t->ZtX,
		   *ZtY = eig_t->ZtY;
	double sigma = eig_t->sigma;
	double loglik;

	polygenic_REML_logLik_eig (
	        n, wXL,
	        X, Y, Z, W, ZtX, ZtY,
            sigma, h,
	        &loglik, &eig_t->res_sigma, eig_t->beta
	);

	eig_t->h = h;

/*	printf("Loglik: %.15e\n", loglik);
	int i;
	for ( i = 0; i < wXL; i++ )
		printf("Beta[%d]: %.15e\n", i, eig_t->beta[i]);
*/
	return loglik;
}

/*
 * MR3 eigendecomposition
 */
void eigen_mr3(int n, double *Phi, double *Z, double *W)
{
    int nb = 192;
    int idummy, nCompPairs, *isuppz, *iwork, info,
        lwork = n * (nb + 6),
        liwork = 10 * n;
    double ddummy = -1.0, *work;

    work = (double *) fgls_malloc ( lwork * sizeof(double) );
    iwork   = (int *) fgls_malloc ( liwork * sizeof(int) );
    isuppz  = (int *) fgls_malloc ( 2 * n * sizeof(int) );

    dsyevr_("V", "A", "L", &n, Phi, &n, 
			&ddummy, &ddummy, &idummy, &idummy, &ddummy, 
			&nCompPairs, W, Z, &n, isuppz, 
			work, &lwork, iwork, &liwork, &info);
    if (info != 0)
    {
        fprintf(stderr, "Error factoring Phi\n");
        exit(-1);
    }

	int i;
	int negs = 0;
	for (i = 0; i < n; i++)
	{
		if (W[i] < 0)
			negs++;
		/*else*/
		/*break;*/
		if (i > 0 && W[i-1] > W[i])
			printf("Not sorted!!!!!!!!!!!!!\n");
	}
	/*printf("# of negative eigenvalues: %d\n", negs);*/
	/*printf("First eigenvalue: %.15e\n", W[0]);*/
	/*printf("Lowest eigen values [-1]: %.15e\n", W[n-1]);*/
	/*printf("Lowest eigen values [-2]: %.15e\n", W[n-2]);*/

	free( work );
	free( iwork );
	free( isuppz );
}

// Assumes leading dimension = n for every matrix
void premultiplyXL( 
		int n, int wXL,
		double *Z, // Eigenvectors
		double *XL,
		double *ZtXL
)
{
	int ldZ, ldXL, ldZtXL;
	
	double ZERO = 0.0,
	       ONE  = 1.0;

	ldZ = ldXL = ldZtXL = n;

	// Z' * X
	dgemm_(
			TRANS, NO_TRANS, 
			&n, &wXL, &n, 
			&ONE, Z, &ldZ, XL, &ldXL, 
			&ZERO, ZtXL, &ldZtXL
	);
}
