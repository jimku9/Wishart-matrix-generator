// ************************************************************************************ //
//  WishPack - Density functions and Generators for Wishart and Inverse Wishart         //
//  Yu-Cheng Ku & Peter Bloomfield, Dept of Stat, North Carolina State University       //
//                                                                                      //
//  ranwish    - Random Wishart matrix generator                                        //
//  raninvwish - Random inverse Wishart matrix generator                                //
//  denwish    - Density function from the Wishart distribution                         //
//  deninvwish - Density function from the inverse Wishart distribution                 //
//                                                                                      //
//  Arguments:                                                                          //
//  df is the degrees of freedom.                                                       //
//  Sc is the scale matrix of the Wishart distribution.                                 //
//  W  is the positive definite matrix to be evaluatated.                               //
//                                                                                      //
//  References:                                                                         //
//  1. Anderson, T. W., (2003), An Introduction to Statitistcal Multivariate Analysis.  //
//  2. Srivastava, M. S., (2003), Singular Wishart and Multivariate Beta Distributions, //
//     Annals of Statistics.                                                            //
// ************************************************************************************ //

// Generating Random Wishart Matrices
ranwish(const df, const Sc)
{
	decl Scatter = 0 * Sc;
	decl p = sizer(Sc);	
	decl A = choleski(Sc);
	if(A == 0) oxrunerror("The scale matrix is not a symmetric and positive definte matrix", 1);

	if(df >= p)
	{
		decl B = 0 * Sc;
		decl i,j;
		for(j=0; j<p; ++j)
		{
			B[j][j] = sqrt( rangamma(1, 1, (df-j)/2, 1/2) );
		}

		for(i=1; i<p; ++i)
		{
			for(j=0; j<i; ++j)
			{
				B[i][j] = rann(1,1);
			}
		}

		Scatter = B * B'; 
	}
	else
	{
		oxwarning("df is less than the dimension of the scale matrix.");
		oxwarning("If df is fractional, it will be replaced by its floor value.");
		oxwarning(1);
		decl intdf = floor(df); 
		decl X = rann(p, intdf);
		Scatter = X * X';
	}
 
	return A * Scatter * A';  
}


// Generating Random Inverse Wishart Matrices
raninvwish(const df, const Sc)
{
	decl p = sizer(Sc);

	if(df >= p)
	{
		return(invert(ranwish(df, invert(Sc))));
	}
	else
	{
		oxrunerror("df is less than the dimension of the scale matrix. This matrix cannot be produced.", 1);
	}
}


// Density function from the Wishart Distribution
denwish(const W, const df, const Sc)
{
	decl pi = 3.141592653589793;
	decl p  = sizer(Sc);
	decl val;

	if(sizer(Sc) != sizec(Sc)) oxrunerror("The scale matrix is not symmetric.", 1);
	if(sizer(W) != sizec(Sc)) oxrunerror("X and the scale matrix are of different dimension.", 1);		
	if(sizer(W) != sizec(W)) oxrunerror("X is not symmetric.", 1);

	if(df >= p)
	{
		decl gamma = prodr( gammafact(.5 * (df + 1 - range(1,p))) );
		decl denom = gamma * (2^(df * p/2)) * pi^(p * (p - 1)/4);
		decl detSc = determinant(Sc);
		decl detW  = determinant(W);
		decl exptr = invert(Sc)*W;
		decl num   = detSc^(-df/2) * detW^((df - p - 1)/2) * exp(-0.5 * trace(exptr));
		val = num/denom;		 
	}
	else
	{
		oxwarning("***************************************************************************");
		oxwarning("* df is less than the dimension of the scale matrix.                      *");
		oxwarning("* If df is fractional, it will be replaced by its floor value.            *");		
		oxwarning("* The density for the singular case is given by Srivastava (2003, p1549). *");
		oxwarning("* Note that the density is now defined on a different spcae.              *");		
		oxwarning("***************************************************************************");		
		decl intdf = floor(df);		
		decl S11   = W[0:(intdf-1)][0:(intdf-1)];
		decl etr   = exp(-0.5 * trace(invert(Sc) * W));
		decl gamma = pi^(intdf * (intdf - 1)/4) * prodr( gammafact(.5 * (intdf + 1 - range(1,intdf))) );
		decl denom = gamma * (determinant(Sc)^(intdf/2));
		decl num   = pi^(intdf * (intdf - p)/2) * 2^(-(intdf * p)/2);
		val = (num/denom) * (determinant(S11)^((intdf - p - 1)/2)) * etr;
	}

	return val;  
}


// Density function from the Inverse Wishart Distribution
deninvwish(const W, const df, const Sc)
{
	decl pi = 3.141592653589793;
	decl p  = sizer(Sc);
	decl val;

	if(sizer(Sc) != sizec(Sc)) oxrunerror("The scale matrix is not symmetric.", 1);
	if(sizer(W) != sizec(Sc)) oxrunerror("X and the scale matrix are of different dimension.", 1);		
	if(sizer(W) != sizec(W)) oxrunerror("X is not symmetric.", 1);

	if(df >= p)
	{
		decl gamma = prodr( gammafact(.5 * (df + 1 - range(1,p))) );
		decl denom = gamma * (2^(df * p/2)) * pi^(p * (p - 1)/4);
		decl detSc = determinant(Sc);
		decl detW  = determinant(W);
		decl exptr = Sc * invert(W);
		decl num   = detSc^(df/2) * detW^(-(df + p + 1)/2) * exp(-0.5 * trace(exptr));
		val = num/denom;		 
	}
	else
	{
		oxrunerror("df is less than the dimension of Sc. The pdf is not defined.", 1);
	}

	return val;  
}
