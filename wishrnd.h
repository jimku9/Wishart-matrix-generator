// ************************************************************************************ //
//  Wishrnd.h  - Generators for Wishart and Inverse Wishart                             //
//  Yu-Cheng Ku & Peter Bloomfield, Dept of Stat, North Carolina State University       //
//                                                                                      //
//  ranwish    - Random Wishart matrix generator                                        //
//  raninvwish - Random inverse Wishart matrix generator                                //
//                                                                                      //
//  Arguments:                                                                          //
//  df is the degrees of freedom.                                                       //
//  Sc is the scale matrix of the Wishart distribution.                                 //
//  W  is the positive definite matrix to be evaluatated.                               //
//                                                                                      //
//  References:                                                                         //
//  1. Anderson, T. W., (2003), An Introduction to Statitistcal Multivariate Analysis.  //
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