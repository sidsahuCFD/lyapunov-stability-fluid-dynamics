double** Eigen2std(Eigen::MatrixXd U_D, const int N, const int N_p, int q_n, int r_m)
{
   int i,n(N_p),m(N);
   double ** a = new double*[q_n];
   double ** r = new double*[n];
   for(i = 0; i < n; i++) 
   {
        a[i] = new double[m];
   }
   for(; i < q_n; i++) 
   {
        a[i] = new double[m];
   }    
    
   
   for(int i=0;i<N;i++)
   {
   	for(int j=0;j<N_p;j++)
   	{
		a[j][i] = U_D(i,j);
	}		
   }	
return a;
}
