Eigen::MatrixXd std2Eigen(double** a, const int N, const int N_p)
{
   Eigen::MatrixXd Q(Eigen::MatrixXd::Random(N,N_p));
   
   for(int i=0;i<N;i++)
   {
   	for(int j=0;j<N_p;j++)
   	{
		Q(i,j) = a[j][i];
	}
		
   }	
return Q;
}
