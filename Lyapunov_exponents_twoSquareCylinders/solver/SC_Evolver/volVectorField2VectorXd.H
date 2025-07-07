 Eigen::VectorXd volVectorField2VectorXd(volVectorField U_mat,const int N)
 {
   Eigen::VectorXd U_vec(N); 
   int k(0);
   
   // for loop that converts a volVectorField of 3 components into a vector
   
   forAll (U_mat, id)
{
    U_vec(k) = U_mat[id].component(0);
    U_vec(k+1) = U_mat[id].component(1);
    k=k+2;    
}
return(U_vec);
}
