clear all;

R_Seed = 949477711;
R_Stream = RandStream.create('mt19937ar', 'Seed', R_Seed);
RandStream.setGlobalStream (R_Stream);
W = zeros ( 1 , 1 );

M_X =  10 ;
S_X =  1 ;
M_E =  0;
S_E = 1;
 
N_mom = 10000000;
N_sam  = 10000000;

X_mom = mvnrnd (M_X , S_X * S_X' , N_mom );
E_mom = mvnrnd (M_E , S_E * S_E' , N_mom );
X_mom = X_mom';
E_mom = E_mom';
Theta  = X_mom ( 1 , : ).^3.0;
Y_mom = Theta + E_mom;
x_inv = nthroot( Y_mom ( 1 , : ) , 3 ); 
Z_mom = [ X_mom ; Y_mom ; x_inv ];

K_ZZ = cov (Z_mom');
M_Y = mean ( Z_mom , 2 );
P = K_ZZ ( 1 : 1 , 2 : 3 ) * pinv ( K_ZZ ( 2 : 3 , 2 : 3 ) )

X_mom = mvnrnd (M_X , S_X * S_X' , N_sam );
E_mom = mvnrnd (M_E , S_E * S_E' , N_sam );
X_mom = X_mom';
E_mom = E_mom';
Theta = X_mom ( 1 , : ).^3.0;
Y_mom = Theta + E_mom ;
x_inv = nthroot( Y_mom ( 1 , : ) , 3 );  
Y_mom = [ Y_mom ; x_inv];

X_hat = M_X * ones ( 1 , N_sam ) + P * ( Y_mom - M_Y ( 2 : 3 , 1 ) * ones ( 1 , N_sam ) );
Err = X_hat - X_mom;
M_Err_CM = mean ( Err , 2 )
K_Err_CM = cov ( Err' )
K_Err_th_CM = S_X^2 - P * K_ZZ ( 2 : 3 , 2 : 3 ) * P'

P = K_ZZ ( 1 : 1 , 2 : 2 ) * pinv ( K_ZZ ( 2 : 2 , 2 : 2 ) );
X_hat = M_X * ones ( 1 , N_sam ) + P * ( Y_mom (  1 : 1 , : ) - M_Y ( 2 : 2 , 1 ) * ones ( 1 , N_sam ) );
Err = X_hat - X_mom;
M_Err_Lin = mean ( Err , 2 )
K_Err_Lin = cov ( Err' )
K_Err_th_Lin = S_X^2 - P * K_ZZ ( 2 : 2 , 2 : 2 ) * P'
