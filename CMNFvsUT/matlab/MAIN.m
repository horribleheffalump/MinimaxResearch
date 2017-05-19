clear all;

R_Seed = 949477711;
R_Stream = RandStream.create('mt19937ar', 'Seed', R_Seed);
RandStream.setGlobalStream (R_Stream);
W = zeros ( 2 , 1 );

M_X = [ 30 ; 40 ];
S_X = [ 30 0 ; 0 30 ];
M_E = [ 0 ; 0 ];
S_E = [ 5 * 3.1415926 / 180 0 ; 0 30 ];
 
N_mom = 10000000;
N_sam  = 10000000;

X_mom = mvnrnd (M_X , S_X * S_X' , N_mom );
E_mom = mvnrnd (M_E , S_E * S_E' , N_mom );
X_mom = X_mom';
E_mom = E_mom';
[ Theta , Rho ] = cart2pol ( X_mom ( 1 , : ) , X_mom ( 2 , : ) );
Y_mom = [ Theta ; Rho ] + E_mom;
[ x1_inv , x2_inv ] = pol2cart ( Y_mom ( 1 , : ) , Y_mom ( 2 , : ) ); 
Z_mom = [ X_mom ; Y_mom ; x1_inv ; x2_inv ];

K_ZZ = cov (Z_mom');
M_Y = mean ( Z_mom , 2 );
P = K_ZZ ( 1 : 2 , 3 : 6 ) * pinv ( K_ZZ ( 3 : 6 , 3 : 6 ) )

X_mom = mvnrnd (M_X , S_X * S_X' , N_sam );
E_mom = mvnrnd (M_E , S_E * S_E' , N_sam );
X_mom = X_mom';
E_mom = E_mom';
[ Theta , Rho ] = cart2pol ( X_mom ( 1 , : ) , X_mom ( 2 , : ) );
Y_mom = [ Theta ; Rho ] + E_mom ;
[ x1_inv , x2_inv ] = pol2cart ( Y_mom ( 1 , : ) , Y_mom ( 2 , : ) ); 
Y_mom = [ Y_mom ; x1_inv ; x2_inv ];

X_hat = M_X * ones ( 1 , N_sam ) + P * ( Y_mom - M_Y ( 3 : 6 , 1 ) * ones ( 1 , N_sam ) );
Err = X_hat - X_mom;
M_Err_CM = mean ( Err , 2 )
K_Err_CM = cov ( Err' )
K_Err_th_CM = S_X^2 - P * K_ZZ ( 3 : 6 , 3 : 6 ) * P'

P = K_ZZ ( 1 : 2 , 3 : 4 ) * pinv ( K_ZZ ( 3 : 4 , 3 : 4 ) );
X_hat = M_X * ones ( 1 , N_sam ) + P * ( Y_mom (  1 : 2 , : ) - M_Y ( 3 : 4 , 1 ) * ones ( 1 , N_sam ) );
Err = X_hat - X_mom;
M_Err_Lin = mean ( Err , 2 )
K_Err_Lin = cov ( Err' )
K_Err_th_Lin = S_X^2 - P * K_ZZ ( 3 : 4 , 3 : 4 ) * P'

Alpha_0 = 0.33;
%Alpha_0 = 0.45;
Alpha = ( 1 - Alpha_0 ) / 4 * ones ( 1 , 5 );
Alpha ( 1 , 1 ) = Alpha_0;
X_UT = zeros ( 2 , 5 );
X_UT ( : , 1 ) = M_X;
SQ_X = sqrtm ( S_X * S_X' );
X_UT ( : , 2 ) = M_X - sqrt ( 2 / ( 1 - Alpha_0 ) ) * SQ_X ( : , 1 );
X_UT ( : , 3 ) = M_X - sqrt ( 2 / ( 1 - Alpha_0 ) ) * SQ_X ( : , 2 );
X_UT ( : , 4 ) = M_X + sqrt ( 2 / ( 1 - Alpha_0 ) ) * SQ_X ( : , 1 );
X_UT ( : , 5 ) = M_X + sqrt ( 2 / ( 1 - Alpha_0 ) ) * SQ_X ( : , 2 );
[ Theta_UT , Rho_UT ] = cart2pol ( X_UT ( 1 , : ) , X_UT ( 2 , : ) );
Y_UT = [ Theta_UT ; Rho_UT ];
%M_UT = mean ( Y_UT * diag ( Alpha ) , 2 )
M_UT = zeros ( 2 , 1 );
for i = 1 : 5
    M_UT = M_UT + Alpha ( 1 , i ) * Y_UT ( : , i );
end;
M_UT;
Z_UT = [ X_UT ; Y_UT ]; 
K_UT = zeros ( 4 , 4 );
MZ_UT = [ M_X ; M_UT ];
        for k = 1 : 5
            K_UT =  K_UT +  Alpha ( 1 , k ) * ( Z_UT ( : , k ) - MZ_UT ) * ( Z_UT ( : , k ) - MZ_UT )';
        end;    
P_UT = K_UT ( 1 : 2 , 3 : 4 ) * pinv ( K_UT ( 3 : 4 , 3 : 4 ) + S_E * S_E' );

X_hat_UT = M_X * ones ( 1 , N_sam ) + P_UT * ( Y_mom ( 1:2 , : ) - M_UT * ones ( 1 , N_sam ) );
Err_UT = X_hat_UT - X_mom;
M_Err_UT = mean ( Err_UT , 2 )
K_Err_UT = cov ( Err_UT' )
K_Err_th_UT = S_X^2 - P_UT * ( K_UT ( 3 : 4 , 3 : 4 ) + S_E * S_E' ) * P_UT'

figure(1)
plot ( x1_inv , x2_inv )
figure(2)
plot ( X_mom ( 1 , : ) , X_mom ( 2 , : ) );
