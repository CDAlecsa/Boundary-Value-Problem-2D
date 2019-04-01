%% Main

    clear ;
    close all ; 
    clc ;

    
%% Loading Coefficients

    load phiExpNmod10000 ;
    load wavenumberExp0Nmod10000 ;
    load wavenumberExp1Nmod10000 ;
     
    
%% Generate Numerical Solutions for different values of the parameters

    Nmod =[10^2, 10^3, 10^4];
    varK = [0.1, 1, 2, 4, 6, 8, 10];

    fileID = fopen('resultsEXP2D_solMan.txt','w');

    for i=1:3
        for j = 1:7 
            tic
            [errorL2,dx,dy] = BVP_2D_DN_Exp(Nmod(i),varK(j),phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000);
            T = toc;
            fprintf(fileID,'Nmod = %d \t varK = %f \t Err = %2.2e dx = %2.2e dy = %2.2e  time = %f\n', Nmod(i),varK(j),errorL2,dx,dy,T);
            fprintf('Nmod = %d \t varK = %f \t Err = %2.2e dx = %2.2e dy = %2.2e  time = %f\n', Nmod(i),varK(j),errorL2,dx,dy,T);
        end
    end
    
    fclose(fileID);