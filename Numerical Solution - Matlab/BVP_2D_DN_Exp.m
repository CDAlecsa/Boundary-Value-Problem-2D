%%   Non-Homogenuous Stationary BVP with Variable Coefficients and Dirichlet-Neumann BC (Case 2D)
%%   ============================================================================================

function [errorL2,dx,dy] = BVP_2D_DN_Exp(Nmod,varK,phiExpNmod10000,wavenumberExp0Nmod10000,wavenumberExp1Nmod10000)
   

%%   Variable Initialization

        KMean = 15;
             
        phi = phiExpNmod10000(1:Nmod) ;
        C1 = wavenumberExp0Nmod10000(1:Nmod) ;
        C2 = wavenumberExp1Nmod10000(1:Nmod) ;
        
    
%%   Grid Initialization
 

        Nx = 10^3 + 1 ;
        Ny = 501 ;
        
        a = 0 ;
        b = 20 ;
        
        c = 0 ;
        d = 10 ;
        
        dx = (b-a)/(Nx-1) ;
        dy = (d-c)/(Ny-1) ;
        
        x = a + ( (1:Nx)-1 )*dx ;
        y = c + ( (1:Ny)-1 )*dy ;
        
        [X,Y] = meshgrid(x,y) ;
        
        
%%   Exact Solution 
  
        solE = @(u,v) 1 + sin(2*u + v) ;
        uE = solE(X,Y);
        
        
 %%   Matrix Initializations  
        
       
        f = func(x,y,Nmod,KMean,varK,C1,C2,phi,Nx,Ny) ;
        CM1 = K(x + dx/2,y,Nmod,KMean,varK,C1,C2,phi,Nx,Ny) ;
        CM2 = K(x,[y(1) - dy, y] + dy/2,Nmod,KMean,varK,C1,C2,phi,Nx,Ny+1) ;
        
        
        A = CM1(1:end-2,:)*(1/dx^2) ;
        D = CM1(2:end-1,:)*(1/dx^2) ;
        
        B = CM2(2:end-1,1:end-1)*(1/dy^2) ;
        E = CM2(2:end-1,2:end)*(1/dy^2) ;
        
        C = (-1) * (A+B+D+E) ; 
        
        
%%   The LHS of the Numerical Solution      

       Diag0 = C(:) ;
       
       auxLW = A(2:end,:) ;
       auxLW = [ auxLW ; zeros(1,Ny) ] ;
       DiagLW = auxLW(:) ;
       DiagLW = DiagLW(1:end-1) ;
       
       auxUP = D(1:end-1,:) ; 
       auxUP = [ auxUP ; zeros(1,Ny) ] ;
       DiagUP = auxUP(:) ;
       DiagUP = DiagUP(1:end-1) ;
       
       DiagUP2 = E(:,1:end-1) ;
       DiagUP2 = DiagUP2(:) + [ B(:,1) ; zeros(((Nx-2)*Ny-(Nx-2))-length(B(:,1)),1)] ;
       
       DiagLW2 = B(:,2:end);
       DiagLW2 = DiagLW2(:) + [zeros(((Nx-2)*Ny-(Nx-2))-length(E(:,end)),1) ; E(:,end)] ;
 
       LHS = sparse((Nx-2)*Ny,(Nx-2)*Ny) ;
       LHS = spdiags(Diag0,0,LHS) ;                   % main diagonal
       LHS = spdiags(DiagLW,-1,LHS) ;                 % lower diagonal
       LHS = spdiags(DiagLW2,-(Nx-2),LHS) ;           % -(Nx-2) diagonal    
       
       DiagUP = [0 DiagUP']';                         % upper diagonal
       LHS = spdiags(DiagUP,1,LHS) ;
       
       DiagUP2 = [zeros(Nx-2,1) ; DiagUP2] ;          % (Nx-2) upper diagonal
       LHS = spdiags(DiagUP2,Nx-2,LHS) ;
       

%%   Boundary Conditions
        
        BCXL = @(v) 1 + sin(2*a + v) ;
        BCXR = @(v) 1 + sin(2*b + v) ;
        
        BCYB = @(u) cos(2*u + c) ; 
        BCYU = @(u) cos(2*u + d) ;

        uN = zeros(Nx,Ny) ;
        uN(1,:) = BCXL(y) ;
        uN(Nx,:) = BCXR(y) ;

        

%%   The RHS of the Matrix Solution
            
         vRHS_F = f(2:end-1,:) ;
         vRHS_F = vRHS_F(:) ;
         
         
         auxRHS_A = [ A(1,:).*BCXL(y)  ; zeros(Nx-3,length(A(1,:))) ] ;
         auxRHS_A = auxRHS_A(:) ; 
         vRHS_A = (-1) *auxRHS_A ;
         
         auxRHS_D = [ zeros(Nx-3,length(D(end,:))) ; D(end,:).*BCXR(y) ] ;
         vRHS_D = (-1) * auxRHS_D(:) ;
         
         
        vRHS_Start = [2*dy * (BCYB(x(2:Nx-1)) .* B(:,1)')  zeros(1,(Nx-2)*Ny-(Nx-2))]';
        vRHS_End = [zeros(1,(Nx-2)*Ny-(Nx-2))  (-2)*dy * (BCYU(x(2:Nx-1)) .* E(:,Ny)')]';
        RHS = vRHS_F + vRHS_A + vRHS_D + vRHS_Start + vRHS_End ;


%%   Numerical Solution

        rez = (LHS\RHS)' ;
        rez = reshape(rez,Nx-2,Ny) ; 
        uN(2:Nx-1,:) = rez ; 
         
        
 %%   Optional Plots between the Exact Solution and the Numerical Solution
 
%     figure ;
%         mesh(X,Y,uE) ;
%         xlabel('x') ;
%         ylabel('y') ;
%         title('Exact Solution') ;
%     figure ;
%         mesh(X,Y,uN') ;
%         xlabel('x') ;
%         ylabel('y') ;
%         title('Numerical Solution') ;
%     figure ;
%         mesh(X,Y,abs(uE-uN')) ;
%         xlabel('x') ;
%         ylabel('y') ;
%         title('Absolute Error') ;
        
      errorL2 = ( dx * dy )^(1/2) * norm(abs(uE-uN')) ;
      
end