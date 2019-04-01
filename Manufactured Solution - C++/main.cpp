//  compute K values on a grid; lnK is log-normally distributed with exponential correlation
//  generated with function "calc_K_2018-exp-lambda.c"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PI 3.141592653589

int i,j,imod,XF,YF;  
double x,y,dx,dy,coeff,phase,vr;
double Nmod,K_MEAN,varK,K,S1,S2,S3,F;
double C[10000][2],phi[10000];
FILE *fK,*fF, *file1, *file2, *file3;

//double calc_K(double x, double y, double Nmod, double K_MEAN, double varK, double C[10000][2], double phi, double K)
void calc_K()
{
  coeff = sqrt(varK*2/Nmod) ;
  vr=0;
  for (imod=0; imod<Nmod; imod++) {
    phase=2*PI*(C[imod][0]*x+C[imod][1]*y)+phi[imod];
    vr+=coeff*cos(phase);
  } 
  K= K_MEAN*exp(-varK/2.0)* exp(vr);
}

void calc_F()
{
  S1=0; S2=0; S3=0;
  for (imod=0; imod<Nmod; imod++) {
    S1+=(-2*PI)*C[imod][1]*sin(phi[imod] + (2*PI)*(C[imod][1]*x + C[imod][2]*y)) ;                    
    S2+=cos(phi[imod] + (2*PI)*(C[imod][1]*x + C[imod][2]*y)) ;                    
    S3+=(-2*PI)*C[imod][2]*sin(phi[imod] + (2*PI)*(C[imod][1]*x + C[imod][2]*y)) ;
  }       
  F = 2*K_MEAN*exp(-varK/2)*sqrt(varK*2/Nmod)*S1*exp(sqrt(varK*2/Nmod)*S2)*cos(2*x+y)
    - 5*K_MEAN*exp(-varK/2)*exp(sqrt(varK*2/Nmod)*S2)*sin(2*x+y)
    + K_MEAN*exp(-varK/2)*sqrt(varK*2/Nmod)*S3*exp(sqrt(varK*2/Nmod)*S2)*cos(2*x+y);
}  
  
void write_KF()
{
  printf("dx= %10.5f dy= %10.5f XF= %10d YF= %10d\n",dx,dy,XF,YF);  
  for (i=0; i<XF; i++) {
    for (j=0; j<YF; j++) {
      x=i*dx; y=j*dy;
      calc_K();
      fprintf(fK,"%le %le %le \n",x,y,K);
      calc_F();
      fprintf(fF,"%le %le %le \n",x,y,F);

    }
    printf("x-section i= %5d\n",i);
  }
}

main()
{
  fK=fopen("KF/k","w");
  fF=fopen("KF/f","w");
  file1 = fopen("wavenumberExp0Nmod10000", "r");
  file2 = fopen("wavenumberExp1Nmod10000", "r");
  file3 = fopen("phiExpNmod10000", "r");

  K_MEAN = 14.268;//100.0;
  Nmod = 64;//00;
  varK = 0.1;

  XF=100;
  YF=100;
  dx=0.1;
  dy=0.1;

  for (int j = 0; j<Nmod; j++)
    {
      fscanf(file1, "%le \n", &C[j][0]);
      fscanf(file2, "%le \n", &C[j][1]);
      fscanf(file3, "%le \n", &phi[j]);
    }

  write_KF();

  fclose(file1);
  fclose(file2);
  fclose(file3);
  fclose(fK);
  fclose(fF);   
}
