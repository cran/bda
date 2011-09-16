/* Compile under R 
 R CMD SHLIB kernel.c -o kernel.so
*/

#include <R.h>
#include <Rmath.h>
#define INLINE

void kspvalue(double *pv){
  double ksp=0.,t;
  int i;
  t = pv[0];
  for(i=1;i<100;i++){
    ksp += exp(-2.*pow(i*t,2.));
    i++;
    ksp -= exp(-2.*pow(i*t,2.));
  }
  pv[0] = 2.*ksp;
}

double g0(double a, double c){
  //first derivative of logarithm Gamma(a)
  return(-digamma(a)+log(a)+c);
}

double g1(double a){
  return(-trigamma(a)+1./a);
}

void FitGamma(double *x, int *n, double *para){
  double ex1=0.,ex2=0., elx1=0., lex1=0.;
  double c=0.;
  int i;
  for(i=0;i<n[0];i++){
    ex1 += x[i];
    elx1 += log(x[i]);
  }
  ex1 = ex1/n[0];
  elx1 = elx1/n[0];
  lex1 = log(ex1);
  for(i=0;i<n[0];i++){
    ex2 += pow(x[i]-ex1,2.0);
  }
  ex2 = ex2/(n[0]-1.);
  c = elx1 - lex1;
  
  double x0=0.,x1, eps=0.00000001;
  int Iter=20;
  i=0;
  x1 = pow(ex1,2.0)/ex2; //mme of alpha
  //    x1 = ex2/ex1; //mme of beta
  while(i<Iter && fabs(x1-x0)>eps){
    x0 = x1;
    x1 = x1 - g0(x1,c)/g1(x1);
    i=i+1;
  }
  x0 = x1/ex1;
    
  para[0] = x1;
  para[1] = x0;
}

void FitBeta(double *x, int *n, double *para){
  double ex1=0.,ex2=0., elx1=0., elx2=0.;
  int i;
  for(i=0;i<n[0];i++){
    ex1 += x[i];
    elx1 += log(x[i]);
    elx2 += log(1.-x[i]);
  }
  ex1 = ex1/n[0];
  elx1 = elx1/n[0];
  elx2 = elx2/n[0];
  for(i=0;i<n[0];i++){
    ex2 += pow(x[i]-ex1,2.0);
  }
  ex2 = ex2/(n[0]-1.);
  
  double a0=0.,a1=0., b0=0., b1=0., eps=0.00000001;
  int Iter=20;
  double F1, F2, a,b,c;
  i=0;
  a1 = ex1*(ex1*(1.-ex1)/ex2-1.);
  b1 = (1-ex1)*(ex1*(1.-ex1)/ex2-1.);
  while(i<Iter && fabs(a1-a0)>eps && fabs(b1-b0)>eps){
    a0 = a1; b0=b1;
    b = -trigamma(a1+b1);
    a = trigamma(a1) +b;
    c = trigamma(b1)+b;
    F1 = digamma(a1)-digamma(a1+b1)-elx1;
    F2 = digamma(b1)-digamma(a1+b1)-elx2;
    a1 = a1 - (c*F1-b*F2)/(a*c-pow(b,2.));
    b1 = b1 - (-b*F1+a*F2)/(a*c-pow(b,2.));
    i=i+1;
  }
  para[0] = a1;
  para[1] = b1;
}


void FitWeibull(double *x, int *n, double *para){
  double mulx=0.;
  int i,j;
  double x0=1.,x1=1.5, eps=0.00000001;
  double g0x, g1x, a1,a2,a3,a4,c1,c2;
  int Iter=20;
  i=0;
  while(i<Iter && fabs(x1-x0)>eps){
    x0 = x1;
    a1=0.;a2=0.;a3=0.;a4=0.;
    c1=0.;c2=0.;
    for(j=0;j<n[0];j++){
      c1 = pow(x[j],x1);
      c2 = log(x[j]);
      a1 += c2*c1;
      a2 += c1;
      a3 += c2;
      a4 += pow(c2,2.)*c1;
    }
    g0x= 1./x1-a1/a2+a3/n[0];
    g1x = -1./pow(x1,2.) - (a4*a2 - pow(a1,2.))/pow(a2,2.);
    x1 = x1 - g0x/g1x;
    i=i+1;
  }

  x0 = pow(a2/n[0],1./x1);
    
  para[0] = x1;
  para[1] = x0;
}

