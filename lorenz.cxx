#include <iostream>
#include <cmath>
#include <fstream>
using namespace std; 

// f(p(x,y.z))
void kutta(double *f, double *p, const double a, const double b, 
	   const double c)
{
  f[0]= a*(p[1]-p[0]);
  f[1]=p[0]*(b-p[2])-p[1];
  f[2]=p[0]*p[1]-c*p[2];
}

int main(){

const double dt=0.01;
const double T = 100;
const int N = T/dt;
const double a=10.0;
const double b=28.0;
const double c=8.0/3.0;
double k1[3];
double k2[3];
double k3[3];
double k4[3];
double p[3];
p[0]=1;
p[1]=1;
p[2]=1;

double ytemp [3];

for(int i=0;i<N;i++){

  kutta(k1,p,a,b,c);  // k1 = f(y_n)
    
    ytemp[0]=p[0]+dt*0.5*k1[0];
    ytemp[1]=p[1]+dt*0.5*k1[1];
    ytemp[2]=p[2]+dt*0.5*k1[2];
  
  kutta(k2,p,a,b,c); // k2 = f(y_n+dt*0.5*k1)
    
    ytemp[0]=p[0]+dt*0.5*(k1[0]+k2[0]);
    ytemp[1]=p[1]+dt*0.5*(k1[1]+k2[1]);
    ytemp[2]=p[2]+dt*0.5*(k1[2]+k2[2]);
    
  kutta(k3,p,a,b,c); // k3 = f(y_n+dt*0.5*(k1+k2))
    
    ytemp[0]=p[0]+dt*0.5*(k1[0]+k2[0]+2*k3[0]);
    ytemp[1]=p[1]+dt*0.5*(k1[1]+k2[1]+2*k3[1]);
    ytemp[2]=p[2]+dt*0.5*(k1[2]+k2[2]+2*k3[2]);  
 
  
  kutta(k4,p,a,b,c); //k4 =f(y_n+dt*0.5(k1+k2+2k3))

  p[0]+= dt/6.0*(k1[0]+k2[0]+2*k3[0]);
  p[1]+= dt/6.0*(k1[1]+k2[1]+2*k3[1]);
  p[2]+= dt/6.0*(k1[2]+k2[2]+2*k3[2]);
   cout << i*dt <<'\t'<< p[0] << '\t' << p[1] << '\t' << p[2] <<endl;
  
}

 return 0;
}

