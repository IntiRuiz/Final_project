//implementacion metodo de runge kutta 
#include <iostream>
#include <iomanip>
#include <cmath>
//declaracion de constantes
const long double M = 0.012277471;//relacion entre las masas
const long double T = 17.0652165601579625588917206249;//periodo
const long double delta= T/6000;

long double f0(long double t, long double x0, long double x1, long double y0, long double y1);
long double f1(long double t, long double x0, long double x1, long double y0, long double y1);
long double g0(long double t, long double y0, long double y1, long double x0, long double x1);
long double g1(long double t, long double y0, long double y1, long double x0, long double x1);
void rk4(long double t, long double h, long double  &x0, long double &x1, long double &y0, long double &y1);

int main(void)
{
  long double x, vx, y, vy, time;
  x=0.994;
  vx=0.0;
  y=0.0;
  vy=-2.0015851063790825224053786224;
  for(time=0.0; time<=T; time+=delta)
    {
      std::cout<<x<<"\t"<<vx<<std::endl;
      rk4(time, delta, x, vx, y, vy);
    }
  return 0;
}

long double f0(long double t, long double x0, long double x1, long double y0, long double y1)
{
  return x1;
}

long double f1(long double t, long double x0, long double x1, long double y0, long double y1)
{
  long double a, b, c, r, s;
  a=1-M;
  b=x0+M;
  c=x0-1+M;
  r=std::sqrt(b*b + y0*y0);
  s=std::sqrt(c*c + y0*y0);
  return x0 + 2.0*y1 -(a*b/pow(r,3))-(M*c/pow(s,3));
}

long double g0(long double t, long double y0, long double y1, long double x0, long double x1)
{
  return y1;
}

long double g1(long double t, long double y0, long double y1, long double x0, long double x1)
{
  long double a, b, c, r, s;
  a=1-M;
  b=x0+M;
  c=x0-1+M;
  r=std::sqrt(b*b + y0*y0);
  s=std::sqrt(c*c + y0*y0);
  return y0-2.0*x1-((a*y0)/pow(r,3))-(M*y0/pow(s,3));
}												  

void rk4(long double t, long double h, long double &x0, long double &x1, long double &y0, long double &y1)
{
  long double k10, k11, k20, k21, k30, k31, k40, k41;
  long double p10, p11, p20, p21, p30, p31, p40, p41;

  k10=h*f0(t, x0, x1, y0, y1);
  k11=h*f1(t, x0, x1, y0, y1);
  k20=h*f0(t+ h/2, x0+ k10/2, x1+ k11/2, y0, y1);
  k21=h*f1(t+ h/2, x0+ k10/2, x1+ k11/2, y0, y1);
  k30=h*f0(t+ h/2, x0+ k20/2, x1+ k21/2, y0, y1);
  k31=h*f1(t+ h/2, x0+ k20/2, x1+ k21/2, y0, y1);
  k40=h*f0(t+h, x0+k30, x1+k31, y0, y1);
  k41=h*f1(t+h, x0+k30, x1+k31, y0, y1);

  x0= x0+ (1.0/6.0)*(k10+(2*k20)+(2*k30)+k40);
  x1= x1+ (1.0/6.0)*(k11+(2*k21)+(2*k31)+k41);

  p10=h*g0(t, y0, y1, x0, x1);
  p11=h*g1(t, y0, y1, x0, x1);
  p20=h*g0(t+ h/2, y0+ p10/2, y1+ p11/2, x0, x1);
  p21=h*g1(t+ h/2, y0+ p10/2, y1+ p11/2, x0, x1);
  p30=h*g0(t+ h/2, y0+ p20/2, y1+ p21/2, x0, x1);
  p31=h*g1(t+ h/2, y0+ p20/2, y1+ p21/2, x0, x1);
  p40=h*g0(t+h, y0+p30, y1+p31, x0, x1);
  p41=h*g1(t+h, y0+p30, y1+p31, x0, x1);

  y0= y0+ (1.0/6.0)*(p10+(2*p20)+(2*p30)+p40);
  y1= y1+ (1.0/6.0)*(p11+(2*p21)+(2*p31)+p41);
}


													    
												    
  
  
