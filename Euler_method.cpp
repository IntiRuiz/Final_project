//implementacion del metodo de euler para resolver el sistema de ecuaciones
#include <iostream>
#include <cmath>

//declaracion de constantes
const double M = 0.012277471;//relacion entre las masas
const double T = 17.0652165601579625588917206249;//periodo


//declaracion de funciones
double f1 (double p, double u, double x, double y);
double g1 (double p, double u, double x, double y);
double g2 (double p, double u, double x, double y);
double f2 (double p, double u, double x, double y);
void euler (double h, double &p0, double &u0, double &y0, double &x0);

int main (void)
{
  double x0, u0, y0, p0, h;
  x0 = 0.994;
  u0 = 0.0;
  y0 = 0.0;
  p0 = -2.0015851063790825224053786224;
  h = T/24000;
  
  for (int ii = 1;ii<=24000;ii++){
    euler (h,p0,u0,y0,x0);
    std::cout<<x0<<"\t"<<u0<<std::endl;}
  
  return 0;
}

double f1 (double p, double u, double x, double y)
{return p;}

double g1 (double p, double u, double x, double y)
{return u;}

double g2 (double p, double u, double x, double y)
{
 long double a, b, c, r, s;
  a=1-M;
  b=x+M;
  c=x-1+M;
  r=std::sqrt(b*b + y*y);
  s=std::sqrt(c*c + y*y);
  return x + 2.0*p -(a*b/pow(r,3))-(M*c/pow(s,3));
}


double f2 (double p, double u, double x, double y)
{
long double a, b, c, r, s;
  a=1-M;
  b=x+M;
  c=x-1+M;
  r=std::sqrt(b*b + y*y);
  s=std::sqrt(c*c + y*y);
  return y-2.0*u-((a*y)/pow(r,3))-(M*y/pow(s,3));
}

void euler (double h, double &p0, double &u0, double &y0, double &x0)
{
  x0 = x0 + (h*g1(p0,u0,x0,y0));
  y0 = y0 + (h*f1(p0,u0,x0,y0));
  p0 = p0 + (h*f2(p0,u0,x0,y0));
  u0 = u0 + (h*g2(p0,u0,x0,y0));
 
}
