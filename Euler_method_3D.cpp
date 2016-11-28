//implementacion del metodo de euler para resolver el sistema de ecuaciones
#include <iostream>
#include <cmath>

//declaracion de constantes
const double M = 0.012277471;//relacion entre las masas
const double T = 17.0652165601579625588917206249;//periodo
const double h = T/240000;

//declaracion de funciones
double f1 (double p, double u, double x, double y, double z, double q);
double g1 (double p, double u, double x, double y, double z, double q);
double g2 (double p, double u, double x, double y, double z, double q);
double f2 (double p, double u, double x, double y, double z, double q);
double h1 (double p, double u, double x, double y, double z, double q);
double h2 (double p, double u, double x, double y, double z, double q);
void euler (double h, double &p0, double &u0, double &y0, double &x0, double &z0, double &q0);

int main (void)
{
  double x0, u0, y0, p0, z0,q0;
  x0 = 0.994;
  u0 = 0.0;
  y0 = 0.0;
  p0 = -2.0015851063790825224053786224;
  z0 = 0.0001;
  q0 = -0.000001;
  
  for (int ii = 0;ii<=240000;ii ++){
   
    std::cout<<x0<<"\t"<<y0<<'\t'<<z0<<std::endl;
    euler (h,p0,u0,y0,x0,z0,q0);}
  return 0;
}

double f1 (double p, double u, double x, double y, double z, double q)
{return p;}

double g1 (double p, double u, double x, double y, double z, double q)
{return u;}

double g2 (double p, double u, double x, double y, double z, double q)
{
 long double a, b, c, r, s;
  a=1-M;
  b=x+M;
  c=x-1+M;
  r=std::sqrt(b*b + y*y+ z*z);
  s=std::sqrt(c*c + y*y+ z*z);
  return x + 2.0*p -a*b/pow(r,3)-M*c/pow(s,3);
}


double f2 (double p, double u, double x, double y, double z, double q)
{
long double a, b, c, r, s;
  a=1-M;
  b=x+M;
  c=x-1+M;
  r=std::sqrt(b*b + y*y+ z*z);
  s=std::sqrt(c*c + y*y+ z*z);
  return y-2.0*u-(a*y)/pow(r,3)-M*y/pow(s,3);
}

double h1 (double p, double u, double x, double y, double z, double q)
{
  return q;
}
double h2 (double p, double u, double x, double y, double z, double q)
{
long double a, b, c, r, s;
  a=1-M;
  b=x+M;
  c=x-1+M;
  r=std::sqrt(b*b + y*y+ z*z);
  s=std::sqrt(c*c + y*y+ z*z);
  return -(a*z/pow(r,3))-(M*z/pow(s,3));
}

void euler (double h, double &p0, double &u0, double &y0, double &x0, double &z0, double &q0)
{
  x0 += h*g1(p0,u0,x0,y0,z0,q0);
  y0 += h*f1(p0,u0,x0,y0,z0,q0);
  z0 += h*h1(p0,u0,x0,y0,z0,q0);
  p0 += h*f2(p0,u0,x0,y0,z0,q0);
  u0 += h*g2(p0,u0,x0,y0,z0,q0);
  q0 += h*h2(p0,u0,x0,y0,z0,q0);
  
}
