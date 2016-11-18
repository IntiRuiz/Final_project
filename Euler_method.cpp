//implementacion del metodo de euler para resolver el sistema de ecuaciones
#include <iostream>
#include <cmath>

//declaracion de constantes
const double miu = 0.012277471;//relacion entre las masas
const double T = 17.0652165601579625588917206249;//periodo

/*el sistema de ecuaciones a resolver es:
f1(p,u,x,y) = dy/dt = p
g1(p,u,x,y) = dx/dt = u
g2(p,u,x,y) = du/dt = x-2p-(((1-miu)*(x+miu))/(pow (sqrt(pow((x+miu),2)+pow(y,2)),3)))-((miu*(x-1+miu))/(pow (sqrt(pow((x-1+miu),2)+pow(y,2)),3)))
f2(p,u,x,y) = dp/dt = y-2u-(((1-miu)*y)/(pow (sqrt(pow((x+miu),2)+pow(y,2)),3)))-((miu*y)/(pow (sqrt(pow((x-1+miu),2)+pow(y,2)),3)))
*/

//declaracion de funciones
double f1 (double p);
double g1 (double u);
double g2 (double p, double u, double x, double y);
double f2 (double p, double u, double x, double y);
void euler (double h, double &p0, double &u0, double &y0, double &x0);

int main (void)
{
  double p0,u0,x0,y0,h;
  h = T/24000;
  p0 = -2.0015851063790825224053786224;
  u0 = 0;
  x0 = 0.994;
  y0 = 0;

    for (int ii=1; ii<= 24000; ii++)
      {euler (h, p0, u0, y0, x0);
	std::cout<<(h*ii)<<'\t'<<p0<<'\t'<<(h*ii)<<'\t'<<(h*ii)<<'\t'<<'\t'<<x0<<'\n';}
  
  return 0;
}

double f1 (double p)
{return p;}

double g1 (double u)
{return u;}

double g2 (double p, double u, double x, double y)
{return  (x-(2*p))-(((1-miu)*(x+miu))/(std::pow (std::sqrt(std::pow((x+miu),2)+std::pow(y,2)),3)))-((miu*(x-1+miu))/(std::pow (std::sqrt(std::pow((x-1+miu),2)+std::pow(y,2)),3)));}

double f2 (double p, double u, double x, double y)
{return (y-(2*u))-(((1-miu)*y)/(std::pow (std::sqrt(std::pow((x+miu),2)+std::pow(y,2)),3)))-((miu*y)/(std::pow (std::sqrt(std::pow((x-1+miu),2)+std::pow(y,2)),3)));}

void euler (double h, double &p0, double &u0, double &y0, double &x0)
{
  p0 = p0 + (h*f1(p0));
  u0 = u0 + (h*g1(u0));
  x0 = x0 + (h*g2(p0,u0,x0,y0));
  y0 = y0 + (h*f2(p0,u0,x0,y0));
}
