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
void dormand(long double t, long double h, long double  &x0, long double &x1, long double &y0, long double &y1);

int main(void)
{
  long double x, vx, y, vy, time;
  x=0.994;
  vx=0.0;
  y=0.0;
  vy=-2.0015851063790825224053786224;
  for(time=0.0; time<=1; time+=delta)
    {
      std::cout<<time<<"\t"<<x<<"\t"<<vx<<std::endl;
      dormand(time, delta, x, vx, y, vy);
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

void dormand(long double t, long double h, long double &x0, long double &x1, long double &y0, long double &y1)
{
  long double k10, k11, k20, k21, k30, k31, k40, k41, k50, k51, k60, k61, k70, k71;
  long double p10, p11, p20, p21, p30, p31, p40, p41, p50, p51, p60, p61, p70, p71;

  k10=h*f0(t, x0, x1, y0, y1);
  k11=h*f1(t, x0, x1, y0, y1);
  k20=h*f0(t+ h/5, x0+ k10/5, x1+ k11/5, y0, y1);
  k21=h*f1(t+ h/5, x0+ k10/5, x1+ k11/5, y0, y1);
  k30=h*f0(t+ (3/10)*h, x0+ (3/40)*k10+ (9/40)*k20, x1+ (3/40)*k11+ (9/40)*k21, y0, y1);
  k31=h*f1(t+ (3/10)*h, x0+ (3/40)*k10+ (9/40)*k20, x1+ (3/40)*k11+ (9/40)*k21, y0, y1);
  k40=h*f0(t+ (4/5)*h, x0+ (44/45)*k10- (56/15)*k20+ (32/9)*k30, x1+(44/45)*k11- (56/15)*k21+ (32/9)*k31, y0, y1);
  k41=h*f1(t+ (4/5)*h, x0+ (44/45)*k10- (56/15)*k20+ (32/9)*k30, x1+(44/45)*k11- (56/15)*k21+ (32/9)*k31, y0, y1);
  k50=h*f0(t+ (8/9)*h, x0+ (19372/6561)*k10- (25360/2187)*k20+ (64448/6561)*k30 -(212/729)*k40, x1+ (19372/6561)*k11- (25360/2187)*k21+ (64448/6561)*k31-(212/729)*k41, y0, y1);
  k51=h*f1(t+ (8/9)*h, x0+ (19372/6561)*k10- (25360/2187)*k20+ (64448/6561)*k30 -(212/729)*k40, x1+ (19372/6561)*k11- (25360/2187)*k21+ (64448/6561)*k31-(212/729)*k41, y0, y1);
  k60=h*f0(t+ h, x0+ (9017/3168)*k10- (355/33)*k20 -(46732/5247)*k30+ (49/176)*k40- (5103/18656)*k50, x1+ (9017/3168)*k11- (355/33)*k21 -(46732/5247)*k31+ (49/176)*k41- (5103/18656)*k51, y0, y1);
  k61=h*f1(t+h, x0+ (9017/3168)*k10- (355/33)*k20 -(46732/5247)*k30+ (49/176)*k40- (5103/18656)*k50, x1+ (9017/3168)*k11- (355/33)*k21 -(46732/5247)*k31+ (49/176)*k41- (5103/18656)*k51, y0, y1);
    
  x0= x0+ (35/384)*k10+ (500/1113)*k30+ (125/192)*k40- (2187/6784)*k50+ (11/84)*k60;
  x1= x1+ (35/384)*k11+ (500/1113)*k31+ (125/192)*k41- (2187/6784)*k51+ (11/84)*k61;

  p10=h*g0(t, y0, y1, x0, x1);
  p11=h*g1(t, y0, y1, x0, x1);
  p20=h*g0(t+ h/5, y0+ p10/5, y1+ p11/5, x0, x1);
  p21=h*g1(t+ h/5, y0+ p10/5, y1+ p11/5, x0, x1);
  p30=h*g0(t+ (3/10)*h, y0+ (3/40)*p10+ (9/40)*p20, y1+ (3/40)*p11+ (9/40)*p21, x0, x1);
  p31=h*g1(t+ (3/10)*h, y0+ (3/40)*p10+ (9/40)*p20, y1+ (3/40)*p11+ (9/40)*p21, x0, x1);
  p40=h*g0(t+ (4/5)*h, y0+ (44/45)*p10- (56/15)*p20+ (32/9)*p30, y1+(44/45)*p11- (56/15)*p21+ (32/9)*p31, x0, x1);
  p41=h*g1(t+ (4/5)*h, y0+ (44/45)*p10- (56/15)*p20+ (32/9)*p30, y1+(44/45)*p11- (56/15)*p21+ (32/9)*p31, x0, x1);
  p50=h*g0(t+ (8/9)*h, y0+ (19372/6561)*p10- (25360/2187)*p20+ (64448/6561)*p30 -(212/729)*p40, y1+ (19372/6561)*p11- (25360/2187)*p21+ (64448/6561)*p31-(212/729)*p41, x0, x1);
  p51=h*g1(t+ (8/9)*h, y0+ (19372/6561)*p10- (25360/2187)*p20+ (64448/6561)*p30 -(212/729)*p40, y1+ (19372/6561)*p11- (25360/2187)*p21+ (64448/6561)*p31-(212/729)*p41, x0, x1);
  p60=h*g0(t+h, y0+ (9017/3168)*p10- (355/33)*p20 -(46732/5247)*p30+ (49/176)*p40- (5103/18656)*p50, y1+ (9017/3168)*p11- (355/33)*p21 -(46732/5247)*p31+ (49/176)*p41- (5103/18656)*p51, x0, x1);
  p61=h*g1(t+h, y0+ (9017/3168)*p10- (355/33)*p20 -(46732/5247)*p30+ (49/176)*p40- (5103/18656)*p50, y1+ (9017/3168)*p11- (355/33)*p21-(46732/5247)*p31+ (49/176)*p41- (5103/18656)*p51, x0, x1);
  
  y0= y0+ (35/384)*p10+ (500/1113)*p30+ (125/192)*p40- (2187/6784)*p50+ (11/84)*p60;
  y1= y1+ (35/384)*p11+ (500/1113)*p31+ (125/192)*p41- (2187/6784)*p51+ (11/84)*p61;
}




