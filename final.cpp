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
  for(time=0.0; time<=(10*T); time+=delta)
    {
      std::cout<<x<<"\t"<<vx<<std::endl;
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
  k30=h*f0(t+ 3*(h/10), x0+ 3*(k10/40)+ 9*(k20/40), x1+ 3*(k11/40)+ 9*(k21/40), y0, y1);
  k31=h*f1(t+ 3*(h/10), x0+ 3*(k10/40)+ 9*(k20/40), x1+ 3*(k11/40)+ 9*(k21/40), y0, y1);
  k40=h*f0(t+ 4*(h/5), x0+ 44*(k10/45)- 56*(k20/15)+ 32*(k30/9), x1+ 44*(k11/45)- 56*(k21/15)+ 32*(k31/9), y0, y1);
  k41=h*f1(t+ 4*(h/5), x0+ 44*(k10/45)- 56*(k20/15)+ 32*(k30/9), x1+ 44*(k11/45)- 56*(k21/15)+ 32*(k31/9), y0, y1);
  k50=h*f0(t+ 8*(h/9), x0+ 19372*(k10/6561)- 25360*(k20/2187)+ 64448*(k30/6561)- 212*(k40/729), x1+ 19372*(k11/6561)- 25360*(k21/2187)+ 64448*(k31/6561)- 212*(k41/729), y0, y1);
  k51=h*f1(t+ 8*(h/9), x0+ 19372*(k10/6561)- 25360*(k20/2187)+ 64448*(k30/6561)- 212*(k40/729), x1+ 19372*(k11/6561)- 25360*(k21/2187)+ 64448*(k31/6561)- 212*(k41/729), y0, y1);
  k60=h*f0(t+ h, x0+ 9017*(k10/3168)- 355*(k20/33)- 46732*(k30/5247)+ 49*(k40/176)- 5103*(k50/18656), x1+ 9017*(k11/3168)- 355*(k21/33)- 46732*(k31/5247)+ 49*(k41/176)- 5103*(k51/18656), y0, y1);
  k61=h*f1(t+h, x0+ 9017*(k10/3168)- 355*(k20/33)- 46732*(k30/5247)+ 49*(k40/176)- 5103*(k50/18656), x1+ 9017*(k11/3168)- 355*(k21/33)- 46732*(k31/5247)+ 49*(k41/176)- 5103*(k51/18656), y0, y1);
    
  x0= x0+ 35*(k10/384)+ 500*(k30/1113)+ 125*(k40/192)- 2187*(k50/6784)+ 11*(k60/84);
  x1= x1+ 35*(k11/384)+ 500*(k31/1113)+ 125*(k41/192)- 2187*(k50/6784)+ 11*(k61/84);

  p10=h*g0(t, y0, y1, x0, x1);
  p11=h*g1(t, y0, y1, x0, x1);
  p20=h*g0(t+ h/5, y0+ p10/5, y1+ p11/5, x0, x1);
  p21=h*g1(t+ h/5, y0+ p10/5, y1+ p11/5, x0, x1);
  p30=h*g0(t+ 3*(h/10), y0+ 3*(p10/40)+ 9*(p20/40), y1+ 3*(p11/40)+ 9*(p21/40), x0, x1);
  p31=h*g1(t+ 3*(h/10), y0+ 3*(p10/40)+ 9*(p20/40), y1+ 3*(p11/40)+ 9*(p21/40), x0, x1);
  p40=h*g0(t+ 4*(h/5), y0+ 44*(p10/45)- 56*(p20/15)+ 32*(p30/9), y1+ 44*(p11/45)- 56*(p21/15)+ 32*(p31/9), x0, x1);
  p41=h*g1(t+ 4*(h/5), y0+ 44*(p10/45)- 56*(p20/15)+ 32*(p30/9), y1+ 44*(p11/45)- 56*(p21/15)+ 32*(p31/9), x0, x1);
  p50=h*g0(t+ 8*(h/9), y0+ 19372*(p10/6561)- 25360*(p20/2187)+ 64448*(p30/6561)- 212*(p40/729), y1+ 19372*(p11/6561)- 25360*(p21/2187)+ 64448*(p31/6561)- 212*(p41/729), x0, x1);
  p51=h*g1(t+ 8*(h/9), y0+ 19372*(p10/6561)- 25360*(p20/2187)+ 64448*(p30/6561)- 212*(p40/729), y1+ 19372*(p11/6561)- 25360*(p21/2187)+ 64448*(p31/6561)- 212*(p41/729), x0, x1);
  p60=h*g0(t+h, y0+ 9017*(p10/3168)- 355*(p20/33)- 46732*(p30/5247)+ 49*(p40/176)- 5103*(p50/18656), y1+ 9017*(p11/3168)- 355*(p21/33)- 46732*(p31/5247)+ 49*(p41/176)- 5103*(p51/18656), x0, x1);
  p61=h*g1(t+h, y0+ 9017*(p10/3168)- 355*(p20/33)- 46732*(p30/5247)+ 49*(p40/176)- 5103*(p50/18656), y1+ 9017*(p11/3168)- 355*(p21/33)- 46732*(p31/5247)+ 49*(p41/176)- 5103*(p51/18656), x0, x1);
  
  y0= y0+ 35*(p10/384)+ 500*(p30/1113)+ 125*(p40/192)- 2187*(p50/6784)+ 11*(p60/84);
  y1= y1+ 35*(p11/384)+ 500*(p31/1113)+ 125*(p41/192)- 2187*(p51/6784)+ 11*(p61/84);
}




