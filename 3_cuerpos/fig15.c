/*Este código imprime 4 columnas de datos*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void leap_frog (double h, double *q1, double *p1, double *q3, double *p3,double eps);
void RK4(double h, double *q1, double *p1, double *q3, double *p3, double eps);
double der_q1(double p1);
double der_p1(double q1, double eps);
double der_q3(double p3);
double der_p3(double q1, double q3, double eps);

int main ()
{
  /* Se crean las variables del hamiltoniano para resolver con Runge Kutta*/
  double q1_RK, q3_RK, p1_RK, p3_RK;
  
  /* Se crean las variables del hamiltoniano para resolver con el Integrador Simplectico*/
  double q1_LF, q3_LF, p1_LF, p3_LF;
  
  /* Distancia vertical entre los cuerpos masivos*/
  double eps=1.0;

  /* Tiempo total*/
  double T=2800;
  /* Paso de tiempo */
  double h=0.006;
  /* Numero de pasos*/
  int n=(int)(T/h);
  
  // Elementos aleatorios
  long seed;
  seed=n;
  srand48( seed );
  
  /* Condiciones iniciales*/
  q1_RK=q1_LF=0.425;
  q3_RK=q3_LF=2*drand48()-1;
  p1_RK=p1_LF=0.0;
  p3_RK=p3_LF=2*drand48()-1;
  double t=0.0;
    
  int i;
  
  
  for (i=0;i<n;i++)
    {
      printf("%.4f %.4f %.4f %.4f %.4f\n",q3_RK,p3_RK,q3_LF,p3_LF,t);
      RK4(h,&q1_RK,&p1_RK,&q3_RK,&p3_RK,eps);
      leap_frog(h,&q1_LF,&p1_LF,&q3_LF,&p3_LF,eps);
      t+=h;
    }
  return 0;
}

double der_q1(double p1)
{
  return p1;
}

double der_p1(double q1, double eps)
{
  return -2*q1*pow((4*pow(q1,2)+pow(eps,2)),-1.5);
}

double der_q3(double p3)
{
  return p3;
}

double der_p3(double q1, double q3, double eps)
{
  return (q1-q3)*pow((pow((q1-q3),2)+pow(eps,2)/4),-1.5) - (q1+q3)*pow((pow((q1+q3),2)+pow(eps,2)/4),-1.5);
}

void RK4(double h, double *q1, double *p1, double *q3, double *p3, double eps)
{
  /* Constantes RK4 para primera ecuacion diferencial*/
  double kq1, kq2, kq3, kq4;
  double kp1, kp2, kp3, kp4;
  /* Constantes RK4 para segunda ecuacion diferencial*/
  double lq1, lq2, lq3, lq4;
  double lp1, lp2, lp3, lp4;

  double nq1, nq3, np1, np3;
  nq1=*q1;
  nq3=*q3;
  np1=*p1;
  np3=*p3;

  kq1=h*der_q1(np1);
  kp1=h*der_p1(nq1,eps);
  lq1=h*der_q3(np3);
  lp1=h*der_p3(nq1, nq3, eps);

  kq2=h*der_q1(np1+kq1*0.5);
  kp2=h*der_p1(nq1+kp1*0.5,eps);
  lq2=h*der_q3(np3+lq1*0.5);
  lp2=h*der_p3(nq1+lp1*0.5, nq3+lp1*0.5, eps);

  kq3=h*der_q1(np1+kq2*0.5);
  kp3=h*der_p1(nq1+kp2*0.5,eps);
  lq3=h*der_q3(np3+lq2*0.5);
  lp3=h*der_p3(nq1+lp2*0.5, nq3+lp2*0.5, eps);

  kq4=h*der_q1(np1+kq3);
  kp4=h*der_p1(nq1+kp3,eps);
  lq4=h*der_q3(np3+lq3);
  lp4=h*der_p3(nq1+lp3, nq3+lp3, eps);

  nq1+= (kq1+2*kq2+2*kq3+kq4)/6.0;
  nq3+= (lq1+2*lq2+2*lq3+lq4)/6.0;
  np1+= (kp1+2*kp2+2*kp3+kp4)/6.0;
  np3+= (lp1+2*lp2+2*lp3+lp4)/6.0;

  *q1=nq1;
  *q3=nq3;
  *p1=np1;
  *p3=np3;
}

void leap_frog (double h, double *q1, double *p1, double *q3, double *p3, double eps)
{
  double a0=-1.702414384;
  double a1=1.351207192;
  double nq1, nq3, np1, np3;
  nq1=*q1;
  nq3=*q3;
  np1=*p1;
  np3=*p3;

    /*Coordenadas 1*/
  //drift
  np1+=0.5*h*a1*der_p1(nq1,eps);
  //kick
  nq1+=a1*h*der_q1(np1);
  //drift+drift
  np1+=0.5*h*der_p1(nq1,eps)*(a1+a0);
  //kick
  nq1+=a0*h*der_q1(np1);
  //drift+drift
  np1+=0.5*h*der_p1(nq1,eps)*(a0+a1);
  //kick
  nq1+=a1*h*der_q1(np1);
  //drift
  np1+=0.5*h*a1*der_p1(nq1,eps);

  /*Coordenadas 3*/
  //drift
  np3+=0.5*h*a1*der_p3(nq1,nq3,eps);
  //kick
  nq3+=a1*h*der_q3(np3);
  //drift+drift
  np3+=0.5*h*der_p3(nq1,nq3,eps)*(a1+a0);
  //kick
  nq3+=a0*h*der_q3(np3);
  //drift+drift
  np3+=0.5*h*der_p3(nq1,nq3,eps)*(a0+a1);
  //kick
  nq3+=a1*h*der_q3(np3);
  //drift
  np3+=0.5*h*a1*der_p3(nq1,nq3,eps);

  *q1=nq1;
  *q3=nq3;
  *p1=np1;
  *p3=np3;
}
