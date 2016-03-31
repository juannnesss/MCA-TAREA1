#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "riemann.h"


/*
 Referencias y proceso
 Riemann solvers and numerical methods for fluid dynamics
 Eleuterio F. Toro
 Chapter 11
 Computational Gasdynamics 
 Culbert Laney
 Chapter 5
 */

void Riemann(double *U1, double *U4, double *F)
{
    const double gamma = 1.4;
    double rr = U1[0];
    double ur = U1[1] / rr;
    double pr = (U1[2] - rr * ur * ur / 2.0) * (gamma - 1.0);
    double hr = (U1[2] + pr)/rr;
    double rl = U4[0];
    double ul = U4[1] / rl;
    double pl = (U4[2] - rl * ul * ul / 2.0) * (gamma - 1.0);
    double hl = (U4[2] + pl)/rl;
    
    
    
    // 1) Calculo promedios Roe 11.60
    
    double u_pro = (sqrt(rl)*ul+sqrt(rr)*ur)/(sqrt(rl)+sqrt(rr));
    double r_pro = sqrt(rr*rl);
    double h_pro = (sqrt(rl)*hl+sqrt(rr)*hr)/(sqrt(rl)+sqrt(rr));
    double a_pro = sqrt(fabs((gamma-1)*(h_pro-0.5*u_pro*u_pro)));
    
    // 2) Calcular los valores propios 11.58
   
    
    double lambda1=((u_pro-a_pro)<0.0)?(u_pro-a_pro):0.0 ;
    double lambda2=(u_pro<0)?u_pro:0.0 ;
    double lambda3=((u_pro+a_pro)<0.0)?(u_pro+a_pro):0.0 ;
    double XX=0.0;
    
    
    

    // 3) delta

    double dr = rl - rr;
    double dp = pl - pr;
    double du = ul - ur;
    
    // 4) wave strengs
    
    double dv1=dr-dp/(pow(a_pro,2));
    double dv2=du+dp/(r_pro*a_pro);
    double dv3=du-dp/(r_pro*a_pro);
    
    //5) flujos
   
    
    F[0] = rr*ur+lambda2*dv1+(0.5*r_pro/a_pro)*(lambda3*dv2-lambda1*dv3);
    
    
    F[1] = rr*(pow(ur,2)) +pr+u_pro*lambda2*dv1+(0.5*r_pro/a_pro)*((u_pro+a_pro)*lambda3*dv2-(u_pro-a_pro)*lambda1*dv3);
    
    
    
    F[2] = rr*hr*ur+(pow(u_pro,2))*lambda2*dv1*0.5+(0.5*r_pro/a_pro)*((h_pro + a_pro*u_pro)*lambda3*dv2 -(h_pro - a_pro*u_pro)*lambda1*dv3);
    

   

    
   
} 












