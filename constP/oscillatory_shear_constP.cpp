#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <float.h>
using namespace std;
#define Np 1024
#define Nn 200
#define dim 2
#define a1 1.0
#define a2 1.4
#define Pt 1.e-3
#define gamma0 1.e-7
#define omega 1.e-5
#define time_period 2.0*M_PI/omega
#define time_max 100.0*time_period

void update(double (*x_update)[dim],double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x_update[i][j]=x[i][j];
}

void calc_disp_max(double *disp_max,double (*x)[dim],double (*x_update)[dim],double L,double gamma){
  double dx,dy,dy_temp;
  double disp;
  for(int i=0;i<Np;i++){
    dx=x[i][0]-x_update[i][0];
    dy=x[i][1]-x_update[i][1];
    dy_temp=dy;
    dy-=L*floor((dy+0.5*L)/L);
    dx-=gamma*L*floor((dy_temp+0.5*L)/L);
    dx-=L*floor((dx+0.5*L)/L);
    disp = dx*dx+dy*dy;
    if(disp > *disp_max)
      *disp_max =disp;
  }
}

void output_shear1(int ns,double P,double G1,double G2){
  char filename[128];
  ofstream file;
  sprintf(filename,"time_shear.dat");
  file.open(filename, ios::app);
  file << fixed << setprecision(10) << ns << " " << P << " " << G1 << " " << G2 << endl;
  file.close();
}

double seed_change(int ensemble_num){
  srand(ensemble_num);
  return 0;
}

double unif_rand(double left, double right)
{
  return left + (right - left) * rand() / RAND_MAX;
}

double gaussian_rand(void)
{
  static double iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  if (iset == 0) {
    do {
      v1 = unif_rand(-1, 1);
      v2 = unif_rand(-1, 1);
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    iset = 0.50;
    return v2 * fac;
  }
  else {
    iset = 0;
    return gset;
  }
}

void ini_coord_rand(double (*x)[dim],double L){
  for(int i = 0 ; i < Np; i++){
    x[i][0] = L*unif_rand(0.,1.);
    x[i][1] = L*unif_rand(0.,1.);
  }
}

void set_diameter(double *a){
  int p=0;
  for(int i=0;i<Np;i++){
    if(p==0){
        a[i] = a1;
        p = p+1;
    }
    else{
        a[i] = a2;
        p = p-1;
    }
  }
}

void ini_array(double* x) {
  for (int j = 0; j < Np; j++){
    x[j] = 0.0;
  }
}

void ini_matrix(double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x[i][j]=0.0;
}

void gamma_update(double *gamma,double d_gamma,double dt){
  *gamma += d_gamma*dt;
  if(*gamma>1.0){
    *gamma = 0.0;
  }
}

int f(int i,int M)
{
  int k;
  
  k=i;
  
  if(k<0)
    k+=M;
  if(k>=M)
    k-=M;
  
  return k;
}

void cell_list(int (*list)[Nn],double (*x)[dim],int M,double L,double Rcell,double gamma){
  int i, j, k;
  int nx, ny;
  int l, m;
  double dx, dy, dy_temp,r;
  int(*map)[Np] = new int[M * M][Np];

  for (j = 0; j < M; j++) {
    for (i = 0; i < M; i++) {
      map[i + M * j][0] = 0;
    }
  }

  for (i = 0; i < Np; i++) {
    nx = f((int)((x[i][0]-gamma*x[i][1]) * M / L), M); 
    ny = f((int)(x[i][1] * M / L), M); 
    for (m = ny - 1; m <= ny + 1; m++) { 
      for (l = nx - 1; l <= nx + 1; l++) {
        map[f(l, M) + M * f(m, M)][map[f(l, M) + M * f(m, M)][0] + 1] = i;
        map[f(l, M) + M * f(m, M)][0]++;
        if (map[f(l, M) + M * f(m, M)][0] >= 499){
          printf("map error\n");
          delete[] map;
        }
      }
    }
  }

  for (i = 0; i < Np; i++) {
    list[i][0] = 0;
    nx = f((int)((x[i][0]-gamma*x[i][1]) * M / L), M); 
    ny = f((int)(x[i][1] * M / L), M);
    for (k = 1; k <= (map[nx + M * ny][0]); k++) { 
      j = map[nx + M * ny][k]; 
      if (j > i) {
        dx = x[i][0] - x[j][0];
        dy = x[i][1] - x[j][1];
        dy_temp = dy;
        dy-=L*floor((dy+0.5*L)/L);
        dx-=gamma*L*floor((dy_temp+0.5*L)/L);
        dx-=L*floor((dx+0.5*L)/L);
        r = dx * dx + dy * dy;
        if (r < Rcell * Rcell) {
          list[i][0]++;
          list[i][list[i][0]] = j;
        }
      }
    }
  }
  delete[] map;
}

void calc_force_hs(double (*x)[dim],double L,double* a,double (*f)[dim], double *U,int(*list)[Nn],double* P,double gamma,double *txy){
  int i, j, k;
  double r, t, dUr, dx, dy, aij, cut,dy_temp;
  ini_matrix(f);
  *P=0.0;
  *U = 0.0;
  *txy = 0.0;
  for (i = 0; i < Np; i++)
    {
      for (j = 1; j <= list[i][0]; j++)
        {
          dx = x[i][0] - x[list[i][j]][0];
          dy = x[i][1] - x[list[i][j]][1];
          dy_temp = dy;
          dy-=L*floor((dy+0.5*L)/L);
          dx-=gamma*L*floor((dy_temp+0.5*L)/L);
          dx-=L*floor((dx+0.5*L)/L); //periodic boundary condition of the distance between particles
          aij = (a[i] + a[list[i][j]]) / 2.0;
          r = sqrt(dx * dx + dy * dy);
          cut = aij;
          if (r < cut) {
            t = r / aij;
            dUr = -(1.0 - t) / aij; // analytical calculation of the 1'st derivative
            *P -= (dx * dx / r * dUr/L/L + dy *dy / r *dUr/L/L)/dim;
            *txy += (dx*dy/r*dUr/L/L);
            *U += (1.0 - t) * (1.0 - t) / 2.0/ double(Np);
          }
          else {
            dUr = 0.0;
            continue;
          }
          f[i][0] -= 1.0 * dUr * dx / r;
          f[i][1] -= 1.0 * dUr * dy / r;
          f[list[i][j]][0] += 1.0 * dUr * dx / r;
          f[list[i][j]][1] += 1.0 * dUr * dy / r;
        }
    }
} //calculation of the force and the stress

void control_P(double *L,double (*x)[dim],double dt,double *P,int M,double RCHK,double Rcell){
  int j,k;
  double zv=1.e-3,l,dL;
  l = *L;
  *L = sqrt((*L)*(*L)+dt*(*P-Pt)/zv); //systemsizeの変更, Pe is the target pressure
  dL = *L-l;

  for(j=0;j<Np;j++){
    for(k=0;k<dim;k++)
    x[j][k] *= (1.0+dL/(*L));
  } //affine transformation of coordinate
  
  M = (int)(*L / RCHK);
  Rcell = *L/M; //the length of the region

}

void p_boundary(double (*x)[dim],double L){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x[i][j]-=L*floor(x[i][j]/L);
}

void scale_sliding_blick(double (*x)[dim],double gamma,double L){
  for(int i=0;i<Np;i++){
    if(x[i][0]<gamma*(x[i][1]))
      x[i][0]+=L;
    if(x[i][0]>L+gamma*(x[i][1]))
      x[i][0]-=L;

    if(x[i][1]<0.0){
      x[i][1]+=L;
      x[i][0]+=gamma*L;
    }
    if(x[i][1]>L){
      x[i][1]-=L;
      x[i][0]-=gamma*L;
    }
  }
}

void eq_motion(double (*x)[dim],double dt,double (*f)[dim],double L){
  int j,k;
  for(j=0;j<Np;j++){
    for(k=0;k<dim;k++){
      x[j][k] += f[j][k]*dt;
    }
  }
  p_boundary(x,L);
}

void eq_motion_shear(double (*x)[dim],double (*v)[dim],double dt,double (*f)[dim],double L,double gamma,double d_gamma){
  for(int i=0;i<Np;i++){
    v[i][0] = d_gamma*x[i][1]+f[i][0];
    v[i][1] = f[i][1];
    for(int j=0;j<dim;j++){
      x[i][j] += v[i][j]*dt;
    }
  }
  scale_sliding_blick(x,gamma,L);
}

void auto_list_update(double *disp_max,double (*x)[dim],double (*x_update)[dim],int (*list)[Nn],int M,double L,double Rcell,double gamma){
  calc_disp_max(&(*disp_max),x,x_update,L,gamma);
  if(*disp_max>0.05*0.05){
    cell_list(list,x,M,L,Rcell,gamma);
    *disp_max=0.0;
    update(x_update,x);
  }
}

void calc_shear(double *txy,double t,double dt,double *G1,double *G2){
  *G1 += (*txy)*sin(omega*t)*dt/M_PI/gamma0*omega;
  *G2 += (*txy)*cos(omega*t)*dt/M_PI/gamma0*omega;
}

void output_shear2(double t,double gamma,double txy,double P){
  char filename2[128];
  ofstream file2;
  sprintf(filename2,"gamma_sigma.dat");
  file2.open(filename2,ios::app);
  file2 << t << " " << gamma << " " << P << " " << txy << endl;
  file2.close();
}

int main()
{
  double x[Np][dim],x_update[Np][dim],v[Np][dim],f[Np][dim],a[Np];
  int list[Np][Nn];
  double txy=0.0,txy0=0.0,P=0.0,U=0.0,dt=0.01,phi0=0.86,gap=0.0,t,gamma=0.0,d_gamma,G1=0.0,G2=0.0;
  double L = sqrt(M_PI*Np/2.0*(a1*a1+a2*a2)/4.0/phi0),disp_max=0.0;
  int count=0,ns=0,cn=0;
  char filename[128];
  ofstream file;
  
  double RCHK = 3.0;
  int M = (int)(L / RCHK);
  double Rcell = L/M;

  set_diameter(a);
  ini_coord_rand(x,L);
  ini_matrix(v);
  cell_list(list,x,M,L,Rcell,0.0);
  calc_force_hs(x,L,a,f,&U,list,&P,0.0,&txy);

  for(;;){
    t += dt;
    control_P(&L,x,dt,&P,M,RCHK,Rcell); //systemsizeの変更
    auto_list_update(&disp_max,x,x_update,list,M,L,Rcell,0.0);
    calc_force_hs(x,L,a,f,&U,list,&P,0.0,&txy);
    eq_motion(x,dt,f,L); //EoM

    gap = P-Pt;
    if(gap<0.0){
      gap = -gap;
    }
    if(gap<1.e-7){
      break;
    }
  }

  t=0.0;
  dt = 0.1;
  txy0 = txy;
  for(t=dt;t<time_max;t+=dt){
    count++;
    d_gamma = gamma0*omega*cos(omega*t);
    gamma_update(&gamma,d_gamma,dt);
    control_P(&L,x,dt,&P,M,RCHK,Rcell); //systemsizeの変更
    auto_list_update(&disp_max,x,x_update,list,M,L,Rcell,gamma);
    calc_force_hs(x,L,a,f,&U,list,&P,gamma,&txy);
    eq_motion_shear(x,v,dt,f,L,gamma,d_gamma);
    calc_shear(&txy,t,dt,&G1,&G2);
    if(count==int(time_period/dt)){
      ns++;
      output_shear1(ns,P,G1,G2);
      G1 = 0.0;
      G2 = 0.0;
      count = 0;
    }
  }

  t=0.0;
  for(t=dt;t<time_period;t+=dt){
    cn++;
    d_gamma = gamma0*omega*cos(omega*t);
    gamma_update(&gamma,d_gamma,dt);
    control_P(&L,x,dt,&P,M,RCHK,Rcell); //systemsizeの変更
    auto_list_update(&disp_max,x,x_update,list,M,L,Rcell,gamma);
    calc_force_hs(x,L,a,f,&U,list,&P,gamma,&txy);
    eq_motion_shear(x,v,dt,f,L,gamma,d_gamma);
    calc_shear(&txy,t,dt,&G1,&G2);
    if(cn==int(time_period/100.0/dt)){
      output_shear2(t,gamma,txy,P);
      cn=0;
    }
  }

  return 0;
}
