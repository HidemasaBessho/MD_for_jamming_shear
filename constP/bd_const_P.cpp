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

void update(double (*x_update)[dim],double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x_update[i][j]=x[i][j];
}

void calc_disp_max(double *disp_max,double (*x)[dim],double (*x_update)[dim],double L){
  double dx,dy;
  double disp;
  for(int i=0;i<Np;i++){
    dx=x[i][0]-x_update[i][0];
    dy=x[i][1]-x_update[i][1];
    dx-=L*floor((dx+0.5*L)/L);
    dy-=L*floor((dy+0.5*L)/L);
    disp = dx*dx+dy*dy;
    if(disp > *disp_max)
      *disp_max =disp;
  }
}

void output_stress(double avP,double t,double L){
  char filename[128];
  ofstream file;
  sprintf(filename,"FIRE_stess_N%d_P%.4f.dat",Np,Pe);
  file.open(filename, ios::app);
  file << fixed << setprecision(10) << t << " " << avP << " " << L << endl;
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

void cell_list(int (*list)[Nn],double (*x)[dim],int M,double L,double Rcell) {
  int i, j, k;
  int nx, ny;
  int l, m;
  double dx, dy, r;
  int(*map)[Np] = new int[M * M][Np];

  for (j = 0; j < M; j++) {
    for (i = 0; i < M; i++) {
      map[i + M * j][0] = 0;
    }
  }

  for (i = 0; i < Np; i++) {
    nx = f((int)(x[i][0] * M / L), M); 
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
    nx = f((int)(x[i][0] * M / L), M); 
    ny = f((int)(x[i][1] * M / L), M);
    for (k = 1; k <= (map[nx + M * ny][0]); k++) { 
      j = map[nx + M * ny][k]; 
      if (j > i) {
        dx = x[i][0] - x[j][0];
        dy = x[i][1] - x[j][1];
        if (dy < (-L / 2.0)) {
          dy += L;
        }
        else if (dy > (L / 2.0)) {
          dy -= L;
        } //calculate y-direction first!
        if (dx < (-L / 2.0))
          dx += L;
        else if (dx > (L / 2.0))
          dx -= L;
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

void calc_force_hs(double (*x)[dim],double L,double* a,double (*f)[dim], double *U,int(*list)[Nn],double* P){
  int i, j, k;
  double r, t, dUr, dx, dy, aij, cut;
  ini_matrix(f);
  *P=0.0;
  *U = 0.0;
  for (i = 0; i < Np; i++)
    {
      for (j = 1; j <= list[i][0]; j++)
        {
          dx = x[i][0] - x[list[i][j]][0];
          dy = x[i][1] - x[list[i][j]][1];
          if (dy > 0.5 * L) {
            dy -= L;
          }
          else if (dy < -0.5 * L) {
            dy += L;
          } //calculate y-direction first!
          if (dx > 0.5 * L)
            dx -= L;
          else if (dx < -0.5 * L)
            dx += L; //periodic boundary condition of the distance between particles
          aij = (a[i] + a[list[i][j]]) / 2.0;
          r = sqrt(dx * dx + dy * dy);
          cut = aij;
          if (r < cut) {
            t = r / aij;
            dUr = -(1.0 - t) / aij; // analytical calculation of the 1'st derivative
            *P -= (dx * dx / r * dUr/L/L + dy *dy / r *dUr/L/L)/dim;
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

void size_change(double *L,double (*x)[dim],double dt,double *P,int M,double RCHK,double Rcell){
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

void eq_motion(double (*x)[dim],double dt,double (*f)[dim],double L){
  int j,k;
  for(j=0;j<Np;j++){
    for(k=0;k<dim;k++){
      x[j][k] += f[j][k]*dt;
    }
  }
  p_boundary(x,L);
}

void auto_list_update(double *disp_max,double (*x)[dim],double (*x_update)[dim],int (*list)[Nn],int M,double L,double Rcell){
  calc_disp_max(&(*disp_max),x,x_update,L);
  if(*disp_max>0.05*0.05){
    cell_list(list,x,M,L,Rcell);
    *disp_max=0.0;
    update(x_update,x);
  }
}

int main()
{
  double x[Np][dim],x_update[Np][dim],v[Np][dim],f[Np][dim],a[Np];
  int list[Np][Nn];
  double txy=0.0,P=0.0,U=0.0,dt=0.01,phi0=0.86,gap=0.0,t;
  double L = sqrt(M_PI*Np/2.0*(a1*a1+a2*a2)/4.0/phi0),disp_max=0.0,sampling_time=5.0*dt;
  char filename[128];
  ofstream file;
  
  double RCHK = 3.0;
  int M = (int)(L / RCHK);
  double Rcell = L/M;

  set_diameter(a);
  ini_coord_rand(x,L);
  ini_matrix(v);
  cell_list(list,x,M,L,Rcell);
  calc_force_hs(x,L,a,f,&U,list,&P);

  //output_stress(txy,Np,Pe,U,P,t,L);
  sampling_time=5.*dt;
  for(;;){
    t += dt;
    size_change(&L,x,dt,&P,M,RCHK,Rcell); //systemsizeの変更
    auto_list_update(&disp_max,x,x_update,list,M,L,Rcell);
    calc_force_hs(x,L,a,f,&U,list,&P);
    eq_motion(x,dt,f,L); //EoM

    //logarithmic sampling
    if(int(t/dt) == int(sampling_time/dt)){
      calc_force_hs(x,L,a,f,&U,list,&P);
      output_stress(P,t,L);
      sampling_time*=pow(10.,0.1);
      sampling_time=int(sampling_time/dt)*dt;
    }
    gap = P-Pe;
    if(gap<0.0){
      gap = -gap;
    }
    if(gap<1.e-7){
      output_stress(P,t,L);
      break;
    }
  }
  return 0;
}
