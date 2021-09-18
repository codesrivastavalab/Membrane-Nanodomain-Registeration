#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// System Variables

  //Dimensions of each stacked sheet
  int m = 100;
  int n = 100;

  // Full set of parameters

  double epsilon = 0.3;                                                         // 0.3 kcal/mol : Temperature scaling constant
  double Temp[1] = {2.048};                                                     // Scaled Temperature
  double Vpaa[5] = {4.0, 7.0, 10.0, 13.0, 16.0};                                // -----------------------------------------------------
  double Vpbb[5] = {4.0, 7.0, 10.0, 13.0, 16.0};                                // ------- In-plane enthalpic strength constants -------
  double Vpab[5] = {4.0, 7.0, 10.0, 13.0, 16.0};                                // ---------------------------------------------------//
  double Viaa[3] = {8.0, 20.0, 32.0};                                           // -----------------------------------------------------
  double Viab[3] = {8.0, 20.0, 32.0};                                           // ----- Interleaflet enthalpic strength constants -----
  double Vibb[3] = {8.0, 20.0, 32.0};                                           // ---------------------------------------------------//
  double VSp[1] = {0.06};                                                       // In-plane entropic term strength constants
  double VSi[4] = {0.005, 0.02, 0.035, 0.050};                                  // Interleaflet entropic term strength constants

  int population[2] = {5000,5000};                                              // number of lipids of each species in order : DPPC, D34
  double salpha[2] = {0.3565843, 0.2400412};                                    // average site variable of species in same order
  //double length[2] = {16, 22};                                                  // length of tail in the same order (currently unused)
  double unsat[2] = {0.0, 1.3181818};                                           // sum(pos)/numlength for each species in the same order

  //Monte Carlo simulation parameters
  int steps = 1000000000;                                                       // steps per set of variables

  //Output parameters
  int interval = 50000;                                                         // print the system configuration every 'interval' steps

  //Current Parameters
    // Note : these values are updated from the above full set of parameters as the runs progress

    double kT = 0.0;                                                            // Variable holding current scaled temperature
    double Vplane[2][2] = {                                                     // Array holding current values of Vplane_(alpha,alpha')
      {0.0, 0.0},
      {0.0, 0.0}
    };
    double Vinter[2][2] = {                                                     // Array holding current values of Vinter_(alpha,alpha')
      {0.0, 0.0},
      {0.0, 0.0}
    };
    double VSplane = 0.0;                                                       // Variable holding current in-plane entropic strength constant
    double VSinter = 0.0;                                                       // Variable holding current interleaflet entropic strength constant

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Function Declarations

/* unif. dist. from 0.0 to 1.0 */
double rand_double();

/* discrete integer unif. dist. from a to b */
int uniform_distribution(int rangeLow, int rangeHigh);

/* calculates total energy of system */
double energy_calc(int** upper, int** lower);

/* calculates difference in energy upon switching (x1,y1) with (x2,y2) */
double mc_move(int** upper, int** lower, int x1, int y1, int x2, int y2, int flag);

/* calculates local energy of a chosen site */
double local_energy(int** upper, int** lower, int x1, int y1, int* x1n, int* y1n, int flag);

/* calculates Ni, phii */
void entropic_term(double* count, double* frac, int* x1n, int* y1n, int x1, int y1, int** upper, int** lower, int leaflet, int flag);

/* gives PBC conforming neighbours of input lattice site */
void pbc_neighbours(int x1, int y1, int* x, int* y);

/* prints the upper and lower leaflets when called */
void print_config(int** upper, int** lower, FILE* fp);

/* sets global variables based on input */
void set_global_params(int tem, int Spl, int Sin, int Vi00, int Vi11, int Vi10, int Vp00, int Vp11, int Vp10);


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Main

void main(){

  // Variable Declarations
  int i, j, t, type, x1, y1, x2, y2, r, quot;
  double energy, delta;

  int tem = 0, Spl = 0, Sin = 1, Vi00 = 0, Vi11 = 2, Vi10 = 2, Vp00 = 4, Vp11 = 2, Vp10 = 0;
  // Setting parameters for run
  set_global_params(tem, Spl, Sin, Vi00, Vi11, Vi10, Vp00, Vp11, Vp10);

  energy = 0.0; delta = 0.0;

  // Config Output File Initialization
  char buffer[32];                                                    // The filename buffer.
  snprintf(buffer,sizeof(char)*32,"Energy-%d-%d-%d-%d-%d-%d-%d-%d-%d.txt",tem,Spl,Sin,Vi00,Vi11,Vi10,Vp00,Vp11,Vp10);
  FILE *fp;
  fp = fopen(buffer, "w+");
  snprintf(buffer,sizeof(char)*32,"%d-%d-%d-%d-%d-%d-%d-%d-%d.txt",tem,Spl,Sin,Vi00,Vi11,Vi10,Vp00,Vp11,Vp10);
  FILE *cf;
  cf = fopen(buffer, "w+");

  // Memory Allocation for lattice layers
  int** upper = (int**)malloc(m*sizeof(int*));                        // 2D array for upper leaflet
  for(i=0; i<m; i++){
    upper[i] = (int*)malloc(n*sizeof(int));
  }
  int** lower = (int**)malloc(m*sizeof(int*));                        // 2D array for lower leaflet
  for(i=0; i<m; i++){
    lower[i] = (int*)malloc(n*sizeof(int));
  }

  // resets the random number seed (otherwise the same random numbers are generated in each run)
  srand((long)time(NULL));

  // System Initialization
  int count[2] = {0,0};                      // Both upper and lower leaflet are randomly populated with lipids of each species,
  for(i=0;i<m;i++){                          // so that the filled lattice has the required populations of each species
    for(j=0;j<n;j++){
      type = uniform_distribution(0,1);
      if(count[type]>=population[type]){
        upper[i][j] = 1 - type;
        count[1-type] = count[1-type] + 1;
      }
      else{
        upper[i][j] = type;
        count[type] = count[type] + 1;
      }
    }
  }
  count[0] = 0; count[1] = 0;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      type = uniform_distribution(0,1);
      if(count[type]>=population[type]){
        lower[i][j] = 1 - type;
        count[1-type] = count[1-type] + 1;
      }
      else{
        lower[i][j] = type;
        count[type] = count[type] + 1;
      }
    }
  }

  // Initial total energy calculation
  energy = energy_calc(upper,lower);

  // MC loop
  t = 0;
  while(t<steps){
    if(t%interval==0){                                                // prints system config and energy every 'interval' steps
      fprintf(fp,"%d %lf\n",t,energy);
      fprintf(cf,"Move %d, Energy = %lf\n",t,energy);
      print_config(upper, lower, cf);
    }
    //trial move for upper leaflet
    x1 = uniform_distribution(0,m-1);
    y1 = uniform_distribution(0,n-1);
    r=0;
    while(r!=1){                                                      // selects a random site distinct from the first one
      x2 = uniform_distribution(0, m-1);
      y2 = uniform_distribution(0, n-1);
      if(x2==x1 && y2==y1){r=0;}
      else if(upper[x1][y1]==upper[x2][y2]){r=0;}
      else{r=1;}
    }
    delta = mc_move(upper, lower, x1, y1, x2, y2, 1);                 // the last argument is a flag specifying the leaflet of the chosen site
    energy = energy + delta;                                          // 0 = lower, 1 = upper

    //trial move for lower leaflet
    x1 = uniform_distribution(0,m-1);
    y1 = uniform_distribution(0,n-1);
    r=0;
    while(r!=1){                                                      // selects a random site distinct from the first point
      x2 = uniform_distribution(0, m-1);
      y2 = uniform_distribution(0, n-1);
      if(x2==x1 && y2==y1){r=0;}
      else if(lower[x1][y1]==lower[x2][y2]){r=0;}
      else{r=1;}
    }
    delta = mc_move(upper, lower, x1, y1, x2, y2, 0);
    energy = energy + delta;

    if(t%5000==0){                                                // corrector to prevent drift (if required)
      energy = energy_calc(upper,lower);
    }
    t++;

    // Percent Completion Indicator for testing
    if(t%10000000==0){
      quot = t/10000000;
      if(quot==100){
        printf("%d%%\n\n",quot);
      }
      else{
        printf("%d%% ",quot);
      }
    }
  }

  // Wrap Up
  fclose(fp);                                                         // closes text file used for config output
  for (i = 0; i < m; i++){free(upper[i]); free(lower[i]);}            // freeing allocated 2D arrays for upper and lower leaflet config
  free(upper); free(lower);
// End of runs with all required parameter combinations
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Function Definitions

double rand_double()
{
   return rand()/(double)RAND_MAX;
}

int uniform_distribution(int rangeLow, int rangeHigh)
{
  int myRand = (int)rand();
  int range = rangeHigh - rangeLow + 1; //+1 makes it [rangeLow, rangeHigh], inclusive.
  int myRand_scaled = (myRand % range) + rangeLow;
  return myRand_scaled;
}

double energy_calc(int** upper,int** lower)
{
  int i, j, k, tu, tl, t, x, y;
  int in[8], jn[8];
  double temp = 0.0;
  double count[2], frac[2];
  double H, Hp, Hi;
  H = 0.0;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      temp = 0.0; Hp = 0.0; Hi = 0.0;
      tu = upper[i][j];
      tl = lower[i][j];
      pbc_neighbours(i, j, in, jn);         // in and jn array indices - 0 :: right, 1 :: down, 2 :: left, 3 :: up,
                                            // 4 :: bottom-right, 5 :: bottom-left, 6 :: top-left, 7 :: top-right.
      for(k=0;k<4;k++){
        x = in[k]; y = jn[k];
        t = lower[x][y];
        temp = temp + Vplane[tl][t]*salpha[tl]*salpha[t];
        t = upper[x][y];
        temp = temp + Vplane[tu][t]*salpha[tu]*salpha[t];
      }
      Hp = Hp - 0.5*temp;
      temp = 0.0;
      count[0] = 0.0; count[1] = 0.0; frac[0] = 0.0; frac[1] = 0.0;
      entropic_term(count, frac, in, jn, i, j, upper, lower, 0, 0);             // the last argument means : 0 = in-plane term, 1 = interleaflet term
      temp = temp + VSplane*kT*((1 + unsat[0])*count[0]*log(frac[0]) + (1 + unsat[1])*count[1]*log(frac[1]))/4.0;
      count[0] = 0.0; count[1] = 0.0; frac[0] = 0.0; frac[1] = 0.0;
      entropic_term(count, frac, in, jn, i, j, upper, lower, 1, 0);             // the last argument means : 0 = in-plane term, 1 = interleaflet term
      temp = temp + VSplane*kT*((1 + unsat[0])*count[0]*log(frac[0]) + (1 + unsat[1])*count[1]*log(frac[1]))/4.0;
      Hp = Hp + temp;

      Hi = Hi - Vinter[tu][tl]*salpha[tu]*salpha[tl];
      temp = 0.0;
      count[0] = 0.0; count[1] = 0.0; frac[0] = 0.0; frac[1] = 0.0;
      entropic_term(count, frac, in, jn, i, j, upper, lower, 0, 1);
      temp = temp + VSinter*kT*((1 + unsat[0])*count[0]*log(frac[0]) + (1 + unsat[1])*count[1]*log(frac[1]));
      count[0] = 0.0; count[1] = 0.0; frac[0] = 0.0; frac[1] = 0.0;
      entropic_term(count, frac, in, jn, i, j, upper, lower, 1, 1);
      temp = temp + VSinter*kT*((1 + unsat[0])*count[0]*log(frac[0]) + (1 + unsat[1])*count[1]*log(frac[1]));
      Hi = Hi + temp;

      H = H + Hp + Hi;
    }
  }
  H = H/epsilon;
  return H;
}

double mc_move(int** upper, int** lower, int x1, int y1, int x2, int y2, int flag)
{
  int tu1, tl1, tu2, tl2, t;
  double temp, Ei, Ef, r;
  double count[2], frac[2];
  int x1n[8], y1n[8], x2n[8], y2n[8];
  pbc_neighbours(x1, y1, x1n, y1n);          // xin and yin array indices - 0 :: right, 1 :: down, 2 :: left, 3 :: up,
  pbc_neighbours(x2, y2, x2n, y2n);          // 4 :: bottom-right, 5 :: bottom-left, 6 :: top-left, 7 :: top-right.

  Ei = 0.0; Ef = 0.0;
  tu1 = upper[x1][y1]; tl1 = lower[x1][y1];
  tu2 = upper[x2][y2]; tl2 = lower[x2][y2];

  // switching attempt for lower leaflet
  // Function call to calculate local energy of initial config
  Ei = Ei + local_energy(upper, lower, x1, y1, x1n, y1n, flag);
  Ei = Ei + local_energy(upper, lower, x2, y2, x2n, y2n, flag);

  // switching the lipids at sites (x1,y1) and (x2,y2) in the correct leaflet
  if(flag==0){
    lower[x2][y2] = tl1; lower[x1][y1] = tl2;
    tl1 = lower[x1][y1]; tl2 = lower[x2][y2];
  }
  if(flag==1){
    upper[x2][y2] = tu1; upper[x1][y1] = tu2;
    tu1 = upper[x1][y1]; tu2 = upper[x2][y2];
  }

  // Function call to calculate local energy of final config
  Ef = Ef + local_energy(upper, lower, x1, y1, x1n, y1n, flag);
  Ef = Ef + local_energy(upper, lower, x2, y2, x2n, y2n, flag);

  // Metropolis acceptance
  r = rand_double();
  temp = exp(-(Ef-Ei)/(epsilon*kT));
  if(r<=temp){
    return (Ef - Ei)/epsilon;
  }
  else{ // switch back the interchanged sites
    if(flag==0){
      lower[x2][y2] = tl1; lower[x1][y1] = tl2;
    }
    if(flag==1){
      upper[x2][y2] = tu1; upper[x1][y1] = tu2;
    }
    return 0.0;
  }
}

double local_energy(int** upper, int** lower, int x1, int y1, int* x1n, int* y1n, int flag){
  int i, j, x, y, u1, l1, xt, yt, t;
  int  xn[8], yn[8];
  double temp, Et;
  double count[2], frac[2];

  Et = 0.0;
  for(i=0;i<9;i++){
    if(i==8){
      x = x1; y = y1;                                            // starting calculation of energy due to (x1,y1)
    }
    else{
      x = x1n[i]; y = y1n[i];                                    // starting calculation of energy due to neighbours of (x1,y1)
    }
    pbc_neighbours(x, y, xn, yn);
    u1 = upper[x][y]; l1 = lower[x][y];

    temp = 0.0;
    // planar and interleaflet enthalpic interaction is only required for the site (x1,y1) and not neighbours
    if(i==8){
      for(j=0;j<4;j++){
        xt = xn[j]; yt = yn[j];
        if(flag==0){
          t = lower[xt][yt];
          temp = temp + Vplane[l1][t]*salpha[l1]*salpha[t];
        }
        if(flag==1){
          t = upper[xt][yt];
          temp = temp + Vplane[u1][t]*salpha[u1]*salpha[t];
        }
      }
      temp = temp + Vinter[u1][l1]*salpha[u1]*salpha[l1];
    }
    // planar and interleaflet entropic interaction calculation for (x,y) and all 8 neighbours
    count[0] = 0.0; count[1] = 0.0; frac[0] = 0.0; frac[1] = 0.0;
    entropic_term(count, frac, xn, yn, x, y, upper, lower, flag, 0);        // the last argument means : 0 = in-plane term, 1 = interleaflet term
    temp = temp - VSplane*kT*((1 + unsat[0])*count[0]*log(frac[0]) + (1 + unsat[1])*count[1]*log(frac[1]))/4.0;
    count[0] = 0.0; count[1] = 0.0; frac[0] = 0.0; frac[1] = 0.0;
    entropic_term(count, frac, xn, yn, x, y, upper, lower, flag, 1);
    temp = temp - VSinter*kT*((1 + unsat[0])*count[0]*log(frac[0]) + (1 + unsat[1])*count[1]*log(frac[1]));
    // interleaflet entropic interaction of neighbours of (x,y) in the opposite leaflet
    if(i<=3 || i==8){
      count[0] = 0.0; count[1] = 0.0; frac[0] = 0.0; frac[1] = 0.0;
      int f = 1 - flag;
      entropic_term(count, frac, xn, yn, x, y, upper, lower, f, 1);
      temp = temp - VSinter*kT*((1 + unsat[0])*count[0]*log(frac[0]) + (1 + unsat[1])*count[1]*log(frac[1]));
    }
    Et = Et - temp;
  }
  return Et;
}

void entropic_term(double* count, double* frac, int* x1n, int* y1n, int x1, int y1, int** upper, int** lower, int leaflet, int flag)
{
  int i, x, y;
  if(leaflet == 0){
    if(lower[x1][y1] == 0){
      count[0] = count[0] + 1.0;
    }
    if(lower[x1][y1] == 1){
      count[1] = count[1] + 1.0;
    }
    if(flag == 0){
      for(i=0;i<8;i++){
        x = x1n[i]; y = y1n[i];
        if(lower[x][y] == 0){
          count[0] = count[0] + 1.0;
        }
        else{
          count[1] = count[1] + 1.0;
        }
      }
    }
    if(flag == 1){
      if(upper[x1][y1] == 0){
        count[0] = count[0] + 1.0;
      }
      if(upper[x1][y1] == 1){
        count[1] = count[1] + 1.0;
      }
      for(i=0;i<4;i++){
        x = x1n[i]; y = y1n[i];
        if(upper[x][y] == 0){
          count[0] = count[0] + 1.0;
        }
        else{
          count[1] = count[1] + 1.0;
        }
      }
    }
  }
  if(leaflet == 1){
    if(upper[x1][y1] == 0){
      count[0] = count[0] + 1.0;
    }
    if(upper[x1][y1] == 1){
      count[1] = count[1] + 1.0;
    }
    if(flag == 0){
      for(i=0;i<8;i++){
        x = x1n[i]; y = y1n[i];
        if(upper[x][y] == 0){
          count[0] = count[0] + 1.0;
        }
        else{
          count[1] = count[1] + 1.0;
        }
      }
    }
    if(flag == 1){
      if(lower[x1][y1] == 0){
        count[0] = count[0] + 1.0;
      }
      if(lower[x1][y1] == 1){
        count[1] = count[1] + 1.0;
      }
      for(i=0;i<4;i++){
        x = x1n[i]; y = y1n[i];
        if(lower[x][y] == 0){
          count[0] = count[0] + 1.0;
        }
        else{
          count[1] = count[1] + 1.0;
        }
      }
    }
  }
  frac[0] = count[0]/(count[0] + count[1]);
  frac[1] = count[1]/(count[0] + count[1]);
  if(count[0]==0.0){                                                            // doesn't change the energy. This is simply to
    frac[0] = 1.0;                                                              // prevent undefined logarithm calculations
  }
  if(count[1]==0.0){
    frac[1] = 1.0;
  }
}

void pbc_neighbours(int x1, int y1, int* x, int* y)
{
  int i;

  // co-ordinates of neighbours of (x,y)
  x[0] = x1;     y[0] = y1 + 1;
  x[1] = x1 + 1; y[1] = y1;
  x[2] = x1;     y[2] = y1 - 1;
  x[3] = x1 - 1; y[3] = y1;
  x[4] = x1 + 1; y[4] = y1 + 1;
  x[5] = x1 + 1; y[5] = y1 - 1;
  x[6] = x1 - 1; y[6] = y1 - 1;
  x[7] = x1 - 1; y[7] = y1 + 1;

  // implementing PBC
  for(i=0;i<8;i++){
    if(x[i]==m){x[i]=0;}
    else if(x[i]==-1){x[i]=m-1;}
    if(y[i]==n){y[i]=0;}
    else if(y[i]==-1){y[i]=n-1;}
  }
}

void print_config(int** upper, int** lower, FILE* fp)
{
  int i, j, size;
  size = (2*n-14)/2;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      fprintf(fp,"%d",upper[i][j]);
    }
    fprintf(fp," ");
    for(j=0;j<n;j++){
      fprintf(fp,"%d",lower[i][j]);
    }
    fprintf(fp,"\n");
  }
}

void set_global_params(int tem, int Spl, int Sin, int Vi00, int Vi11, int Vi10, int Vp00, int Vp11, int Vp10)
{
  kT = Temp[tem];

  Vinter[0][0] = Viaa[Vi00];
  Vinter[0][1] = Viab[Vi10];
  Vinter[1][0] = Viab[Vi10];
  Vinter[1][1] = Vibb[Vi11];

  Vplane[0][0] = Vpaa[Vp00];
  Vplane[1][1] = Vpbb[Vp11];
  Vplane[0][1] = Vpab[Vp10];
  Vplane[1][0] = Vpab[Vp10];

  VSplane = VSp[Spl];
  VSinter = VSi[Sin];
}
