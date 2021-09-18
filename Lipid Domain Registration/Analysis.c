#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// System Variables

  //Dimensions of each stacked sheet
  int m = 100;
  int n = 100;

  // Full set of parameters

  double epsilon = 0.3;                                                         // 0.3 kcal/mol : Temperature scaling constant
  double kT = 0.0;                                                              // Variable holding current scaled temperature
  double Vplane[2][2] = {                                                       // Array holding current values of Vplane_(alpha,alpha')
    {0.0, 0.0},
    {0.0, 0.0}
  };
  double Vinter[2][2] = {                                                       // Array holding current values of Vinter_(alpha,alpha')
    {0.0, 0.0},
    {0.0, 0.0}
  };
  double VSplane = 0.0;                                                         // Variable holding current in-plane entropic strength constant
  double VSinter = 0.0;                                                         // Variable holding current interleaflet entropic strength constant

  int p0 = 5000;
  int p1 = 5000;

  int num_configs = 2000;                                                       // number of configurations printed for each run
  int PScutoff0[2] = {1000, 4000};
  int PScutoff1[2] = {1000, 4000};
  double regcutoff[4] = {0.0333, 0.10, 0.2333, 0.30};

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Function Declarations

/* unif. dist. from 0.0 to 1.0 */
double rand_double();

/* discrete integer unif. dist. from a to b */
int uniform_distribution(int rangeLow, int rangeHigh);

/* reads the parameters from the first line of a config file */
void read_params(FILE* fp);

/* counts # of islands(domains) in each leaflet, stores population distribution by size */
void countIslands(int matrix[][n], int* pop0, int* pop1);

/* the depth-first search algorithm used to recursively find size of domains */
void DFS(int matrix[][n], int i, int j, int visit[][n], int* size, int species);

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Main

void main(){

  // Variable Declarations
  int i, j, k, t, type, x1, y1, x2, y2, r, quot, ucat0, ucat2, lcat0, lcat2, reg, preg, ureg, pareg, areg;
  double uscore0, lscore0, uscore1, lscore1, uscoredr, lscoredr, scoredr, avscoredr;
  char buff[255], skipline[255], gridskip[255], params[255];

  // Output File Initializations
  FILE *PSstate;
  PSstate = fopen("PSstate.txt", "w+");     // For phase-separation classification output
  FILE *regstate;
  regstate = fopen("Regstate.txt", "w+");   // For domain registration classification output

  // Outer Loops to change parameters between runs
  int tem, Spl, Sin, Vi00, Vi11, Vi10, Vp00, Vp11, Vp10;
  for(tem = 0; tem<1; tem++){                                                   // increments through the entries in the temperature array
    for(Spl = 0; Spl<1; Spl++){                                                 // increments through the entries in the VSp array
      for(Sin = 0; Sin<4; Sin++){                                               // increments through the entries in the VSi array
        for(Vi11 = 0; Vi11<3; Vi11++){                                          // increments through the entries in the Viaa array
          for(Vi00 = 0; Vi00<3; Vi00++){                                        // increments through the entries in the Vibb array
            for(Vi10 = 0; Vi10<3; Vi10++){                                      // increments through the entries in the Viab array
              for(Vp00 = 0; Vp00<5; Vp00++){                                    // increments through the entries in the Vpaa array
                for(Vp11 = 0; Vp11<=Vp00; Vp11++){                              // increments through the entries in the Vpbb array such that Vpbb <= Vpaa
                  for(Vp10 = 0; Vp10<=Vp11; Vp10++){                            // increments through the entries in the Vpab array such that Vpab <= Vpbb

                    // resets the random number seed (otherwise the same random numbers are generated in each run)
                    srand((long)time(NULL));

                    // Config Input File input initialization
                    char buffer[32];                                            // The filename buffer.
                    snprintf(buffer,sizeof(char)*32,"%d-%d-%d-%d-%d-%d-%d-%d-%d.txt",tem,Spl,Sin,Vi00,Vi11,Vi10,Vp00,Vp11,Vp10);
                    FILE *fp;
                    fp = fopen(buffer, "r");

                    //Reading the parameters from the first line of the file
                    read_params(fp);
                    fgets(params, 255, (FILE*)fp);

                    // Memory Allocation for lattice layers
                    int upper[m][n];                                            // 2D array for upper leaflet
                    int lower[m][n];                                            // 2D array for lower leaflet
                    int upop0[p0], upop1[p1], lpop0[p0], lpop1[p1];             // arrays for storing population distribution of both species in each leaflet

                    double uscore[3] = {0.0, 0.0, 0.0}; double lscore[3] = {0.0, 0.0, 0.0}; double score[3] = {0.0, 0.0, 0.0};
                    reg = 0; preg = 0; ureg = 0; pareg = 0; areg = 0; avscoredr = 0.0;

                    for(t=0;t<num_configs;t++){
                      // skip line where move number and energy are printed before each config
                      fgets(skipline, 255, (FILE*)fp);
                      // skip first 1800 configs
                      if(t<1800){
                        for(i=0;i<m;i++){
                          fgets(gridskip, 255, (FILE*)fp);
                        }
                      }
                      else{
                        // Read config for both upper and lower leaflet to respective 2D arrays
                        for(i=0;i<m;i++){
                          fgets(buff, 255, (FILE*)fp);
                          for(j=0;j<=2*n;j++){
                            if(j<n){
                              upper[i][j] = buff[j] - '0';
                            }
                            else if(j>n){
                              lower[i][j-n-1] = buff[j] - '0';
                            }
                          }
                        }

                        // Phase Separation Characterisation

                        // Counting size of domains in each leaflets and obtaining size distribution

                        memset(upop0, 0, sizeof(upop0));
                        memset(upop1, 0, sizeof(upop1));
                        memset(lpop0, 0, sizeof(lpop0));
                        memset(lpop1, 0, sizeof(lpop1));
                        countIslands(upper,upop0,upop1);
                        countIslands(lower,lpop0,lpop1);

                        ucat0 = 0; ucat2 = 0; lcat0 = 0; lcat2 = 0;
                        // Check if fully phase separated
                        uscore0 = 0; lscore0 = 0; uscore1 = 0; lscore1 = 0;
                        for(k=1600;k<p0;k++){
                          uscore0 = uscore0 + upop0[k]*(k+1);
                          lscore0 = lscore0 + lpop0[k]*(k+1);
                        }
                        for(k=1600;k<p1;k++){
                          uscore1 = uscore1 + upop1[k]*(k+1);
                          lscore1 = lscore1 + lpop1[k]*(k+1);
                        }
                        if(uscore0>=PScutoff0[1] && uscore1>=PScutoff1[1]){uscore[2] = uscore[2] + 1.0; ucat2 = 1;}
                        if(lscore0>=PScutoff0[1] && lscore1>=PScutoff1[1]){lscore[2] = lscore[2] + 1.0; lcat2 = 1;}

                        // Check if not phase separated, and consists mostly of small local aggregations
                        uscore0 = 0; lscore0 = 0; uscore1 = 0; lscore1 = 0;
                        for(k=0;k<200;k++){
                          uscore0 = uscore0 + upop0[k]*(k+1);
                          lscore0 = lscore0 + lpop0[k]*(k+1);
                        }
                        for(k=0;k<200;k++){
                          uscore1 = uscore1 + upop1[k]*(k+1);
                          lscore1 = lscore1 + lpop1[k]*(k+1);
                        }
                        if(uscore0>=PScutoff0[1] || uscore1>=PScutoff1[1]){uscore[0] = uscore[0] + 1.0; ucat0 = 1;}
                        if(lscore0>=PScutoff0[1] || lscore1>=PScutoff1[1]){lscore[0] = lscore[0] + 1.0; lcat0 = 1;}

                        // If none of the above, it is partially phase separated
                        if(ucat2==0 && ucat0==0){uscore[1] = uscore[1] + 1.0;}
                        if(lcat2==0 && lcat0==0){lscore[1] = lscore[1] + 1.0;}

                        // Domain Registration Characterisation

                        scoredr = 0.0;
                        double norm = (p0/2) + p1;    // Normalization factor for pdfs
                        // KL Divergence calculation for each config
                        for(i=0;i<m;i++){
                          for(j=0;j<n;j++){
                            uscoredr = 0.5/norm; lscoredr = 0.5/norm;
                            if(upper[i][j]==1){
                              uscoredr = 1.0/norm;
                            }
                            if(lower[i][j]==1){
                              lscoredr = 1.0/norm;
                            }
                            scoredr = scoredr + uscoredr*log(uscoredr/lscoredr)/log(2);
                          }
                        }
                        // Domain registration classification
                        if(scoredr<=regcutoff[0]){
                          reg = reg + 1;
                        }
                        else if(scoredr<=regcutoff[1] && scoredr>regcutoff[0]){
                          preg = preg + 1;
                        }
                        else if(scoredr<=regcutoff[2] && scoredr>regcutoff[1]){
                          ureg = ureg + 1;
                        }
                        else if(scoredr<=regcutoff[3] && scoredr>regcutoff[2]){
                          pareg = pareg + 1;
                        }
                        else if(scoredr>regcutoff[3]){
                          areg = areg + 1;
                        }
                        avscoredr = avscoredr + scoredr;
                      }

                    }
                    avscoredr = avscoredr/(num_configs-200.0);

                    // Phase Separation Output

                    if(uscore[2]>=uscore[1] && uscore[2]>uscore[0]){score[2] = score[2] + 0.5;}
                    if(lscore[2]>=lscore[1] && lscore[2]>lscore[0]){score[2] = score[2] + 0.5;}
                    if(uscore[1]>uscore[2] && uscore[1]>uscore[0]){score[1] = score[1] + 0.5;}
                    if(lscore[1]>lscore[2] && lscore[1]>lscore[0]){score[1] = score[1] + 0.5;}
                    if(uscore[0]>=uscore[1] && uscore[0]>=uscore[2]){score[0] = score[0] + 0.5;}
                    if(lscore[0]>=lscore[1] && lscore[0]>=lscore[2]){score[0] = score[0] + 0.5;}

                    // Classifying the run as phase separated, partially phase separated or not phase separated
                    if(score[2]>=score[1] && score[2]>score[0]){
                      fprintf(PSstate,"%lf %lf %lf %lf %lf %lf %lf %lf %lf 6\n",kT,VSplane,VSinter,Vinter[0][0],Vinter[1][1],Vinter[0][1],Vplane[0][0],Vplane[1][1],Vplane[0][1]);
                    }
                    else if(score[1]>=score[2] && score[1]>=score[0]){
                      fprintf(PSstate,"%lf %lf %lf %lf %lf %lf %lf %lf %lf 5\n",kT,VSplane,VSinter,Vinter[0][0],Vinter[1][1],Vinter[0][1],Vplane[0][0],Vplane[1][1],Vplane[0][1]);
                    }
                    else if(score[0]>=score[1] && score[0]>score[2]){
                      fprintf(PSstate,"%lf %lf %lf %lf %lf %lf %lf %lf %lf 4\n",kT,VSplane,VSinter,Vinter[0][0],Vinter[1][1],Vinter[0][1],Vplane[0][0],Vplane[1][1],Vplane[0][1]);
                    }

                    // Domain Registration Output

                    // Classifying run as converging to registered, partially registered or unregistered
                    if(reg>=preg && reg>ureg && reg>pareg && reg>areg){
                      fprintf(regstate,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf 6\n",kT,VSplane,VSinter,Vinter[0][0],Vinter[1][1],Vinter[0][1],Vplane[0][0],Vplane[1][1],Vplane[0][1],avscoredr);
                    }
                    else if(preg>=reg && preg>=ureg && preg>pareg && preg>areg){
                      fprintf(regstate,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf 5\n",kT,VSplane,VSinter,Vinter[0][0],Vinter[1][1],Vinter[0][1],Vplane[0][0],Vplane[1][1],Vplane[0][1],avscoredr);
                    }
                    else if(ureg>=preg && ureg>reg && ureg>=pareg && ureg>areg){
                      fprintf(regstate,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf 4\n",kT,VSplane,VSinter,Vinter[0][0],Vinter[1][1],Vinter[0][1],Vplane[0][0],Vplane[1][1],Vplane[0][1],avscoredr);
                    }
                    else if(pareg>reg && pareg>preg && pareg>=ureg && pareg>=areg){
                      fprintf(regstate,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf 3\n",kT,VSplane,VSinter,Vinter[0][0],Vinter[1][1],Vinter[0][1],Vplane[0][0],Vplane[1][1],Vplane[0][1],avscoredr);
                    }
                    else if(areg>reg && areg>preg && areg>ureg && areg>=pareg){
                      fprintf(regstate,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf 2\n",kT,VSplane,VSinter,Vinter[0][0],Vinter[1][1],Vinter[0][1],Vplane[0][0],Vplane[1][1],Vplane[0][1],avscoredr);
                    }

                    // Clean up
                    fclose(fp);                                                 // closes text file used for config output

                  }
                }
              }
            }
          }
        }
      }
    }
  }
// End of runs with all required parameter combinations
}

double rand_double()
{
   return rand()/(double)RAND_MAX;
}

int uniform_distribution(int rangeLow, int rangeHigh)
{
  int myRand = (int)rand();
  int range = rangeHigh - rangeLow + 1; // +1 makes it [rangeLow, rangeHigh], inclusive.
  int myRand_scaled = (myRand % range) + rangeLow;
  return myRand_scaled;
}

void read_params(FILE* fp)
{
  fseek(fp, 32, SEEK_SET);
  fscanf(fp, "%lf",&kT);
  fseek(fp, 11, SEEK_CUR);
  fscanf(fp, "%lf",&VSplane);
  fseek(fp, 11, SEEK_CUR);
  fscanf(fp, "%lf",&VSinter);
  fseek(fp, 17, SEEK_CUR);
  fscanf(fp, "%lf",&Vinter[0][0]);
  fseek(fp, 17, SEEK_CUR);
  fscanf(fp, "%lf",&Vinter[1][1]);
  fseek(fp, 17, SEEK_CUR);
  fscanf(fp, "%lf",&Vinter[0][1]);
  Vinter[1][0] = Vinter[0][1];
  fseek(fp, 17, SEEK_CUR);
  fscanf(fp, "%lf",&Vplane[0][0]);
  fseek(fp, 17, SEEK_CUR);
  fscanf(fp, "%lf",&Vplane[1][1]);
  fseek(fp, 17, SEEK_CUR);
  fscanf(fp, "%lf",&Vplane[0][1]);
  Vplane[1][0] = Vplane[0][1];
  fseek(fp, 0, SEEK_SET);

}

void countIslands(int matrix[][n], int* pop0, int* pop1)
{
  int count, i, j, s;
  int size[1];
  int visit[m][n];                                                             // Array storing visit state of points on the lattice

  // Carrying out the count for species 1
  for(i=0;i<m;i++){                                                             // initializing all to unvisited (value of 0 means unvisited, 1 means visited)
    for(j=0;j<n;j++){
      visit[i][j] = 0;
    }
  }
  count = 0;                                                                    // Count of different islands
  for (i = 0; i < m; ++i){
    for (j = 0; j < n; ++j){
      size[0] = 0;                                                              // Holds and updates to current size of an island during DFS
      if (matrix[i][j]==visit[i][j]+1)                                          // Check if cell is visited, if not it is a separate domain
      {
        DFS(matrix, i, j, visit, size, 1);                                      // Recursively visits all connected points in the island using the DFS algorithm
        ++count;
        s = size[0] - 1;
        pop1[s] = pop1[s] + 1;                                                  // Adds size of island to the size distribution
      }
    }
  }

  // Carrying out the count for species 0
  for(i=0;i<m;i++){                                                             // initializing all to unvisited (value of 1 means unvisited, 0 means visited)
    for(j=0;j<n;j++){
      visit[i][j] = 1;
    }
  }
  count = 0;                                                                    // Count of different islands
  for (i = 0; i < m; ++i){
    for (j = 0; j < n; ++j){
      size[0] = 0;                                                              // Holds and updates to current size of an island during DFS
      if (matrix[i][j]==visit[i][j]-1)                                          // Check if cell is visited, if not it is a separate domain
      {
        DFS(matrix, i, j, visit, size, 0);                                      // Recursively visits all connected points in the island using the DFS algorithm
        ++count;
        s = size[0] - 1;
        pop0[s] = pop0[s] + 1;                                            // Adds size of island to the size distribution
      }
    }
  }
}

void DFS(int matrix[][n], int i, int j, int visit[][n], int* size, int species)
{
  int x, y, k, temp;
  if(species==0){temp = -1;}
  if(species==1){temp = 1;}

  int rowNbr[] = { -1, 0, 0, 1 };
  int colNbr[] = { 0, -1, 1, 0 };

  // Mark this cell as visited
  visit[i][j] = species;

  size[0] = size[0] + 1;
  // Recur for all connected neighbours
  for (k = 0; k < 4; ++k){
    x = i + rowNbr[k];
    y = j + colNbr[k];
    if(x==m){x=0;}
    if(x==-1){x=m-1;};
    if(y==n){y=0;}
    if(y==-1){y=n-1;}
    if (matrix[x][y]==visit[x][y]+temp){
      DFS(matrix, x, y, visit, size, species);
    }
  }
}
