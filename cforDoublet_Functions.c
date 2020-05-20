#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// use "gcc -fPIC -shared -o cforDoublet_Functions.so cforDoublet_Functions.c" to compile and create .so file
// use "gcc -shared -o cforDoublet_Functions.so -fPIC -fopenmp cforDoublet_Functions.c" to compile and create .so file

void correctDiffpos(double *diffpos, const double *COMs, int number_of_pairs, int timesteps, double Dm, const int *indice_pairs, const int* dim){
    int k, t;
    int i, j;
    double correct_diffpos_square;
    double dx, dy, dz;
    // #pragma omp parallel shared(diffpos,COMs,number_of_pairs,timesteps,Dm,indice_pairs,dim) private(i,j,k,t,correct_diffpos_square,dx,dy,dz) 
    // #pragma omp for  schedule(static)
    for (k = 0; k < number_of_pairs; k ++){
        for (t = 1; t < timesteps; t ++){
            i = indice_pairs[2*k];
            j = indice_pairs[2*k + 1];

            // Deal with periodic boundary condition
            dx = fabs(COMs[i*timesteps*3 + t*3] - (COMs[j*timesteps*3 + t*3]));
            dy = fabs(COMs[i*timesteps*3 + t*3 + 1] - (COMs[j*timesteps*3 + t*3 + 1]));
            dz = fabs(COMs[i*timesteps*3 + t*3 + 2] - (COMs[j*timesteps*3 + t*3 + 2]));
            if ((dx > dim[0]/2) || (dz > dim[2]/2)){
                correct_diffpos_square = dy*dy;
                if (dx > dim[0]/2) dx = dim[0] - dx;
                if (dz > dim[2]/2) dz = dim[2] - dz;
                correct_diffpos_square = dx*dx + dy*dy + dz*dz;
                diffpos[k*timesteps + (t)] = sqrt(correct_diffpos_square);
            }    
        }
    }
}

void calcDFGeneral(int *doublet_or_not, const double *diffpos, int period, int timesteps, int number_of_pairs, double Dm, double criteria_Dm){
    int i, j, t;
    double mean;

    // #pragma omp parallel shared(diffpos,doublet_or_not,number_of_pairs,timesteps,Dm,criteria_Dm,period) private(i,t,mean) 
    // #pragma omp for  schedule(static)
    for (i = 0; i < number_of_pairs; i++){
        // Initialiaze the mean here
        mean = 0;
        for (j = 0; j < period; j++){
            mean += diffpos[i*timesteps + j] / period;
        }
        if (mean < criteria_Dm*Dm){
            doublet_or_not[i*(timesteps - period)] = 1;
        }

        // Run over t
        for (t = 1; t < (timesteps - period); t++){
            mean = (mean - diffpos[i*timesteps + (t-1)]/period + diffpos[i*timesteps + (t + period - 1)]/period );
            if (mean < criteria_Dm*Dm){
                doublet_or_not[i*(timesteps - period) + t] = 1;
            }
        }
    }
}

void calcDFStrict(int *doublet_or_not, const double *diffpos, int period, int timesteps, int number_of_pairs, double Dm, double criteria_Dm){
    int i, j, t, found;
    // #pragma omp parallel shared(diffpos,doublet_or_not,number_of_pairs,timesteps,Dm,criteria_Dm,period) private(i,j,t,found) 
    // #pragma omp for  schedule(static)
    for (i = 0; i < number_of_pairs; i ++){
        // Initialize here
        t = 0;
        found = 0;
        for (j = t + period - 1; j >= t; j --){
            // Exceed the criteria for r
            if (diffpos[i*timesteps + j] > criteria_Dm*Dm){
                t = j+1;
                found = 1;
                break;
            }
        }
        if (!found){
            doublet_or_not[i*(timesteps-period) + t] = 1;
            t += 1;
        }

        // Run over t
        while (t < (timesteps - period)){
            // this pair forms doublet at time t-1
            if (doublet_or_not[i*(timesteps-period)+t-1] == 1){
                // this doublet survives
                if (diffpos[i*timesteps + t + period - 1] < criteria_Dm*Dm){
                    doublet_or_not[i*(timesteps-period) + t] = 1;
                    t += 1;
                }
                // this doublet breaks
                else{
                    t += period;
                }
            }
            
            // this pair doesn't form doublet at time t-1, i.e. diffpos[i][t-1] > criteria_Dm*Dm
            else{
                found = 0;
                for (j = t + period - 1; j >= t; j --){
                    // Exceed the criteria for r
                    if (diffpos[i*timesteps + j] > criteria_Dm*Dm){
                        t = j+1;
                        found = 1;
                        break;
                    }
                }

                // Form a doublet
                if (!found){
                    doublet_or_not[i*(timesteps-period) + t] = 1;
                    t += 1;
                }
            }
        }
    }
}