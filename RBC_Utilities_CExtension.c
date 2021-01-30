#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// use "gcc -fPIC -shared -o RBC_Utilities_CExtension.so RBC_Utilities_CExtension.c" to compile and create .so file
// use "gcc -shared -o RBC_Utilities_CExtension.so -fPIC -fopenmp RBC_Utilities_CExtension.c" to compile and create .so file

void correctDiffpos(double *diffpos, const double *COMs, int number_of_pairs, int timesteps, double Dm, const int *indice_pairs, const int* dim){
    int k, t;
    int i, j;
    double correct_diffpos_square;
    double dx, dy, dz;
    // #pragma omp parallel shared(diffpos,COMs,number_of_pairs,timesteps,Dm,indice_pairs,dim) private(i,j,k,t,correct_diffpos_square,dx,dy,dz) 
    // #pragma omp for  schedule(static)
    for (k = 0; k < number_of_pairs; k ++){
        i = indice_pairs[2*k];
        j = indice_pairs[2*k + 1];
        for (t = 0; t < timesteps; t ++){ 
            dx = fabs(COMs[i*timesteps*3 + t*3] - (COMs[j*timesteps*3 + t*3]));
            dy = fabs(COMs[i*timesteps*3 + t*3 + 1] - (COMs[j*timesteps*3 + t*3 + 1]));
            dz = fabs(COMs[i*timesteps*3 + t*3 + 2] - (COMs[j*timesteps*3 + t*3 + 2]));

            // Deal with periodic boundary condition
            if (dx > dim[0]/2) dx = dim[0] - dx;
            if (dz > dim[2]/2) dz = dim[2] - dz;
            
            correct_diffpos_square = dx*dx + dy*dy + dz*dz;
            diffpos[k*timesteps + (t)] = sqrt(correct_diffpos_square); 
        }
    }
}

int extreme(const double *array, int t, const int k, const int timesteps){
    /*
    * Function:  extreme 
    * --------------------
    *  Distinguish whether current point is min, max, or neither.
    * 
    *  returns: 2 if is min, 1 if is max, 0 otherwise.
    */
    t += (k*timesteps);
    if (((array[t] < array[t-1]) && (array[t] < array[t+1])) &&
    ((array[t] < array[t-3]) && (array[t] < array[t+3]))){
        return 2;
    }
    else if (((array[t] > array[t-1]) && (array[t] > array[t+1])) &&
    ((array[t] > array[t-3]) && (array[t] > array[t+3]))){
        return 1;
    }
    else return 0;
}

void calcDF(int *doublet_or_not, const double *diffpos, int period, int timesteps, int number_of_pairs, double Dm, double criteria_Dm,
const double *COMs, const int *indice_pairs, int *state_series, int *end_time, const int* dim){
    int k, t, inner_t, i, j, inner_st;
    double mean, sum_x;
    double dl, sum_l;
    double prev_min, cur_min, max, range; // for excluding close uncoupled partlces
    double dz, sum_dz, mean_dz; double mean_dz_threshold = 5; // for excluding close uncoupled partlces
    double min_distance_1 = 11.5; double min_distance_2 = 7.5; // criteria for min in each period
    double min_range = 3;
    int prev_min_t, pprev_min_t;
    // The value of the following two criteria may need further consideration
    double velocity_sliding_threshold = 0.5*Dm/period;
    double velocity_kayaking_threshold = 1.6/period;
    int state, prev_state;

    // #pragma omp parallel shared(diffpos,doublet_or_not,number_of_pairs,timesteps,Dm,criteria_Dm,period) private(i,t,mean) 
    // #pragma omp for  schedule(static)
    for (k = 0; k < number_of_pairs; k++){
        i = indice_pairs[2*k]; j = indice_pairs[2*k+1]; // particle id

        // Initilaization
        pprev_min_t = prev_min_t = 0;
        prev_state = state = 4;
        // run over t
        for (t = 3; t < (timesteps - period); t ++){
            // min 
            if (extreme(diffpos, t, k, timesteps) == 2){
                state = 4;
                // is sliding state
                if (((fabs(COMs[i*timesteps*3 + t*3] - COMs[i*timesteps*3 + pprev_min_t*3])/(t-prev_min_t)) > velocity_sliding_threshold)
                && ((fabs(COMs[j*timesteps*3 + t*3] - COMs[j*timesteps*3 + pprev_min_t*3])/(t-prev_min_t)) > velocity_sliding_threshold)){
                    state = 3;
                }
                // is not sliding state
                else{
                    // calc the mean
                    sum_x = 0;
                    for (inner_t = prev_min_t; inner_t <= t; inner_t++){
                        sum_x += diffpos[k*timesteps + inner_t];
                    }
                    mean = sum_x/(t-prev_min_t+1);

                    // calc range
                    prev_min = diffpos[k*timesteps + prev_min_t];
                    cur_min = diffpos[k*timesteps + t];
                    range = ((max - cur_min) + (max - prev_min))/2;

                    // calc path of length
                    sum_l = 0;
                    for (inner_t = prev_min_t; inner_t <= t; inner_t++){
                        dl = fabs(COMs[i*timesteps*3 + inner_t*3]-COMs[i*timesteps*3 + (inner_t-1)*3]);
                        if (dl > dim[0]) dl = (dim[0]-dl); // cross the boundary
                        sum_l += dl;
                    }

                    // calc dz
                    sum_dz = 0;
                    for (inner_t = prev_min_t; inner_t <= t; inner_t++){
                        dz = fabs(COMs[i*timesteps*3 + inner_t*3 + 2]-COMs[j*timesteps*3 + inner_t*3 + 2]);
                        if (dz > dim[2]) dz = (dim[2]-dz); // cross the boundary
                        sum_dz += dz;
                    }
                    mean_dz = sum_dz/(t-prev_min_t+1);

                    // recognize if is doublet state
                    if ((mean < criteria_Dm*Dm) && (max < (Dm+1)) && (mean_dz < mean_dz_threshold)){
                        if ((prev_min < min_distance_2) && (cur_min < min_distance_2)){
                            state = 1;
                        }
                        else{
                            // have relative motion (may use speed to exclude long cycle in kayaking state)
                            if ((prev_min < min_distance_1) && (cur_min < min_distance_1) && (range > min_range)){
                                state = 1;
                            }
                            else{
                                if ((sum_l/(t-prev_min_t+1)) < velocity_kayaking_threshold) state = 2;
                            }
                        }
                    }
                    // recognize if is two-kayaking-singlet state
                    else{
                        if ((sum_l/(t-prev_min_t+1)) < velocity_kayaking_threshold) state = 2;
                    }
                }
                //state series
                if (prev_state == 3) inner_st = pprev_min_t;
                else inner_st = prev_min_t;
                for (inner_t = inner_st; inner_t <= t; inner_t++){
                    state_series[inner_t] = state;
                }
                // doublet fraction
                if (state == 1){
                    for (inner_t = inner_st; inner_t <= t; inner_t++){
                        doublet_or_not[k*(timesteps - period) + inner_t] = 1;
                    }
                }
                pprev_min_t = prev_min_t;
                prev_min_t = t;
                prev_state = state;
            }

            // max
            else if (extreme(diffpos, t, k, timesteps) == 1){
                max = diffpos[k*timesteps + t];
            }
        }
        end_time[0] = prev_min_t;
    }
}