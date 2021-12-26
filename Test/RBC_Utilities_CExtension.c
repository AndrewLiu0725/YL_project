#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(x, y) (((x) > (y)) ? (x) : (y))
#define DEBUG 0

// use "gcc -fPIC -shared -o RBC_Utilities_CExtension.so RBC_Utilities_CExtension.c" to compile and create .so file

void correctDiffpos(double *distance, const double *COMs, int number_of_pairs, int timesteps, double Dm, const int *indice_pairs, const int* dim){
    int k, t;
    int i, j;
    double dx, dy, dz;

    for (k = 0; k < number_of_pairs; k++){
        i = indice_pairs[2*k];
        j = indice_pairs[2*k + 1];
        for (t = 0; t < timesteps; t++){ 
            dx = fabs(COMs[i*timesteps*3 + t*3] - (COMs[j*timesteps*3 + t*3]));
            dy = fabs(COMs[i*timesteps*3 + t*3 + 1] - (COMs[j*timesteps*3 + t*3 + 1]));
            dz = fabs(COMs[i*timesteps*3 + t*3 + 2] - (COMs[j*timesteps*3 + t*3 + 2]));

            // Deal with periodic boundary condition
            if (dx > dim[0]/2) dx = dim[0] - dx;
            if (dz > dim[2]/2) dz = dim[2] - dz;
            
            distance[k*timesteps + (t)] = sqrt(dx*dx + dy*dy + dz*dz); 
        }
    }
}

// this extremum recognition function may be improved
int extremum(const double *array, int t, const int k, const int timesteps){
    /*
    * Function:  extremum
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

void calcDF(int *doublet_or_not, const double *distance, const double *uncorrected_distance, const double *COMs, const double *COMs_NB,
const int *indice_pairs, const int* dim, int period, int timesteps, int number_of_pairs, double Dm, double criteria_Dm, int *end_time){
    int k, i, j; // indices related to particles
    int t, inner_t, prev_min_distance_t;
    int current_end_time = timesteps;
    double mean_d, sum_d; // inter-particle distance

    double prev_min_distance, cur_min_distance, max_distance;
    double max_max_distance = Dm*4/3; // allow relative motion
    double max_min_distance = Dm*3/4; // whereas need to be close enough
    
    // exclude close uncoupled partlces (mainly happens in the two-cell system)
    double dz, sum_dz, mean_dz; 
    double max_dz = Dm/3; 

    // exclude those transition states (two particles nearly touch each other but do not form or break)
    int caution; double caution_percentage;
    double caution_distance = Dm; double max_caution_percentage = 0.1;

    //double max_min_distance_difference = Dm*0.5;

    // exclude the condition that two particles keep moving along the streamline inside a small box
    // for the two-cell system, when the volume fraction is high, this case may satisfy all conditions above
    double prev_relative_x_NB, current_relative_x_NB, relative_x_motion;
    double max_relative_x_motion = Dm;

    // run over each pair of particles
    for (k = 0; k < number_of_pairs; k++){

        i = indice_pairs[2*k]; j = indice_pairs[2*k+1]; // two particles' indices
        prev_min_distance_t = 0;

        // run over t
        for (t = 3; t < (timesteps - 3); t ++){
            // min distance
            if (extremum(distance, t, k, timesteps) == 2){
                // calc the mean inter-particle distance (d) and caution percentage
                sum_d = 0; caution = 0;
                for (inner_t = prev_min_distance_t; inner_t <= t; inner_t++){
                    sum_d += distance[k*timesteps + inner_t];
                    if (distance[k*timesteps + inner_t] > caution_distance){
                        caution += 1;
                    }
                }
                mean_d = sum_d/(t-prev_min_distance_t+1);
                caution_percentage = (double)caution/(t-prev_min_distance_t+1);

                // record the min distances for this time window
                prev_min_distance = distance[k*timesteps + prev_min_distance_t];
                cur_min_distance = distance[k*timesteps + t];

                // calc the mean inter-particle distance in z-axis (dz)
                sum_dz = 0;
                for (inner_t = prev_min_distance_t; inner_t <= t; inner_t++){
                    dz = fabs(COMs[i*timesteps*3 + inner_t*3 + 2]-COMs[j*timesteps*3 + inner_t*3 + 2]);
                    if (dz > dim[2]) dz = (dim[2]-dz); // cross the boundary
                    sum_dz += dz;
                }
                mean_dz = sum_dz/(t-prev_min_distance_t+1);

                
                prev_relative_x_NB = fabs(COMs_NB[i*timesteps*3 + prev_min_distance_t*3] - COMs_NB[j*timesteps*3 + prev_min_distance_t*3]);
                current_relative_x_NB = fabs(COMs_NB[i*timesteps*3 + t*3] - COMs_NB[j*timesteps*3 + t*3]);
                relative_x_motion = fabs(prev_relative_x_NB-current_relative_x_NB);
                // printf("[%d,%d]: |%lf-%lf|=%lf\n", prev_min_distance_t, t, prev_relative_x_NB, current_relative_x_NB, fabs(prev_relative_x_NB-current_relative_x_NB));
                
                if (DEBUG){
                    if ((t > 350) && (t < 450)){
                        printf("[%d,%d]:\n", prev_min_distance_t, t);
                        printf("caution_percentage = %lf\n", caution_percentage);
                        printf("relative_x = %lf\n", relative_x_motion);
                    }
                }

                // determine if this time window is in doublet state
                if ((mean_d < criteria_Dm*Dm) && (mean_dz < max_dz) && (relative_x_motion < max_relative_x_motion) &&
                (caution_percentage < max_caution_percentage) &&
                (max_distance < max_max_distance) && (max(prev_min_distance, cur_min_distance) < max_min_distance)){
                    for (inner_t = prev_min_distance_t; inner_t <= t; inner_t++){
                        doublet_or_not[k*(timesteps - period) + inner_t] = 1;
                    }
                }
                prev_min_distance_t = t; // update the prev_min_distance_t 
            }

            // max distance (peak)
            else if (extremum(distance, t, k, timesteps) == 1){
                max_distance = distance[k*timesteps + t];
            }
        }
        current_end_time = min(current_end_time, prev_min_distance_t);
    }
    end_time[0] = current_end_time;
}