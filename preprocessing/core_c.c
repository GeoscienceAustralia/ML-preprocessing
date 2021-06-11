#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int dx[8] = {1, 1, 1, 0, -1, -1, -1, 0 };
int dy[8] = { -1, 0, 1, 1, 1, 0, -1, -1 };

void process_cell(float *data,
                 float *result,
                 int *mask,
                 float nodataval,
                 float tolerance,
                 int m,
                 int n,
                 int target,
                 int current,
                 int mask_value)
{
    int i, ti, tj, ci, cj, ni, nj;

    ti = target / n;
    tj = target % n;
    ci = current / n;
    cj = current % n;

    mask[ci*n+cj] = mask_value;
    result[ti*n+tj] += 1;

    int nis[8] = {0,0,0,0,0,0,0,0};
    int njs[8] = {0,0,0,0,0,0,0,0};

    for (i=0; i<8; i++)
    {
        nis[i] = ci + dx[i];
        njs[i] = cj + dy[i];
    }

    for (i=0; i<8; i++)
    {
        ni = nis[i];
        nj = njs[i];

        if(ni<0 || ni>m-1 || nj<0 || nj>n-1)
        {
            continue;
        }
        else
        {
            if(data[ni*n+nj] == nodataval)
            {
                mask[ni*n+nj] = mask_value;
            }
            else if(data[ni*n+nj] + tolerance >= data[ci*n+cj])
            {
                if(mask[ni*n+nj] != mask_value)
                {
                    process_cell(data, result, mask, nodataval, tolerance,
                                 m, n, target, ni*n+nj, mask_value);
                }
            }
        }
    }
}

void process_cell_non_recursive( float *data,
                                 float *result,
                                 int *mask,
                                 int *queue,
                                 float nodataval,
                                 float tolerance,
                                 int m,
                                 int n,
                                 int target,
                                 int mask_value )
{
    int i, ci, cj, ni, nj, current, neighbour;
    int qidx = 0;
    
    queue[qidx++] = target;
    do
    {
        current = queue[--qidx];
        ci = current / n;
        cj = current % n;

        if(mask[current] != mask_value)
        {
            mask[current] = mask_value;
            result[target] += 1;
        }

        for (i=0; i<8; i++)
        {
            ni = ci + dx[i];
            nj = cj + dy[i];
            
            neighbour = ni*n+nj;

            if(ni<0 || ni>m-1 || nj<0 || nj>n-1)
            {
                continue;
            }
            else
            {
                if(data[neighbour] == nodataval)
                {
                    mask[neighbour] = mask_value;
                }
                else if(data[neighbour] + tolerance >= data[current])
                {
                    if(mask[neighbour] != mask_value)
                    {
                        queue[qidx++] = neighbour;
                    }
                }
            }
        }
    }while(qidx);
}

void compute_upness(float *data, float *result,
                      float nodataval,
                      float tolerance,
                      int m, int n, int nproc)
{
    int NTHREADS = nproc;
    size_t array_size_bytes = (size_t)m * (size_t)n * sizeof(float);
    int **mask = (int**)malloc(sizeof(int*)*NTHREADS);
    int **queue = (int**)malloc(sizeof(int*)*NTHREADS);
    
    memset(result, 0, sizeof(float)*m*n);
    
    for(int i=0; i<NTHREADS; i++)
    {
        printf("%ld\n", array_size_bytes);
        mask[i] = (int*)malloc((size_t)m * (size_t)n * sizeof(int));
        queue[i] = (int*)malloc((size_t)m * (size_t)n * sizeof(int));

        memset(mask[i], 0, (size_t)m * (size_t)n * sizeof(int));
        memset(queue[i], 0, (size_t)m * (size_t)n * sizeof(int));
    }

    int brk=0;
    while(brk){
        brk = brk + 1 - 1;
    };

    omp_set_num_threads(NTHREADS);

#pragma omp parallel for schedule(dynamic, 5)
    for (size_t i=0; i<m; i++)
    {
        int threadid = omp_get_thread_num();
        int *local_mask = mask[threadid];
        int *local_queue = queue[threadid];
        int mask_value = 0;

        printf("%d : %d\n", threadid, i);

        for (size_t j=0; j<n; j++)
        {
            size_t idx = i*n+j;

            if(data[idx] == nodataval)
            {
                continue;
            }
            else
            {
                result[idx] = 0;
                mask_value = idx+1;
                
                // original algorithm
                //mask_value = 1;
                //memset(local_mask, 0, array_size_bytes);
                
#if 0
                process_cell(data, result, local_mask,
                             nodataval, tolerance,
                             m, n, idx, idx, mask_value);
#else
                process_cell_non_recursive(data, result, local_mask, local_queue,
                             nodataval, tolerance,
                             m, n, idx, mask_value);
#endif                
            }
        }
    }

    for(int i=0; i<NTHREADS; i++)
    {
        free(mask[i]);
        free(queue[i]);
    }

    free(mask);
    free(queue);

    printf("Done..\n");
}

