#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "pcg_basic.h"

int flag;      // for debugging
clock_t begin; // start time of the main process, used if flag = 1 for debugging
clock_t end;
double time_spent;

int n;   //dim
int nc = 1;  //crossover
pcg32_random_t rng;

float rand_float()
{
    // return (float)rand() / RAND_MAX;    //C rand seems somewhat bad for the purpose
    // use external rng instead, but remember to seed it!
    return ldexp(pcg32_random_r(&rng), -32);
}

int rand_bool(float p){ //p is success probability
    if (rand_float()>p){
        return 0;
    }
    else{
        return 1;
    }
}

void matinit(int mode, int **M, int d, float p){
    switch (mode)
    {
    case 1:
        for(int i=0;i<d;i++){
            M[i] = (int *)malloc(d * sizeof(int));
            for(int j=0;j<d;j++){
                M[i][j] = rand_bool(p);
            }
        }
        break;

    case 2:
        for(int i=0;i<d;i++){
            M[i] = (int *)malloc(d * sizeof(int));
            for(int j=0;j<d;j++){
                M[i][j] = rand_bool(p) * (j % 3 - 1);
            }
        }
        break;

    default:
        for(int i=0;i<d;i++){
            M[i] = (int *)malloc(d * sizeof(int));
            for(int j=0;j<d;j++){
                M[i][j] = 0;
            }
        }
        break;
    }
}

// A, B is dxd square matrices,
// multiply AB to get C which is dxd,
// C is passed by reference
void matmult(int **A, int **B, int **C, int d){
    // traditional
    // for(int i = 0; i < d; i++){
    //     for(int j=0;j<d;j++){
    //         for(int k=0;k<d;k++){
    //             C[i][j] += A[i][k] * B[k][j]; // definition of matrix mult
    //         }
    //     }
    // }

    // optimized - changed looping order
    for (int i = 0; i < d; i++) {
        for (int k = 0; k < d; k++) {
            for (int j = 0; j < d; j++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void matprint(int **M, int d){
    for(int i = 0; i < d; i++){
        for(int j=0;j<d;j++){
            printf("%d ", M[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void matdiagprint(int **M, int d){
    for(int i = 0; i < d; i++){
        printf("%d\n", M[i][i]);
    }
    printf("\n");
}

void mathalve(int **M, int **A, int **B, int **C, int **D, int d){
    for(int i=0;i<d;i++){
        for(int j=0;j<d;j++){
            A[i][j] = M[i][j];
            B[i][j] = M[i][j+d];
            C[i][j] = M[i+d][j];
            D[i][j] = M[i+d][j+d];
        }
    }
}

// the reverse of mathalve, combine ABCD into M 
void matcombine(int **M, int **A, int **B, int **C, int **D, int d, int odd){
    if(odd==0){
        for(int i=0;i<d;i++){
        for(int j=0;j<d;j++){
            M[i][j] = A[i][j];
            M[i][j+d] = B[i][j];
            M[i+d][j] = C[i][j];
            M[i+d][j+d] = D[i][j];
        }
        }
    }
    else{
        for(int i=0;i<d;i++){
            for(int j=0;j<d;j++){
                M[i][j] = A[i][j];
            }
        }
        for(int i=0;i<d-1;i++){
            for(int j=0;j<d;j++){
                M[i+d][j] = C[i][j];
            }
        }
        for(int i=0;i<d;i++){
            for(int j=0;j<d-1;j++){
                M[i][j+d] = B[i][j];
            }
        }
        for(int i=0;i<d-1;i++){
            for(int j=0;j<d-1;j++){
                M[i+d][j+d] = D[i][j];
            }
        }
    }
    
}

void matadd(int **A, int **B, int **C, int d){
    for(int i=0;i<d;i++){
        for(int j=0;j<d;j++){
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

void matsub(int **A, int **B, int **C, int d){
    for(int i=0;i<d;i++){
        for(int j=0;j<d;j++){
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

int matdiagsum(int **M, int n){
    int total = 0;
    for(int i=0;i<n;i++){
        total += M[i][i];
    }
    return total;
}

void matcopy(int **A, int **B, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            B[i][j] = A[i][j];
        }
    }
}

void straussen_mult(int **A, int **B, int **C, int d){
    int augmented = 0;
    if(d%2==1){
        augmented = 1;
        d += 1;
        A = (int**)realloc(A, d* sizeof(int *));
        B = (int**)realloc(B, d* sizeof(int *));
        C = (int**)realloc(C, d* sizeof(int *));
        for(int row=0; row<d-1; row++){
            A[row] = (int*)realloc(A[row], d* sizeof(int));
            B[row] = (int*)realloc(B[row], d* sizeof(int));
            C[row] = (int*)realloc(C[row], d* sizeof(int));
            A[row][d-1] = 0;
            B[row][d-1] = 0;
            C[row][d-1] = 0;
        }
        A[d-1] = (int*)malloc(d* sizeof(int));
        B[d-1] = (int*)malloc(d* sizeof(int));
        C[d-1] = (int*)malloc(d* sizeof(int));
        for(int col=0; col<d; col++){
            A[d-1][col] = 0;
            B[d-1][col] = 0;
            C[d-1][col] = 0;
        }
    }
    //do normal mult if below nc (crossover point)
    if(nc>=d){
        matmult(A,B,C,d);
    }
    else{ // above crossover, split
        d/=2;
        int **A2 = (int **)malloc(d * sizeof(int *));
        int **B2 = (int **)malloc(d * sizeof(int *));
        int **C2 = (int **)malloc(d * sizeof(int *));
        int **D = (int **)malloc(d * sizeof(int *));
        int **E = (int **)malloc(d * sizeof(int *));
        int **F = (int **)malloc(d * sizeof(int *));
        int **G = (int **)malloc(d * sizeof(int *));
        int **H = (int **)malloc(d * sizeof(int *));
        int **P1 = (int **)malloc(d * sizeof(int *));
        int **P2 = (int **)malloc(d * sizeof(int *));
        int **P3 = (int **)malloc(d * sizeof(int *));
        int **P4 = (int **)malloc(d * sizeof(int *));
        int **P5 = (int **)malloc(d * sizeof(int *));
        int **P6 = (int **)malloc(d * sizeof(int *));
        int **P7 = (int **)malloc(d * sizeof(int *));
        int **I = (int **)malloc(d * sizeof(int *));
        int **J = (int **)malloc(d * sizeof(int *));
        int **K = (int **)malloc(d * sizeof(int *));
        int **L = (int **)malloc(d * sizeof(int *));
        int **X = (int **)malloc(d * sizeof(int *));
        int **Y = (int **)malloc(d * sizeof(int *));
        for(int i=0; i<d; i++){
            A2[i] = (int *)calloc(d, sizeof(int));
            B2[i] = (int *)calloc(d, sizeof(int));
            C2[i] = (int *)calloc(d, sizeof(int));
            D[i] = (int *)calloc(d, sizeof(int));
            E[i] = (int *)calloc(d, sizeof(int));
            F[i] = (int *)calloc(d, sizeof(int));
            G[i] = (int *)calloc(d, sizeof(int));
            H[i] = (int *)calloc(d, sizeof(int));
            P1[i] = (int *)calloc(d, sizeof(int));
            P2[i] = (int *)calloc(d, sizeof(int));
            P3[i] = (int *)calloc(d, sizeof(int));
            P4[i] = (int *)calloc(d, sizeof(int));
            P5[i] = (int *)calloc(d, sizeof(int));
            P6[i] = (int *)calloc(d, sizeof(int));
            P7[i] = (int *)calloc(d, sizeof(int));
            I[i] = (int *)calloc(d, sizeof(int));
            J[i] = (int *)calloc(d, sizeof(int));
            K[i] = (int *)calloc(d, sizeof(int));
            L[i] = (int *)calloc(d, sizeof(int));
            X[i] = (int *)calloc(d, sizeof(int));
            Y[i] = (int *)calloc(d, sizeof(int));
        }
        //split
        mathalve(A,A2,B2,C2,D,d);
        mathalve(B,E,F,G,H,d);

        //compute
        matsub(F,H,X,d);
        straussen_mult(A2,X,P1,d);
        matadd(A2,B2,X,d);
        straussen_mult(X,H,P2,d);
        matadd(C2,D,X,d);
        straussen_mult(X,E,P3,d);
        matsub(G,E,X,d);
        straussen_mult(D,X,P4,d);
        matadd(A2,D,X,d);
        matadd(E,H,Y,d);
        straussen_mult(X,Y,P5,d);
        matsub(B2,D,X,d);
        matadd(G,H,Y,d);
        straussen_mult(X,Y,P6,d);
        matsub(A2,C2,X,d);
        matadd(E,F,Y,d);
        straussen_mult(X,Y,P7,d);

        //finally
        matadd(P5,P4,I,d);
        matsub(I,P2,I,d);
        matadd(I,P6,I,d);

        matadd(P1,P2,J,d);

        matadd(P3,P4,K,d);

        matadd(P5,P1,L,d);
        matsub(L,P3,L,d);
        matsub(L,P7,L,d);

        //combine into C!!
        matcombine(C, I, J, K, L, d, augmented);

        //free memory
        free(A2);
        free(B2);
        free(C2);
        free(D);
        free(E);
        free(F);
        free(G);
        free(H);
        free(P1);
        free(P2);
        free(P3);
        free(P4);
        free(P5);
        free(P6);
        free(P7);
        free(I);
        free(J);
        free(K);
        free(L);
    }
}

void read_file(char filename[], int **A, int **B) {
    FILE* file = fopen(filename, "r");

    int num = 0;
    for(int i=0;i<n;i++){
        A[i] = (int *)malloc(n * sizeof(int));
        for(int j=0;j<n;j++){
            int ret = fscanf(file, "%d", &num);
            A[i][j] = num;
        }
    }
    
    for(int i=0;i<n;i++){
        B[i] = (int *)malloc(n * sizeof(int));
        for(int j=0;j<n;j++){
            int ret = fscanf(file, "%d", &num);
            B[i][j] = num;
        }
    }
    
    fclose(file);
}

int main(int argc, char *argv[])
{
    flag = atoi(argv[1]);
    n = atoi(argv[2]);

    if(flag == 0){
        // more random seed based on external entropy - the time and some program addresses
        pcg32_srandom_r(&rng, time(NULL) ^ (intptr_t)&printf,
                        (intptr_t)&scanf);
    } else{
        pcg32_srandom_r(&rng, 50u, 124u); // fixed seed 
    }
    
    int **A = (int **)malloc(n * sizeof(int *));
    int **B = (int **)malloc(n * sizeof(int *));
    int **C = (int **)malloc(n * sizeof(int *));

    switch(flag)
    {
        case 0:
        {
            read_file(argv[3], A, B);
            matinit(0, C, n, 0);
            nc=80;
            straussen_mult(A,B,C,n);
            matdiagprint(C, n);
            break;
        }
        case 1:
        {
            float param1 = atof(argv[3]);
            int step = (int)param1;
            matinit(2, A, n, 0.5);
            matinit(2, B, n, 0.5);   
            printf("Now scanning the running time of Strassen's algo with n = %d and step size = %d.\n", n, step);
            double min_time = 1e20;
            int min_nc = 0;
            for(nc=2; nc<=n/2+1; nc=nc+step){
                matinit(0, C, n, 0);
                begin = clock();
                straussen_mult(A,B,C,n);
                end = clock();
                time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
                if (time_spent < min_time) {
                    min_time = time_spent;
                    min_nc = nc;
                }
                printf("Total time at nc=%d: %lfs.\n", nc, time_spent);
            }
            printf("Min time at nc=%d: %lfs.\n", min_nc, min_time);
            break;
        }
        case 2:
        {
            matinit(1, A, n, 0.5);
            matinit(1, B, n, 0.5);   
            printf("Now calculating matrix product using conventional multiplication.\n");
            matinit(0, C, n, 0);
            begin = clock();
            matmult(A,B,C,n);
            end = clock();
            time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
            printf("Total time at nc=%d: %lfs.\n", nc, time_spent);
            break;
        }
        case 3:
        {
            nc = 80; // crossover point
            float p = atof(argv[3]);
            int **D = (int **)malloc(n * sizeof(int *));
            matinit(1, A, n, p);
            matinit(0, B, n, 0);
            matcopy(A, B, n);
            matinit(0, C, n, 0);
            matinit(0, D, n, 0);
            straussen_mult(A,B,C,n);
            straussen_mult(A,C,D,n);
            printf("Number of triangles in the random graph with p=%3f: %d.\n", p, matdiagsum(D,n)/6);
            break;
        }
        default:
            printf("Meow! You have provided an invalid flag.\n");
            break;
    }
    
    free(A);
    free(B);
    free(C);
}


