#include "get_divs.h"
// test the divgence caculation
int main(){
    double  div;
    int i, region_length=3;
    // this is a special region
    /*
    double S[9][3];
    for(i = 0; i<9; i++){
        C[i][0] = 0.4;
        C[i][1] = 0;
        C[i][2] = 1;
    }
    */

    // sample of non-turbulence
    double A[9][3] = {{0.4, 0,2},{0.4,0,1},{0.4,0,2},{0.4,0,1},{0.4,0,0},{0.4,0,1},{0.4,0,2},{0.4,0,1},{0.4,0,2}};

    // a sample turbulence
    double B[9][3] = {{0.2,0.2, 2},{0.4, 0 ,1}, {0.2, -0.2, 2},{0, 0.4,1}, {0 ,0 ,0 }, {0, -0.4, 1 },{-0.2 , 0.2, 2}, {-0.4, 0, 1}, {-0.2, -0.2, 2}};

    // also a turbulence
    double C[9][3];
    for(i = 0; i < 9 ; i++){
        C[i][0] = -B[i][0];
        C[i][1] = -B[i][1];
        C[i][2] = B[i][2];
    }

    // another pattern
    double D[9][3] = {{0.2,0.2, 2},{0.4, 0 ,1}, {0.2, -0.2, 2},{0.4, 0,1}, {0.4 ,0 ,0 }, {0.4, 0, 1 },{0.2 , -0.2, 2}, {0.4, 0, 1}, {0.2, 0.2, 2}};




    // this is a 3*3 matrix with each has three features
    /*
    for(i = 0; i< 9;i++){
        A[i][0] = i*0.5;
        A[i][1] = i*0.5;
        A[i][2] = i*0.5;
    }
    */

    int k =4 ;
    //div = get_divs(B, B, 3, 2);
    printf("all matrix are 3*3, use %d-nearest neighbours\n\n",k);
    printf("A and A:\n");
    printf("non-turbulence itself is %lf\n\n", get_divs(A, A, 3,k));

    printf("B and B:\n");
    printf("turbulence itself is %lf\n\n", get_divs(B, B, 3,k));

    printf("A and B:\n");
    printf("the divergence of non-turbulence vs turbulence  is %lf\n\n", get_divs(A, B, 3, k));

    printf("B and A:\n");
    printf("turbulence vs non-turbulence is %lf\n\n", get_divs(B, A, 3, k));
    

    // how about two turbulence in oppsite direction?
    printf("B and C:\n");
    printf(" two turbulence in opposite direction, div is %lf\n\n", get_divs(B, C, 3, k));

    printf("C and C:\n");
    printf("C itself, div is %lf\n\n", get_divs(C, C, 3, k));

    printf("B and D:\n");
    printf("turbulence vs certain pattern, div is %lf\n\n", get_divs(B, D, 3, k));
    return 0;
}






