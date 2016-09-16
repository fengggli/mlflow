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

    int k =4;
    //div = get_divs(B, B, 3, 2);

    int div_funcs[]={0,1,2};
    int num_div_funcs = 3;
    double a_a, a_b, a_c, a_d,b_a, b_b, b_c, b_d,c_a, c_b, c_c, c_d, d_a, d_b, d_c, d_d;

    for(i = 0; i < num_div_funcs; i++ ){
        // printf the current divergence function
        // slightly different from 'divergence' and kernel: divergence can be used to construct kernel
        if(i == 0){
            printf("*****linear kernel********\n");
        }else if(i == 1){
            printf("*****L2 divergence********\n");
        }
        else if( i = 2){
            printf("*****KL divergence********\n");
        }
        else{
            printf("div function not found!\n");
            exit(-1);
        }

        
        a_a = get_divs(A, A, 3, k, i); 
        printf("all matrix are 3*3, use %d-nearest neighbours\n\n",k);
        printf("A and A:\n");
        printf("non-turbulence itself is %lf\n\n", a_a );

        a_b = get_divs(A, B, 3,k, i);
        printf("A and B:\n");
        printf("the divergence of non-turbulence vs turbulence  is %lf\n\n", a_b);


        printf("A and C\n");
        a_c= get_divs(A, C, 3, k ,i);

        printf("A and D\n");
        a_d= get_divs(A, D, 3, k ,i);

        b_a = get_divs(B, A, 3,k, i);
        printf("B and A:\n");
        printf("turbulence vs non-turbulence is %lf\n\n", b_a);

        b_b = get_divs(B, B, 3,k, i);
        printf("B and B:\n");
        printf("turbulence itself is %lf\n\n",b_b );


        

        // how about two turbulence in oppsite direction?
        b_c= get_divs(B, C, 3,k, i);
        printf("B and C:\n");
        printf(" two turbulence in opposite direction, div is %lf\n\n", b_c);


        b_d= get_divs(B, D, 3,k, i);
        printf("B and D:\n");
        printf("turbulence vs certain pattern, div is %lf\n\n", b_d);

        printf("C and A\n");
        c_a= get_divs(C, A, 3, k ,i);

        printf("C and B\n");
        c_b= get_divs(C, B, 3, k ,i);

        c_c= get_divs(C, C, 3,k, i);
        printf("C and C:\n");
        printf("C itself, div is %lf\n\n", c_c);

        printf("C and D:\n");
        c_d= get_divs(C, D, 3, k ,i);


        printf("D and A:\n");
        d_a= get_divs(D, A, 3, k ,i);

        printf("D and B:\n");
        d_b= get_divs(D, B, 3, k ,i);

        printf("D and C:\n");
        d_c= get_divs(D, C, 3, k ,i);

        printf("D and D:\n");
        d_d= get_divs(D, D, 3, k ,i);
        printf("summary\n");
        printf("%.4f\t%.4f\t%.4f\t%.4f\t\n", a_a, a_b, a_c, a_d);
        printf("%.4f\t%.4f\t%.4f\t%.4f\t\n", b_a, b_b, b_c, b_d);
        printf("%.4f\t%.4f\t%.4f\t%.4f\t\n", c_a, c_b, c_c, c_d);
        printf("%.4f\t%.4f\t%.4f\t%.4f\t\n", d_a, d_b, d_c, d_d);

    }
    return 0;
}






