#include "get_divs.h"
// test the divgence caculation
int main(){
    float  div;
    int i, region_length=2;
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
    float A[27] = 
    {  
        0.4,0,2,
        0.4,0,1,
        0.4,0,2,
        0.4,0,1,
        0.4,0,0, 
        0.4,0,1,
        0.4,0,2, 
        0.4,0,1,
        0.4,0,2};

    // a sample turbulence
    float B[27] = {
        0.2,0.2, 2,
        0.4, 0 ,1,
        0.2, -0.2, 2, 
        0, 0.4,1,
        0 ,0 ,0,
        0, -0.4, 1,
        -0.2 , 0.2, 2,
        -0.4, 0, 1,
        -0.2, -0.2, 2};

    // also a turbulence
    float C[27];
    for(i = 0; i < 27 ; i++){
        if(i%3!=2){
            C[i] = -B[i];
        }
        else{ 
            C[i] = B[i];
        }
    }

    // another pattern
    float D[27] = {
        0.2,0.2, 2,
        0.4, 0 ,1,
        0.2, -0.2, 2,
        0.4, 0,1,
        0.4 ,0 ,0,
        0.4, 0, 1,
        0.2 , -0.2, 2,
        0.4, 0, 1,
        0.2, 0.2, 2};


    // this is a 3*3 matrix with each has three features
    /*
    for(i = 0; i< 9;i++){
        A[i][0] = i*0.5;
        A[i][1] = i*0.5;
        A[i][2] = i*0.5;
    }
    */

    int k = 4;
    //div = get_divs(B, B, 3, 2);

    int div_funcs[]={0,1,2};
    int num_div_funcs = 3;
    float a_a, a_b, a_c, a_d,b_a, b_b, b_c, b_d,c_a, c_b, c_c, c_d, d_a, d_b, d_c, d_d;

    for(i = 0; i < num_div_funcs; i++ ){
        // printf the current divergence function
        // slightly different from 'divergence' and kernel: divergence can be used to construct kernel
        if(i == 0){
            printf("*****linear kernel********\n");
        }else if(i == 1){
            printf("*****L2 divergence********\n");
        }
        else if( i == 2){
            printf("*****KL divergence********\n");
        }
        else{
            printf("div function not found!\n");
            exit(-1);
        }

        
        printf("all matrix are %dx%d, use %d-nearest neighbours\n\n",(region_length+1), (region_length +1), k);
        //printf("A and A:\n");
        a_a = get_divs(A, A, region_length, k, i); 
        //printf("non-turbulence itself is %lf\n\n", a_a );


        //printf("A and B:\n");
        a_b = get_divs(A, B, region_length ,k, i);
        //printf("the divergence of non-turbulence vs turbulence  is %lf\n\n", a_b);


        //printf("A and C\n");
        a_c= get_divs(A, C, region_length, k ,i);

        //printf("A and D\n");
        a_d= get_divs(A, D, region_length, k ,i);

        //printf("B and A:\n");
        b_a = get_divs(B, A, region_length,k, i);
        //printf("turbulence vs non-turbulence is %lf\n\n", b_a);

        //printf("B and B:\n");
        b_b = get_divs(B, B, region_length,k, i);
        //printf("turbulence itself is %lf\n\n",b_b );


        

        // how about two turbulence in oppsite direction?
        //printf("B and C:\n");
        b_c= get_divs(B, C, region_length,k, i);
        //printf(" two turbulence in opposite direction, div is %lf\n\n", b_c);


        //printf("B and D:\n");
        b_d= get_divs(B, D, region_length,k, i);
        //printf("turbulence vs certain pattern, div is %lf\n\n", b_d);

        //printf("C and A\n");
        c_a= get_divs(C, A, region_length, k ,i);

        //printf("C and B\n");
        c_b= get_divs(C, B, region_length, k ,i);

        //printf("C and C:\n");
        c_c= get_divs(C, C, region_length,k, i);
        //printf("C itself, div is %lf\n\n", c_c);

        //printf("C and D:\n");
        c_d= get_divs(C, D, region_length, k ,i);

        //printf("D and A:\n");
        d_a= get_divs(D, A, region_length, k ,i);

        //printf("D and B:\n");
        d_b= get_divs(D, B, region_length, k ,i);

        //printf("D and C:\n");
        d_c= get_divs(D, C, region_length, k ,i);

        //printf("D and D:\n");
        d_d= get_divs(D, D, region_length, k ,i);
        printf("summary\n");
        printf("%.4f\t%.4f\t%.4f\t%.4f\t\n", a_a, a_b, a_c, a_d);
        printf("%.4f\t%.4f\t%.4f\t%.4f\t\n", b_a, b_b, b_c, b_d);
        printf("%.4f\t%.4f\t%.4f\t%.4f\t\n", c_a, c_b, c_c, c_d);
        printf("%.4f\t%.4f\t%.4f\t%.4f\t\n", d_a, d_b, d_c, d_d);

    }
    return 0;
}






