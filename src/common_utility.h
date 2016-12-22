// some common functios 
#ifndef COMMON_UTILITY_H
#define COMMON_UTILITY_H

#include <sys/time.h>
#include <stdio.h>

/*******************************
 * for sequential use
 *******************************/

// get current time
double get_cur_time(); 


/*******************************
 * parallel use
 *******************************/
void my_message(char *msg, int rank);

#endif
