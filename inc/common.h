#ifndef COMMON_H 
#define COMMON_H
#include <sys/time.h>

typedef unsigned int vertex_id; 
typedef unsigned int counter;

inline void print_nowtime()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long long mslong = (long long) tp.tv_sec * 1000L + tp.tv_usec / 1000; //get current timestamp in milliseconds
    std::cout << "nowtime: " << mslong << std::endl;
}
#endif
