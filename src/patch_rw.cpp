#include<iostream>
#include "common.h"
#include <cstdint>

//#define DISABLE_SCHEDULER

struct vertex_val_t{
    uint32_t cur;
    uint32_t inwalks;

    vertex_val_t(vertex_id id){
        if (id % 50 == 0)
            inwalks = 100;
        else
            inwalks = 0;
        cur = 0;
    }

    vertex_val_t(){}
};

typedef uint32_t msg_t;

#include "engine.hpp"

inline uint32_t fastrand(uint32_t &seed) {
    seed = 214013 * seed + 2531011;
    return seed;
}

template<typename edge_t, typename vertex_val_t, typename msg_t>
inline bool gen_msg(graphzx::adjlst_inmem<edge_t, vertex_val_t> &adj, const int &iter, vertex_id &tid, msg_t &msg){
    uint32_t seed = tid >> 4 + iter;
    msg = 0;
    for (uint32_t i = 0; i < adj.val.inwalks; i++) {
        if (fastrand(seed) % adj.num == 0) msg++;
    }
    return msg != 0;
}

template<typename edge_t, typename vertex_val_t, typename msg_t>
inline void process_vertex(graphzx::adjlst<edge_t, vertex_val_t> &adj, 
                           const int &iter,
                           msg_iter<edge_t, vertex_val_t, msg_t> &it){
    uint32_t agg = 0, msg;
    while (it(msg)) agg += msg;
    //singles.fetch_add(1);
    adj.val.cur += agg;
    adj.val.inwalks = agg;
    if (agg == 0) return;
    for (size_t i = 0; i < adj.num; i++) {
       // edges.fetch_add(1);
       if (!vertices[adj.edges[i]].scheduled) vertices[adj.edges[i]].scheduled = true;
    }
    need_more_interation();
    
}

int main(int argc, char *argv[]){
    init_parameters(argc, argv);
    max_iter = 5;
    process();    
}

