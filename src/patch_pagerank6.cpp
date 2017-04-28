#include<iostream>
#include "common.h"

#define DISABLE_SCHEDULER

struct vertex_val_t{
    float cur;

    vertex_val_t(vertex_id id){
        cur = 1.0 ;
    }
    vertex_val_t(){}
};

typedef float msg_t;

#include "engine.hpp"


template<typename edge_t, typename vertex_val_t, typename msg_t>
inline bool gen_msg(graphzx::adjlst_inmem<edge_t, vertex_val_t> &adj, const int &iter, vertex_id &tid, msg_t &msg){
    msg = adj.val.cur / adj.num;
    return true;
}

template<typename edge_t, typename vertex_val_t, typename msg_t>
inline void process_vertex(graphzx::adjlst<edge_t, vertex_val_t> &adj, 
                           const int &iter,
                           msg_iter<edge_t, vertex_val_t, msg_t> &it){
    float agg = 0.0, msg;
    while (it(msg)) agg += msg;
    //singles.fetch_add(1);
    float newval = 0.15 + 0.85 * agg;
//    cout << newval << ' ' << adj.val.cur << endl;
    swap(adj.val.cur, newval);
//    if (ended) ended = false;
    if (fabs(adj.val.cur - newval) < 0.01) {
        return;
    }
    need_more_interation();
    
}

int main(int argc, char *argv[]){
    init_parameters(argc, argv);
    process();    
}

