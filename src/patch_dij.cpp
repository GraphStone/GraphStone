#include<iostream>
#include "common.h"
#include <climits>

//#define DISABLE_SCHEDULER
extern vertex_id start;

struct vertex_val_t{
    vertex_id cur;

    vertex_val_t(vertex_id id){
        if (id == start) cur = 0;
        else cur = UINT_MAX - 2;
    }
    vertex_val_t(){}
};

typedef vertex_id msg_t;

#include "engine.hpp"


template<typename edge_t, typename vertex_val_t, typename msg_t>
inline bool gen_msg(graphzx::adjlst_inmem<edge_t, vertex_val_t> &adj, const int &iter, vertex_id &tid, msg_t &msg){
    msg =  adj.val.cur;
    return true;
}

template<typename edge_t, typename vertex_val_t, typename msg_t>
inline void process_vertex(graphzx::adjlst<edge_t, vertex_val_t> &adj, 
                           const int &iter,
                           msg_iter<edge_t, vertex_val_t, msg_t> &it){
//    singles.fetch_add(1);
    vertex_id minimum = UINT_MAX - 2, msg;
    while (it(msg)) minimum = min(minimum, msg);
    if (minimum + 1 >= adj.val.cur && iter != 0) return;
    if (minimum + 1 < adj.val.cur) adj.val.cur = minimum + 1;
    for (size_t i = 0; i < adj.num; i++) {
//        edges2.fetch_add(1);
        schedule_vertex(adj.edges[i]);
    }
    need_more_interation();
}

int main(int argc, char *argv[]){
    init_parameters(argc, argv);
    all_active = false;
    cout << "sizeof(adjlst_inmem) = " << sizeof(adjlst_inmem<edge_t, vertex_val_t>) << endl;
    process();    
}

