#ifndef IO_TASK_HPP
#define IO_TASK_HPP

#include "common.h"
#include "boost/thread/mutex.hpp"

namespace graphzx {
#define QD 64

#ifndef DISABLE_SCHEDULER
    template<typename edge_t, typename vertex_val_t>
    struct adjlst_inmem   {
        bool scheduled;
//        vertex_id fid;
        counter num; //the number of out edges
        vertex_val_t val;
//        counter indegs;
//        edge_t *inedges;
//        edge_t *edges; //out edges
    };
#else
    template<typename edge_t, typename vertex_val_t>
    struct adjlst_inmem   {
//        vertex_id fid;
        counter num; //the number of out edges
        vertex_val_t val;
    };
#endif

    template<typename edge_t> 
    struct adjlst_assist {
        counter indegs = 0;
        edge_t *inedges = nullptr;
        edge_t *edges = nullptr;
    };

    template<typename edge_t, typename vertex_val_t>
    struct adjlst : adjlst_inmem<edge_t, vertex_val_t> {
//    struct adjlst  {
//        vertex_id fid;
//        counter num; //the number of out edges
//        vertex_val_t val;
        vertex_id fid;
        counter indegs;
        edge_t *inedges;
        edge_t *edges; //out edges

   public:
        adjlst<edge_t, vertex_val_t>(adjlst_inmem<edge_t, vertex_val_t> &adj_inmem){
            adjlst_inmem<edge_t, vertex_val_t>::fid = adj_inmem.fid;
            adjlst_inmem<edge_t, vertex_val_t>::num = adj_inmem.num;
            adjlst_inmem<edge_t, vertex_val_t>::val = adj_inmem.num;
//            fid = adj_inmem.fid;
//            num = adj_inmem.num;
//            val = adj_inmem.num;
        }

        //adjlst<edge_t, vertex_val_t>() : indegs(0), inedges(nullptr){
        adjlst<edge_t, vertex_val_t>() {
        }

    public:
        
        void tell_self(){
//            std::cout << "adj: " << fid << " " << num << std::endl;
//            for(int i=0; i<num; i++)
//                std::cout << edges[i] << " ";
//            std::cout << std::endl;
            //std::cout << "inadj: " << adjlst_inmem<edge_t, vertex_val_t>::fid << " " << indegs << std::endl;
            //for(int i=0; i<indegs; i++)
            //    std::cout << inedges[i] << " ";
            std::cout << std::endl;
        }
    };
    
}
#endif

