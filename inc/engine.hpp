#include<iostream>
#include<fstream>
#include "io_task.hpp"
#include "common.h"
#include "GraphProperty.hpp"
#include <sys/stat.h>
#include <vector>
#include <atomic>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <boost/program_options.hpp>
#include <omp.h>

using namespace boost::program_options;

using namespace graphzx;

size_t win_size = 1000000;

size_t tracker_win_id = 0;

typedef uint32_t edge_t;
std::atomic<long long> edges(0);
std::atomic<long long> edges2(0);

volatile bool ended = false;
std::atomic<long long> singles(0);

bool do_random_shuffle = false;
bool all_active = true; //used for initialize scheduler, all vertices are activated
vertex_id start = 0; // used for initialize scheduler, only vertex start is activated

//bool *scheduled;

adjlst_inmem<edge_t, vertex_val_t> *vertices; //store all vertices' values and accumlated messages
adjlst<edge_t, vertex_val_t> *patch; //store all vertices' values and accumlated messages
//adjlst_inmem<edge_t, vertex_val_t> *patch; //store all vertices' values and accumlated messages
adjlst_assist<edge_t> *assists; //add support for edges offsets to be used for shuffling

vector<edge_t *> bufs;
vector<edge_t *> inbufs;
vector<size_t *> offs;
vector<size_t *> inoffs;

//vertices *vmaps; //map old vertex id to new id

string fpath;
int max_iter = 100;
GraphProperty *gp;

inline void need_more_interation() {
    if (ended) ended = false;
//        std::cout << "ended " << ended << std::endl;
}

inline void schedule_vertex(vertex_id &v) {
#ifndef DISABLE_SCHEDULER
    if (!vertices[v].scheduled) vertices[v].scheduled = true;
#endif
}

size_t getFilesize(const std::string& filename) {
    struct stat st;
    assert(stat(filename.c_str(), &st) == 0); 
    return st.st_size;   
}

void tell_vertices() {
    for (size_t i = 0; i < gp->vertices_num; i++) {
//        vertices[i].tell_self();
        std::cout << i << ' '
                  << vertices[i].val.cur << ' '
                  //<< vertices[i].scheduled 
                  << std::endl;
    }

}

vertex_id map2newid(vertex_id id) {
    vertex_id res = id % gp->shuffle_size * gp->nrow;
    res += id / gp->shuffle_size;
    return res;
//    return gp->vertices_num128 - 1 - res;
}

vertex_id map2oldid(vertex_id id) {
//    id = gp->vertices_num128 - 1 - id;
    vertex_id res = id % gp->nrow * gp->shuffle_size;
    res += id / gp->nrow;
    return res;
}


void *my_mmap(const char *fpath) {
    int fd = open(fpath, O_RDONLY);
    if (fd == -1) cout << fpath << " " << fd << ' '  << strerror(errno) << endl;
    assert (fd > 0);
    int flag = PROT_READ;
//    if (do_random_shuffle) flag |= PROT_WRITE;
    void *res = mmap(NULL, getFilesize(fpath), flag,
                       MAP_PRIVATE, fd, 0);
    if (res == MAP_FAILED) cout << "mmap error "  << strerror(errno) << endl;
    assert(res != MAP_FAILED);
    return res;
}
    
void init_per_par(int nth_par){
    string base_path = gp->get_par_path(nth_par) + string("bases.idx");
    string edge_path = gp->get_par_path(nth_par) + string("output.graph");
    offs[nth_par] = (size_t*)my_mmap(base_path.c_str());
    bufs[nth_par] = (edge_t*)my_mmap(edge_path.c_str());

    size_t block_size = 256;
    //size_t block_size = 1024;
    size_t startid, lastid;
    startid = gp->nvertices_per_partition * nth_par;
    if (nth_par != gp->par_num -1)
        lastid = gp->nvertices_per_partition * (nth_par + 1);
    else
        lastid = gp->raw_vertices_num;

    size_t *bases = offs[nth_par];
    size_t edge_base = (size_t)bufs[nth_par];
    #pragma omp parallel for num_threads(40)
    for(size_t id = startid; id < lastid; id += block_size) {
        for (size_t i = 0; i < block_size && id + i < lastid; i++) {
            size_t vid = id + i;
            size_t start = edge_base, end = bases[vid-startid] + edge_base;
            if (vid - startid != 0) start = bases[vid-startid-1] + edge_base;
            vertices[vid].num = (end - start) / sizeof(edge_t);
        }
    }
}

void init_per_par_inedges(int nth_par){
    string base_path = gp->get_par_path(nth_par) + string("inbases.idx");
    string edge_path = gp->get_par_path(nth_par) + string("inoutput.graph");
    inoffs[nth_par] = (size_t*)my_mmap(base_path.c_str());
    inbufs[nth_par] = (edge_t*)my_mmap(edge_path.c_str());
}

void init_global(string _fpath){
    gp = new GraphProperty(_fpath);
    //vertices = new adjlst_inmem<edge_t, vertex_val_t>[gp->vertices_num];
    vertices = (adjlst_inmem<edge_t, vertex_val_t> *) malloc (
                sizeof(adjlst_inmem<edge_t, vertex_val_t>) *
                gp->vertices_num);
    patch = new adjlst<edge_t, vertex_val_t>[win_size];
    //patch = new adjlst_inmem<edge_t, vertex_val_t>[win_size];
    //scheduled = new bool[gp->vertices_num];
    #pragma omp parallel for
    for (size_t i = 0; i < gp->vertices_num; i++) {
        vertices[i].val = vertex_val_t(i);
#ifndef DISABLE_SCHEDULER
        vertices[i].scheduled = all_active;
#endif
    }

    for (int i = 0; i < gp->par_num; i++) {
        bufs.push_back(nullptr);
        inbufs.push_back(nullptr);
        offs.push_back(nullptr);
        inoffs.push_back(nullptr);
    }
    cout << "start loading" << endl;
    #pragma omp parallel for
    for (int i = 0; i < gp->par_num; i++) {
        init_per_par(i);
        init_per_par_inedges(i);
    }
    cout << "loading down" << endl;
}



template<typename edge_t, typename vertex_val_t, typename msg_t>
inline bool gen_msg(graphzx::adjlst_inmem<edge_t, vertex_val_t> &adj, int &iter, vertex_id &tid, msg_t &msg);

template<typename edge_t, typename vertex_val_t, typename msg_t>
struct msg_iter{
    size_t i = 0;
    graphzx::adjlst<edge_t, vertex_val_t> &adj;
    const int &iter;

    msg_iter<edge_t, vertex_val_t, msg_t>(graphzx::adjlst<edge_t, vertex_val_t> &_adj, const int _iter) :
    adj(_adj), iter(_iter) {}
    
    inline bool operator() (msg_t &msg) {
        while (i < adj.indegs) {
            if(gen_msg(vertices[adj.inedges[i]], iter, adj.fid, msg)) {
                i++;
                return true;
            }
            else i++; 
        }
        return false;
    }
};

template<typename edge_t, typename vertex_val_t, typename msg_t>
inline void process_vertex(graphzx::adjlst<edge_t, vertex_val_t> &adj,
                           const int &iter,
                           msg_iter<edge_t, vertex_val_t, msg_t> &it);

void iter1win(const int &iter, const size_t  win_id){
    tracker_win_id = win_id;
    //cout << "start processing in-mem: " << iter \
         << '\t' << win_id << endl;
    //size_t block_size = 200;
    size_t block_size = 256;
//    size_t block_size = 128;
    //size_t block_size = 1024;
    size_t lastid = min(win_id+win_size, (size_t)gp->vertices_num);
    #pragma omp parallel for schedule(static)
    for(size_t bid = win_id; bid < lastid; bid += block_size) {
        for (size_t id = bid; id < bid + block_size && id < lastid; id++) {
            vertex_id off_id = id-win_id;
#ifndef DISABLE_SCHEDULER
            if (!vertices[id].scheduled) { 
                patch[off_id].scheduled = false;
                continue;
            }
#endif
 //           patch[off_id] = (adjlst<edge_t, vertex_val_t>)vertices[id];
 //           if (vertices[id].fid == 101) assert(0 != vertices[id].scheduled);
 //           if (patch[off_id].fid == 101) assert(0 != patch[off_id].scheduled);
 //           vertices[id].scheduled = false;
//            continue;
            //patch[off_id].fid = vertices[id].fid;
            patch[off_id].fid = id;
            patch[off_id].val = vertices[id].val;
            patch[off_id].num = vertices[id].num;
#ifndef DISABLE_SCHEDULER
            patch[off_id].scheduled = vertices[id].scheduled;
            vertices[id].scheduled = false;
#endif
            size_t par = id / gp->nvertices_per_partition;
            size_t ith_id = id % gp->nvertices_per_partition;
            if (ith_id == 0) {
                patch[off_id].inedges = inbufs[par];
                patch[off_id].indegs = inoffs[par][0]/sizeof(edge_t);
                //patch[off_id].indegs = inoffs[par][0] >> 2;
#ifndef DISABLE_SCHEDULER
                patch[off_id].edges = bufs[par];
#endif
            }
            else {
                patch[off_id].inedges = (edge_t*)((size_t)inbufs[par] + (size_t)inoffs[par][ith_id - 1]);
                patch[off_id].indegs = (inoffs[par][ith_id] - inoffs[par][ith_id - 1])
                                       / sizeof(edge_t);
#ifndef DISABLE_SCHEDULER
                patch[off_id].edges = (edge_t*)((size_t)bufs[par] + (size_t)offs[par][ith_id - 1]);
#endif
            }
        }
    }
    //memcpy(patch, vertices + win_id, sizeof(adjlst<edge_t, vertex_val_t>)*(lastid-win_id));  
    lastid -= win_id;
//if (false) {

#ifdef DISABLE_SCHEDULER
    #pragma omp parallel for schedule(dynamic)
    for(size_t bid = 0; bid < lastid; bid += block_size) {
        for (size_t id = bid; id < bid + block_size && id < lastid; id++) {
             msg_iter<edge_t, vertex_val_t, msg_t> it(patch[id], iter);
             process_vertex(patch[id], iter, it);
        }
    }
#else
    #pragma omp parallel for schedule(dynamic)
    for(size_t bid = 0; bid < lastid; bid += block_size) {
        for (size_t id = bid; id < bid + block_size && id < lastid; id++) {
            //if (patch[id].fid == 101) assert(0 != patch[id].scheduled);
            if (patch[id].scheduled) {
                 //singles.fetch_add(1);
                 msg_iter<edge_t, vertex_val_t, msg_t> it(patch[id], iter);
                 process_vertex(patch[id], iter, it);
            }
        }
    }
#endif
    lastid += win_id;

    #pragma omp parallel for schedule(static)
    for(size_t bid = win_id; bid < lastid; bid += block_size) {
        for (size_t id = bid; id < bid + block_size && id < lastid; id++) {
#ifndef DISABLE_SCHEDULER
            if (patch[id-win_id].scheduled)
#endif
            vertices[id].val = patch[id-win_id].val;
        }
    }
}

void iterate(int iter){
    cout << "start processing in-mem: " << iter 
         << '\t' << gp->vertices_num << endl;
    size_t n0 = 0;
    size_t sum = 0;
#pragma omp parallel for reduction(+:sum) 
    for (size_t id = 0; id < gp->vertices_num; id++) {
//        if (vertices[id].scheduled) sum++;
    }
    cout << "sum = " << sum << endl;
    for (size_t id = 0; id < gp->vertices_num; id+=win_size) {
        iter1win(iter, id);
    }
     sum = 0;

//if (false)
#pragma omp parallel for reduction(+:n0,sum) 
    for (size_t id = 0; id < gp->vertices_num; id++) {
        if (vertices[id].val.cur == 0) n0++;
        sum += vertices[id].val.cur;
    }
    cout << "n0 = " << n0 << endl;
    cout << "sum = " << sum << endl;

}

/*
void process(string fpath, 
             int max_iter, 
             bool _all_active = false, 
             vertex_id _start = 0
            )
*/
void process()
{
//    do_random_shuffle = false;
    cout << " start = " << start << endl;
    //init_global(fpath.c_str(), do_random_shuffle);
    init_global(fpath.c_str());
#ifndef DISABLE_SCHEDULER
    if (do_random_shuffle) vertices[map2newid(start)].scheduled = true;
    else vertices[start].scheduled = true;
#endif
    std::cout << "win_size = " << win_size << std::endl;
    print_nowtime();
    for (int iter = 0; iter < max_iter && !ended; iter++) {
        ended = true;
//        ended = false;
//        __sync_synchronize();
        singles.store(0);
        iterate(iter);
 //       __sync_synchronize();
        std::cout << edges.load() << std::endl;
//        tell_vertices();
        if (do_random_shuffle) {
            std::cout << map2newid(9766) << " " << vertices[map2newid(9766)].val.cur << std::endl;
            std::cout << map2newid(9765) << " " << vertices[map2newid(9765)].val.cur << std::endl;
        } else {
            std::cout << vertices[9766].val.cur << std::endl;
            std::cout << vertices[9765].val.cur << std::endl;
        }
        std::cout << vertices[9766].val.cur << std::endl;
        std::cout << vertices[9765].val.cur << std::endl;
        std::cout << vertices[0].val.cur << std::endl;
        std::cout << vertices[1].val.cur << std::endl;
        std::cout << vertices[2].val.cur << std::endl;
        std::cout << "ended " << ended << std::endl;
/*
*/
        std::cout << "singles = "<< singles.load() << std::endl;
        std::cout << "edges = "<< edges.load() << std::endl;
        std::cout << "edges2 = "<< edges2.load() << std::endl;
    }
    print_nowtime();

    atomic<size_t> singles(0);
    std::cout << "signles = "<< singles.load() << std::endl;
    std::cout << edges.load() << std::endl;
}

void init_parameters(int argc, char *argv[]) {
    try
    {
        int nthreads = 0;
        options_description desc{"Options"};
        desc.add_options()
            ("help,h", "Help screen")
            ("file,f", value<string>(&fpath)->default_value(""), "original file path")
            ("winsize,w", value<size_t>(&win_size)->default_value(1000000), "The size of window")
            ("maxiter,m", value<int>(&max_iter)->default_value(100), "max number of interations")
            ("all_active,a", value<bool>(&all_active)->default_value(true), "whether activate all vertices at first iteration")
            ("start,s", value<vertex_id>(&start)->default_value(101), "if not all vertices activated, the start point")
            ("random_shuffle,r", value<bool>(&do_random_shuffle)->default_value(false), "whether do the random shuffling")
            ("nthreads,n", value<int>(&nthreads)->default_value(0), "whether do the random shuffling");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        if (nthreads) omp_set_num_threads(nthreads);
        
        #define get_val(val)  cout << #val" = " << val << endl
        get_val(fpath);
        get_val(win_size);
        get_val(max_iter);
        get_val(all_active);
        get_val(start);
        get_val(do_random_shuffle);
        
        assert(win_size != 0);
    }
    catch (const error &ex)
    {
        std::cerr << ex.what() << '\n';
    }
//    exit(0);
}

