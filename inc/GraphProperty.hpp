#ifndef GRAPHPROPERTY_HPP
#define	GRAPHPROPERTY_HPP

#include<iostream>
#include<libconfig.h++>
#include <sstream>

using namespace libconfig;
using namespace std;

const string val_suffix = "val";
const string ops_suffix = "val.ops";

class GraphProperty {
public:
    libconfig::Config cfg;
    string config_path;
    string root_dir;
    string fpath;
    unsigned long long active_vertices_num; //the number of vertices that has a out-edge
    unsigned long long vertices_num; //number of all vertices
    unsigned long long raw_vertices_num; //number of vertices without adding extras
    unsigned long long nvertices_per_partition;
    unsigned long long par_num; //the number of partitions
    unsigned long long active_par_num; //the number of partitions that have a vertex whose degree > 1
    bool converged;
    double diff; // record the  variations for each iteration 
    size_t shuffle_size;
    size_t nrow;

    
    GraphProperty(string _fpath, const size_t _shuffle_size = 0) 
          : shuffle_size(_shuffle_size) 
    {
        fpath = _fpath;
        root_dir = fpath+".nodos_dir/";
        config_path = root_dir+"properties.cfg";
        init();
        active_par_num = (active_vertices_num-1)/nvertices_per_partition +1;
 //       if (_shuffle_size != 0) {
        if (shuffle_size != 0) {
            vertices_num = (raw_vertices_num-1)/shuffle_size * shuffle_size + shuffle_size;
            nrow = vertices_num / shuffle_size;
        }
    }
    
//    GraphProperty(string _fpath) : GraphProperty(_fpath, _shuffle_size) 
//    {
//    }

    virtual ~GraphProperty(){}
    
private:
    void init() {
        cfg.readFile(config_path.c_str());
        const Setting& root = cfg.getRoot();
        const Setting& graph = root["graph"];
        graph.lookupValue("active_vertices_num", active_vertices_num);
        graph.lookupValue("vertices_num", vertices_num);
        raw_vertices_num = vertices_num;
        graph.lookupValue("partition_num", par_num);
        graph.lookupValue("shuffle_size", (long long&)shuffle_size);
        graph.lookupValue("nvertices_per_partition", nvertices_per_partition);
    }

public:

    string get_vals_path(int nth_par) {
        stringstream ss;
        ss << nth_par;
        string valspath;
        ss >> valspath;
        valspath = root_dir + valspath + std::string("/") + val_suffix;
        return valspath;
    }

    string get_ops_path(int nth_par) {
        stringstream ss;
        ss << nth_par;
        string opspath;
        ss >> opspath;
        opspath = root_dir + opspath + std::string("/") + ops_suffix;
        return opspath;
    }
    
    string get_par_path(int nth_par){
        stringstream ss;
        ss << nth_par;
        string parpath;
        ss >> parpath;
        parpath = root_dir + parpath + std::string("/");
        return parpath;        
    }
    
    size_t get_vertices_num(int nth_par){
        if(nth_par < par_num-1)
            return nvertices_per_partition;
        else
            return vertices_num % nvertices_per_partition;
    }
};

//int main(){
//    GraphProperty gp(string("/home/ubu/Downloads/com-lj.ungraph.txt.1000"));
//    cout<<gp.vertices_num<<endl;
//    
//}

#endif	/* GRAPHPROPERTY_HPP */


