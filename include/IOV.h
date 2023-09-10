/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef _IOV_H
#define _IOV_H

#include "vec.h"
#include "parameter.h"
#include "class.h"

#include <fstream>
#include <iostream>
#include <dirent.h>
#include <unistd.h>
#include <cstring>
#include <ctime>
#include <tuple>
#include <cstdio>

using namespace std;

namespace readV{
    void read_organ_txt(Organ*,int);

    void oneVTK(Organ*,string,string);
    void oneCell(Organ*,string);
    void oneLine(Organ*,string);
    void oneCellLine(Organ*);

    //void allVTK(Organ*,string);

    void final_VTK(Organ*,string);
    vector<Vertex*> xy_txt_to_vertexXY(string);
    vector<double> read_vd(string);
    vector<Vertex*> read_vv(string);

    Geometrics_analysis_class read_geo_output_txt(string);
}

namespace output{
    void VTK(Organ*);
    void geo(Organ*);
    void division(Organ*);
    void batch(Batch*);
    void batch_final_analysis(vector<double>,vector<double>,vector<double>,string);
    void simulation_log(vector<time_t>, vector<time_t>);
    void geo_initial(void);
}

namespace file_process{

int getFileNum(const string&);

}

namespace cout_fout_debug{
    void cout_vector_vertex(vector<Vertex*>);
    void fout_vector_vertex(vector<Vertex*>,string);
    void cout_vector_double(vector<double>);
    void cout_vector_int(vector<int>);
    void cout_vector_string(vector<string>);
    void cout_vector_line(vector<Line*>);
    void fout_vector_double(vector<double>,string);
    void fout_vector_int(vector<int>,string);
    void fout_pair_double(vector<pair<double,double>>,string);
    void fout_tuple_double(vector<tuple<double,double,double>>,string);
}

#endif