/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#ifndef VORGAN_H
#define VORGAN_H

#include <vector>
#include "vec.h"
#include "class.h"

class Vertex;
class Line;
class DivisionRecord;
class Cell;
class Geometrics;
class ForceCheck;

class Organ{
 public:
    int step;
    std::vector<Vertex *> p_v;
    std::vector<Line *>p_l;
    std::vector<Cell *>p_c;
    std::vector<DivisionRecord *> d_r;

    std::vector<ForceCheck* >forceCheck;
    double initial_cell_number;

    //geometrics analysis
    //1. epidermal identity
    int N_inner_cell;
    int N_epi_cell;
    std::vector<int> surface_vertex;
    std::vector<int> surface_line;
    //2. organ center
    _vec<double> center;
    //3. organ cell layer
    int cell_layer_number;
    //4. organ area
    double area; //the summed area for all cells
    double epiArea; //the summed area for all epidermal cells
    double inArea;  //the summed area for all inner cells
    double epiArea_averaged;
    double inArea_averaged;
    double area_averaged;
    //5. organ perimeter
    double perimeter;
    double perimeter_averaged;
    //6. organ circularity
    double circularity; 
    //7. similarity index
    double similarity_index;

    //analysis of area control division frequency ratio between epidermal cells and inner cells
    double av_in_division_frequency_modifier;
    double av_epi_division_frequency_modifier;
    double F_ratio; //divsion frequency ratio between epidermal cells and inner cells
    

    //more geometrics
    double organ_potential_energy;
    double y_min_v;
    double x_min_v;
    ~Organ();

    //debugs
    void print_basics(void){
        std::cout<<"Cell number: "<<p_c.size()<<"; line number: "<<p_l.size()<<"; vertex number: "<<p_v.size()<<std::endl;
    }
};



#endif