/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#ifndef VVERTEX_H
#define VVERTEX_H

#include "vec.h"
#include <vector>
#include <iostream>
#include "class.h"


class Vertex;
class Line;
class Cell;
class Organ;

class Vertex{
 public:
    std::vector<int> li; 
    std::vector<int> ci;
    _vec<double> loc; //Cartesian coordinates
    _vec<double> frc_area;
    _vec<double> frc_edge;
    
    double r,theta; //Polar coordinates
    unsigned int occurrenceInCell;
    bool IsSurface;

    int vi;
    std::vector<int> li_array;
    int vi_array;
    //destructor
   
    void print_Cartesian(void);
    void print_Polar(void);
    double distance_from_vertex(Vertex);
    void Cartesian_to_Polar(void);
    void Polar_to_Cartesian(void);
    bool collinear_points(Vertex,Vertex);
    bool same_vertex(Vertex);

};



#endif