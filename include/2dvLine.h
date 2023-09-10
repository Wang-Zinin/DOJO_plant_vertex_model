/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#ifndef VLINE_H
#define VLINE_H

#include <vector>
#include <iostream>
#include "class.h"

class Vertex;
class Organ;
class Line;
class Cell;

class Line{
 public:
    int vi[2];
    int li;
    std::vector<int> ci;
    Vertex d1;
    Vertex d2;
    bool IsOutermost;

    double length;
    double edgeForce;
    double frc_edge;
    int ordered_array;

    double slope, intercept;
    double A,B,C; //general form of a line: Ax+By=C

    double calc_length(Organ);
    double calc_length(void);
    void calc_slope_intercept(Organ);
    void calc_slope_intercept(void);
    void calc_general_ABC(Organ);
    double distance_from_point(Vertex);
    void set_endpoints(Vertex,Vertex);
    void set_endpoints(Vertex*,Vertex*);
    bool same_segment(Line);
    void print_slope_intercept(void);

    
    //~Line();
};



#endif