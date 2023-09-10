/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#ifndef VCELL_H
#define VCELL_H

#include <vector>
#include <iostream>
#include "class.h"

class Vertex;
class Line;
class Cell;
class Organ;

class Cell{
 public:
    std::vector<int> li;
    std::vector<int> vi;

    int n_edges;

    bool IsEpidermal;

    double area;
    double areaForce;

    int cellDivisionCount;
    double cellTime;

    double axisTheta;
    _vec<double> axis;
    double surfaceVertex[2];
    _vec<double> center;

    double outermostLength;
    
    double perimeter;

    int layer=-1;
    

    

    double regularity;
    double area_R;

    double y_rank;
    double area_rank;
    double frequency_modifier=1;

    int tag=0;

    double Gaussian_modifier=1;
    double area_modifier=1;

    //~Cell();

};


#endif