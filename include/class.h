/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef _CLASS_H
#define _CLASS_H

#include <iostream>
#include <vector>
#include "vec.h"
#include "2dvVertex.h"

constexpr int DEGREE_ACCURACY =2;
constexpr double EPS_geo = 1e-5;
using namespace std;
//list of class
class Organ;
class Vertex;
class Line;
class Cell;
class DivisionRecord;
class Geometrics;
class VertexCheck;
class CellCheck;
class Circle;
class Batch;
class parameterList;
class Intersection_relationship;
class Distance_point_line;
class Ordered_boundary;
class Quadratic_Solution;

//Quadratic_Solution Quadratic_Equation_Solve_(double,double,double);

class parameterList{
 public:
    std::vector<double> parameter_1;
    std::vector<double> parameter_2;
};

class Batch{
 public:
    std::vector<double> organ_perimeter;
    std::vector<double> organ_area;
    std::vector<double> averaged_perimeter;
    std::vector<double> averaged_inner_area;
    std::vector<double> averaged_epi_area;
    std::vector<int> N_in;
    std::vector<int> N_epi;
    std::vector<double> overlap_area;
    std::vector<double> real_area;
    std::vector<double> overlap_index;
    std::vector<double> averaged_regularity_in;
    std::vector<double> averaged_regularity_epi;
};

class Geometrics_analysis_class{
 public:
    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> value;
};


class DivisionRecord{
 public:
    int time;
    int cidx;
    bool IsEpidermal;
    double axisTheta; 

    double center_x;
    double center_y;

    int division_count;
    double av_in_frequency_modifier;
    double av_epi_frequency_modifier;
};

class Circle{
 public:
    _vec<double> center;
    double radius;
};

class Quadratic_Solution{
 public:
    double delta;
    double x1;
    double x2;
};

class Intersection_relationship{
 public:
    int Relationship;
    vector<Vertex> cross_points;
};

class Distance_point_line{
 public:
    double distance;
    Vertex Closest_Point;
    double t;
};

class Ordered_boundary{
 public:
    vector<int> li;
    vector<Vertex> vi;
};

#endif
