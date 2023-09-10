/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef _DIVISION_H
#define _DIVISION_H

#include "class.h"
#include "vec.h"
#include "wangMath.h"
#include "parameter.h"
#include "geo2dv.h"
#include "force.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <random>
#include <math.h> 


namespace division_frequency{
    void no_control(Organ*);
    void area_control(Organ*);
    void balance_control(Organ*);
    double calcGaussian(double,double,double);
    void Gaussian_control(Organ*,double,double);
}

namespace division_direction{
    _vec<double> angles(Organ*,int);
    _vec<double> random(Organ*,int);
    _vec<double> epi_anticlinal(Organ*,int);
    _vec<double> epi_periclinal(Organ*,int);
    _vec<double> constant_0(Organ*,int);
    _vec<double> constant_90(Organ*,int);
    _vec<double> mochizuki_bias(Organ*,int,double,double);
}

namespace division{
    void cell_time_initialization(Organ*);
    void cellTimeAdd(Organ *p_g, int increase_index);
    void Global(Organ*);
    void One(Organ*, int);
    void Record(Organ*,int); //have issue in recording frequency_modifier: the new cell's frequency_modifier is set to be 1 by default
    void tag_control(Organ*);
}

#endif