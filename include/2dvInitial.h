/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#ifndef VINITIAL_H
#define VINITIAL_H

#include "2dvOrgan.h"
#include "class.h"
#include "force.h"
#include "geo2dv.h"
#include "division.h"
#include "IOV.h"

extern string real_organ_contour_imagej_txt;

namespace initialization{
    void organ(Organ*);
}


#endif