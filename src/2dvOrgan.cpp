/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#include "../include/2dvOrgan.h"
#include "../include/2dvCell.h"
#include "../include/2dvLine.h"
#include "../include/2dvVertex.h"
#include "../include/class.h"

Organ::~Organ(){
        for(Vertex* ptr_v : p_v){
                delete ptr_v;
                ptr_v = nullptr;
        }
        for(Line* ptr_l : p_l) {
                delete ptr_l;
                ptr_l=nullptr;
        }
        for(Cell* ptr_c : p_c) {
                delete ptr_c;
                ptr_c=nullptr;
        }
        for(DivisionRecord* ptr_dr : d_r){
                delete ptr_dr;
                ptr_dr=nullptr;
        }
}

