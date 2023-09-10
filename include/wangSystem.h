/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#ifndef WANGSYSTEM_H
#define WANGSYSTEM_H
#include <stdio.h>
#include <math.h>

#include <ctime>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
 
extern std::vector<time_t> initial_time_v, terminal_time_v; 

namespace wangSystem{
    time_t initial_time(void);
    time_t terminal_time(void);
    void printMemoryUsage(void);
}



#endif