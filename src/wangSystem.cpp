/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#include "../include/wangSystem.h"

std::vector<time_t> initial_time_v, terminal_time_v; 

namespace wangSystem{

    time_t initial_time(void){
        time_t curr_time;
        tm * curr_tm;
        char initial_time_string [100];
        time(&curr_time);
        curr_tm = localtime(&curr_time);
        strftime(initial_time_string, 100, "Simulation starts at %T, %Y %B %d", curr_tm);
        cout << initial_time_string << endl;
        initial_time_v.push_back(curr_time);
        return curr_time;
        //reference: https://www.programiz.com/cpp-programming/library-function/ctime/strftime
    }

    time_t terminal_time(void){
        time_t curr_time;
        tm * curr_tm;
        char terminal_time_string [100];
        time(&curr_time);
        curr_tm = localtime(&curr_time);
        strftime(terminal_time_string, 100, "Simulation terminates at %T, %Y %B %d", curr_tm);
        cout << terminal_time_string << endl;
        terminal_time_v.push_back(curr_time);
        return curr_time;
        //reference: https://www.programiz.com/cpp-programming/library-function/ctime/strftime
    }

    void printMemoryUsage() {
    std::ifstream statusFile("/proc/self/status");
    std::string line;

    while (std::getline(statusFile, line)) {
        if (line.find("VmSize") != std::string::npos) {
                std::cout << line << std::endl;
            }
        }
    }
}