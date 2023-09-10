/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>
#include "class.h"
#include "IOV.h"
#include "vec.h"
#include <fstream>
#include <iostream>
using namespace std;

//debug_mode
constexpr const char *debug_force_record ="OFF"; //when "ON", force will be recorded for each T_vtkoutput

//potential_energy_mode 
constexpr const char *potential_energy_mode = "simple"; //1. simple, suggested by Mochizuki-sensei, mostly used; 2. L_std

//list of parameters
extern int repeat_time;
extern int checkend_tmp;

//parameters about force
extern double sigma_L;
extern double sigma_O;
extern double kappa_S;
extern double S_std;
extern double eta;
extern double L_std;

//parameters about time and cell division
extern double standard_cell_period_length;
extern double F_apparent;
extern double F_modifier; //modify cell_division_threshold
extern int T_division_check;
extern int T_vtkoutput;
extern int T_cell_time;
extern double delta_time ; //originally set as 5e-3
extern int step_end ;
extern int end_cell_number ; //end when cell number reaches this number 

//"area_control"
extern double area_control_lower_limit;
extern double area_control_slope;

//"balance_control"
extern double BF;

//"Gaussian_control"
extern double gau_sigma;
extern double gau_mu;

//"mochizuki_bias" parameters
extern double mochizuki_bias_phi;
extern double mochizuki_bias_beta;

//list of files for input and output
extern string parameterInputFile;
extern string parameterRecordFile;
extern string initialOrganFile;
extern string modeFile;

//list of modes used
extern string major_mode;
extern string minor_mode;

extern string division_control;
extern string in_division_direction;
extern string epi_division_direction;
extern string temporal_mode;

extern string mechanics_mode;


namespace parameter{
    void read_mode(string);
    void read(string);
    void record();
    void batchRead(parameterList*);
}

namespace termination{
    bool checkEnd(Organ*);
}

#endif