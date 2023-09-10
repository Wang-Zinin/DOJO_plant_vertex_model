/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <vector>
#include <random>
#include <time.h>
#include <cstring>
#include <assert.h>
#include <ctime>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>

#include "../include/class.h"
#include "../include/division.h"
#include "../include/vec.h"
#include "../include/parameter.h"
#include "../include/force.h"
#include "../include/geo2dv.h"
#include "../include/IOV.h"
#include "../include/wangSystem.h"
#include "../include/wangMath.h"


#include "../include/2dvOrgan.h"
#include "../include/2dvCell.h"
#include "../include/2dvLine.h"
#include "../include/2dvVertex.h"
#include "../include/2dvInitial.h"


//Cell division frequency control:
//1. No control "no"; 2. "area"; 3. "balance"; 4. "Gaussian_area"; 



//Cell division angles control:
//1. "random"; 2."anticlinal" (only for epidermal cells); 3. "periclinal" (only for epidermal cells); 4. "constant_0"; 5. "constant_90"
//6. "mochizuki_bias"; 

//Cell division frequency and angles combination control: 
//1. "biregion_frequency_angles_position"


using namespace std;

//string major_mode = "test"; 
//Major Mode: analysis, simualtion, test, or preparation

//string minor_mode = "single"; 

//minor mode for simulation
//single: only for a single simulation; 
//repeat: repeat the simulation with same parameter for several times; 
//batch: do simulation with defined chaning parameters

//minor mode for analysis
//single : single time frame; 
//time_lapse: all time frames of a single simulaiton;
//final: the final time frame of a single simulation; 
//batch_final: the final time frame of batch simulations


int main(){

parameter::read_mode(modeFile);


//recording time
vector<time_t> initial_time, terminal_time;
time_t initial_time_first = wangSystem::initial_time();
initial_time.push_back (initial_time_first);

cout<<"******************| Mode Selection| ********************************"<<endl;
cout<<"Major mode: "<<major_mode<<endl;
if(major_mode=="simulation"||major_mode=="analysis"||major_mode=="experiment"||major_mode=="plot")
cout<<"Minor mode: "<<minor_mode<<endl;
cout<<"********************************************************************"<<endl;


if(major_mode=="simulation"){

if(minor_mode=="single"){

//initialization
//parameter::read_and_record();
parameter::read(parameterInputFile);
parameter::record();

//exit(1);

Organ *p_g = new Organ;
readV::read_organ_txt(p_g,0);
organ_geo::epidermal_identity(p_g);

force::forceShapeInitiation(p_g,200000);
division::cell_time_initialization(p_g);

//simulation
cout<<"*********************| Simulation Start |***************************"<<endl;
for(int step=0; step<step_end; step++){
    
    p_g->step=step;
    force::calcForceMotion(p_g);

    //check cell division
    if(step%T_division_check==0){
        cout<<"At time step "<<step<<endl;
        //cout<<"Start geometrics calculation ---- ";
        geo::calcGeometrics(p_g);
        //cout<<"Output ---- ";
        output::geo(p_g);
        //cout<<"Finished "<<endl;

        //cout<<"Start division event judgement ---- ";
        division::Global(p_g);
        //cout<<"Division events output ---- ";
        output::division(p_g);
        //cout<<"Finished "<<endl;

        cout<<"********************************************************************"<<endl;
    }

    //cell time add
    if(step%T_cell_time==0){
        division::cellTimeAdd(p_g,1);
    }

    //output VTK
    if(step%T_vtkoutput==0){
        output::VTK(p_g);            
        if(termination::checkEnd(p_g)==1){
                    goto End_One_Simulation1;
        }
    }
}
End_One_Simulation1:;
    cout<<"********************************************************************"<<endl;

}

else if(minor_mode=="repeat"){

int repeat_time = 3;
parameter::read(parameterInputFile);
parameter::record();
initialOrganFile= "../"+initialOrganFile;

for(int repeat_i=0;repeat_i<repeat_time;repeat_i++){

char dir_file[100];
sprintf(dir_file,"%d/",repeat_i);
mkdir(dir_file,0777);
chdir(dir_file);

//initialization
Organ *p_g = new Organ;
readV::read_organ_txt(p_g,0);
organ_geo::epidermal_identity(p_g);

force::forceShapeInitiation(p_g,200000);
division::cell_time_initialization(p_g);

//simulation
cout<<"*********************| Simulation Start |***************************"<<endl;
for(int step=0; step<step_end; step++){
    p_g->step=step;
    force::calcForceMotion(p_g);

    //check cell division
    if(step%T_division_check==0){
        cout<<"At time step "<<step<<endl;
            //cout<<"Calculate geometrics information ---- ";
            geo::calcGeometrics(p_g);
            //cout<<"Output ---- ";
            output::geo(p_g);
            //cout<<"Finished "<<endl;

            //cout<<"Start division events judgement ---- "<<endl;
            division::Global(p_g);
            //cout<<"Division events output ---- ";            
            output::division(p_g);
            //cout<<"Finished "<<endl;
        cout<<"********************************************************************"<<endl;
    }

    //cell time add
    if(step%T_cell_time==0){
        division::cellTimeAdd(p_g,1);
    }

    //output VTK
    if(step%T_vtkoutput==0){
        output::VTK(p_g);            
        if(termination::checkEnd(p_g)==1){
                    goto End_One_Simulation2;
                }
    }
}
End_One_Simulation2:;
    cout<<"********************************************************************"<<endl;
chdir("../");
}

}

else{
    cout<<"Fatal error: no minor mode selected ! (major mode is simulation)"<<endl;
    exit(-1);
}

}

else{
    cout<<"Fatal error: no major mode selected !"<<endl;
    exit(-1);
}


//record termination time
time_t terminal_time_last = wangSystem::terminal_time();
terminal_time.push_back(terminal_time_last);

output::simulation_log(initial_time,terminal_time);

return 0;
}