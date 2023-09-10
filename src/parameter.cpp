/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#include "../include/parameter.h"
#include "../include/2dvOrgan.h"
#include "../include/2dvCell.h"
#include "../include/2dvLine.h"
#include "../include/2dvVertex.h"

string initialOrganFile= "../input/61cells.txt";
string parameterInputFile = "../input/parameter.txt";
string modeFile = "../input/mode.txt";
string parameterRecordFile = "parameterRecord.txt";
string parameterListFile = "batchParameter.txt";

double kappa_S;
double S_std;
double delta_time;
double eta;
double standard_cell_period_length;
double F_apparent;
int T_division_check;
int T_vtkoutput;
int T_cell_time;
int end_cell_number;
int step_end;
double sigma_L;
double sigma_O;
int checkend_tmp=0;

//"balance_control"
double BF;

//"area_control"
double area_control_lower_limit;
double area_control_slope;

double gau_mu;
double gau_sigma;

double temporal_Gaussian_a;
double temporal_Gaussian_b;

//"temporal_Gaussian_linear_mu" parameters
double temporal_Gaussian_initial_mu;
double temporal_Gaussian_terminal_mu;

//"temporal_Gaussian_linear_sigma" parameters
double temporal_Gaussian_initial_sigma;
double temporal_Gaussian_terminal_sigma;

// "Gaussian_xy" parameters
double gau_mu_x;
double gau_sigma_x;
double gau_mu_y;
double gau_sigma_y;

//"uniform_to_Gaussian"
double uniform_to_Gaussian_transition_time;

int temporal_t1;
int temporal_t2;

//"mochizuki_bias" parameters
double mochizuki_bias_phi;
double mochizuki_bias_beta;

//"mochizuki_bias_asymmetrical" parameters
double mochizuki_bias_beta_left;
double mochizuki_bias_beta_right;
double mochizuki_bias_phi_left;
double mochizuki_bias_phi_right;

//"mochizuki_bias_apical_basal" parameters
double mochizuki_bias_beta_apical;
double mochizuki_bias_beta_basal;
double mochizuki_bias_phi_apical;
double mochizuki_bias_phi_basal;
double angle_bias_y_boundary;

//"biregion_angles_position" parameters
double biregion_angles_position_apical_gaussian_beta;
double biregion_angles_position_basal_gaussian_beta;
double biregion_angles_position_apical_gaussian_phi;
double biregion_angles_position_basal_gaussian_phi;
double biregion_angles_position_y_boundary;

//"biregion_angles_identity" parameters
double biregion_angles_identity_apical_gaussian_beta;
double biregion_angles_identity_basal_gaussian_beta;
double biregion_angles_identity_apical_gaussian_phi;
double biregion_angles_identity_basal_gaussian_phi;
double biregion_angles_identity_y_boundary;

//"biregion_frequency_position" parameters
double biregion_frequency_position_y_boundary;
double biregion_frequency_position_relative_frequency;

//"biregion_frequency_identity" parameters
double biregion_identity_y_boundary;
double biregion_identity_relative_frequency;

//"biregion_frequency_angles_position" parameters
double biregion_frequency_angles_position_apical_gaussian_beta;
double biregion_frequency_angles_position_basal_gaussian_beta;
double biregion_frequency_angles_position_apical_gaussian_phi;
double biregion_frequency_angles_position_basal_gaussian_phi;
double biregion_frequency_angles_position_relative_frequency;
double biregion_frequency_angles_position_y_boundary;

//"temporal_angle_bias" parameters
double temporal_bias_beta_initial;
double temporal_bias_beta_terminal;
double temporal_bias_phi_initial;
double temporal_bias_phi_terminal;

//"Gaussian_bias" parameters
double Gaussian_bias_phi;
double Gaussian_bias_beta;

//"Gaussian_bias_continuous_Gaussian" parameters
double Gaussian_bias_beta_A;
double Gaussian_bias_beta_mu;
double Gaussian_bias_beta_sigma;
//extern double Gaussian_bias_phi;

//"temporal_biregion_angles" parameters
double y_boundary_initial;
double y_boundary_change;
double y_boundary_terminal;
double temporal_biregion_angles_apical_gaussian_beta;
double temporal_biregion_angles_basal_gaussian_beta;
double temporal_biregion_angles_apical_gaussian_phi;
double temporal_biregion_angles_basal_gaussian_phi;

//"arrest_front" parameters
double y_arrest_front;
double t_arrest_front;
double k_arrest_front;

//"simga_O_spatial" parameters
double sigma_O_spatial_max;
double sigma_O_spatial_min;
double sigma_O_spatial_y_boundary;

//"temporal_angle_bias_Gaussian"
double bias_sigma_initial;
double bias_sigma_terminal;
double bias_mu_initial;
double bias_mu_terminal;
int bias_t_init;
int bias_t_term;

//meristem_position_relative_constant
double meristem_position_relative_constant_y;

//list of modes used
string major_mode;
string minor_mode;

string division_control;
string in_division_direction;
string epi_division_direction;
string mechanics_mode;

string temporal_mode;
string division_mode;


double F_modifier;

using namespace std;

double L_std;

namespace parameter{
void read_mode(string modeFile){
    ifstream fin(modeFile, ios::in);
    if(!fin.is_open()){
        {
            std::cout<<"Error: missing mode.txt ("<<modeFile<<")"<<endl;
            modeFile = "../"+modeFile;
            ifstream fin(modeFile, ios::in);
        }
        
    }
    if(!fin.is_open()){
        std::cout<<"Error: missing mode.txt ("<<modeFile<<")"<<endl;
        exit(1);
    }
    string mode_name;
    fin>>mode_name;
    fin>>major_mode;
    fin>>mode_name;
    fin>>minor_mode;
    fin.close();
}

void read(string parameterFile){
    std::cout<<"**************| Reading parameter settings ("<<parameterFile<<") |**************"<<endl;
    ifstream fin(parameterFile, ios::in);
    if(!fin.is_open()){
        string parameterFile_1 = "../"+parameterFile;
        ifstream fin(parameterFile_1, ios::in);
        if(!fin.is_open()){
            string parameterFile_2 = "../"+parameterFile_1;
            ifstream fin(parameterFile_2, ios::in);
            if(!fin.is_open()){
                std::cout<<"Error: missing parameter.txt ("<<parameterFile<<","<<parameterFile_1<<","<<parameterFile_2<<")"<<endl;
                exit(1);
            }
        }
    }
    string division_mode_tmp;
    string mechanics_mode_tmp;
    fin>>division_mode_tmp;
    fin>>division_control;
    fin>>division_mode_tmp;
    fin>>in_division_direction;
    fin>>division_mode_tmp;
    fin>>epi_division_direction;
    fin>>mechanics_mode_tmp;
    fin>>mechanics_mode;

    std::cout<<"**************| Simulation Mode |**************"<<endl;
    std::cout<<"Division Frequency Control: "<<division_control<<endl;
    std::cout<<"Inner Cell Division Direction: "<<in_division_direction<<endl;
    std::cout<<"Epidermal Cell Division Direction: "<<epi_division_direction<<endl;
    std::cout<<"Mechanics Mode: "<<mechanics_mode<<endl;
    std::cout<<"**************| General Parameters |**************"<<endl;
    vector<string> parameter_name;
    vector<double> parameter_value;

    //string parameter_name_tmp;
    //fin>>parameter_name_tmp;
    //std::cout<<parameter_name_tmp<<endl;
    //std::cout<<"Reading parameter.txt"<<endl;
    while(fin){
        string parameter_name_tmp;
        double parameter_value_tmp;
        fin>>parameter_name_tmp;
        fin>>parameter_value_tmp;
        parameter_name.push_back(parameter_name_tmp);
        parameter_value.push_back(parameter_value_tmp);
        //std::cout<<parameter_name_tmp<<" "<<parameter_value_tmp<<endl;
    }
    std::cout<<"End of reading parameter.txt"<<endl;

    for(int i=0;i<parameter_name.size();i++){
        //general parameters
        if(parameter_name[i]=="sigma_L"){
            sigma_L = parameter_value[i];
        }
        else if(parameter_name[i]=="sigma_O"){
            sigma_O = parameter_value[i];
        }
        else if(parameter_name[i]=="kappa_S"){
            kappa_S = parameter_value[i];
        }
        else if(parameter_name[i]=="S_std"){
            S_std = parameter_value[i];
        }
        else if(parameter_name[i]=="delta_time"){
            delta_time = parameter_value[i];
        }
        else if(parameter_name[i]=="eta"){
            eta = parameter_value[i];
        }
        else if(parameter_name[i]=="standard_cell_period_length"){
            standard_cell_period_length = parameter_value[i];
        }
        else if(parameter_name[i]=="F_apparent"){
            F_apparent = parameter_value[i];
        }
        else if(parameter_name[i]=="T_division_check"){
            T_division_check = parameter_value[i];
        }
        else if(parameter_name[i]=="T_vtkoutput"){
            T_vtkoutput = parameter_value[i];
        }
        else if(parameter_name[i]=="T_cell_time"){
            T_cell_time = parameter_value[i];
        }
        else if(parameter_name[i]=="F_apparent"){
            F_apparent = parameter_value[i];
        }
        else if(parameter_name[i]=="end_cell_number"){
            end_cell_number = parameter_value[i];
        }
        else if(parameter_name[i]=="step_end"){
            step_end = parameter_value[i];
        }
    }

//********************* Division Frequency Control: Parameter Reading *********************
    std::cout<<"**************| Division Control Parameters |**************"<<endl;
    if(division_control!="no"&&division_control!="balance"&&division_control!="Gaussian_balance"){
        bool is_area_control_lower_limit_set = 0;
        bool is_area_control_slope_set = 0;

        for(int i=0;i<parameter_name.size();i++){
            if(parameter_name[i]=="lower_limit"){
                area_control_lower_limit=parameter_value[i];
                is_area_control_lower_limit_set = 1;
            }
            if(parameter_name[i]=="slope"){
                area_control_slope=parameter_value[i];
                is_area_control_slope_set = 1;
            }
        }
        if(is_area_control_lower_limit_set==0){
            std::cout<<"Error: area control parameters unset (area_control_lower_limit)"<<endl;
            exit(1);
        }
        if(is_area_control_slope_set==0){
            std::cout<<"Error: area control parameters unset (area control slope)"<<endl;
            exit(1);
        }
    }

    if(division_control=="balance"||division_control=="Gaussian_balance"){
        bool is_BF_set = 0;
            for(int i=0;i<parameter_name.size();i++){
                if(parameter_name[i]=="BF"){
                    BF = parameter_value[i];
                    is_BF_set = 1;
                }
            }
            if(is_BF_set==0){
                std::cout<<"Error: balance control parameter unset (BF)"<<endl;
                exit(1);
            }
    }

    if(division_control=="Gaussian_area"||division_control=="Gaussian_balance"||division_control=="uniform_to_Gaussian"){
        bool is_gau_mu_set = 0;
        bool is_gau_sigma_set = 0;

        for(int i=0;i<parameter_name.size();i++){
            if(parameter_name[i]=="gau_mu"){
                gau_mu = parameter_value[i];
                is_gau_mu_set = 1;
            }
            else if(parameter_name[i]=="gau_sigma"){
                gau_sigma = parameter_value[i];
                is_gau_sigma_set = 1;
            }
        }

        if(is_gau_mu_set==0){
            std::cout<<"Error: Gaussian division frequency control parameter unset (gau_mu)"<<endl;
            exit(1);
        }
        if(is_gau_sigma_set==0){
            std::cout<<"Error: Gaussian division frequency control parameter unset (gau_sigma)"<<endl;
            exit(1);
        }

    }

    /*
    if(division_control==" "){
        bool is_set =0;
        for(int i=0;i<parameter_name.size();i++){
            if(parameter_name[i] == " "){
                 = parameter_value[i];
                is_set = 1;
            }
        }
        if(is_set==0){
            std::cout<<"Error: control parameter unset () "<<endl;
            exit(1);
        }
    }
    */

//********************* Division Direction Control: Parameter Reading *********************
    if(in_division_direction=="random"){

    }
    else if(in_division_direction=="mochizuki_bias"){
        bool is_mochizuki_bias_beta_set = 0;
        bool is_mochizuki_bias_phi_set = 0;
        for(int i=0;i<parameter_name.size();i++){
            if(parameter_name[i] == "mochizuki_bias_beta"){
                mochizuki_bias_beta = parameter_value[i];
                is_mochizuki_bias_beta_set = 1;
            }
            else if(parameter_name[i] == "mochizuki_bias_phi"){
                mochizuki_bias_phi = parameter_value[i];
                mochizuki_bias_phi = mochizuki_bias_phi * M_PI;
                is_mochizuki_bias_phi_set = 1;
            }
        }

        if(is_mochizuki_bias_beta_set==0){
            std::cout<<"Error: mochizuki_bias control parameter unset (mochizuki_bias_beta)"<<endl;
            exit(1);
        }
        if(is_mochizuki_bias_phi_set==0){
            std::cout<<"Error: mochizuki_bias control parameter unset (mochizuki_phi_beta)"<<endl;
            exit(1);
        }
    }

    else{
            std::cout<<"Fatal Error: No inner cell division angle control is selected"<<endl;
            std::cout<<"Current in_division_direction is "<<in_division_direction<<endl;
            exit(1);
        }
   /*
    std::cout<<"Current parameter list"<<endl;
    for(int i=0;i<parameter_name.size();i++){
        std::cout<<parameter_name[i]<<": "<<parameter_value[i]<<endl;
    }
    */
    //std::cout<<"*****************************| End of reading parameter settings |**************************"<<endl;
    //exit(1);

//********************* Mechanics Mode: Parameter Reading *********************

    if(mechanics_mode=="simple"){

    }
    else if(mechanics_mode=="L_std"){
        bool is_L_std_set =0;
        for(int i=0;i<parameter_name.size();i++){
            if(parameter_name[i] == "L_std"){
                L_std = parameter_value[i];
                is_L_std_set = 1;
            }
        }
        if(is_L_std_set==0){
            std::cout<<"Error: Mechanics mode L_std parameter unset (L_std) "<<endl;
            exit(1);
        }
    }
    else{
        std::cout<<"Fatal Error: No mechanics mode is selected"<<endl;
        std::cout<<"Current mechanics mode is "<<mechanics_mode<<endl;
        exit(1);
    }
    

}

void record(){
    //print modes and parameter on terminal
    std::cout<<"Division control mode : "<<division_control<<endl; 
    std::cout<<"in_division_direction mode : "<<in_division_direction<<endl;
    std::cout<<"epi_division_direction mode : "<<epi_division_direction<<endl;
    std::cout<<"sigma_L"<<" : "<<sigma_L<<"   ";
    std::cout<<"sigma_O"<<" : "<<sigma_O<<"   ";
    
//********************* Division Frequency Control: Parameter std::cout *********************
    if(division_control=="balance"||division_control =="Gaussian_balance"){
        std::cout<<"BF"<<" : "<<BF<<"    ";
    }
    
    else if(division_control=="no"){

    }
    
    else{
        std::cout<<"area_control_lower_limit"<<" : "<<area_control_lower_limit<<"   ";
        std::cout<<"area_control_slope"<<" : "<<area_control_slope<<"   ";
    }

    if(division_control=="Gaussian_area"||division_control =="Gaussian_balance"||division_control=="uniform_to_Gaussian"){
        std::cout<<"gau_mu"<<" : "<<gau_mu<<"   ";
        std::cout<<"gau_sigma"<<" : "<<gau_sigma<<"   ";
    }
    

    std::cout<<endl;

//********************* Division Direction Control: Parameter std::cout *********************
    if(in_division_direction=="mochizuki_bias"){
        std::cout<<"mochizuki_bias_beta"<<" : "<<mochizuki_bias_beta<<"    ";
        std::cout<<"mochizuki_bias_phi"<<" : "<<mochizuki_bias_phi<<" ";
    }

    std::cout<<endl;
//********************* Mechanics Mode: Parameter std::cout *********************
//"simga_O_spatial" parameters

    //record parameters
    ofstream fout(parameterRecordFile);
    fout<<"Information of modes and parameters used in this simulation"<<endl;
    //print mode information on terminal
    fout<<"in_division_direction : "<<in_division_direction<<"    ";
    fout<<"epi_division_direction : "<<epi_division_direction<<"    ";
    fout<<"division_control : "<<division_control<<"    ";
    //fout<<"debug_force_mode : "<<debug_force_mode<<"    ";
    //fout<<"potential_energy_equation : "<<potential_energy_equation<<"    ";
    //fout<<"repetition_mode : "<<repetition_mode<<"    ";
    fout<<endl;
    fout<<endl;
    fout<<"General parameter information: "<<endl;

    fout<<"sigma_L"<<" : "<<sigma_L<<"   ";
    fout<<"sigma_O"<<" : "<<sigma_O<<"   ";
    fout<<"kappa_S"<<" : "<<kappa_S<<"   ";
    fout<<"S_std"<<" : "<<S_std<<"   ";
    fout<<"delta_time"<<" : "<<delta_time<<"   ";
    fout<<"eta"<<" : "<<eta<<"   ";
    fout<<"standard_cell_period_length"<<" : "<<standard_cell_period_length<<"   ";
    fout<<"F_apparent"<<" : "<<F_apparent<<"   ";
    fout<<"T_division_check"<<" : "<<T_division_check<<"   ";
    fout<<"T_vtkoutput"<<" : "<<T_vtkoutput<<"   ";
    fout<<"T_cell_time"<<" : "<<T_cell_time<<"   ";
    fout<<"end_cell_number"<<" : "<<end_cell_number<<"   ";
    fout<<"step_end"<<" : "<<step_end<<"   ";     

    fout<<endl;
    fout<<endl;
    fout<<"Cell division frequency and cell division direction information: "<<endl;

//********************* Division Frequency Control: Parameter Record *********************
    if(division_control=="balance"||division_control=="Gaussian_balance"){
        fout<<"BF"<<" : "<<BF<<endl;
    }
    
    else if(division_control=="no"){

    }
    
    else{
        fout<<"area_control_lower_limit"<<" : "<<area_control_lower_limit<<endl;
        fout<<"area_control_slope"<<" : "<<area_control_slope<<endl;
    }

    if(division_control=="Gaussian_area"||division_control=="Gaussian_balance"){
        fout<<"gau_mu"<<" : "<<gau_mu<<endl;
        fout<<"gau_sigma"<<" : "<<gau_sigma<<endl;
    }
    

    fout<<endl;

//********************* Division Direction Control: Parameter Record *********************
    if(in_division_direction=="mochizuki_bias"){
        fout<<"mochizuki_bias_beta"<<" : "<<mochizuki_bias_beta<<endl;
        fout<<"mochizuki_bias_phi"<<" : "<<mochizuki_bias_phi<<endl;
    }
    

    /*
    else if(in_division_direction==" "){
        fout<<" "<<" : "<< <<endl;
    }
    */
    fout<<endl;
    fout<<"End of parameter record."<<endl;
    std::cout<<"********************************************************************"<<endl;

}

void batchRead(parameterList* pl){
    std::cout<<"Information of parameter list: "<<endl;
    std::ifstream fin(parameterListFile, std::ios::in);
    if(!fin.is_open()){
        std::cout<<"Error: missing batchParameter.txt"<<endl;
        exit(1);
    }
    while(!fin.eof()){
        double parameter_1_tmp, parameter_2_tmp;
        fin>>parameter_1_tmp;
        fin>>parameter_2_tmp;
        pl->parameter_1.push_back(parameter_1_tmp);
        pl->parameter_2.push_back(parameter_2_tmp);
    }
    std::cout<<"parameter List"<<endl;
    for(int i=0;i<pl->parameter_1.size();i++){
        std::cout<<"parameter_1 "<<pl->parameter_1[i]<<" parameter_2 "<<pl->parameter_2[i]<<endl;
    }
}

}

namespace termination{
    bool checkEnd(Organ *p_g){
        bool checkEnd_tmp =0;
        if(division_control=="arrest_front"){
                double arrest_front_position_tmp;
            if((p_g->p_c.size()-p_g->initial_cell_number)<t_arrest_front){
                arrest_front_position_tmp = y_arrest_front;
            }
            else{
                arrest_front_position_tmp = y_arrest_front+k_arrest_front*(p_g->p_c.size()-p_g->initial_cell_number-t_arrest_front);
            }

            if(arrest_front_position_tmp<=0){
                std::cout<<"The arrest front reached the base. There will be no more cell division. The simulation will go to end"<<endl;
                //force::forceShapeInitiation(p_g,200000);
                output::VTK(p_g);
                checkEnd_tmp=1;
            }
        }
        
        if(p_g->p_c.size()>=end_cell_number){
            checkend_tmp ++;
            output::VTK(p_g);
            if(checkend_tmp>10){
                std::cout<<"End of the Simulation"<<std::endl;
                output::VTK(p_g);
                checkEnd_tmp = 1;
            }
            
        }
        return checkEnd_tmp;
    }
}