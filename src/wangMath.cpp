/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#define _USE_MATH_DEFINES
#include "../include/wangMath.h"

using namespace std;
namespace wangMath{
    double mochizuki_bias_sampling_pdf(double x, double bias_beta, double bias_phi){
        return 1+bias_beta*cos(2*x-bias_phi);
    }  

    //requires the sampling_pdf
    double mochizuki_bias_single_random_sampling(double bias_beta, double bias_phi){
        double x_range_min=0;
        double x_range_max=M_PI;
        double y_range_max=1+bias_beta;
        std::random_device rnd;
        std::mt19937 mt(rnd());
        std::uniform_real_distribution<> rand_x(x_range_min,x_range_max);
        std::uniform_real_distribution<> rand_y(0.0,y_range_max);
        double x_tmp=0,y_tmp=0,loop_index=0;
        do{
            x_tmp=rand_x(mt);
            y_tmp=rand_y(mt);
            if(y_tmp<mochizuki_bias_sampling_pdf(x_tmp, bias_beta, bias_phi))
                loop_index=1;
        }
        while(loop_index==0);
        return x_tmp;   
    }

    //using quadratic formula to calculate quadratic equation
    Quadratic_Solution Quadratic_Equation_Solve(double A_tmp,double B_tmp,double C_tmp)
    {
        //quadratic equation: 𝐴𝑥^2+𝐵𝑥+𝐶=0
        //The solution is 𝑥_1,2=(−𝐵±√(𝐵^2−4𝐴𝐶))/2𝐴, 𝑖𝑓 √(𝐵^2−4𝐴𝐶)≥0
        Quadratic_Solution Result_Quadratic;
        Result_Quadratic.delta = B_tmp*B_tmp - 4*A_tmp*C_tmp;

        if(Result_Quadratic.delta>0){
            Result_Quadratic.x1 = (-B_tmp+sqrt(Result_Quadratic.delta))/(2*A_tmp);
            Result_Quadratic.x2 = (-B_tmp-sqrt(Result_Quadratic.delta))/(2*A_tmp);
        }
        else if(Result_Quadratic.delta==0){
            Result_Quadratic.x1 = -B_tmp/(2*A_tmp);
        }
        else{

        }

        return Result_Quadratic;
    }

}