/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#include "../include/2dvInitial.h"


string real_organ_contour_imagej_txt = "../input/obovate_petal.txt";

namespace initialization{
    void organ(Organ* p_g){
        

        organ_geo::epidermal_identity(p_g);
        cout<<"similarity_calculation_required: "<<similarity_calculation_required<<endl;
        if(similarity_calculation_required==1){
            real_organ_contour_processed_for_similarity_index=boundary_geo::read_and_process_real_organ_contour_imagej(real_organ_contour_imagej_txt);
        }
        force::forceShapeInitiation(p_g,200000);
        division::cell_time_initialization(p_g);
        p_g->step=0;
        output::VTK(p_g);
        output::geo_initial();
        geo::calcGeometrics(p_g);
        output::geo(p_g);    
    }
}