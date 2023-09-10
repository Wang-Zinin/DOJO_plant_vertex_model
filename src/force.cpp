/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#include "../include/force.h"
#include "../include/2dvOrgan.h"
#include "../include/2dvCell.h"
#include "../include/2dvLine.h"
#include "../include/2dvVertex.h"

using namespace std;

namespace force{

//calculate edge elastic force and area elastic force for each vertices and execute the vertex motion

//F_i = -sigma_L * L_k (vector) / L_k (scalar)
void line_elastic_force(Organ* p_g){
    //calculate the elastic force for each edge
    
    //has two potential energy modes for elastic forces of edges: 1. simple; 2. L_std
    if(mechanics_mode=="simple"){
        for(int li=0;li<(int)p_g->p_l.size();li++){
            _vec<double> relative = p_g->p_v[p_g->p_l[li]->vi[0]]->loc-p_g->p_v[p_g->p_l[li]->vi[1]]->loc;
            double dist = relative.norm();
            p_g->p_l[li]->length = dist;
                
            _vec<double> frc;
            if(p_g->p_l[li]->IsOutermost==0){
                frc = (-1.0)*sigma_L* relative/dist;
            }
            else{
                frc = (-1.0)*sigma_O* relative/dist;
            }
                
            p_g->p_v[p_g->p_l[li]->vi[0]]->frc_edge += frc;
            p_g->p_v[p_g->p_l[li]->vi[1]]->frc_edge -= frc;
            p_g->p_l[li]->edgeForce = frc.norm();

        }
    }
    
    else if(mechanics_mode == "L_std"){
        for(int li=0;li<(int)p_g->p_l.size();li++){
            _vec<double> relative = p_g->p_v[p_g->p_l[li]->vi[0]]->loc-p_g->p_v[p_g->p_l[li]->vi[1]]->loc;
            double dist = relative.norm();
            p_g->p_l[li]->length = dist;
            _vec<double> frc;
            if(p_g->p_l[li]->IsOutermost==0){
                frc = (-1.0)*sigma_L* (dist -L_std)*relative/dist;
            }
            else{
                frc = (-1.0)*sigma_O* (dist -L_std)*relative/dist;
            }
                
            p_g->p_v[p_g->p_l[li]->vi[0]]->frc_edge += frc;
            p_g->p_v[p_g->p_l[li]->vi[1]]->frc_edge -= frc;
            p_g->p_l[li]->edgeForce = frc.norm();
        }
    }
    
    else{
        cout<<"Fatal Error: No mechanics mode is selected ! Current mechanics mode is "<<mechanics_mode<<endl;
        cout<<"End of simulation"<<endl;
        exit(1);
    }
}

//requires the assumption that vertex indices inside a cell should be arranged in anticlockwise direction
void cell_elastic_force(Organ *p_g){
    
    for(int ci=0;ci<(int)p_g->p_c.size();ci++){
        p_g->p_c[ci]->area = cell_geo::cell_area(p_g,ci);
        Cell *cp = p_g->p_c[ci];

        //calculate delta S/ delta x_i
        for (int j = 0; j < (int)cp->vi.size(); j++) {
        Vertex *vp[3];
        vp[1] = p_g->p_v[cp->vi[j]];

        if (j == 0) {
            vp[0] = p_g->p_v[cp->vi[cp->vi.size() - 1]];
        }
        else {
            vp[0] = p_g->p_v[cp->vi[j - 1]];
        }

        if (j == (int)cp->vi.size() - 1) {
            vp[2] = p_g->p_v[cp->vi[0]];
        }
        else {
            vp[2] = p_g->p_v[cp->vi[j + 1]];
        }

        _vec<double> s_grad = _vec<double>(0.0, 0.0, 0.0);
        s_grad.x =  (vp[2]->loc.y - vp[0]->loc.y);
        s_grad.y =  (vp[0]->loc.x - vp[2]->loc.x);

        _vec<double> frc_tmp = (-1.0) * kappa_S * (p_g->p_c[ci]->area - S_std) * s_grad;
          vp[1]->frc_area += frc_tmp;
        }
        
        
    }
}

void force_reset(Organ *p_g){
    //force reset to zero
    for(int vi=0; vi<(int)p_g->p_v.size();vi++){
        p_g->p_v[vi]->frc_edge = _vec<double>(0.0,0.0,0.0);
        p_g->p_v[vi]->frc_area = _vec<double>(0.0,0.0,0.0);
    }
}

void calcForceMotion(Organ* p_g){
    
    force::line_elastic_force(p_g);
    force::cell_elastic_force(p_g);
    
    //the elastic forces drive vertices to move
    for(int vi=0; vi<(int)p_g->p_v.size();vi++){
        p_g->p_v[vi]->loc +=(p_g->p_v[vi]->frc_area+p_g->p_v[vi]->frc_edge)*delta_time;
    }
    for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(p_g->p_c[ci]->area>10.0 || p_g->p_c[ci]->area<0.1){
                std::cout<<"The area of cell "<<ci<<" is abnormal: "<<p_g->p_c[ci]->area<<std::endl;
            }
    }
    
    force::force_reset(p_g);
    
}

void forceShapeInitiation(Organ *p_g, int initiation_time){
    for(int pre_step=0;pre_step<initiation_time;pre_step++){
        force::calcForceMotion(p_g);
    }
}

}