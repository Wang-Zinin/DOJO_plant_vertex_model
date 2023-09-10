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

    
    void Vertex::print_Cartesian(void){
        std::cout<<"("<<loc.x<<","<<loc.y<<")"<<std::endl;
    }

    void Vertex::print_Polar(void){
        std::cout<<"("<<r<<","<<theta<<")"<<std::endl;
    }

    double Vertex::distance_from_vertex(Vertex v2){
        return sqrt((loc.x-v2.loc.x)*(loc.x-v2.loc.x)+(loc.y-v2.loc.y)*(loc.y-v2.loc.y));
    }

    void Vertex::Cartesian_to_Polar(void){
        r = sqrt(loc.x*loc.x+loc.y*loc.y);
        theta = atan(loc.y/loc.x); 
    }

    void Vertex::Polar_to_Cartesian(void){
        loc.x = r*cos(theta);
        loc.y = r*sin(theta);
    }

    bool Vertex::collinear_points(Vertex v2,Vertex v3){
        //if the slopes of l12 and l13 are the same and so does the intercepts of l12 and l13, then point p1,p2,p3 is collinear, otherwise not collinear
        double slope_12 = (loc.y-v2.loc.y)/(loc.x-v2.loc.x);
        double slope_23 = (v2.loc.y-v3.loc.y)/(v2.loc.x-v3.loc.x);

        double intercept_12 = (loc.x*v2.loc.y-v2.loc.x*loc.y)/(loc.x-v2.loc.x);
        double intercept_23 = (v2.loc.x*v3.loc.y-v3.loc.x*v2.loc.y)/(v2.loc.x-v3.loc.x);

        if(abs(slope_12-slope_23)<EPS_geo&&abs(intercept_12-intercept_23)<EPS_geo)
        {
            return 1;
        }
        else{
            return 0;
        }
    }
    
    bool Vertex::same_vertex(Vertex v2){
        if(abs(distance_from_vertex(v2))<EPS_geo){
            return 1;
        }
        else{
            return 0;
        }
    }