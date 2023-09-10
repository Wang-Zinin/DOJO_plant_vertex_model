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


double Line::calc_length(Organ p_g){
        length = p_g.p_v[vi[0]]->distance_from_vertex(*p_g.p_v[vi[1]]);
        return length;
    }

double Line::calc_length(void){
        length = d1.distance_from_vertex(d2);
        return length;
    }

void Line::calc_slope_intercept(Organ p_g){
        slope = (p_g.p_v[vi[0]]->loc.y-p_g.p_v[vi[1]]->loc.y)/(p_g.p_v[vi[0]]->loc.x-p_g.p_v[vi[1]]->loc.x);
        intercept = (p_g.p_v[vi[0]]->loc.x*p_g.p_v[vi[1]]->loc.y-p_g.p_v[vi[1]]->loc.x*p_g.p_v[vi[0]]->loc.y)/(p_g.p_v[vi[0]]->loc.x-p_g.p_v[vi[1]]->loc.x);
    }
    
void Line::calc_slope_intercept(void){
        slope = (d1.loc.y-d2.loc.y)/(d1.loc.x-d2.loc.x);
        intercept = d2.loc.y - slope*d2.loc.x;
    }

double Line::distance_from_point(Vertex p1){
        //reference 1: http://www.csharphelper.com/howtos/howto_point_segment_distance.html
        //reference 2: https://www.youtube.com/watch?v=egmZJU-1zPU

        Vertex closest;
        double dx = d2.loc.x-d1.loc.x;
        double dy = d2.loc.y-d1.loc.y;
        if((dx==0)&(dy==0)){
            //It's a point not a line segment
            closest = d1;
            dx = p1.loc.x - d1.loc.x;
            dy = p1.loc.y - d1.loc.y;
            return sqrt(dx*dx+dy*dy);
        }
        
        //Calculate the t that minimizes the distance
        double t = ((p1.loc.x-d1.loc.x)*dx + (p1.loc.y-d1.loc.y)*dy)/(dx*dx+dy*dy);
        
        //See if this represents one of the segment's endpoint or a point in the middle
        if(t<0){
            closest = d1;
        }
        else if(t>1){
            closest = d2;
        }
        else{
            closest.loc = _vec<double>{d1.loc.x+  t*dx, d1.loc.y + t*dy, 0.0}; 
            if(p1.same_vertex(closest)==1){
                return 0;
            }
        }

        return p1.distance_from_vertex(closest);
    }

    void Line::set_endpoints(Vertex v1, Vertex v2){
        d1 = v1;
        d2 = v2;
        calc_slope_intercept();
        calc_length();
    }
    void Line::set_endpoints(Vertex* v1, Vertex* v2){
        Vertex p1;
        Vertex p2;
        p1.loc = v1->loc;
        p2.loc = v2->loc;
        d1 = p1;
        d2 = p2;
        calc_slope_intercept();
        calc_length();
    }

    bool Line::same_segment(Line ls2){
        if(length==ls2.length){
            if((d1.same_vertex(ls2.d1)&&d2.same_vertex(ls2.d2))||(d2.same_vertex(ls2.d1)&&d1.same_vertex(ls2.d2))){
                return true;
            }
            else{
                return false;
            }
        }else{
            return false;
        }
    }

    void Line::print_slope_intercept(void){
        std::cout<<"y="<<slope<<"x+"<<intercept<<std::endl;
    }

    void Line::calc_general_ABC(Organ p_g){
        A = p_g.p_v[vi[1]]->loc.y - p_g.p_v[vi[0]]->loc.y;
        B = p_g.p_v[vi[0]]->loc.x - p_g.p_v[vi[1]]->loc.x;
        C = A*p_g.p_v[vi[0]]->loc.x + B*p_g.p_v[vi[0]]->loc.y;
    }