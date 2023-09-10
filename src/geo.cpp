/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#include "../include/geo2dv.h"
#include "../include/2dvOrgan.h"
#include "../include/2dvCell.h"
#include "../include/2dvLine.h"
#include "../include/2dvVertex.h"

bool similarity_calculation_required=0;
vector<Vertex*> real_organ_contour_processed_for_similarity_index;

namespace vertex_geo{
    //calculate the distance between vertex v1 and vertex v2
    double vertex_distance(Vertex* v1, Vertex* v2){
        double vertex_distance_tmp=0;
        vertex_distance_tmp = sqrt((v1->loc.y-v2->loc.y)*(v1->loc.y-v2->loc.y)+(v1->loc.x-v2->loc.x)*(v1->loc.x-v2->loc.x));
        return vertex_distance_tmp;
    }

    double vertex_distance(Vertex v1, Vertex v2){
        double vertex_distance_tmp=0;
        vertex_distance_tmp = sqrt((v1.loc.y-v2.loc.y)*(v1.loc.y-v2.loc.y)+(v1.loc.x-v2.loc.x)*(v1.loc.x-v2.loc.x));
        return vertex_distance_tmp;
    }

    //judge the relationship between vertex v1 and vertex v2; if return=0, they are different points; if return=1, they are so close that be seemed as the same point 
    bool vertex_relationship(Vertex* v1, Vertex* v2){
        bool index_tmp = 0;
        double epsilon_tmp =  10e-6;

        if(vertex_geo::vertex_distance(v1,v2)<epsilon_tmp){
            index_tmp = 1;
        }
        else{}

        return index_tmp;
    }

}

namespace line_geo{
    //calculate the length of line li
    double line_length(Organ* p_g, int li){
        p_g->p_l[li]->length = (p_g->p_v[p_g->p_l[li]->vi[0]]->loc - p_g->p_v[p_g->p_l[li]->vi[1]]->loc).norm();
        return p_g->p_l[li]->length;
    }


    //calculate the intersection of a line defined with slope and intercept, with cell wall li. if they are parallel _vec.z = 1.0; if they simply have no intersection _vec.z = -1.0
    _vec<double> line_cell_wall_intersection(double slope, double intercept, Organ* p_g, int li){
        _vec<double> intersection_result;
        //std::cout<<"li "<<li<<endl;
        double cell_wall_slope = (p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y-p_g->p_v[p_g->p_l[li]->vi[1]]->loc.y)/(p_g->p_v[p_g->p_l[li]->vi[0]]->loc.x-p_g->p_v[p_g->p_l[li]->vi[1]]->loc.x);
        double cell_wall_intercept = p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y-cell_wall_slope*p_g->p_v[p_g->p_l[li]->vi[0]]->loc.x;
        //std::cout<<"cell wall slope: "<<cell_wall_slope<<" , cell_wall_intercept: "<<cell_wall_intercept<<endl;
        
        //now we have the line slope, line intercept, cell wall slope and cell wall intercept: y=a1x+b1, y=a2x+b2
        // => there are three possibilities for two lines: they could be parallel, intersected, or identical (ignore identical situation)
        //parallel: a1=a2
        //std::cout<<slope<<" "<<cell_wall_slope<<endl;
        if(abs(slope-cell_wall_slope)<1e-9){
            intersection_result = _vec<double> {0.0,0.0,1.0};
            //std::cout<<"parallel"<<endl;
        }
        //intersected
        else{
                double intersection_x = (cell_wall_intercept - intercept)/(slope - cell_wall_slope);
                double intersection_y = cell_wall_slope*intersection_x + cell_wall_intercept;
                //std::cout<<"intersection_x "<<intersection_x<<", intersection_y "<<intersection_y<<endl;
                double segment_min_y = std::min(p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y,p_g->p_v[p_g->p_l[li]->vi[1]]->loc.y);
                double segment_max_y = std::max(p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y,p_g->p_v[p_g->p_l[li]->vi[1]]->loc.y); 
                double segment_min_x = std::min(p_g->p_v[p_g->p_l[li]->vi[0]]->loc.x,p_g->p_v[p_g->p_l[li]->vi[1]]->loc.x);
                double segment_max_x = std::max(p_g->p_v[p_g->p_l[li]->vi[0]]->loc.x,p_g->p_v[p_g->p_l[li]->vi[1]]->loc.x);
                
                intersection_result = _vec<double> {intersection_x,intersection_y,0.0};
                if(intersection_y <= segment_max_y && intersection_y >= segment_min_y &&intersection_x >= segment_min_x&&intersection_x <= segment_max_x){
                        //this intersection is correct
                        //std::cout<<"li "<<li<<"; endpoint 1 for cell wall "<<p_g->p_v[p_g->p_l[li]->vi[0]]->loc.x<<","<<p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y<<"; endpoint 2 for cell wall "<<p_g->p_v[p_g->p_l[li]->vi[1]]->loc.x<<","<<p_g->p_v[p_g->p_l[li]->vi[1]]->loc.y<<"; Line slope: "<<slope<<", line intercept "<<intercept<<";cell wall slope: "<<cell_wall_slope<<" , cell_wall_intercept: "<<cell_wall_intercept<<"; intersection_x "<<intersection_x<<", intersection_y "<<intersection_y<<endl;
                        //std::cout<<" "<<std::endl;
                }
                else{
                    intersection_result.z=-1.0;
                    //this intersection is out of the line segment
                }    
    }
    return intersection_result;
}

    Distance_point_line distance_line_segment_to_vertex(Line li, Vertex p1){
        //reference 1: http://www.csharphelper.com/howtos/howto_point_segment_distance.html
        //reference 2: https://www.youtube.com/watch?v=egmZJU-1zPU
        Distance_point_line Result;
        Vertex closest_point;
        Vertex d1 = li.d1;
        Vertex d2 = li.d2;
        double dx = d2.loc.x - d1.loc.x;
        double dy = d2.loc.y - d1.loc.y;
        if((dx==0)&&(dy==0)){
            //It's a point not a line segment
            closest_point = d1;
            dx = p1.loc.x - d1.loc.x;
            dy = p1.loc.y - d1.loc.y;
            Result.distance=sqrt(dx*dx+dy*dy);
            Result.Closest_Point=d1;
            return Result;
        }

        // Calculate the t that minimizes the distance
        double t = ((p1.loc.x - d1.loc.x)*dx + (p1.loc.y-d1.loc.y)*dy)/(dx*dx+dy*dy);
        Result.t =t;
        //std::cout<<"t "<<t<<endl;
        // See if this represents one of the segment's ends points or a point in the middle
        if(t<0){
            //std::cout<<"the closest point should be the end point d1"<<endl;
            closest_point = d1;
        }
        else if(t>1){
            //std::cout<<"the closest point should be the end point d2"<<endl;
            closest_point = d2;
        }
        else{
            //std::cout<<"the closest point is between d1 and d2"<<endl;
            closest_point.loc = _vec<double>{d1.loc.x + t*dx, d1.loc.y + t*dy,0.0};
            
        }
        //std::cout<<"Closest";
        //closest.print_Cartesian();
        Result.Closest_Point = closest_point;
        Result.distance = p1.distance_from_vertex(closest_point);
        return Result;
    }

    //calculate line segment distance to vertex
    Distance_point_line distance_line_segment_to_vertex(Organ* p_g,int li,Vertex p1){
        //reference 1: http://www.csharphelper.com/howtos/howto_point_segment_distance.html
        //reference 2: https://www.youtube.com/watch?v=egmZJU-1zPU
        Distance_point_line Result;
        
        Vertex closest_point;
        Vertex d1 = *p_g->p_v[p_g->p_l[li]->vi[0]];
        Vertex d2 = *p_g->p_v[p_g->p_l[li]->vi[1]];

        double dx = d2.loc.x - d1.loc.x;
        double dy = d2.loc.y - d1.loc.y;
        if((dx==0)&&(dy==0)){
            //It's a point not a line segment
            closest_point = d1;
            dx = p1.loc.x - d1.loc.x;
            dy = p1.loc.y - d1.loc.y;
            Result.distance=sqrt(dx*dx+dy*dy);
            Result.Closest_Point=d1;
            return Result;
        }

        // Calculate the t that minimizes the distance
        double t = ((p1.loc.x - d1.loc.x)*dx + (p1.loc.y-d1.loc.y)*dy)/(dx*dx+dy*dy);
        Result.t =t;
        //std::cout<<"t "<<t<<endl;
        // See if this represents one of the segment's ends points or a point in the middle
        if(t<0){
            //std::cout<<"the closest point should be the end point d1"<<endl;
            closest_point = d1;
        }
        else if(t>1){
            //std::cout<<"the closest point should be the end point d2"<<endl;
            closest_point = d2;
        }
        else{
            //std::cout<<"the closest point is between d1 and d2"<<endl;
            closest_point.loc = _vec<double>{d1.loc.x + t*dx, d1.loc.y + t*dy,0.0};
            
        }
        //std::cout<<"Closest";
        //closest.print_Cartesian();
        Result.Closest_Point = closest_point;
        Result.distance = p1.distance_from_vertex(closest_point);
        return Result;
    }
    
    //=0, the vertex is outside of the line; =1, the vertex is inside the line
    bool line_vertex_relationship(Line l1, Vertex p1){
        //check if (p1.loc.x,p1.loc.y) satisfy the l1 eqution y = slope*x+intercept 
        double y_p1 = p1.loc.x*l1.slope+l1.intercept;
        if(abs(y_p1-p1.loc.y)<EPS_geo){
            //the vertex is inside the line
            return true;
        }else{
            //the vertex is outside the line
            return false;
        }

    }

    //=0, the vertex is outside of the linesegment; =1, the vertex is inside the line segment
    bool line_segment_vertex_relationship(Line ls1, Vertex p1){
        //check if the relationship between Vertex p1 and Line 
        if(line_vertex_relationship(ls1,p1)==0){
            return false;
        }
        else{
            //the Vertex p1 is inside the line, but maybe outside of the line segment
            //check if the (p1.x,p1.y) is inside the two endpoints
            if(p1.loc.x >= min(ls1.d1.loc.x, ls1.d2.loc.x) && p1.loc.x <= max(ls1.d1.loc.x, ls1.d2.loc.x) &&
               p1.loc.y >= min(ls1.d1.loc.y, ls1.d2.loc.y) && p1.loc.y <= max(ls1.d1.loc.y, ls1.d2.loc.y)){
                return true;
               }else{
                return false;
               }
        }
    } 

    //=0, intersected; =1, parallel; =2 perpendicular; =3, identical
    pair<int,Vertex> lines_relationship(Line l1,Line l2)
    {   
        int relationship;
        Vertex intersection;
        //in the case of inifinely large slope x=1, and x=2
        if(l1.d1.loc.x==l1.d2.loc.x&&l2.d1.loc.x==l2.d2.loc.x){
            if(l1.d1.loc.x==l2.d1.loc.x){
                relationship=3;
            }else{
                relationship=1;
            }
        }
        else if(l1.d1.loc.x==l1.d2.loc.x){
            if(l2.d1.loc.y==l2.d2.loc.y){
                relationship=2;
                intersection.loc.x=l1.d1.loc.x;
                intersection.loc.y=l2.d1.loc.y;
            }else{
                relationship=0;
                l2.calc_slope_intercept();
                intersection.loc.x=l1.d1.loc.x;
                intersection.loc.y=intersection.loc.x*l2.slope+l2.intercept;
            }
        }
        else if(l2.d1.loc.x==l2.d2.loc.x){
            if(l1.d1.loc.y==l1.d2.loc.y){
                relationship=2;
                intersection.loc.x=l2.d1.loc.x;
                intersection.loc.y=l1.d1.loc.y;
            }else{
                relationship=0;
                l1.calc_slope_intercept();
                intersection.loc.x=l2.d1.loc.x;
                intersection.loc.y=intersection.loc.x*l1.slope+l1.intercept;
            }
        }
        else
        {
            l1.calc_slope_intercept();
            l2.calc_slope_intercept();
            if(l1.slope==l2.slope){
                if(l1.intercept==l2.intercept){
                    relationship=3;
                }
                else{
                    relationship=1;
                }
            }
            else if(l1.slope*l2.slope==-1){
                relationship = 2;
                intersection.loc.x = (l2.intercept-l1.intercept)/(l1.slope-l2.slope);
                intersection.loc.y = intersection.loc.x*l1.slope + l1.intercept;
            }
            else{
                relationship=0;
                intersection.loc.x = (l2.intercept-l1.intercept)/(l1.slope-l2.slope);
                intersection.loc.y = intersection.loc.x*l1.slope + l1.intercept;
            }
        }

        return make_pair(relationship,intersection);
    }

    //=0, intersected; =1, parallel; =2, perpendicular; =3, identical; =4, overlapping; =5, touching; =6, disjoint
    pair<int,Vertex> segments_relationship(Line ls1, Line ls2)
    {   
        int segments_rela;
        Vertex segments_intersection;
        pair<int, Vertex> lines_rela = line_geo::lines_relationship(ls1,ls2);
        
        //if the two lines are intersected
        if(lines_rela.first==0){
            //judge if the intersection still existed in the segments 
            if(line_segment_vertex_relationship(ls1,lines_rela.second)){
                //the intersection existed in the segments 
                if(lines_rela.second.same_vertex(ls1.d1)||lines_rela.second.same_vertex(ls1.d2)){
                    //the intersection is one of the endpoint => touching
                    segments_rela = 5;
                    segments_intersection = lines_rela.second;
                }else{
                    segments_rela = 0;
                    segments_intersection = lines_rela.second;
                }
            }else{
                //the intersection existed outside the segments
                segments_rela = 6;
            }
        } 
        else if(lines_rela.first==1){
        //if the two lines are parallel
            segments_rela = 1;
        }
        else if(lines_rela.first==2){
            //judge if the perpendicular intersection still existed in the segments 
            if(line_segment_vertex_relationship(ls1,lines_rela.second)){
                //the intersection existed in the segments 
                if(lines_rela.second.same_vertex(ls1.d1)||lines_rela.second.same_vertex(ls1.d2)){
                    //the intersection is one of the endpoint => touching
                    segments_rela = 5;
                    segments_intersection = lines_rela.second;
                }else{
                    //the perpendicular intersection still exists
                    segments_rela = 2;
                }
            }else{
                //the intersection existed outside the segments
                segments_rela = 6;
            }
        }
        else{
            //line segment is the same
            if(ls1.same_segment(ls2)){
                segments_rela = 3;
            }
            else if((ls1.d1.loc.x>max(ls2.d1.loc.x,ls2.d2.loc.x))&&(ls1.d2.loc.x>max(ls2.d1.loc.x,ls2.d2.loc.x))||
                   (ls1.d1.loc.x<min(ls2.d1.loc.x,ls2.d2.loc.x))&&(ls1.d2.loc.x<min(ls2.d1.loc.x,ls2.d2.loc.x))){
                segments_rela = 1;
            }
            else if(ls1.d1.same_vertex(ls2.d1)||ls1.d1.same_vertex(ls2.d2)){
                segments_rela = 5;
                segments_intersection = ls1.d1;
            }
            else if(ls1.d2.same_vertex(ls2.d1)||ls1.d2.same_vertex(ls2.d2)){
                segments_rela = 5;
                segments_intersection = ls1.d2;
            }
            else{
                segments_rela = 4;
            }

        }
   
        return make_pair(segments_rela,segments_intersection);
    }
}

namespace cell_geo{

    //calculate the center of cell ci in p_g
    _vec<double> cell_center(Organ* p_g, int ci){
        _vec<double> center_tmp = _vec<double>{0.0,0.0,0.0};
        for(int vi=0;vi<(int)p_g->p_c[ci]->vi.size();vi++){
            center_tmp = center_tmp + p_g->p_v[p_g->p_c[ci]->vi[vi]]->loc;
        }
        center_tmp = center_tmp/(double)p_g->p_c[ci]->vi.size();
        p_g->p_c[ci]->center = center_tmp;
        return center_tmp;
    }

    //calculate the area of cell ci in p_g
    //assumptions: vertex indices are sorted in anticlockwise direction
    double cell_area(Organ* p_g, int ci){
        double i_area =0;
        Cell *cp = p_g->p_c[ci];

        //i番目の細胞の面積を計算。点が反時計回りに格納されていることを前提にしている。
        for (int j = 0; j < (int)cp->vi.size(); j++) {
        Vertex *vp[2];
        vp[0] = p_g->p_v[cp->vi[j]];
        if (j != (int)cp->vi.size() - 1) {
            vp[1] = p_g->p_v[cp->vi[j + 1]];
        }
        else if (j == (int)cp->vi.size() - 1) {
            vp[1] = p_g->p_v[cp->vi[0]];
        }
        else {
            std::cout << "Bug.Area" << std::endl;
            exit(0);
        }
        i_area += 0.5 * (vp[0]->loc.x * vp[1]->loc.y - vp[1]->loc.x * vp[0]->loc.y);
        }
        p_g->p_c[ci]->area=i_area;
        return i_area;
    }

    //calculate the perimeter of cell ci in p_g
    double cell_perimeter(Organ* p_g, int ci){

        double perimeter_tmp=0;
        for(int li=0;li<p_g->p_c[ci]->li.size();li++){
            double length_tmp = (p_g->p_v[p_g->p_l[p_g->p_c[ci]->li[li]]->vi[0]]->loc-p_g->p_v[p_g->p_l[p_g->p_c[ci]->li[li]]->vi[1]]->loc).norm();
            perimeter_tmp += length_tmp;
        }
        p_g->p_c[ci]->perimeter = perimeter_tmp;
        //std::cout<<"Cell "<<ci<<" perimeter: "<<perimeter_tmp<<endl;
        return perimeter_tmp;
    }

    //calculate the regularity of cell ci in p_g
    //regularity defines the similarity between a defined polygon (eg., cell ci) and the regular polygon that has the same perimeter
    //regularity ranges from (0,1.0], when the polygon is a regular polygon, regularity = 1.0;
    //regularity = area (cell ci)/area (the regularity polygon that has the same perimeter) 
    //reference: https://doi.org/10.1016/j.cad.2012.07.012
    double cell_regularity(Organ* p_g, int ci){
        double area_P = p_g->p_c[ci]->area;
        double perimeter_P = p_g->p_c[ci]->perimeter;
        double area_R = perimeter_P*perimeter_P/(4.0*(double)p_g->p_c[ci]->li.size()*tan(3.145926/(double)p_g->p_c[ci]->li.size()));
        p_g->p_c[ci]->regularity = area_P/area_R;
        //std::cout<<area_R<<" "<<perimeter_P<<" "<<p_g->p_c[ci]->li.size()<<" "<<tan(3.145926/(double)p_g->p_c[ci]->li.size())<<endl;
        //std::cout<<"Cell "<<ci<<" regularity "<<p_g->p_c[ci]->regularity<<endl;
        return p_g->p_c[ci]->regularity;
    }

    //sort lines in counterclockwise direction, based on vertices in counterclockwise direction
    vector<int> cell_counterclock_line(Organ* p_g, int ci){
        std::vector<int> anticlockwise_lidx;
        Cell *cp = p_g->p_c[ci];
        int cp_vsize = cp->vi.size();
        for(int i = 0; i<cp_vsize;++i){
            for(int lidx: cp->li) {
            if(p_g->p_l[lidx]->vi[0] == cp->vi[i] && p_g->p_l[lidx]->vi[1]==cp->vi[(i+1)%cp_vsize]){
                anticlockwise_lidx.push_back(lidx);
                break;
            }else if(p_g->p_l[lidx]->vi[0] == cp->vi[(i+1)%cp_vsize] && p_g->p_l[lidx]->vi[1] == cp->vi[i]){
                std::swap(p_g->p_l[lidx]->vi[0], p_g->p_l[lidx]->vi[1]);
                anticlockwise_lidx.push_back(lidx);
            }
            }
        }
        return anticlockwise_lidx;
    }

}

namespace organ_geo{

    //calculate the length of all lines within an organ
    void organ_line_length(Organ* p_g){
        for(int li=0;li<(int)p_g->p_l.size();li++){
            line_geo::line_length(p_g,li);
        }
    }

    //calculate the perimeter of all cells within an organ
    double organ_cell_perimeter(Organ* p_g){
        double av_perimeter_tmp=0;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            av_perimeter_tmp +=cell_geo::cell_perimeter(p_g,ci);
        }
        av_perimeter_tmp = av_perimeter_tmp/(double)p_g->p_c.size();
        return av_perimeter_tmp;
    }

    //calculate the center of the organ: sum of the center of all cells
    _vec<double> organ_center(Organ* p_g){
        _vec<double> center_tmp = _vec<double>{0.0,0.0,0.0};
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            center_tmp += cell_geo::cell_center(p_g,ci);
        }
        p_g->center = center_tmp/(double)p_g->p_c.size();
        return center_tmp;
    }

    //calculate the area of the organ: sum of the area of all cells
    double organ_area(Organ* p_g){
        double area_tmp=0.0;
        double epiArea=0.0, inArea=0.0;
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            area_tmp+=cell_geo::cell_area(p_g,ci);
            if(p_g->p_c[ci]->IsEpidermal==1){
                epiArea+=p_g->p_c[ci]->area;
            }
            else{
                inArea+=p_g->p_c[ci]->area;
            }
        }
        p_g->area=area_tmp;
        p_g->area_averaged = p_g->area/(double)p_g->p_c.size();
        p_g->epiArea = epiArea;
        p_g->epiArea_averaged = epiArea/p_g->N_epi_cell;
        p_g->inArea = inArea;
        p_g->inArea_averaged = inArea/p_g->N_inner_cell; 
        return area_tmp;
    }

    //determine the epidermal identity of all cells and lines within an organ: 
    //cell->IsSurface=1, epidermal cell;=0, inner cell; Line->IsOutermost=1, epidermal line;=0, inner line; also count the epidermal cell number and inner cell number
    void epidermal_identity(Organ* p_g){

        vector<int> surface_line;
        vector<int> surface_vertex;

        for(int vi=0; vi<(int)p_g->p_v.size();vi++){
            p_g->p_v[vi]->occurrenceInCell=0;
        }


        //count the occurrence of each vertex in all cells
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            for(int vi=0;vi<(int)p_g->p_c[ci]->vi.size();vi++){
                int tmp = p_g->p_c[ci]->vi[vi];
                p_g->p_v[tmp]->occurrenceInCell++;
            }
        }

        //find surface vertex: if the occurrence of vertex is less than 3, then this vertex is a surface vertex
        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
            if(p_g->p_v[vi]->occurrenceInCell<3){
                p_g->p_v[vi]->IsSurface=1;
                surface_vertex.push_back(vi);
            }
            else{
                p_g->p_v[vi]->IsSurface=0;
            }
        }
        //Debug: check the occurrence of vertex in cell and judgement of surface vertex
        /*
        std::cout<<"Debug: check the occurrence of vertex in cell and judgement of surface vertex"<<std::endl;
        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
            std::cout<<vi<<" "<<p_g->p_v[vi]->occurrenceInCell<<" "<<p_g->p_v[vi]->IsSurface<<std::endl;
        }
        std::cout<<"End Debug: check the occurrence of vertex in cell and judgement of surface vertex"<<std::endl;
        */

    //find outermost edge: if the two vertices connecting the edge are surface vertices, then this edge is an outermost edge
    for(int li=0; li<(int)p_g->p_l.size();li++){
        if(p_g->p_v[p_g->p_l[li]->vi[0]]->IsSurface==1 && p_g->p_v[p_g->p_l[li]->vi[1]]->IsSurface==1){
            p_g->p_l[li]->IsOutermost=1;
            surface_line.push_back(li);
        }
        else{
            p_g->p_l[li]->IsOutermost=0;
        }
    }
    //Debug: check each line's IsOutermost parameter
    /*
    std::cout<<"Debug: check IsOutermost parameter for Line"<<std::endl;
    for(int li=0;li<(int)p_g->p_l.size();li++){
        std::cout<<li<<" "<<p_g->p_l[li]->IsOutermost<<std::endl;
    }
    std::cout<<"End Debug: check IsOutermost parameter for Line"<<std::endl;
        */

        //find epidermal cell: if a cell have an outermost edge, then it is an epidermal cell
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            bool IsEpidermal_tmp=0;
            for(int li=0;li<(int)p_g->p_c[ci]->li.size();li++){
                if(p_g->p_l[p_g->p_c[ci]->li[li]]->IsOutermost==1){
                    IsEpidermal_tmp=1;
                    p_g->p_c[ci]->surfaceVertex[0]=p_g->p_l[p_g->p_c[ci]->li[li]]->vi[0];
                    p_g->p_c[ci]->surfaceVertex[1]=p_g->p_l[p_g->p_c[ci]->li[li]]->vi[1];
                    p_g->p_c[ci]->outermostLength = p_g->p_l[p_g->p_c[ci]->li[li]]->length;
                }
            }
            p_g->p_c[ci]->IsEpidermal=IsEpidermal_tmp;
        }

        int inner_cell_number=0, peripheral_cell_number=0;

        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(p_g->p_c[ci]->IsEpidermal==0){
                inner_cell_number++;
            }
            else{
                peripheral_cell_number++;
            }
        }
        
        //std::cout<<"Currently, we have "<<inner_cell_number<<" inner cells, and "<<peripheral_cell_number<<" peripheral cells."<<std::endl;
        p_g->N_epi_cell = peripheral_cell_number;
        p_g->N_inner_cell = inner_cell_number;
        p_g->surface_line=surface_line;
        p_g->surface_vertex=surface_vertex;
    }

    //calculate the perimeter of the organ: sum of the length of all epidermal line
    double organ_perimeter(Organ* p_g){
        double perimeter_tmp=0.0;
        for(int li=0; li<(int)p_g->p_l.size();li++){
        if(p_g->p_l[li]->IsOutermost==1){
                perimeter_tmp += p_g->p_l[li]->length;
            }
        }
        p_g->perimeter = perimeter_tmp;
        p_g->perimeter_averaged = perimeter_tmp/(double)p_g->p_c.size();
        return perimeter_tmp;
    }

    //calculate the circularity of the organ: circularity quantifies the similarity between the defined polygon (eg., p_g) and a circle
    //circularity ranges from (0,1], when circularity=1, the polygon is a perfect circle; also could be refered as roundness
    //circularity = Perimeter^2/(4pi*Area)
    //reference: https://en.wikipedia.org/wiki/Roundness
    double organ_circularity(Organ* p_g){
        p_g->circularity = 4*3.1415926*p_g->area/(p_g->perimeter*p_g->perimeter);
        return p_g->circularity;
    }

    //calculate the maximum y of the organ
    double organ_maximum_y(Organ* p_g){
        double y_max_tmp=p_g->p_c[0]->center.y;

        for(int i=0; i<(int)p_g->p_c.size();i++){
            if(y_max_tmp<p_g->p_c[i]->center.y){
                y_max_tmp=p_g->p_c[i]->center.y;
            }
        }

        return y_max_tmp;
    }

    //calculate the minimum y of the organ
    double organ_minimum_y(Organ* p_g){
        double y_min_tmp=p_g->p_c[0]->center.y;

        for(int i=0; i<(int)p_g->p_c.size();i++){
            if(y_min_tmp>p_g->p_c[i]->center.y){
                y_min_tmp=p_g->p_c[i]->center.y;
            }
        }

        return y_min_tmp;
    }

    double organ_minimum_y_v(Organ* p_g){
        double y_min_tmp=p_g->p_v[0]->loc.y;

        for(int i=0; i<(int)p_g->p_v.size();i++){
            if(y_min_tmp>p_g->p_v[i]->loc.y){
                y_min_tmp=p_g->p_v[i]->loc.y;
            }
        }
        p_g->y_min_v=y_min_tmp;
        return y_min_tmp;
    }

    double organ_minimum_x_v(Organ* p_g){
        double x_min_tmp=p_g->p_v[0]->loc.x;

        for(int i=0; i<(int)p_g->p_v.size();i++){
            if(x_min_tmp>p_g->p_v[i]->loc.x){
                x_min_tmp=p_g->p_v[i]->loc.x;
            }
        }
        p_g->x_min_v=x_min_tmp;
        return x_min_tmp;
    }

    //order all epidermal vertices (boundary vertices)/ boundary lines into a anticlockwise vector
    Ordered_boundary organ_ordered_anticlockwise_boundary(Organ* p_g){
        //std::cout<<"Start ordering surface vertices"<<endl;
        Ordered_boundary Result;
        vector<Vertex> surface_vertices;
        vector<Vertex> anticlockwise_surface_vertices;
        vector<Line> surface_line;
        vector<int> anticlockwise_surface_lines;
        
    //0. collecting all surface lines into surface_line, all surface vertices into surface_vertices
        for(int li=0;li<(int)p_g->p_l.size();li++){
            if(p_g->p_l[li]->IsOutermost==1){
                surface_line.push_back(*p_g->p_l[li]);
            }
        }

        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                surface_vertices.push_back(*p_g->p_v[vi]);
            }
        }

    //1. the first surface vertex: the surface vertex with lowest y position
        Vertex first_boundary_vertex;
        first_boundary_vertex = surface_vertices[0];

        for(int vi=0;vi<(int)surface_vertices.size();vi++){
            if(first_boundary_vertex.loc.y>surface_vertices[vi].loc.y){
                first_boundary_vertex=surface_vertices[vi];
            }
        }

        anticlockwise_surface_vertices.push_back(first_boundary_vertex);

        //std::cout<<" The first boundary vertices is "<<anticlockwise_surface_vertices[0]->vi<<", and its (x,y) is ("<<anticlockwise_surface_vertices[0]->loc.x<<","<<anticlockwise_surface_vertices[0]->loc.y<<")."<<endl;

    //2. the second boundary vertex: the surface vertex that is the neighbor of the first surface vertex (they share the same lines)
        //                              and this surface vertex should have smaller x values to make the direction anticlockwise
        Vertex v_tmp1;
        Vertex v_tmp2;
        int li_tmp1,li_tmp2;
        bool v_tmp1_found=0;
            for(int li=0; li<(int)surface_line.size(); li++){
                //std::cout<<"li "<<surface_line[li]->li<<" connecting vertex: "<<surface_line[li]->vi[0]<<","<<surface_line[li]->vi[1]<<"; size of surface line"<<surface_line.size()<<endl;
                if(surface_line[li].vi[0]==anticlockwise_surface_vertices[0].vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[surface_line[li]->vi[1]]->loc.x: "<<p_g->p_v[surface_line[li]->vi[1]]->loc.x<<endl;
                    v_tmp1 = *p_g->p_v[surface_line[li].vi[1]];
                    li_tmp1=li;
                    v_tmp1_found=1;
                    //if(p_g->p_v[surface_line[li]->vi[1]]->loc.x<anticlockwise_surface_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    anticlockwise_surface_vertices.push_back(p_g->p_v[surface_line[li]->vi[1]]);
                    //}
                }
                else if(surface_line[li].vi[0]==anticlockwise_surface_vertices[0].vi&&v_tmp1_found==1){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[surface_line[li]->vi[1]]->loc.x: "<<p_g->p_v[surface_line[li]->vi[1]]->loc.x<<endl;
                    v_tmp2 = *p_g->p_v[surface_line[li].vi[1]];
                    li_tmp2=li;
                    break;
                }
                else if(surface_line[li].vi[1]==anticlockwise_surface_vertices[0].vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[surface_line[li]->vi[0]]->loc.x: "<<p_g->p_v[surface_line[li]->vi[0]]->loc.x<<endl;
                    v_tmp1 = *p_g->p_v[surface_line[li].vi[0]];
                    li_tmp1=li;
                    v_tmp1_found=1;
                    //if(p_g->p_v[surface_line[li]->vi[0]]->loc.x<anticlockwise_surface_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    anticlockwise_surface_vertices.push_back(p_g->p_v[surface_line[li]->vi[0]]);
                    //}
                }
                else if(surface_line[li].vi[1]==anticlockwise_surface_vertices[0].vi&&v_tmp1_found==1){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[surface_line[li]->vi[0]]->loc.x: "<<p_g->p_v[surface_line[li]->vi[0]]->loc.x<<endl;
                    v_tmp2 = *p_g->p_v[surface_line[li].vi[0]];
                    li_tmp2=li;
                    break;
                    //if(p_g->p_v[surface_line[li]->vi[0]]->loc.x<anticlockwise_surface_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    anticlockwise_surface_vertices.push_back(p_g->p_v[surface_line[li]->vi[0]]);
                    //}
                }
                else{

                }
            }
            if(v_tmp1.loc.x < v_tmp2.loc.x){
                anticlockwise_surface_vertices.push_back(v_tmp1);
                anticlockwise_surface_lines.push_back(surface_line[li_tmp1].li);
                //std::cout<<"v_tmp1 "<<v_tmp1<<endl;
            }
            else{
                anticlockwise_surface_vertices.push_back(v_tmp2);
                anticlockwise_surface_lines.push_back(surface_line[li_tmp2].li);
                //std::cout<<"v_tmp2 "<<v_tmp2<<endl;
            }
        //std::cout<<" The secondary boundary vertices is "<<anticlockwise_surface_vertices[1]->vi<<", and its (x,y) is ("<<anticlockwise_surface_vertices[1]->loc.x<<","<<anticlockwise_surface_vertices[1]->loc.y<<")."<<endl;

    //3. the tertiary and following boundary vertices could be retrieved using incremental algorithm: 
        for(int vi=1; vi<(int)surface_line.size()-1;vi++){
            for(int li=0;li<(int)surface_line.size();li++){
                if(surface_line[li].vi[0]==anticlockwise_surface_vertices[vi].vi&&surface_line[li].vi[1]!=anticlockwise_surface_vertices[vi-1].vi){
                    anticlockwise_surface_vertices.push_back(*p_g->p_v[surface_line[li].vi[1]]);
                    anticlockwise_surface_lines.push_back(surface_line[li].li);
                }
                else if(surface_line[li].vi[1]==anticlockwise_surface_vertices[vi].vi&&surface_line[li].vi[0]!=anticlockwise_surface_vertices[vi-1].vi){
                    anticlockwise_surface_vertices.push_back(*p_g->p_v[surface_line[li].vi[0]]);
                    anticlockwise_surface_lines.push_back(surface_line[li].li);
                }
            }
        }
        //std::cout<<" The third boundary vertices is "<<anticlockwise_surface_vertices[2]->vi<<", and its (x,y) is ("<<anticlockwise_surface_vertices[2]->loc.x<<","<<anticlockwise_surface_vertices[2]->loc.y<<")."<<endl;
        //cout_fout_debug::cout_vector_vertex(anticlockwise_surface_vertices);
        Result.li = anticlockwise_surface_lines;
        Result.vi = anticlockwise_surface_vertices;
        return Result;
    }

    //find boundary points with equal distance on a close contour along polygon
    vector<Vertex*> organ_boundary_points_along_polygon(Organ* p_g, vector<Vertex*> anticlockwise_surface_vertices, int boundary_points_number){
        vector<Vertex*> boundary_points;
        double boundary_points_average_distance = p_g->perimeter/(double)boundary_points_number;

        //now we have the ordered boundary vertices: the boundary vertex with lowest y is ranked as the first; the other boundary vertices are ordered anticlockwise
    //1. the first boundary point: just pick the first boundary vertex (with lowest y)
        boundary_points.push_back(anticlockwise_surface_vertices[0]);
        anticlockwise_surface_vertices.push_back(anticlockwise_surface_vertices[0]); //link the last boundary vertices with the first boundary vertices
        
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    //2. the i-th boundary point: incremental algorithm

        double required_boundary_points_distance = boundary_points_average_distance;
        int vj = 0; //the current boundary_vertex index
        bool distance_intial_is_boundary_vertex=0;
        double surface_vertex_distance;

    for(int pi=0; pi<boundary_points_number-1; pi++){
        
            boundary_points_distance_larger_than_surface_vertex_distance:;
            if(distance_intial_is_boundary_vertex==0) {
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1]->loc-boundary_points[pi]->loc).norm();
            }
            else{
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1]->loc-anticlockwise_surface_vertices[vj]->loc).norm();
        }


        //judgement: if boundary_point_average_distance >=< surface_vertex_distance
        if(required_boundary_points_distance < surface_vertex_distance){
        // the new boundary point will be generated on the line segment between surface vertex i and surface vertex i+1
            Vertex* boundary_points_new = new Vertex;
            if(distance_intial_is_boundary_vertex==0) {
                boundary_points_new->loc = boundary_points[pi]->loc+(anticlockwise_surface_vertices[vj+1]->loc-boundary_points[pi]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }
            else{
                boundary_points_new->loc = anticlockwise_surface_vertices[vj]->loc+(anticlockwise_surface_vertices[vj+1]->loc-anticlockwise_surface_vertices[vj]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }

            boundary_points.push_back(boundary_points_new);
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<","<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance == surface_vertex_distance){
        // the new boundary point will be exactly the suface vertex i+1
            boundary_points.push_back(anticlockwise_surface_vertices[vj+1]);
            vj++;
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance > surface_vertex_distance){
        // the new boundary point might be generated on the line segment between surface vertex i+1 and surface vertex i+2
            required_boundary_points_distance -= surface_vertex_distance;
            vj++;
            distance_intial_is_boundary_vertex=1;
            goto boundary_points_distance_larger_than_surface_vertex_distance;
        }
        
        required_boundary_points_distance = boundary_points_average_distance;
        distance_intial_is_boundary_vertex=0;
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    }

        //cout_fout_debug::fout_vector_vertex(boundary_points,"50_boundary_points.txt");
        return boundary_points;

    }

    vector<Vertex> organ_boundary_points_along_polygon(Organ* p_g, vector<Vertex> anticlockwise_surface_vertices, int boundary_points_number){
        vector<Vertex> boundary_points;
        double boundary_points_average_distance = p_g->perimeter/(double)boundary_points_number;

        //now we have the ordered boundary vertices: the boundary vertex with lowest y is ranked as the first; the other boundary vertices are ordered anticlockwise
    //1. the first boundary point: just pick the first boundary vertex (with lowest y)
        boundary_points.push_back(anticlockwise_surface_vertices[0]);
        anticlockwise_surface_vertices.push_back(anticlockwise_surface_vertices[0]); //link the last boundary vertices with the first boundary vertices
        
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    //2. the i-th boundary point: incremental algorithm

        double required_boundary_points_distance = boundary_points_average_distance;
        int vj = 0; //the current boundary_vertex index
        bool distance_intial_is_boundary_vertex=0;
        double surface_vertex_distance;

    for(int pi=0; pi<boundary_points_number-1; pi++){
        
            boundary_points_distance_larger_than_surface_vertex_distance:;
            if(distance_intial_is_boundary_vertex==0) {
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1].loc-boundary_points[pi].loc).norm();
            }
            else{
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1].loc-anticlockwise_surface_vertices[vj].loc).norm();
        }


        //judgement: if boundary_point_average_distance >=< surface_vertex_distance
        if(required_boundary_points_distance < surface_vertex_distance){
        // the new boundary point will be generated on the line segment between surface vertex i and surface vertex i+1
            Vertex boundary_points_new;
            if(distance_intial_is_boundary_vertex==0) {
                boundary_points_new.loc = boundary_points[pi].loc+(anticlockwise_surface_vertices[vj+1].loc-boundary_points[pi].loc)*required_boundary_points_distance/surface_vertex_distance;
            }
            else{
                boundary_points_new.loc = anticlockwise_surface_vertices[vj].loc+(anticlockwise_surface_vertices[vj+1].loc-anticlockwise_surface_vertices[vj].loc)*required_boundary_points_distance/surface_vertex_distance;
            }

            boundary_points.push_back(boundary_points_new);
            //std::cout<<pi<<" location: ("<<boundary_points[pi].loc.x<<","<<boundary_points[pi].loc.y<<")."<<"; on line "<<vj<<","<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance == surface_vertex_distance){
        // the new boundary point will be exactly the suface vertex i+1
            boundary_points.push_back(anticlockwise_surface_vertices[vj+1]);
            vj++;
            //std::cout<<pi<<" location: ("<<boundary_points[pi].loc.x<<","<<boundary_points[pi].loc.y<<")."<<"; on line "<<vj<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance > surface_vertex_distance){
        // the new boundary point might be generated on the line segment between surface vertex i+1 and surface vertex i+2
            required_boundary_points_distance -= surface_vertex_distance;
            vj++;
            distance_intial_is_boundary_vertex=1;
            goto boundary_points_distance_larger_than_surface_vertex_distance;
        }
        
        required_boundary_points_distance = boundary_points_average_distance;
        distance_intial_is_boundary_vertex=0;
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    }

        //cout_fout_debug::fout_vector_vertex(boundary_points,"50_boundary_points.txt");
        return boundary_points;

    }


    vector<Vertex*> organ_ordered_boundary_points_finding_pointer(Organ* p_g, int boundary_points_number){
        //std::cout<<"Start boundary points finding"<<endl;
        double boundary_points_average_distance = p_g->perimeter/boundary_points_number; 
        //std::cout<<"Organ perimeter "<<p_g->perimeter<<"; boundary_points: "<<boundary_points_number<<"; the averaged distance between boundary_points: "<<boundary_points_average_distance<<endl;

        vector<Vertex*> boundary_vertices;
        vector<Vertex*> boundary_points;
        vector<Line*> boundary_line;

            for(int li=0;li<(int)p_g->p_l.size();li++){
                if(p_g->p_l[li]->IsOutermost==1){
                    boundary_line.push_back(p_g->p_l[li]);
                }
            }

        //we need first to arrange the surface vertices into an anticlockwisely ordered array

        //the first boundary vertex: the surface vertex with lowest y position
        Vertex* first_boundary_vertex = new Vertex;

        for(int vi=0;vi<(int)p_g->p_v.size(); vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                first_boundary_vertex = p_g->p_v[vi];
                break;
            }
        }

        for(int vi=0;vi<(int)p_g->p_v.size(); vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                if(first_boundary_vertex->loc.y>p_g->p_v[vi]->loc.y){
                    first_boundary_vertex = p_g->p_v[vi];
                }
            }
        }
        
        boundary_vertices.push_back(first_boundary_vertex);

        //std::cout<<" The first boundary vertices is "<<boundary_vertices[0]->vi<<", and its (x,y) is ("<<boundary_vertices[0]->loc.x<<","<<boundary_vertices[0]->loc.y<<")."<<endl;

    //the second boundary vertex: the surface vertex that is the neighbor of the first surface vertex (they share the same lines)
        //                              and this surface vertex should have smaller x values to make the direction anticlockwise
        Vertex* v_tmp1 = new Vertex;
        Vertex* v_tmp2 = new Vertex;
        bool v_tmp1_found=0;
            for(int li=0; li<(int)boundary_line.size(); li++){
                if(boundary_line[li]->vi[0]==boundary_vertices[0]->vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[boundary_line[li]->vi[1]]->loc.x: "<<p_g->p_v[boundary_line[li]->vi[1]]->loc.x<<endl;
                    v_tmp1 = p_g->p_v[boundary_line[li]->vi[1]];
                    v_tmp1_found=1;
                    //if(p_g->p_v[boundary_line[li]->vi[1]]->loc.x<boundary_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[1]]);
                    //}
                }
                else if(boundary_line[li]->vi[0]==boundary_vertices[0]->vi&&v_tmp1_found==1){
                    v_tmp2 = p_g->p_v[boundary_line[li]->vi[1]];
                    break;
                }
                else if(boundary_line[li]->vi[1]==boundary_vertices[0]->vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[boundary_line[li]->vi[0]]->loc.x: "<<p_g->p_v[boundary_line[li]->vi[0]]->loc.x<<endl;
                    v_tmp1 = p_g->p_v[boundary_line[li]->vi[0]];
                    v_tmp1_found=1;
                    //if(p_g->p_v[boundary_line[li]->vi[0]]->loc.x<boundary_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[0]]);
                    //}
                }
                else if(boundary_line[li]->vi[1]==boundary_vertices[0]->vi&&v_tmp1_found==1){
                    v_tmp2 = p_g->p_v[boundary_line[li]->vi[0]];
                    break;
                }
                else{

                }
            }
            if(v_tmp1->loc.x < v_tmp2->loc.x){
                boundary_vertices.push_back(v_tmp1);
            }
            else{
                boundary_vertices.push_back(v_tmp2);
            }
        //std::cout<<" The secondary boundary vertices is "<<boundary_vertices[1]->vi<<", and its (x,y) is ("<<boundary_vertices[1]->loc.x<<","<<boundary_vertices[1]->loc.y<<")."<<endl;

    //the tertiary and following boundary vertices could be retrieved using incremental algorithm: 
    for(int vi=1; vi<(int)boundary_line.size()-1;vi++){
        for(int li=0;li<(int)boundary_line.size();li++){
                if(boundary_line[li]->vi[0]==boundary_vertices[vi]->vi&&boundary_line[li]->vi[1]!=boundary_vertices[vi-1]->vi){
                    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[1]]);
                }
                else if(boundary_line[li]->vi[1]==boundary_vertices[vi]->vi&&boundary_line[li]->vi[0]!=boundary_vertices[vi-1]->vi){
                    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[0]]);
                }
        }
    }
        //std::cout<<" The third boundary vertices is "<<boundary_vertices[2]->vi<<", and its (x,y) is ("<<boundary_vertices[2]->loc.x<<","<<boundary_vertices[2]->loc.y<<")."<<endl;
        //cout_fout_debug::cout_vector_vertex(boundary_vertices);

    //now we have the ordered boundary vertices: the boundary vertex with lowest y is ranked as the first; the other boundary vertices are ordered anticlockwise
    //the first boundary point: just pick the first boundary vertex (with lowest y)
        boundary_points.push_back(boundary_vertices[0]);
        boundary_vertices.push_back(boundary_vertices[0]); //link the last boundary vertices with the first boundary vertices
        
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    //the second boundary point: the point that is the left neighbor of the first point and have a distance of the boundary_averaged_distance

        double required_boundary_points_distance = boundary_points_average_distance;
        int vj = 0; //the current boundary_vertex index
        bool distance_intial_is_boundary_vertex=0;
        double surface_vertex_distance;

    for(int pi=0; pi<boundary_points_number-1; pi++){
        
        boundary_points_distance_larger_than_surface_vertex_distance:;
        if(distance_intial_is_boundary_vertex==0) {
            surface_vertex_distance = (boundary_vertices[vj+1]->loc-boundary_points[pi]->loc).norm();
        }
        else{
            surface_vertex_distance = (boundary_vertices[vj+1]->loc-boundary_vertices[vj]->loc).norm();
        }


    //judgement: if boundary_point_average_distance >=< surface_vertex_distance
        if(required_boundary_points_distance < surface_vertex_distance){
        // the new boundary point will be generated on the line segment between surface vertex i and surface vertex i+1
            Vertex* boundary_points_new = new Vertex;
            if(distance_intial_is_boundary_vertex==0) {
                boundary_points_new->loc = boundary_points[pi]->loc+(boundary_vertices[vj+1]->loc-boundary_points[pi]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }
            else{
                boundary_points_new->loc = boundary_vertices[vj]->loc+(boundary_vertices[vj+1]->loc-boundary_vertices[vj]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }

            boundary_points.push_back(boundary_points_new);
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<","<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance == surface_vertex_distance){
        // the new boundary point will be exactly the suface vertex i+1
            boundary_points.push_back(boundary_vertices[vj+1]);
            vj++;
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance > surface_vertex_distance){
        // the new boundary point might be generated on the line segment between surface vertex i+1 and surface vertex i+2
            required_boundary_points_distance -= surface_vertex_distance;
            vj++;
            distance_intial_is_boundary_vertex=1;
            goto boundary_points_distance_larger_than_surface_vertex_distance;
        }
        
        required_boundary_points_distance = boundary_points_average_distance;
        distance_intial_is_boundary_vertex=0;
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    }

        //cout_fout_debug::fout_vector_vertex(boundary_points,"50_boundary_points.txt");
        return boundary_points;
    }

    vector<Vertex> organ_ordered_boundary_points_finding(Organ* p_g, int boundary_points_number){
        //std::cout<<"Start boundary points finding"<<endl;
        double boundary_points_average_distance = p_g->perimeter/boundary_points_number; 
        //std::cout<<"Organ perimeter "<<p_g->perimeter<<"; boundary_points: "<<boundary_points_number<<"; the averaged distance between boundary_points: "<<boundary_points_average_distance<<endl;

        vector<Vertex> boundary_vertices;
        vector<Vertex> boundary_points;
        vector<Line> boundary_line;

            for(int li=0;li<(int)p_g->p_l.size();li++){
                if(p_g->p_l[li]->IsOutermost==1){
                    boundary_line.push_back(*p_g->p_l[li]);
                }
            }

        //we need first to arrange the surface vertices into an anticlockwisely ordered array

        //the first boundary vertex: the surface vertex with lowest y position
        Vertex first_boundary_vertex;

        for(int vi=0;vi<(int)p_g->p_v.size(); vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                first_boundary_vertex = *p_g->p_v[vi];
                break;
            }
        }

        for(int vi=0;vi<(int)p_g->p_v.size(); vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                if(first_boundary_vertex.loc.y>p_g->p_v[vi]->loc.y){
                    first_boundary_vertex = *p_g->p_v[vi];
                }
            }
        }
        
        boundary_vertices.push_back(first_boundary_vertex);

        //std::cout<<" The first boundary vertices is "<<boundary_vertices[0]->vi<<", and its (x,y) is ("<<boundary_vertices[0]->loc.x<<","<<boundary_vertices[0]->loc.y<<")."<<endl;

    //the second boundary vertex: the surface vertex that is the neighbor of the first surface vertex (they share the same lines)
        //                              and this surface vertex should have smaller x values to make the direction anticlockwise
        Vertex v_tmp1;
        Vertex v_tmp2;
        bool v_tmp1_found=0;
            for(int li=0; li<(int)boundary_line.size(); li++){
                if(boundary_line[li].vi[0]==boundary_vertices[0].vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[boundary_line[li]->vi[1]]->loc.x: "<<p_g->p_v[boundary_line[li]->vi[1]]->loc.x<<endl;
                    v_tmp1 = *p_g->p_v[boundary_line[li].vi[1]];
                    v_tmp1_found=1;
                    //if(p_g->p_v[boundary_line[li]->vi[1]]->loc.x<boundary_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[1]]);
                    //}
                }
                else if(boundary_line[li].vi[0]==boundary_vertices[0].vi&&v_tmp1_found==1){
                    v_tmp2 = *p_g->p_v[boundary_line[li].vi[1]];
                    break;
                }
                else if(boundary_line[li].vi[1]==boundary_vertices[0].vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[boundary_line[li]->vi[0]]->loc.x: "<<p_g->p_v[boundary_line[li]->vi[0]]->loc.x<<endl;
                    v_tmp1 = *p_g->p_v[boundary_line[li].vi[0]];
                    v_tmp1_found=1;
                    //if(p_g->p_v[boundary_line[li]->vi[0]]->loc.x<boundary_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[0]]);
                    //}
                }
                else if(boundary_line[li].vi[1]==boundary_vertices[0].vi&&v_tmp1_found==1){
                    v_tmp2 = *p_g->p_v[boundary_line[li].vi[0]];
                    break;
                }
                else{

                }
            }
            if(v_tmp1.loc.x < v_tmp2.loc.x){
                boundary_vertices.push_back(v_tmp1);
            }
            else{
                boundary_vertices.push_back(v_tmp2);
            }
        //std::cout<<" The secondary boundary vertices is "<<boundary_vertices[1].vi<<", and its (x,y) is ("<<boundary_vertices[1].loc.x<<","<<boundary_vertices[1].loc.y<<")."<<endl;

    //the tertiary and following boundary vertices could be retrieved using incremental algorithm: 
    for(int vi=1; vi<(int)boundary_line.size()-1;vi++){
        for(int li=0;li<(int)boundary_line.size();li++){
                if(boundary_line[li].vi[0]==boundary_vertices[vi].vi&&boundary_line[li].vi[1]!=boundary_vertices[vi-1].vi){
                    boundary_vertices.push_back(*p_g->p_v[boundary_line[li].vi[1]]);
                }
                else if(boundary_line[li].vi[1]==boundary_vertices[vi].vi&&boundary_line[li].vi[0]!=boundary_vertices[vi-1].vi){
                    boundary_vertices.push_back(*p_g->p_v[boundary_line[li].vi[0]]);
                }
        }
    }
        //std::cout<<" The third boundary vertices is "<<boundary_vertices[2]->vi<<", and its (x,y) is ("<<boundary_vertices[2]->loc.x<<","<<boundary_vertices[2]->loc.y<<")."<<endl;
        //cout_fout_debug::cout_vector_vertex(boundary_vertices);

    //now we have the ordered boundary vertices: the boundary vertex with lowest y is ranked as the first; the other boundary vertices are ordered anticlockwise
    //the first boundary point: just pick the first boundary vertex (with lowest y)
        boundary_points.push_back(boundary_vertices[0]);
        boundary_vertices.push_back(boundary_vertices[0]); //link the last boundary vertices with the first boundary vertices
        
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    //the second boundary point: the point that is the left neighbor of the first point and have a distance of the boundary_averaged_distance

        double required_boundary_points_distance = boundary_points_average_distance;
        int vj = 0; //the current boundary_vertex index
        bool distance_intial_is_boundary_vertex=0;
        double surface_vertex_distance;

    for(int pi=0; pi<boundary_points_number-1; pi++){
        
        boundary_points_distance_larger_than_surface_vertex_distance:;
        if(distance_intial_is_boundary_vertex==0) {
            surface_vertex_distance = (boundary_vertices[vj+1].loc-boundary_points[pi].loc).norm();
        }
        else{
            surface_vertex_distance = (boundary_vertices[vj+1].loc-boundary_vertices[vj].loc).norm();
        }


    //judgement: if boundary_point_average_distance >=< surface_vertex_distance
        if(required_boundary_points_distance < surface_vertex_distance){
        // the new boundary point will be generated on the line segment between surface vertex i and surface vertex i+1
            Vertex boundary_points_new;
            if(distance_intial_is_boundary_vertex==0) {
                boundary_points_new.loc = boundary_points[pi].loc+(boundary_vertices[vj+1].loc-boundary_points[pi].loc)*required_boundary_points_distance/surface_vertex_distance;
            }
            else{
                boundary_points_new.loc = boundary_vertices[vj].loc+(boundary_vertices[vj+1].loc-boundary_vertices[vj].loc)*required_boundary_points_distance/surface_vertex_distance;
            }

            boundary_points.push_back(boundary_points_new);
            //std::cout<<pi<<" location: ("<<boundary_points[pi].loc.x<<","<<boundary_points[pi].loc.y<<")."<<"; on line "<<vj<<","<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance == surface_vertex_distance){
        // the new boundary point will be exactly the suface vertex i+1
            boundary_points.push_back(boundary_vertices[vj+1]);
            vj++;
            //std::cout<<pi<<" location: ("<<boundary_points[pi].loc.x<<","<<boundary_points[pi].loc.y<<")."<<"; on line "<<vj<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance > surface_vertex_distance){
        // the new boundary point might be generated on the line segment between surface vertex i+1 and surface vertex i+2
            required_boundary_points_distance -= surface_vertex_distance;
            vj++;
            distance_intial_is_boundary_vertex=1;
            goto boundary_points_distance_larger_than_surface_vertex_distance;
        }
        
        required_boundary_points_distance = boundary_points_average_distance;
        distance_intial_is_boundary_vertex=0;
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    }

        //cout_fout_debug::fout_vector_vertex(boundary_points,"50_boundary_points.txt");
        return boundary_points;
    }

    void organ_vertex_counterclockwise_sort(Organ* p_g){
        for(int cidx=0; cidx < (int)p_g->p_c.size();++cidx){
            std::vector<int> tmp_vi;
            std::vector<int> tmp_li;

            int prev_vidx = -1, curr_vidx = p_g->p_c[cidx]->vi[0];
            tmp_vi.push_back(curr_vidx);
            while(tmp_li.size() < p_g->p_c[cidx]->li.size()) {
                for(int lidx: p_g->p_c[cidx]->li) {
                    if(p_g->p_l[lidx]->vi[0] == curr_vidx && p_g->p_l[lidx]->vi[1] != prev_vidx) {
                    prev_vidx = curr_vidx;
                    curr_vidx = p_g->p_l[lidx]->vi[1];
                    tmp_li.push_back(lidx);
                    if(curr_vidx != tmp_vi[0]){
                        tmp_vi.push_back(curr_vidx);
                    }
                        }else if(p_g->p_l[lidx]->vi[0] != prev_vidx && p_g->p_l[lidx]->vi[1] == curr_vidx) {
                        prev_vidx = curr_vidx;
                        curr_vidx = p_g->p_l[lidx]->vi[0];
                        tmp_li.push_back(lidx);
                        if(curr_vidx != tmp_vi[0]){
                            tmp_vi.push_back(curr_vidx);
                        }
                    }
                }
            }

            assert(tmp_vi.size() == p_g->p_c[cidx]->vi.size());
            assert(tmp_li.size() == p_g->p_c[cidx]->li.size());

            double area = 0.0;
            for(int i = 0; i < (int)tmp_vi.size(); ++i) {
            _vec<double> r1 = p_g->p_v[tmp_vi[i]]->loc;
            _vec<double> r2 = p_g->p_v[tmp_vi[(i + 1) % tmp_vi.size()]]->loc;
            area += 0.5 * (r1 % r2).z;
            }
            if(area < 0) {
            std::reverse(tmp_vi.begin(), tmp_vi.end());
            std::reverse(tmp_li.begin(), tmp_li.end());
            }

            for(int i = 0; i < (int)tmp_vi.size(); ++i) {
            p_g->p_c[cidx]->vi[i] = tmp_vi[i];
            p_g->p_c[cidx]->li[i] = tmp_li[i];
            }
        }

    }
}

namespace geo{

void calcGeometrics(Organ* p_g){ 
    
    organ_geo::organ_line_length(p_g);
    organ_geo::organ_cell_perimeter(p_g);
    organ_geo::organ_minimum_y_v(p_g);
    organ_geo::organ_minimum_x_v(p_g);

    //1. epidermal identity
    organ_geo::epidermal_identity(p_g);
    //2. organ center
    organ_geo::organ_center(p_g);
    //3. organ area
    organ_geo::organ_area(p_g);
    //4. organ perimeter
    organ_geo::organ_perimeter(p_g);
    //5. organ circularity
    organ_geo::organ_circularity(p_g);
    //6. similarity index
    if(similarity_calculation_required==1){
        boundary_geo::similarity_cal_during_simulation(p_g,real_organ_contour_processed_for_similarity_index);
    }   
}

void basic_VTK_analysis(Organ *p_g){

    //organ_geo::organ_vertex_counterclockwise_sort(p_g);
    organ_geo::organ_line_length(p_g);
    organ_geo::organ_cell_perimeter(p_g);
    organ_geo::organ_center(p_g);
    organ_geo::organ_area(p_g);
    organ_geo::epidermal_identity(p_g);
    organ_geo::organ_perimeter(p_g);
    organ_geo::organ_circularity(p_g);
    organ_geo::organ_minimum_y_v(p_g);
    organ_geo::organ_minimum_x_v(p_g);
    boundary_geo::similarity_cal_during_simulation(p_g,real_organ_contour_processed_for_similarity_index);
}

double distP(double x1, double y1, double x2, double y2){
        return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    }

bool relaP(double x1, double y1, double x2, double y2){
        double sigma = 1e-7;
        if(abs(x1-x2)<sigma&&abs(y1-y2)<sigma){
            return 1;
        }
        else{
            return 0;
        }
    }
    
//relationship between lines, =2, identical, =1, parallel, =0, intersection
int relaL(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2){
        //y=ax+b => a = (y1-y2)/(x1-x2), b = (x1y2-y1x2)/(x1-x2)
        double a1 = (yi1-yi2)/(xi1-xi2);
        double a2 = (yj1-yj2)/(xj1-xj2);
        double b1 = (xi1*yi2-yi1*xi2)/(xi1-xi2);
        double b2 = (xj1*yj2-yj1*xj2)/(xj1-xj2);

        double sigma = 1e-7;
        if(abs(a1-a2)<sigma){
            if(abs(b1-b2)<sigma){
                return 2;
            }
            else{
                return 1;
            }   
        }
        else{
            return 0;
        }

    }
    
//find intersection between lines
pair<double, double> intersectLine(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2){
        //confirm that these two lines intersect
        if(relaL(xi1,yi1, xi2, yi2, xj1, yj1, xj2, yj2)!=0){
            std::cout<<"These two lines do not intersect !"<<endl;
        }

        //y1=ax1+b, y2=ax2+b => a = (y1-y2)/(x1-x2), b = (x1y2-y1x2)/(x1-x2)
        double a1 = (yi1-yi2)/(xi1-xi2);
        double a2 = (yj1-yj2)/(xj1-xj2);
        double b1 = (xi1*yi2-yi1*xi2)/(xi1-xi2);
        double b2 = (xj1*yj2-yj1*xj2)/(xj1-xj2);
        //y=a1x+b1, y=a2x+b2 => x = (b2-b1)/(a1-a2), y = (a1b2-a2b1)/(a1-a2)
        double intersect_x = (b2-b1)/(a1-a2);
        double intersect_y = (a1*b2-a2*b1)/(a1-a2);

        pair<double, double> intersect;
        intersect.first = intersect_x;
        intersect.second = intersect_y;
        return intersect;
    }

//=0, if not intersected, =1 if intersected
int relaRS_1(double xi1, double yi1, double xj1, double yj1, double xj2, double yj2){
        double max_yj = max(yj1,yj2);
        double min_yj = min(yj1,yj2);
        if(yi1>max_yj||yi1<min_yj){
            return 0;
        }
        
        //find intersection between y=yi1(x>xi1) and y=ax+b
        double a = (yj1-yj2)/(xj1-xj2);
        double b = (xj1*yj2-yj1*xj2)/(xj1-xj2);
        //do not consider coincidence
        double intersect_x = (yi1-b)/a;
        double intersect_y = yi1;
        if(intersect_x>xi1){
            return 1;
        }
        else{
            return 0;
        }
    }




}

namespace geo_basic{
 
    double EPS = 1e-9;

    bool isIntersectSegmentLine(pLine &s, pLine &l){
        //line s and line l are not parallel to each other
        if(abs(((s.second - s.first) % (l.second - l.first )).z)< EPS){
            return false;
        }
        
        _vec<double> m1 = (l.second - l.first) % (s.first - l.first)  ;
        _vec<double> m2 = (l.second - l.first) % (s.second - l.first) ;
    
        if( m1.z * m2.z > EPS) return false;
        return true;
    }

    _vec<double> crossPoint(geo_basic::pLine &s, geo_basic::pLine &t){
        _vec<double> sv = s.second - s.first;
        _vec<double> tv = t.second - t.first;
        assert((sv % tv).norm() > geo_basic::EPS);
        double length = (tv % (t.first - s.first)).z / (tv % sv).z;
        return s .first + sv * length;
    }
}

namespace geo_vv{


    vector<Vertex*> after_ImageJ_process(vector<Vertex*> vv){
        vector<Vertex*> vv_processed;

        //1. normalization
        vector<Vertex*> vv_normalized = geo_vv::normalization(vv);
        
        //2. swap from bottom to top 
        
        for(int vi=0; vi<vv_normalized.size(); vi++){
                vv_normalized[vi]->loc.y=1-vv_normalized[vi]->loc.y;
        }
        
        return vv_normalized;
        
    }

    vector<Vertex*> normalization_by_perimeter(vector<Vertex*> vv){
        vector<Vertex*> vv_normalized;

        //calculate y_max and y_min
        double y_max = vv[0]->loc.y;
        double y_min = vv[0]->loc.y;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi]->loc.y<y_min){
                y_min=vv[vi]->loc.y;
            }
            if(vv[vi]->loc.y>y_max){
                y_max=vv[vi]->loc.y;
            }
        }
        //std::cout<<"y_min "<<y_min<<"; y_max "<<y_max<<endl;

        //calculate vv_normalization
        double norm_length = perimeter_vv_boundary(vv);
        double x_center=0;
        for(int vi=0;vi<vv.size();vi++){
            x_center+=vv[vi]->loc.x/(double)vv.size();
        }
        //std::cout<<"x_center: "<<x_center<<endl;
        //std::cout<<"Norm length: "<<norm_length<<endl;
        for(int vi=0;vi<(int)vv.size();vi++){
            Vertex* v_tmp = new Vertex;
            v_tmp->loc.x = (vv[vi]->loc.x-x_center)/norm_length;
            v_tmp->loc.y = (vv[vi]->loc.y-y_min)/norm_length;
            //std::cout<<"y "<<vv[vi]->loc.y<<";y_norm "<<v_tmp->loc.y<<endl;
            vv_normalized.push_back(v_tmp);
        }
        return vv_normalized;
    }

    vector<Vertex*> normalization(vector<Vertex*> vv){
        vector<Vertex*> vv_normalized;

        //calculate y_max and y_min
        double y_max = vv[0]->loc.y;
        double y_min = vv[0]->loc.y;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi]->loc.y<y_min){
                y_min=vv[vi]->loc.y;
            }
            if(vv[vi]->loc.y>y_max){
                y_max=vv[vi]->loc.y;
            }
        }
        //std::cout<<"y_min "<<y_min<<"; y_max "<<y_max<<endl;

        //calculate vv_normalization
        double norm_length = y_max-y_min;
        double x_center=0;
        for(int vi=0;vi<vv.size();vi++){
            x_center+=vv[vi]->loc.x/(double)vv.size();
        }
        //std::cout<<"x_center: "<<x_center<<endl;
        //std::cout<<"Norm length: "<<norm_length<<endl;
        for(int vi=0;vi<(int)vv.size();vi++){
            Vertex* v_tmp = new Vertex;
            v_tmp->loc.x = (vv[vi]->loc.x-x_center)/norm_length;
            v_tmp->loc.y = (vv[vi]->loc.y-y_min)/norm_length;
            //std::cout<<"y "<<vv[vi]->loc.y<<";y_norm "<<v_tmp->loc.y<<endl;
            vv_normalized.push_back(v_tmp);
        }
        return vv_normalized;
    }

    vector<Vertex> normalization(vector<Vertex> vv){
        vector<Vertex> vv_normalized;

        //calculate y_max and y_min
        double y_max = vv[0].loc.y;
        double y_min = vv[0].loc.y;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi].loc.y<y_min){
                y_min=vv[vi].loc.y;
            }
            if(vv[vi].loc.y>y_max){
                y_max=vv[vi].loc.y;
            }
        }
        //std::cout<<"y_min "<<y_min<<"; y_max "<<y_max<<endl;

        //calculate vv_normalization
        double norm_length = y_max-y_min;
        double x_center=0;
        for(int vi=0;vi<vv.size();vi++){
            x_center+=vv[vi].loc.x/(double)vv.size();
        }
        //std::cout<<"x_center: "<<x_center<<endl;
        //std::cout<<"Norm length: "<<norm_length<<endl;
        for(int vi=0;vi<(int)vv.size();vi++){
            Vertex v_tmp;
            v_tmp.loc.x = (vv[vi].loc.x-x_center)/norm_length;
            v_tmp.loc.y = (vv[vi].loc.y-y_min)/norm_length;
            //std::cout<<"y "<<vv[vi].loc.y<<";y_norm "<<v_tmp.loc.y<<endl;
            vv_normalized.push_back(v_tmp);
        }
        return vv_normalized;
    }

    bool equal_Euclidean_distance_except_last(vector<Vertex*> vv)
    {   
        vector<double> Eu_dist;
        for(int vi=0;vi<(int)vv.size()-1;vi++){
            Eu_dist.push_back(vertex_geo::vertex_distance(vv[vi],vv[vi+1]));
        }

        //compare all element with the first element

        for(int vi=0;vi<(int)Eu_dist.size();vi++){
            if(abs(Eu_dist[vi]-Eu_dist[0])>EPS_geo){
                return false;
            }
        }
        return true;
    }

    //the vv should be ordered whether clockwise or anticlockwise
    double area_vv_boundary(vector<Vertex*> vv){
        double area =0.0;
        int n=vv.size();
        for(int i=0; i<n; ++i){
            int next = (i+1)%n;
            area += vv[i]->loc.x*vv[next]->loc.y-vv[next]->loc.x*vv[i]->loc.y;
        }
        return abs(area)/2.0;
    }

    double perimeter_vv_boundary(vector<Vertex*> vv){
        double perimeter =0.0;
        int n = vv.size();
        for(int i=0; i<n; ++i){
            int next=(i+1)%n;
            perimeter += vertex_geo::vertex_distance(vv[i],vv[next]);
            //std::cout<<"perimeter "<<i<<" "<<perimeter<<endl;
        }
        return perimeter;
    }

    vector<Vertex*> organ_boundary_points_along_polygon(vector<Vertex*> anticlockwise_surface_vertices, int boundary_points_number){
        vector<Vertex*> boundary_points;
        double perimeter_vv = geo_vv::perimeter_vv_boundary(anticlockwise_surface_vertices);
        double boundary_points_average_distance = perimeter_vv/boundary_points_number;

        //now we have the ordered boundary vertices: the boundary vertex with lowest y is ranked as the first; the other boundary vertices are ordered anticlockwise
    //1. the first boundary point: just pick the first boundary vertex (with lowest y)
        boundary_points.push_back(anticlockwise_surface_vertices[0]);
        anticlockwise_surface_vertices.push_back(anticlockwise_surface_vertices[0]); //link the last boundary vertices with the first boundary vertices
        
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    //2. the i-th boundary point: incremental algorithm

        double required_boundary_points_distance = boundary_points_average_distance;
        int vj = 0; //the current boundary_vertex index
        bool distance_intial_is_boundary_vertex=0;
        double surface_vertex_distance;

    for(int pi=0; pi<boundary_points_number-1; pi++){
        
            boundary_points_distance_larger_than_surface_vertex_distance:;
            if(distance_intial_is_boundary_vertex==0) {
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1]->loc-boundary_points[pi]->loc).norm();
            }
            else{
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1]->loc-anticlockwise_surface_vertices[vj]->loc).norm();
        }


        //judgement: if boundary_point_average_distance >=< surface_vertex_distance
        if(required_boundary_points_distance < surface_vertex_distance){
        // the new boundary point will be generated on the line segment between surface vertex i and surface vertex i+1
            Vertex* boundary_points_new = new Vertex;
            if(distance_intial_is_boundary_vertex==0) {
                boundary_points_new->loc = boundary_points[pi]->loc+(anticlockwise_surface_vertices[vj+1]->loc-boundary_points[pi]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }
            else{
                boundary_points_new->loc = anticlockwise_surface_vertices[vj]->loc+(anticlockwise_surface_vertices[vj+1]->loc-anticlockwise_surface_vertices[vj]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }

            boundary_points.push_back(boundary_points_new);
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<","<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance == surface_vertex_distance){
        // the new boundary point will be exactly the suface vertex i+1
            boundary_points.push_back(anticlockwise_surface_vertices[vj+1]);
            vj++;
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance > surface_vertex_distance){
        // the new boundary point might be generated on the line segment between surface vertex i+1 and surface vertex i+2
            required_boundary_points_distance -= surface_vertex_distance;
            vj++;
            distance_intial_is_boundary_vertex=1;
            goto boundary_points_distance_larger_than_surface_vertex_distance;
        }
        
        required_boundary_points_distance = boundary_points_average_distance;
        distance_intial_is_boundary_vertex=0;
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    }

        //cout_fout_debug::fout_vector_vertex(boundary_points,"50_boundary_points.txt");
        return boundary_points;

    }

    double vd_minimum(vector<double> vd){
        double min_vd=vd[0];
        for(int i=1; i<vd.size();i++){
            if(min_vd>vd[i]){
                min_vd=vd[i];
                //std::cout<<i<<" "<<min_vd<<endl;
            }
        }
        return min_vd;
    }

    double vd_maximum(vector<double> vd){
        double max_vd=vd[0];
        for(int i=1; i<vd.size();i++){
            if(max_vd<vd[i]){
                max_vd=vd[i];
                //std::cout<<i<<" "<<max_vd<<endl;
            }
        }
        return max_vd;
    }

    double accumulated_negative(vector<double> vd){
        double accumulated=0;
        for(int i=0; i<vd.size(); i++){
            if(vd[i]<0){
                accumulated+=vd[i]/(double)vd.size();
            }
        }
        return accumulated;
    }

    vector<double> vd_averaged(vector<double> vd, int numNeighbors){
        vector<double> output_vd(vd.size());

        int halfNeighbors = (numNeighbors-1)/2;

        for(int i=0;i<vd.size();++i){
            double sum=0.0;

            for(int j=i-halfNeighbors-1;j<i+halfNeighbors;++j){
                if(j>=0 && j<vd.size()){
                    sum+=vd[j];
                }
                else if(j<0){
                    sum+=vd[vd.size()-1+j];
                }
                else if(j>vd.size()){
                    sum+=vd[j-vd.size()];
                }
            }

            output_vd[i] = sum/numNeighbors;
        }

        return output_vd;
    }

    vector<Vertex*> vector_vertex_sampling(vector<Vertex*> vv,double sampling_distance){
        vector<Vertex*> vv_sampled;
        //divide the organ into left side and right side based on the theta 
        //when sampling_distance = 0.01, if theta_0<theta_100, left side, else right side
        vector<double> theta;
        vector<Vertex*> left_side_vertex;
        vector<Vertex*> right_side_vertex;
        //find center of outline vertices
        double center_x=0,center_y=0;
        for(int i=0;i<(int)vv.size();i++){
            center_x = vv[i]->loc.x/(double)vv.size();
            center_y = vv[i]->loc.y/(double)vv.size();
        }
        //std::cout<<"center_x "<<center_x<<" center_y "<<center_y<<endl;
        for(int i=0;i<(int)vv.size();i++){
            double x_tmp=vv[i]->loc.x-center_x;
            double y_tmp=vv[i]->loc.y-center_y;
            double theta_tmp = atan2(y_tmp,x_tmp);
            if(theta_tmp<0){

            }
            theta.push_back(theta_tmp);
        }

        double x_0,x_100,theta_0,theta_100;
        int vi_0,vi_100;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi]->loc.y==0.00){
                x_0 = vv[vi]->loc.x;
                theta_0 = theta[vi];
                vi_0=vi;
            }
            else if(vv[vi]->loc.y==1.00){
                x_100 = vv[vi]->loc.x;
                theta_100 = theta[vi];
                vi_100=vi;
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0){
                left_side_vertex.push_back(vv[vi]);
            }
            else if(vi==vi_100){
                right_side_vertex.push_back(vv[vi]);
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0||vi==vi_100){
                
            }
            else if(theta[vi]>theta_0&&theta[vi]<theta_100){
                left_side_vertex.push_back(vv[vi]);
            }
            else{
                right_side_vertex.push_back(vv[vi]);
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0){
                right_side_vertex.push_back(vv[vi]);
            }
            else if(vi==vi_100){
                left_side_vertex.push_back(vv[vi]);
            }
        }

        //sorting left_side and right_side by y value
        left_side_vertex = sort_vector_vertex_ascend(left_side_vertex);
        right_side_vertex = sort_vector_vertex_descend(right_side_vertex);

        Vertex* v_0 = new Vertex;
        Vertex* v_100 = new Vertex;
        v_0->loc.x = x_0;
        v_0->loc.y =0;
        v_100->loc.x=x_100;
        v_100->loc.y=1.0;

        vv_sampled.push_back(v_0);
        for(int i=1;i<100;i++){
            for(int vi=0;vi<left_side_vertex.size();vi++){
                //std::cout<<(double)i*0.01<<" "<<left_side_vertex[vi]->loc.y<<" "<<left_side_vertex[vi+1]->loc.y<<endl;
                if(left_side_vertex[vi]->loc.y==(double)0.01*i){
                    vv_sampled.push_back(left_side_vertex[vi]);
                }
                else if(left_side_vertex[vi]->loc.y<(double)i*0.01&&left_side_vertex[vi+1]->loc.y>(double)i*0.01){
                    //std::cout<<"fit"<<endl;
                    Vertex* vt_tmp = new Vertex;
                    vt_tmp->loc.y = i*0.01;
                    vt_tmp->loc.x = linear_fitting(left_side_vertex[vi],left_side_vertex[vi+1],i*0.01);
                    //vt_tmp->loc.x = i*0.01;
                    vv_sampled.push_back(vt_tmp);
                }
                else{
                    //std::cout<<"Error in finding fitted simulated outline"<<endl;
                }
            }
        }
        vv_sampled.push_back(v_100);

        for(int i=101;i<200;i++){
            for(int vi=0;vi<right_side_vertex.size();vi++){
                if(right_side_vertex[vi]->loc.y==(2-0.01*i)){
                    vv_sampled.push_back(right_side_vertex[vi]);
                }
                else if(right_side_vertex[vi]->loc.y>(2-i*0.01)&&right_side_vertex[vi+1]->loc.y<(2-i*0.01)){
                    Vertex* vt_tmp = new Vertex;
                    vt_tmp->loc.y = 2-i*0.01;
                    vt_tmp->loc.x = linear_fitting(right_side_vertex[vi],right_side_vertex[vi+1],2-i*0.01);
                    vv_sampled.push_back(vt_tmp);
                }
                else{
                    //std::cout<<"Error in finding fitted simulated outline"<<endl;
                }
            }
        }
        /*
        ofstream fout3("sampled.txt");
        fout3<<"x y"<<endl;
        for(int i=0;i<vv_sampled.size();i++){
            fout3<<vv_sampled[i]->loc.x<<" "<<vv_sampled[i]->loc.y<<endl;
        }
        fout3.close();
        
        
        ofstream fout("sampled_left.txt");
        fout<<"x y"<<endl;
        for(int i=0;i<left_side_vertex.size();i++){
            fout<<left_side_vertex[i]->loc.x<<" "<<left_side_vertex[i]->loc.y<<endl;
        }
        fout.close();

        ofstream fout2("sampled_right.txt");
        fout2<<"x y"<<endl;
        for(int i=0;i<right_side_vertex.size();i++){
            fout2<<right_side_vertex[i]->loc.x<<" "<<right_side_vertex[i]->loc.y<<endl;
        }
        fout2.close();
        */
        return vv_sampled;
    }

    vector<Vertex> vector_vertex_sampling(vector<Vertex> vv,double sampling_distance){
        vector<Vertex> vv_sampled;
        //divide the organ into left side and right side based on the theta 
        //when sampling_distance = 0.01, if theta_0<theta_100, left side, else right side
        vector<double> theta;
        vector<Vertex> left_side_vertex;
        vector<Vertex> right_side_vertex;
        //find center of outline vertices
        double center_x=0,center_y=0;
        for(int i=0;i<(int)vv.size();i++){
            center_x = vv[i].loc.x/(double)vv.size();
            center_y = vv[i].loc.y/(double)vv.size();
        }
        //std::cout<<"center_x "<<center_x<<" center_y "<<center_y<<endl;
        for(int i=0;i<(int)vv.size();i++){
            double x_tmp=vv[i].loc.x-center_x;
            double y_tmp=vv[i].loc.y-center_y;
            double theta_tmp = atan2(y_tmp,x_tmp);
            if(theta_tmp<0){

            }
            theta.push_back(theta_tmp);
        }

        double x_0,x_100,theta_0,theta_100;
        int vi_0,vi_100;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi].loc.y==0.00){
                x_0 = vv[vi].loc.x;
                theta_0 = theta[vi];
                vi_0=vi;
            }
            else if(vv[vi].loc.y==1.00){
                x_100 = vv[vi].loc.x;
                theta_100 = theta[vi];
                vi_100=vi;
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0){
                left_side_vertex.push_back(vv[vi]);
            }
            else if(vi==vi_100){
                right_side_vertex.push_back(vv[vi]);
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0||vi==vi_100){
                
            }
            else if(theta[vi]>theta_0&&theta[vi]<theta_100){
                left_side_vertex.push_back(vv[vi]);
            }
            else{
                right_side_vertex.push_back(vv[vi]);
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0){
                right_side_vertex.push_back(vv[vi]);
            }
            else if(vi==vi_100){
                left_side_vertex.push_back(vv[vi]);
            }
        }

        //sorting left_side and right_side by y value
        left_side_vertex = sort_vector_vertex_ascend(left_side_vertex);
        right_side_vertex = sort_vector_vertex_descend(right_side_vertex);

        Vertex v_0;
        Vertex v_100;
        v_0.loc.x = x_0;
        v_0.loc.y =0;
        v_100.loc.x=x_100;
        v_100.loc.y=1.0;

        vv_sampled.push_back(v_0);
        for(int i=1;i<100;i++){
            for(int vi=0;vi<left_side_vertex.size();vi++){
                //std::cout<<(double)i*0.01<<" "<<left_side_vertex[vi]->loc.y<<" "<<left_side_vertex[vi+1]->loc.y<<endl;
                if(left_side_vertex[vi].loc.y==(double)0.01*i){
                    vv_sampled.push_back(left_side_vertex[vi]);
                }
                else if(left_side_vertex[vi].loc.y<(double)i*0.01&&left_side_vertex[vi+1].loc.y>(double)i*0.01){
                    //std::cout<<"fit"<<endl;
                    Vertex vt_tmp;
                    vt_tmp.loc.y = i*0.01;
                    vt_tmp.loc.x = linear_fitting(left_side_vertex[vi],left_side_vertex[vi+1],i*0.01);
                    //vt_tmp->loc.x = i*0.01;
                    vv_sampled.push_back(vt_tmp);
                }
                else{
                    //std::cout<<"Error in finding fitted simulated outline"<<endl;
                }
            }
        }
        vv_sampled.push_back(v_100);

        for(int i=101;i<200;i++){
            for(int vi=0;vi<right_side_vertex.size();vi++){
                if(right_side_vertex[vi].loc.y==(2-0.01*i)){
                    vv_sampled.push_back(right_side_vertex[vi]);
                }
                else if(right_side_vertex[vi].loc.y>(2-i*0.01)&&right_side_vertex[vi+1].loc.y<(2-i*0.01)){
                    Vertex vt_tmp;
                    vt_tmp.loc.y = 2-i*0.01;
                    vt_tmp.loc.x = linear_fitting(right_side_vertex[vi],right_side_vertex[vi+1],2-i*0.01);
                    vv_sampled.push_back(vt_tmp);
                }
                else{
                    //std::cout<<"Error in finding fitted simulated outline"<<endl;
                }
            }
        }
        /*
        ofstream fout3("sampled.txt");
        fout3<<"x y"<<endl;
        for(int i=0;i<vv_sampled.size();i++){
            fout3<<vv_sampled[i]->loc.x<<" "<<vv_sampled[i]->loc.y<<endl;
        }
        fout3.close();
        
        
        ofstream fout("sampled_left.txt");
        fout<<"x y"<<endl;
        for(int i=0;i<left_side_vertex.size();i++){
            fout<<left_side_vertex[i]->loc.x<<" "<<left_side_vertex[i]->loc.y<<endl;
        }
        fout.close();

        ofstream fout2("sampled_right.txt");
        fout2<<"x y"<<endl;
        for(int i=0;i<right_side_vertex.size();i++){
            fout2<<right_side_vertex[i]->loc.x<<" "<<right_side_vertex[i]->loc.y<<endl;
        }
        fout2.close();
        */
        return vv_sampled;
    }
    
    vector<Vertex*> sort_vector_vertex_ascend(vector<Vertex*> vv){
        vector<Vertex*> sort_vertex;

        vector<double> sort_y;
        for(int vi=0;vi<(int)vv.size();vi++){
            sort_y.push_back(vv[vi]->loc.y);
        }

        sort(sort_y.begin(),sort_y.end());
        for(int vi=0;vi<(int)vv.size();vi++){
            for(int vj=0;vj<(int)vv.size();vj++){
                if(sort_y[vi]==vv[vj]->loc.y){
                    sort_vertex.push_back(vv[vj]);
                }
            }
        }

        return sort_vertex;
    }

    vector<Vertex> sort_vector_vertex_ascend(vector<Vertex> vv){
        vector<Vertex> sort_vertex;

        vector<double> sort_y;
        for(int vi=0;vi<(int)vv.size();vi++){
            sort_y.push_back(vv[vi].loc.y);
        }

        sort(sort_y.begin(),sort_y.end());
        for(int vi=0;vi<(int)vv.size();vi++){
            for(int vj=0;vj<(int)vv.size();vj++){
                if(sort_y[vi]==vv[vj].loc.y){
                    sort_vertex.push_back(vv[vj]);
                }
            }
        }

        return sort_vertex;
    }

    vector<Vertex*> sort_vector_vertex_descend(vector<Vertex*> vv){
        vector<Vertex*> sort_vertex;

        vector<double> sort_y;
        for(int vi=0;vi<(int)vv.size();vi++){
            sort_y.push_back(vv[vi]->loc.y);
        }

        sort(sort_y.begin(),sort_y.end(),comp_descend);
        for(int vi=0;vi<(int)vv.size();vi++){
            for(int vj=0;vj<(int)vv.size();vj++){
                if(sort_y[vi]==vv[vj]->loc.y){
                    sort_vertex.push_back(vv[vj]);
                }
            }
        }

        return sort_vertex;
    }

    vector<Vertex> sort_vector_vertex_descend(vector<Vertex> vv){
        vector<Vertex> sort_vertex;

        vector<double> sort_y;
        for(int vi=0;vi<(int)vv.size();vi++){
            sort_y.push_back(vv[vi].loc.y);
        }

        sort(sort_y.begin(),sort_y.end(),comp_descend);
        for(int vi=0;vi<(int)vv.size();vi++){
            for(int vj=0;vj<(int)vv.size();vj++){
                if(sort_y[vi]==vv[vj].loc.y){
                    sort_vertex.push_back(vv[vj]);
                }
            }
        }

        return sort_vertex;
    }

    bool comp_descend(double y1,double y2){
        return (y1>y2);
    }

    double linear_fitting(Vertex* v1, Vertex* v2, double y_i){
        double x_i;
        if(v1->loc.x == v2->loc.x){
            x_i = v1->loc.x;
        }
        else{
            double a_tmp = (v1->loc.y-v2->loc.y)/(v1->loc.x-v2->loc.x);
            double b_tmp = (v1->loc.x*v2->loc.y-v2->loc.x*v1->loc.y)/(v1->loc.x-v2->loc.x);
            x_i = (y_i-b_tmp)/a_tmp;
        }
        //std::cout<<"a_tmp "<<a_tmp<<" b_tmp "<<b_tmp<<" x_i "<<x_i<<endl;
        return x_i;
    }

    double linear_fitting(Vertex v1, Vertex v2, double y_i){
        double x_i;
        if(v1.loc.x == v2.loc.x){
            x_i = v1.loc.x;
        }
        else{
            double a_tmp = (v1.loc.y-v2.loc.y)/(v1.loc.x-v2.loc.x);
            double b_tmp = (v1.loc.x*v2.loc.y-v2.loc.x*v1.loc.y)/(v1.loc.x-v2.loc.x);
            x_i = (y_i-b_tmp)/a_tmp;
        }
        //std::cout<<"a_tmp "<<a_tmp<<" b_tmp "<<b_tmp<<" x_i "<<x_i<<endl;
        return x_i;
    }

    vector<Vertex*> vv_x_swap(vector<Vertex*> vv){
        vector<Vertex*> vv_swaped;
        for(int i=0;i<vv.size();i++){
            Vertex* v_tmp = new Vertex;
            v_tmp->loc.x = -vv[i]->loc.x;
            v_tmp->loc.y = vv[i]->loc.y;
            vv_swaped.push_back(v_tmp);
        }
        return vv_swaped;
    }

    vector<Vertex> vv_x_swap(vector<Vertex> vv){
        vector<Vertex> vv_swaped;
        for(int i=0;i<vv.size();i++){
            Vertex v_tmp;
            v_tmp.loc.x = -vv[i].loc.x;
            v_tmp.loc.y = vv[i].loc.y;
            vv_swaped.push_back(v_tmp);
        }
        return vv_swaped;
    }

    
}

namespace boundary_geo{

    double similarity_Index_1(vector<Vertex*> vv_simu,vector<Vertex*> vv_real, double sampling_distance){
        //adjust the real y_position
        /*
        for(int vi=0;vi<(int)vv_real.size();vi++){
            vv_real[vi]->loc.y = abs(1-vv_real[vi]->loc.y);
        }
        vector<Vertex*> vv_real_align;
        for(int i=0;i<101;i++){
            vv_real_align.push_back(vv_real[100-i]);
        }
        for(int i=101;i<200;i++){
            vv_real_align.push_back(vv_real[i-100]);
        }
        //output::vv_foutXY(vv_real_align,"real_aligned.txt");
        */
        double error_sum=0;
        for(int vi=0;vi<(int)vv_simu.size();vi++){
            double x_error = abs(vv_simu[vi]->loc.x-vv_real[vi]->loc.x);
            error_sum += x_error;
            //std::cout<<vv_simu[vi]->loc.x<<" "<<vv_real[vi]->loc.x<<" "<<x_error<<endl;
        }
        double error_area = error_sum*sampling_distance;
        return error_area;
    }

    double similarity_Index_1(vector<Vertex> vv_simu,vector<Vertex*> vv_real, double sampling_distance){
        //adjust the real y_position
        /*
        for(int vi=0;vi<(int)vv_real.size();vi++){
            vv_real[vi]->loc.y = abs(1-vv_real[vi]->loc.y);
        }
        vector<Vertex*> vv_real_align;
        for(int i=0;i<101;i++){
            vv_real_align.push_back(vv_real[100-i]);
        }
        for(int i=101;i<200;i++){
            vv_real_align.push_back(vv_real[i-100]);
        }
        //output::vv_foutXY(vv_real_align,"real_aligned.txt");
        */
        double error_sum=0;
        for(int vi=0;vi<(int)vv_simu.size();vi++){
            double x_error = abs(vv_simu[vi].loc.x-vv_real[vi]->loc.x);
            error_sum += x_error;
            //std::cout<<vv_simu[vi]->loc.x<<" "<<vv_real[vi]->loc.x<<" "<<x_error<<endl;
        }
        double error_area = error_sum*sampling_distance;
        return error_area;
    }

    double similarity_Index_2(string real_contour_file, string simulated_contour_file, double sampling_distance){
        vector<Vertex*> real_contour = readV::xy_txt_to_vertexXY(real_contour_file);
        vector<Vertex*> real_contour_normalized = geo_vv::after_ImageJ_process(real_contour);
        vector<Vertex*> real_contour_sampled = geo_vv::vector_vertex_sampling(real_contour_normalized,sampling_distance);
        //gnu_plot::organ_contour_plot(real_contour_normalized);

        
        vector<Vertex*> simulated_contour = readV::read_vv(simulated_contour_file);
        vector<Vertex*> simulated_contour_normalized = geo_vv::normalization(simulated_contour);
        vector<Vertex*> simulated_contour_swapped = geo_vv::vv_x_swap(simulated_contour_normalized);
        vector<Vertex*> simulated_outline_sampled = geo_vv::vector_vertex_sampling(simulated_contour_normalized,sampling_distance);
        vector<Vertex*> simulated_outline_swapped_sampled = geo_vv::vector_vertex_sampling(simulated_contour_swapped,sampling_distance);
        //gnu_plot::organ_contour_plot(simulated_contour_normalized);
        cout_fout_debug::fout_vector_vertex(simulated_outline_sampled,"simulated_petal.txt");

        double similarity_index_tmp1 = boundary_geo::similarity_Index_1(simulated_outline_sampled,real_contour_sampled,sampling_distance);
        double similarity_index_tmp2 = boundary_geo::similarity_Index_1(simulated_outline_swapped_sampled,real_contour_sampled,sampling_distance);
        
        double similarity_index;

        if(similarity_index_tmp1>similarity_index_tmp2){
            similarity_index = similarity_index_tmp2;
        }
        else{
            similarity_index = similarity_index_tmp1;
        }
        std::cout<<"similarity_Index1 "<<similarity_index_tmp1<<" similarity_index2 "<<similarity_index_tmp2 <<endl;
        std::cout<<"similarity_index "<<similarity_index<<endl;

        return similarity_index;
    }

    vector<Vertex*> read_and_process_real_organ_contour_imagej(string real_contour_file){
        double sampling_distance=0.01;
        vector<Vertex*> real_contour = readV::xy_txt_to_vertexXY(real_contour_file);
        vector<Vertex*> real_contour_normalized = geo_vv::after_ImageJ_process(real_contour);
        vector<Vertex*> real_contour_sampled = geo_vv::vector_vertex_sampling(real_contour_normalized,sampling_distance);
        //std::cout<<"real_contour_sampled size: "<<real_contour_sampled.size()<<std::endl;
        //cout_fout_debug::cout_vector_vertex(real_contour_sampled);
        return real_contour_sampled;
    }

    double similarity_cal_during_simulation(Organ* p_g,vector<Vertex*> real_contour_sampled){
        //std::cout<<"real_contour_sampled size: "<<real_contour_sampled.size()<<std::endl;
        //cout_fout_debug::cout_vector_vertex(real_contour_sampled);
        double similarity_index;
        vector<Vertex> simulated_contour;
        double sampling_distance=0.01;
        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                simulated_contour.push_back(*p_g->p_v[vi]);
            }
        }
        vector<Vertex> simulated_contour_normalized = geo_vv::normalization(simulated_contour);
        vector<Vertex> simulated_contour_swapped = geo_vv::vv_x_swap(simulated_contour_normalized);
        vector<Vertex> simulated_outline_sampled = geo_vv::vector_vertex_sampling(simulated_contour_normalized,sampling_distance);
        vector<Vertex> simulated_outline_swapped_sampled = geo_vv::vector_vertex_sampling(simulated_contour_swapped,sampling_distance);
        //gnu_plot::organ_contour_plot(simulated_contour_normalized);
        //cout_fout_debug::fout_vector_vertex(simulated_outline_sampled,"simulated_petal.txt");
        //cout_fout_debug::cout_vector_vertex(simulated_outline_sampled);

        double similarity_index_tmp1 = boundary_geo::similarity_Index_1(simulated_outline_sampled,real_contour_sampled,sampling_distance);
        double similarity_index_tmp2 = boundary_geo::similarity_Index_1(simulated_outline_swapped_sampled,real_contour_sampled,sampling_distance);

        if(similarity_index_tmp1>similarity_index_tmp2){
            similarity_index = similarity_index_tmp2;
        }
        else{
            similarity_index = similarity_index_tmp1;
        }
        //std::cout<<"similarity_Index1 "<<similarity_index_tmp1<<" similarity_index2 "<<similarity_index_tmp2 <<endl;
        std::cout<<"similarity_index "<<similarity_index<<endl;
        p_g->similarity_index=similarity_index;
        return similarity_index;

    }
}
