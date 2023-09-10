/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef GEO2DV_H
#define GEO2DV_H

#include "class.h"
#include "vec.h"
#include "parameter.h"
#include "wangMath.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <cmath> 
#include <assert.h>

extern bool similarity_calculation_required;
extern vector<Vertex*> real_organ_contour_processed_for_similarity_index;

namespace vertex_geo{
    double vertex_distance(Vertex*,Vertex*);
    double vertex_distance(Vertex,Vertex);
    bool vertex_relationship(Vertex*, Vertex*);
}

namespace line_geo{
    //lines
    bool line_vertex_relationship(Line,Vertex);
    pair<int,Vertex>  lines_relationship(const Line, const Line);
    _vec<double> line_cell_wall_intersection(double,double,Organ*,int);
    //segments
    double line_length(Organ*,int);
    Distance_point_line distance_line_segment_to_vertex(Organ*,int,Vertex);
    Distance_point_line distance_line_segment_to_vertex(Line,Vertex);
    
    bool line_segment_vertex_relationship(Line,Vertex);
    pair<int,Vertex> segments_relationship(Line,Line);

    vector<Line*> region_partition_lines_EdU(vector<Vertex*>);
    vector<double> angles_normalization_changing_axis(vector<Vertex*>,vector<Vertex*>);
}

namespace cell_geo{
    _vec<double> cell_center(Organ*, int);
    double cell_area(Organ*,int);
    double cell_perimeter(Organ*,int);
    double cell_regularity(Organ*,int);
    vector<int> cell_counterclock_line(Organ*,int);
}

namespace organ_geo{
    //1. geometric analysis of organ
    //1.1 basic 
    void organ_line_length(Organ*);
    double organ_cell_perimeter(Organ*);
    _vec<double> organ_center(Organ*);
    void epidermal_identity(Organ*);
    //1.2 area and perimeter
    double organ_area(Organ*);
    double organ_perimeter(Organ*);
    double organ_circularity(Organ*);

    double organ_maximum_y(Organ*);
    double organ_minimum_y(Organ*);
    double organ_minimum_y_v(Organ*);
    double organ_minimum_x_v(Organ*);

    //2. boundary analysis of organ
    //2.1 ordered boundary points generation
    vector<Vertex> organ_ordered_boundary_points_finding(Organ*,int);
    vector<Vertex*> organ_ordered_boundary_points_finding_pointer(Organ*,int);
    Ordered_boundary organ_ordered_anticlockwise_boundary(Organ*);
    vector<Vertex*> organ_boundary_points_along_polygon(Organ*,vector<Vertex*>,int);
    vector<Vertex> organ_boundary_points_along_polygon(Organ*,vector<Vertex>,int);
    vector<Vertex> organ_boundary_points_euclidean(Organ*,double);
    //2.2 curvature analysis

    //bool point_in_polygon_ray_casting(Vertex,Organ*);
    vector<_vec<double>> line_polygon_intersection(double,double,Organ*);
    void organ_vertex_counterclockwise_sort(Organ*);

}

namespace geo{
    void calcGeometrics(Organ*);

    void batch(Batch *, Organ *);

    void basic_VTK_analysis(Organ *);

    double distP(double x1, double y1, double x2, double y2);
    bool relaP(double x1, double y1, double x2, double y2);
    int relaL(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2);
    pair<double, double> intersectLine(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2);
    int relaRS(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2);
    //this ray is parallel to x axis and goes to infinite large 
    int relaRS_1(double xi1, double yi1, double xj1, double yj1, double xj2, double yj2);
    int relaS(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2);
}

namespace geo_basic{
    typedef std::pair<_vec<double>, _vec<double>> pLine;
    bool isIntersectSegmentLine(pLine&, pLine&);
    extern double EPS;
    _vec<double> crossPoint(geo_basic::pLine &, geo_basic::pLine &);
}

namespace geo_vv{
    double maximum_distance_between_points_brute_force(vector<Vertex>,Organ*);
    double maxmum_distance_between_points_rotating_caliphers(vector<Vertex>,Organ*);
    vector<Vertex*> after_ImageJ_process(vector<Vertex*>);
    vector<Vertex*> normalization(vector<Vertex*>);
    vector<Vertex> normalization(vector<Vertex>);

    double area_vv_boundary(vector<Vertex*>);
    double perimeter_vv_boundary(vector<Vertex*>);

    double vd_minimum(vector<double>);
    double vd_maximum(vector<double>);
    double accumulated_negative(vector<double>);
    vector<double> vd_averaged(vector<double>,int);

    vector<Vertex*> sort_vector_vertex_ascend(vector<Vertex*>);
    vector<Vertex*> sort_vector_vertex_descend(vector<Vertex*>);
    vector<Vertex> sort_vector_vertex_ascend(vector<Vertex>);
    vector<Vertex> sort_vector_vertex_descend(vector<Vertex>);
    bool comp_descend(double,double);
    double linear_fitting(Vertex*,Vertex*,double);
    double linear_fitting(Vertex,Vertex,double);
    vector<Vertex*> vv_x_swap(vector<Vertex*>);
    vector<Vertex> vv_x_swap(vector<Vertex>);
    vector<Vertex*> vector_vertex_sampling(vector<Vertex*> vv,double sampling_distance);
    vector<Vertex> vector_vertex_sampling(vector<Vertex> vv,double sampling_distance);
}

namespace boundary_geo{
    double similarity_Index_1(vector<Vertex*>,vector<Vertex*>,double);
    double similarity_Index_1(vector<Vertex>,vector<Vertex*>,double);
    double similarity_Index_2(string,string,double);
    vector<Vertex*> read_and_process_real_organ_contour_imagej(string);
    double similarity_cal_during_simulation(Organ*,vector<Vertex*>);
}

#endif