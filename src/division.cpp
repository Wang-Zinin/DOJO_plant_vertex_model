/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#define _USE_MATH_DEFINES
#include "../include/division.h"
#include "../include/2dvOrgan.h"
#include "../include/2dvCell.h"
#include "../include/2dvLine.h"
#include "../include/2dvVertex.h"


//const double M_PI =3.1415926;

//biased angles division parameters


namespace division_frequency{

    void no_control(Organ* p_g){
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->frequency_modifier=1.0;
        }
    }

    void balance_control(Organ* p_g){
        //for temporal pattern control
        double BF_tmp = BF;
        F_modifier = BF_tmp + (BF_tmp-1)*(double)p_g->N_epi_cell/(double)p_g->N_inner_cell;
        
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(p_g->p_c[ci]->IsEpidermal==1){
                p_g->p_c[ci]->frequency_modifier=p_g->p_c[ci]->frequency_modifier*F_modifier;
            }
        }
        std::cout<<"BF is "<< BF <<" and current F_modifier is "<<F_modifier<<std::endl;
    }

    void area_control(Organ* p_g){
        
        //sort area for each cell from the smallest to largest

        //preparation for sorting
        std::vector<double>cell_area_sort;
        std::vector<double>cell_area_rank;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            cell_area_sort.push_back(p_g->p_c[ci]->area);
            cell_area_rank.push_back(ci);
        }
        
        
        //get the area rank for each cell
        for(int ci=0;ci<(int)p_g->p_c.size()-1;ci++){
            for(int cj=0;cj<(int)p_g->p_c.size()-ci-1;cj++){
                if(cell_area_sort[cj]>cell_area_sort[cj+1]){
                    //change the elements
                    double sort_temp = cell_area_sort[cj+1];
                    cell_area_sort[cj+1] = cell_area_sort[cj];
                    cell_area_sort[cj] = sort_temp;
                    //change the rank
                    double rank_temp = cell_area_rank[cj+1];
                    cell_area_rank[cj+1] = cell_area_rank[cj];
                    cell_area_rank[cj] = rank_temp;
                }
            }
        }

        for(int i=0; i<(int)cell_area_sort.size();i++){
            p_g->p_c[cell_area_rank[i]]->area_rank=i;
        }

        //assign division frequency modifier to each cell : smallest 30% can not divide; modifier = 5 * (rank - 30%);
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double ci_rank_temp = p_g->p_c[ci]->area_rank/(double)p_g->p_c.size();
            if(ci_rank_temp <=area_control_lower_limit){
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*100;
                p_g->p_c[ci]->area_modifier = 100;
            }
            else{
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*1.0/(area_control_slope*(ci_rank_temp-area_control_lower_limit));
                p_g->p_c[ci]->area_modifier = 1.0/(area_control_slope*(ci_rank_temp-area_control_lower_limit));
            }
            //std::cout<<"p_g->p_c[ci]->area_modifier "<<p_g->p_c[ci]->area_modifier<<endl;
        }

        //for debug
        /*
        std::cout<<"area_control_lower_limit "<<area_control_lower_limit<<endl;
        std::cout<<"area_control_slope "<<area_control_slope<<endl;
        std::cout<<"ci ; cell area ; area_rank; area_modifier; division_frequency"<<std::endl;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            std::cout<<ci<<" ; "<<p_g->p_c[ci]->area<<" ; "<<p_g->p_c[ci]->area_rank<<" ; "<<p_g->p_c[ci]->area_modifier<<" ; "<<1/p_g->p_c[ci]->frequency_modifier<<std::endl;
        } 
        */
        
    }
    
    double calcGaussian(double x, double mu, double sigma){
        return exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
    }
    
    void Gaussian_control(Organ *p_g, double mu,double sigma){
        organ_geo::organ_center(p_g);
        double y_min=p_g->p_c[0]->center.y, y_max=p_g->p_c[0]->center.y;

        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(y_min>p_g->p_c[ci]->center.y){
                y_min = p_g->p_c[ci]->center.y;
            }
            
            if(y_max<p_g->p_c[ci]->center.y){
                y_max = p_g->p_c[ci]->center.y;
            }
        }


        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double y_relative = (p_g->p_c[ci]->center.y-y_min)/(y_max-y_min);
            if(y_relative<0.01){
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*100;
                p_g->p_c[ci]->Gaussian_modifier = 100;
            }
            else{
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*1.0/calcGaussian(y_relative,mu,sigma);
                p_g->p_c[ci]->Gaussian_modifier = 1.0/calcGaussian(y_relative,mu,sigma);
            }
        }
    }
    

}

namespace division_direction{
    _vec<double> angles(Organ* p_g, int ci){
        organ_geo::organ_center(p_g);
        _vec<double> angles_tmp;
        if(p_g->p_c[ci]->IsEpidermal==0){
            if(in_division_direction=="random"){
                division_direction::random(p_g, ci);
            }
            else if(in_division_direction=="constant_0"){
                division_direction::constant_0(p_g,ci);
            }
            else if(in_division_direction=="constant_90"){
                division_direction::constant_90(p_g,ci);
            }
            else if(in_division_direction=="mochizuki_bias"){
                division_direction::mochizuki_bias(p_g,ci,mochizuki_bias_beta,mochizuki_bias_phi);
            }
            
            else{
                std::cout<<"Fatal Error: no inner cell division direction control selected"<<endl;
                exit(-1);
            }
        }
        else if(p_g->p_c[ci]->IsEpidermal==1){
            if(epi_division_direction=="random"){
                division_direction::random(p_g,ci);
            }
            else if(epi_division_direction=="anticlinal"){
                division_direction::epi_anticlinal(p_g,ci);
                //std::cout<<"anticlinal division of cell "<<ci<<endl;
            }
            else if(epi_division_direction=="periclinal"){
                division_direction::epi_periclinal(p_g,ci);
            }
            else if(epi_division_direction=="constant_0"){
                division_direction::constant_0(p_g,ci);
            }
            else if(epi_division_direction=="constant_90"){
                division_direction::constant_90(p_g,ci);
            }
            else{
                std::cout<<"Fatal Error: no epidermal cell division direction control selected"<<endl;
                exit(-1);
            }
        }
        return p_g->p_c[ci]->axis;
    }
    
    //direction: random, anticlinal,periclinal, constant_0 (horizontal), constant_90 (vertical), and yin_distribution, mochizuki bias
    _vec<double> random(Organ* p_g, int ci){
        random_device rnd;
        mt19937 mt(rnd());
        uniform_real_distribution<> rand_axis(0.0,M_PI);
        double axisTheta=rand_axis(mt);
        p_g->p_c[ci]->axisTheta = axisTheta;
        p_g->p_c[ci]->axis = _vec<double>(cos(axisTheta),sin(axisTheta),0.0);

        return p_g->p_c[ci]->axis;
    }

    _vec<double> epi_anticlinal(Organ* p_g, int ci){
        if(p_g->p_c[ci]->IsEpidermal==0){
            std::cout<<"Fatal Error: inner cell is using epidermal anticlinal direction control"<<endl;
            exit(-1);
        }
        _vec<double> center_outermostEdge = (p_g->p_v[p_g->p_c[ci]->surfaceVertex[0]]->loc + p_g->p_v[p_g->p_c[ci]->surfaceVertex[1]]->loc)/2.0;
        p_g->p_c[ci]->axis = center_outermostEdge - p_g->p_c[ci]->center;
        p_g->p_c[ci]->axisTheta = atan(p_g->p_c[ci]->axis.y/p_g->p_c[ci]->axis.x);
        return p_g->p_c[ci]->axis;
    }

    _vec<double> epi_periclinal(Organ* p_g, int ci){
        if(p_g->p_c[ci]->IsEpidermal==0){
            std::cout<<"Fatal Error: inner cell is using epidermal periclinal direction control"<<endl;
            exit(-1);
        }
        p_g->p_c[ci]->axis = p_g->p_v[p_g->p_c[ci]->surfaceVertex[0]]->loc - p_g->p_v[p_g->p_c[ci]->surfaceVertex[1]]->loc;
        p_g->p_c[ci]->axisTheta = atan(p_g->p_c[ci]->axis.y/p_g->p_c[ci]->axis.x);
        return p_g->p_c[ci]->axis;
    }

    _vec<double> constant_0(Organ* p_g,int ci){
        p_g->p_c[ci]->axis = _vec<double>(0.0,1.0,0.0);
        p_g->p_c[ci]->axisTheta = 0.0;
        return p_g->p_c[ci]->axis;
    }

    _vec<double> constant_90(Organ* p_g, int ci){
        p_g->p_c[ci]->axis = _vec<double>(1.0,0.0,0.0);
        p_g->p_c[ci]->axisTheta = M_PI/2.0;
        return p_g->p_c[ci]->axis;
    }

    _vec<double> mochizuki_bias(Organ* p_g, int ci, double bias_beta, double bias_phi){
        double axisTheta = wangMath::mochizuki_bias_single_random_sampling(bias_beta, bias_phi);
        p_g->p_c[ci]->axisTheta = axisTheta;
        p_g->p_c[ci]->axis = _vec<double>(cos(axisTheta),sin(axisTheta),0.0);
        return p_g->p_c[ci]->axis;
    }    
}

namespace division{

    template<typename T>
    void findAndErase(std::vector<T> &vector, T search) {
            auto itr = std::find(vector.begin(), vector.end(), search);
            if(itr == vector.end()) {
            std::cerr << "Couldn't find " << search << std::endl;
            return;
            }
            vector.erase(itr);
        }    

    void Global(Organ *p_g){
        //division::Axis(p_g);
        int cell_division_count=0;
        int inner_cell_division_count=0;
        int epi_cell_division_count=0;
        
        //std::cout<<"Start division frequency judgement ---- ";

        if(division_control=="random_picking"){
            std::random_device rnd;
            std::mt19937 mt(rnd());
            std::uniform_int_distribution<> division_random_picking(0,p_g->p_c.size());

            int ci_dividing = division_random_picking(mt);
            if(p_g->p_c.size()>=end_cell_number){
                goto EndDivision2;
            }
            //std::cout<<"Inner cell "<<ci<<" dividing"<<std::endl;
            division::One(p_g,ci_dividing);
            std::cout<<"Cell "<<ci_dividing<<" divided"<<std::endl;
            cell_division_count =1;
            division::Record(p_g,ci_dividing);
        }
        else{
        //1. prepare for division frequency control
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->frequency_modifier = 1.0;
        }        

        //2. cell division frequency control
        if(division_control=="no"){
            division_frequency::no_control(p_g);
        }
        else if(division_control=="balance"){
            division_frequency::balance_control(p_g);
        }
        else if(division_control=="area"){
            division_frequency::area_control(p_g);
        }
        else if(division_control=="Gaussian_area"){
            division_frequency::area_control(p_g);
            division_frequency::Gaussian_control(p_g,gau_mu,gau_sigma);
        }
        else{
            std::cout<<"Fatal Error: no division frequency control selected"<<endl;
            exit(-1);
        }

        //3. cell division determination
        //std::cout<<"Start cell division determination ---- ";
        
        for(int ci=0;ci<p_g->p_c.size();ci++){
            //3.1 check if end (cell number reached to termination condition)
            if(p_g->p_c.size()>=end_cell_number){
                std::cout<<"It is already "<<end_cell_number<<" cells."<<endl;
                goto EndDivision2;
            }
            //3.2 if the cell time is over threshold, then the cell will divide
            if(p_g->p_c[ci]->cellTime>standard_cell_period_length*p_g->p_c[ci]->frequency_modifier){
                random_device rnd;
                mt19937 mt(rnd());
                uniform_real_distribution<> division_judge(0,1);
                double division_judge_tmp = division_judge(mt);
                if(division_judge_tmp<F_apparent){
                    //3.3 cell division angles judgement 

                    //std::cout<<"Cell division direction judgement ---- "<<endl;
                    division_direction::angles(p_g,ci);
                    //3.4 cell division execution
                    division::One(p_g,ci);
                    //3.5 cell division information std::cout and record
                    cell_division_count++;
                    if(p_g->p_c[ci]->IsEpidermal==1){
                        std::cout<<"Epidermal cell "<<ci<<" divided"<<"; division frequency: "<<1.0/p_g->p_c[ci]->frequency_modifier<<"; division angle: "<<p_g->p_c[ci]->axisTheta<<endl;
                        epi_cell_division_count++;
                    }
                    else{
                        std::cout<<"Inner cell "<<ci<<" divided"<<"; division frequency: "<<1.0/p_g->p_c[ci]->frequency_modifier<<"; division angle: "<<p_g->p_c[ci]->axisTheta<<endl;
                        inner_cell_division_count++;
                    }
                    organ_geo::epidermal_identity(p_g);
                    division::Record(p_g,ci);
                }
            }
        }
        if(cell_division_count==0){
            std::cout<<"No cell divided in this step"<<endl;
        }
        else if(cell_division_count==1){
            std::cout<<"1 cell divided in this step"<<endl;
        }
        else{
            std::cout<<cell_division_count<<" cells divided in this step"<<endl;
        }
        std::cout<<"Now, we have "<<(int)p_g->p_c.size()-cell_division_count<<" + "<<cell_division_count<<" cells: inner cells = "<<p_g->N_inner_cell-inner_cell_division_count<<" + "<<inner_cell_division_count<<"; epidermal cells = "<<p_g->N_epi_cell-epi_cell_division_count<<" + "<<epi_cell_division_count<<endl;
        std::cout<<"Organ area: "<<p_g->area<<"; Organ perimeter: "<<p_g->perimeter<<endl;
        }
        EndDivision2:;
    }

    void One(Organ *p_g, int cellula_idx){
    //this function will add two new vertex, three new lines and one new cell, four cells need to be changed
    Cell *cp = p_g->p_c[cellula_idx];
    
    _vec<double> axis = p_g->p_c[cellula_idx]->axis; //from division::Axis
    
    _vec<double> center = cp->center; //from organ_geo::organ_center
    
    geo_basic::pLine division = std::make_pair(center, center+axis);
    std::vector<std::pair<int, _vec<double>>> crosspoint;
    std::vector<int> epoint_idx;
    
    //std::cout<<"Cell center is "<<center.x<<","<<center.y<<";Cell axisTheta is"<<p_g->p_c[cellula_idx]->axisTheta<<";Cell axis is "<<axis.x<<","<<axis.y<<std::endl;

    // sort lines counterclockwise, based on counterclockwise vertex 
    std::vector<int> anticlockwise_lidx;
    anticlockwise_lidx = cell_geo::cell_counterclock_line(p_g,cellula_idx);
    int cp_vsize = cp->vi.size();
    
        //counterclockwise sorted lines will
        for(int lidx: anticlockwise_lidx){
        int tmp_epoint_idx1=p_g->p_l[lidx]->vi[0];
        int tmp_epoint_idx2=p_g->p_l[lidx]->vi[1];

        geo_basic::pLine segment = std::make_pair(p_g->p_v[tmp_epoint_idx1]->loc,p_g->p_v[tmp_epoint_idx2]->loc);
        if(geo_basic::isIntersectSegmentLine(segment,division)){
            epoint_idx.push_back(tmp_epoint_idx1);
            epoint_idx.push_back(tmp_epoint_idx2);

            crosspoint.push_back(std::make_pair(lidx, geo_basic::crossPoint(segment, division)));
        }
        }

        if(crosspoint.size() != 2) {
            //std::cout<<"Fatal bug of cross point"<<std::endl;
            //p_g->p_c[cellula_idx]->cellDivisionCount=100;
            //throw "Cannot divide cellulas";
            return;
        }
            p_g->p_c[cellula_idx]->cellDivisionCount++;
            p_g->p_c[cellula_idx]->cellTime =0.0;
            
        int lim_num1 = -1, lim_num2 = -1;

        for(int i = 0; i < cp_vsize; ++i){
            int vidx = cp->vi[i];
            if(vidx == epoint_idx[1]){
            lim_num1 = i;
            }
            if(vidx == epoint_idx[2]){
            lim_num2 = i-1;
            }
        }
        assert(lim_num1 != -1);
        assert(lim_num2 != -1);

        int v1_idx = p_g->p_v.size();
        int v2_idx = p_g->p_v.size() + 1;

        int l1_idx = p_g->p_l.size();
        int l2_idx = p_g->p_l.size() + 1;
        int l3_idx = p_g->p_l.size() + 2;

        int c1_idx = p_g->p_c.size();

        // vertex v1
        Vertex *v1 = new Vertex;
        v1->li.push_back(crosspoint[0].first);
        v1->li.push_back(l1_idx);
        v1->li.push_back(l3_idx);
        for(int c: p_g->p_l[crosspoint[0].first]->ci) {
            v1->ci.push_back(c);
        }
        v1->ci.push_back(c1_idx);
        v1->loc = crosspoint[0].second;
        v1->loc.z = 0.0;
        v1->vi = p_g->p_v.size();
        //v1->loc[1] = ???

        // vertex v2
        Vertex *v2 = new Vertex;
        v2->li.push_back(crosspoint[1].first);
        v2->li.push_back(l2_idx);
        v2->li.push_back(l3_idx);
        for(int c: p_g->p_l[crosspoint[1].first]->ci) {
            v2->ci.push_back(c);
        }
        v2->ci.push_back(c1_idx);
        v2->loc = crosspoint[1].second;
        v2->loc.z = 0.0;
        v2->vi = p_g->p_v.size()+1;

        //v2->loc[0][1] = ???

        // vertex epoint_idx[1]
        p_g->p_v[epoint_idx[1]]->li.push_back(l1_idx);
        division::findAndErase(p_g->p_v[epoint_idx[1]]->li, crosspoint[0].first);
        p_g->p_v[epoint_idx[1]]->ci.push_back(c1_idx);
        division::findAndErase(p_g->p_v[epoint_idx[1]]->ci, cellula_idx);

        // vertex epoint_idx[2]
        p_g->p_v[epoint_idx[2]]->li.push_back(l2_idx);
        division::findAndErase(p_g->p_v[epoint_idx[2]]->li, crosspoint[1].first);
        p_g->p_v[epoint_idx[2]]->ci.push_back(c1_idx);
        division::findAndErase(p_g->p_v[epoint_idx[2]]->ci, cellula_idx);

        for(int i = lim_num1 + 1; i <= lim_num2; ++i) {
            int vidx = cp->vi[i];
            p_g->p_v[vidx]->ci.push_back(c1_idx);
            division::findAndErase(p_g->p_v[vidx]->ci, cellula_idx);
        }

        // line l1
        Line *l1 = new Line;
        l1->vi[0] = v1_idx;
        l1->vi[1] = epoint_idx[1];
        for(int cidx: p_g->p_l[crosspoint[0].first]->ci) {
            l1->ci.push_back(cidx);
        }
        l1->ci.push_back(c1_idx);
        division::findAndErase(l1->ci, cellula_idx);
        l1->edgeForce = 0.0;
        l1->li = p_g->p_l.size();
        //l1->Balanced_Length = balanced_length;
        //l1->K_LENGTH = K_LENGTH;
        //l1->K2_LENGTH = K2_LENGTH;
        //l1->LENGTH_EQ = LENGTH_EQ;

        // line l2
        Line *l2 = new Line;
        l2->vi[0] = v2_idx;
        l2->vi[1] = epoint_idx[2];
        for(int cidx: p_g->p_l[crosspoint[1].first]->ci) {
            l2->ci.push_back(cidx);
        }
        l2->ci.push_back(c1_idx);
        division::findAndErase(l2->ci, cellula_idx);
        l2->edgeForce = 0.0;
        l2->li = p_g->p_l.size()+1;

        //l2->Balanced_Length = balanced_length;
        //l2->K_LENGTH = K_LENGTH;
        //l2->K2_LENGTH = K2_LENGTH;
        //l2->LENGTH_EQ = LENGTH_EQ;
        // line l3
        Line *l3 = new Line;
        l3->vi[0] = v1_idx;
        l3->vi[1] = v2_idx;
        l3->ci.push_back(cellula_idx);
        l3->ci.push_back(c1_idx);
        l3->edgeForce = 0.0;
        l3->li = p_g->p_l.size()+2;
        //l3->Balanced_Length = balanced_length;
        //l3->K_LENGTH = K_LENGTH;
        //l3->K2_LENGTH = K2_LENGTH;
        //l3->LENGTH_EQ = LENGTH_EQ;

        // line crosspoint[0].first
        p_g->p_l[crosspoint[0].first]->vi[1] = v1_idx;

        // line crosspoint[1].first
        p_g->p_l[crosspoint[1].first]->vi[0] = v2_idx;

        // line lim_num1-th --- lim_num2-th
        for(int i = lim_num1; i <= lim_num2; ++i) {
            int lidx = anticlockwise_lidx[i];
            p_g->p_l[lidx]->ci.push_back(c1_idx);
            division::findAndErase(p_g->p_l[lidx]->ci, cellula_idx);
        }

        // cellula c1
        Cell *c1 = new Cell;
        c1->vi.push_back(v1_idx);
        for(int i = lim_num1; i <= lim_num2 + 1; ++i) {
            int vidx = cp->vi[i];
            c1->vi.push_back(vidx);
        }
        c1->vi.push_back(v2_idx);
        c1->li.push_back(l1_idx);
        for(int i = lim_num1; i <= lim_num2; ++i) {
            int lidx = anticlockwise_lidx[i];
            c1->li.push_back(lidx);
        }
        c1->li.push_back(l2_idx);
        c1->li.push_back(l3_idx);
        c1->center = cp->center;
        
        //c1->Balanced_Area = balanced_area;
        //c1->K_AREA = K_AREA;

        c1->cellDivisionCount = p_g->p_c[cellula_idx]->cellDivisionCount;
        c1->tag = p_g->p_c[cellula_idx]->tag;
        c1->cellTime =0;
        // cellula cellula_idx
        cp->vi.erase(cp->vi.begin() + lim_num1, cp->vi.begin() + lim_num2 + 2);
        cp->vi.push_back(v1_idx);
        cp->vi.push_back(v2_idx);
        for(int i = lim_num1; i <= lim_num2; ++i) {
            //std::cerr << anticlockwise_lidx[i] << std::endl;
            division::findAndErase(cp->li, anticlockwise_lidx[i]);
        }
        cp->li.push_back(l3_idx);

        // cellula adjusting line-crosspoint[0].first (if exist)
        int adj_cidx1 = -1;
        for(int cidx: p_g->p_l[crosspoint[0].first]->ci) {
            if(cidx != cellula_idx) {
            adj_cidx1 = cidx;
            }
        }
        if(adj_cidx1 != -1) {
            p_g->p_c[adj_cidx1]->vi.push_back(v1_idx);
            p_g->p_c[adj_cidx1]->li.push_back(l1_idx);
        }

        // cellula adjusting line-crosspoint[1].first (if exist)
        int adj_cidx2 = -1;
        for(int cidx: p_g->p_l[crosspoint[1].first]->ci) {
            if(cidx != cellula_idx) {
            adj_cidx2 = cidx;
            }
        }
        if(adj_cidx2 != -1) {
            p_g->p_c[adj_cidx2]->vi.push_back(v2_idx);
            p_g->p_c[adj_cidx2]->li.push_back(l2_idx);
        }

        p_g->p_c.push_back(c1);
        p_g->p_l.push_back(l1);
        p_g->p_l.push_back(l2);
        p_g->p_l.push_back(l3);
        p_g->p_v.push_back(v1);
        p_g->p_v.push_back(v2);

        //isConsistent(p_g);
            //std::cout<<"hello 1"<<endl;
        //sortAntiClockwise(p_g);
        organ_geo::organ_vertex_counterclockwise_sort(p_g);
            //std::cout<<"hello 2"<<endl;

    }

    void Record(Organ *p_g, int cidx){
        DivisionRecord *dr = new DivisionRecord;
        dr->time = p_g->step;
        dr->cidx = cidx;
        dr->IsEpidermal = p_g->p_c[cidx]->IsEpidermal;
        dr->axisTheta = p_g->p_c[cidx]->axisTheta; 
        dr->center_x = p_g->p_c[cidx]->center.x;
        dr->center_y = p_g->p_c[cidx]->center.y;
        dr->division_count = p_g->p_c[cidx]->cellDivisionCount;
        //record frequency_modifier
        int inner_cell_number=0, peripheral_cell_number=0;
        double sum_in_frequency_modifier=0, sum_epi_frequency_modifier=0;
        for(int ci=0; ci<p_g->p_c.size(); ci++){
            if(p_g->p_c[ci]->IsEpidermal==0){
                sum_in_frequency_modifier+=p_g->p_c[ci]->frequency_modifier;
                inner_cell_number++;
            }
            else{
                sum_epi_frequency_modifier+=p_g->p_c[ci]->frequency_modifier;
                peripheral_cell_number++;
            }
        }
        double av_in_frequency_modifier = sum_in_frequency_modifier/inner_cell_number;
        double av_epi_frequency_modifier = sum_epi_frequency_modifier/peripheral_cell_number;
        dr->av_in_frequency_modifier = av_in_frequency_modifier;
        dr->av_epi_frequency_modifier = av_epi_frequency_modifier;
        p_g->d_r.push_back(dr);
    }

    void cell_time_initialization(Organ* p_g){
        if(division_control == "random_picking"){
    }
    else{
        if(division_control=="area"){
                division_frequency::area_control(p_g);
            }
            else if(division_control=="no"){
                division_frequency::no_control(p_g);
            }
            else if(division_control=="Gaussian_balance"){
                division_frequency::balance_control(p_g);
                division_frequency::Gaussian_control(p_g,gau_mu,gau_sigma);
            }
            else if(division_control=="Gaussian_area"){
                division_frequency::area_control(p_g);
                division_frequency::Gaussian_control(p_g,gau_mu,gau_sigma);
            }
            else if(division_control=="balance"){
                division_frequency::balance_control(p_g);
            }
            else{
                std::cout<<"Fatal error: no cell division control is selected! (initialization)"<<endl;
                exit(-1);
            }

        random_device rnd;
        mt19937 mt(rnd());
        uniform_real_distribution <> cell_time_init(0,standard_cell_period_length);
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->cellTime = cell_time_init(mt)*p_g->p_c[ci]->frequency_modifier;
        }
    }
    }

    //cell time will be added for each cell
    void cellTimeAdd(Organ *p_g, int increase_index){
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->cellTime = p_g->p_c[ci]->cellTime+increase_index;
        }
    }

}