/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#include "../include/IOV.h"
#include "../include/2dvOrgan.h"
#include "../include/2dvCell.h"
#include "../include/2dvLine.h"
#include "../include/2dvVertex.h"

namespace readV{

    //read_organ_txt 
    //reads organ information (vertex position, line connectivity and cell connectivity) from the designated "initialOrganFile".txt file
    //the "initialOrganFile" requires a specified format: 
    //description_of_file
    //number_of_vertex
    //vertex_index, x,y,z position of vertex vi
    //i x_i y_i z_i
    //number_of_line
    //line index, the two vertex indices connected by line lj
    //j l(j,1) l(j,2)
    //number of hexagon
    //hegaxon index, the six line indices and vertex indices (vertex indices must be in anticlockwise direction)
    //k cl(k,1) cl(k,2) cl(k,3) cl(k,4) cl(k,5) cl(k,6)
    //k cv(k,1) cv(k,2) cv(k,3) cv(k,4) cv(k,5) cv(k,6)
    //number of pentagon
    //pentagon index, the six line indices and vertex indices (vertex indices must be in anticlockwise direction)
    //k cl(k,1) cl(k,2) cl(k,3) cl(k,4) cl(k,5)
    //k cv(k,1) cv(k,2) cv(k,3) cv(k,4) cv(k,5)
    //number of quandrangle
    //quandrangle index, the six line indices and vertex indices (vertex indices must be in anticlockwise direction)
    //k cl(k,1) cl(k,2) cl(k,3) cl(k,4)
    //k cv(k,1) cv(k,2) cv(k,3) cv(k,4)
    void read_organ_txt(Organ* p_g,int fixed_index){
        cout<<"*************|Initial Shape (primordia) Information ("<<initialOrganFile<<") |***************"<<endl;
        //extern string initialOrganFile;
        //cout<<"initialOrganFile: "<<initialOrganFile<<endl;
        ifstream fin(initialOrganFile,ios::in);
        if(!fin.is_open()){
            initialOrganFile = "../"+initialOrganFile;
            ifstream fin(initialOrganFile,ios::in);
        }
        if(!fin.is_open()){
            cout<<"Error: missing "<<initialOrganFile <<" (the initial state file!)"<<endl;        
            exit(1);
        }
        string description;
        int vertex_number,line_number,hexagon_number,pentagon_number,quadrangle_number;
        fin>>description;
        fin>>vertex_number;

        std::cout<<"Description: "<<description<<std::endl;
        std::cout<<"Number of vertices in initial organ: "<<vertex_number<<std::endl;

        for(int i=0; i<vertex_number;i++)
        {   
            Vertex *tmp = new Vertex;
            int vertex_index;
            //fin>>vertex_index;
            fin>>tmp->loc.x;
            fin>>tmp->loc.y;
            fin>>tmp->loc.z;
            tmp->loc.z =0.0;
            tmp->vi=i;

            p_g->p_v.push_back(tmp);
        }


        //Debug: use initialOrgan to check the force amd motion calculation
        /*
        std::cout<<"Debug: use initialOrgan to check the force amd motion calculation"<<std::endl;
        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
        p_g->p_v[vi]->loc[0].x = 1.1*p_g->p_v[vi]->loc[0].x;
        p_g->p_v[vi]->loc[0].y = 1.1*p_g->p_v[vi]->loc[0].y;
        }
        */
        //debug: check vertex input
        /*
        std::cout<<"Debug: check vertex input"<<std::endl;
        for(int i=0; i<vertex_number;i++)
        {
            std::cout<<i<<" "<<p_g->p_v[i]->loc[0].x<<" "<<p_g->p_v[i]->loc[0].y<<std::endl;
        }
        std::cout<<"End Debug: check vertex input" <<std::endl;
        */

        fin>>line_number;
        if(fixed_index==0)
    std::cout<<"Number of lines in initial state: "<<line_number<<std::endl;

    for(int i=0; i<line_number;i++)
    {
        int line_index;
        Line *tmp = new Line;
        
        fin>>line_index;
        fin>>tmp->vi[0];
        fin>>tmp->vi[1];
        tmp->li=i;

        p_g->p_l.push_back(tmp);
    }

        //debug: check line input
        /*
        std::cout<<"Debug: check line input"<<std::endl;
        for(int i=0; i<line_number;i++){
            std::cout<<i<<" "<<p_g->p_l[i]->vi[0]<<" "<<p_g->p_l[i]->vi[1]<<std::endl;
        }
        std::cout<<"End Debug: check line input"<<std::endl;
        */

        for(int i=0; i<line_number;i++){
        p_g->p_v[p_g->p_l[i]->vi[0]]->li.push_back(i);
        p_g->p_v[p_g->p_l[i]->vi[1]]->li.push_back(i);
    }

    //debug: check which lines vertex belongs to
    /*
    std::cout<<"Debug: check which lines vertex belongs to."<<std::endl;
    for(int i=0; i<vertex_number; i++){
        std::cout<<"For vertex "<<i<<", it belongs to line ";
        for(int j=0; j<p_g->p_v[i]->li.size();j++){
            std::cout<<p_g->p_v[i]->li[j]<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<"End Debug: check which lines vertex belongs to."<<std::endl;
        */


        fin>>hexagon_number;
        
        
        for(int i=0; i<hexagon_number; i++){
            Cell *tmp = new Cell;
            int hexagon_index,line_index,vertex_index;
            fin>>hexagon_index;
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);

            fin>>hexagon_index;
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            
            tmp->n_edges=6;
            p_g->p_c.push_back(tmp);       
        }

        fin>>pentagon_number;
        
        for(int i=0; i<pentagon_number; i++){
            Cell *tmp = new Cell;
            int pentagon_index, line_index, vertex_index;
            fin>>pentagon_index;
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);

            fin>>pentagon_index;
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);

            tmp->n_edges=5;
            p_g->p_c.push_back(tmp);
        }
        fin>>quadrangle_number;
        int cell_number = hexagon_number+pentagon_number+quadrangle_number;
        
        if(fixed_index==0){
            std::cout<<"Number of hexagons in initial state : "<<hexagon_number<<std::endl;
            std::cout<<"Number of pentagons in initial state : "<<pentagon_number<<std::endl;
            std::cout<<"Number of quadrangles in initial state : "<<quadrangle_number<<std::endl;
            std::cout<<"Number of all cells in initial state : "<<cell_number<<std::endl;
        }
        for(int i=0; i<quadrangle_number; i++){
            Cell *tmp = new Cell;
            int quadrangle_index,  line_index, vertex_index;
            fin>>quadrangle_index;
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            

            fin>>quadrangle_index;
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            
            tmp->n_edges=4;
            p_g->p_c.push_back(tmp);
        }

        //Debug: check the vertex index and line index in each cell
        /*
        std::cout<<"Debug: check the vertex index and line index in each cell"<<std::endl;
        for(int i=0; i<p_g->p_c.size(); i++){
            std::cout<<"Cell "<<i<<" has line ";
            for(int j=0; j<p_g->p_c[i]->n_edges; j++){
                std::cout<<p_g->p_c[i]->li[j]<<" ";
            }
            std::cout<<std::endl;
            
            std::cout<<"Cell "<<i<<" has vertex ";
            for(int j=0; j<p_g->p_c[i]->n_edges; j++){
                std::cout<<p_g->p_c[i]->vi[j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<"End Debug: check the vertex index and line index in each cell"<<std::endl;
        */

        for(int i=0; i<cell_number; i++){
            for(int j=0; j<p_g->p_c[i]->n_edges;j++){
                p_g->p_l[p_g->p_c[i]->li[j]]->ci.push_back(i);
                p_g->p_v[p_g->p_c[i]->vi[j]]->ci.push_back(i);
            }
        }

        //initialize the cellDivisionCount
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->cellDivisionCount=0;
        }

        //debug: check which cells vertex belongs to
        /*
        std::cout<<"Debug: check which cells each vertex belongs to."<<std::endl;
        for(int i=0; i<vertex_number; i++){
        std::cout<<"For vertex "<<i<<", it belongs to cell ";
        for(int j=0; j<p_g->p_v[i]->ci.size();j++){
            std::cout<<p_g->p_v[i]->ci[j]<<" ";
        }
        std::cout<<std::endl;
        }
        std::cout<<"End Debug: check which cells each vertex belongs to."<<std::endl;
        */

        //debug: check which cells line belongs to
        /*
        std::cout<<"Debug: check which cells each line belongs to."<<std::endl;
        for(int i=0; i<line_number; i++){
        std::cout<<"For line "<<i<<", it belongs to cell ";
        for(int j=0; j<p_g->p_l[i]->ci.size();j++){
            std::cout<<p_g->p_l[i]->ci[j]<<" ";
        }
        std::cout<<std::endl;
        }
        std::cout<<"End Debug: check which cells each line belongs to."<<std::endl;
        */
       p_g->initial_cell_number = p_g->p_c.size();
       fin.close();
       cout<<"********************************************************************"<<endl;
    }

//read vertex position and line-vertex, cell-vertex connectivity from cell.vtk and line.vtk
    void oneCell(Organ* p_g,string FileName){
        ifstream fin(FileName,ios::in);
        if(!fin.is_open()){
            cout<<"Error: missing cell.vtk ("<<FileName<<")"<<endl;
            exit(1);
        }
        char parameter_name[100];
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        int vertex_number;
        fin>>vertex_number;
        fin>>parameter_name;

        //read vertex position information
        for(int i=0; i<vertex_number;i++){
            Vertex *v_tmp = new Vertex;
            fin>>v_tmp->loc.x;
            fin>>v_tmp->loc.y;
            fin>>v_tmp->loc.z;
            p_g->p_v.push_back(v_tmp);
        }
        //debug: vertex position information
        /*
        for(int i=0; i<vertex_number;i++){
            std::cout<<"Vertex "<< i<<" position: ("<<p_g->p_v[i]->loc.x<<","<<p_g->p_v[i]->loc.y<<","<<p_g->p_v[i]->loc.z<<")."<<std::endl;
        }
        */
    
        //read cell-vertex connection informtation
        fin>>parameter_name;
        int cell_number;
        fin>>cell_number;
        fin>>parameter_name;
        for(int ci=0;ci<cell_number;ci++){
            Cell *c_tmp = new Cell;
            int cell_vertex_number;
            fin>>cell_vertex_number;
            for(int vi=0;vi<cell_vertex_number;vi++){
                int vertex_index;
                fin>>vertex_index;
                c_tmp->vi.push_back(vertex_index);
            }
            p_g->p_c.push_back(c_tmp);
        }
        fin.close();
    }

    void oneLine(Organ* p_g,string FileName){
        ifstream fin(FileName,ios::in);
            if(!fin.is_open()){
            cout<<"Error: missing Line.vtk ("<<FileName<<")"<<endl;
                return;
            }
            char parameter_name[100];
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            int vertex_number;
            fin>>vertex_number;
            fin>>parameter_name;

            //read vertex position information
            for(int i=0; i<vertex_number;i++){
                Vertex *v_tmp = new Vertex;
                fin>>v_tmp->loc.x;
                fin>>v_tmp->loc.y;
                fin>>v_tmp->loc.z;
                //p_g->p_v.push_back(v_tmp);
            }

            //read line-vertex connectivity information
            fin>>parameter_name;
            int line_number;
            fin>>line_number;
            fin>>parameter_name;
            for(int li=0;li<line_number;li++){
                int line_index;
                fin>>line_index;
                Line *l_tmp = new Line;
                fin>>l_tmp->vi[0];
                fin>>l_tmp->vi[1];
                p_g->p_l.push_back(l_tmp);
            }

            //debug: line-vertex connectivity information
            /*
            for(int li=0;li<(int)p_g->p_l.size();li++){
                std::cout<<li<<" "<<p_g->p_l[li]->vi[0]<<" "<<p_g->p_l[li]->vi[1]<<std::endl;
            }
            */
            fin.close();
    }

    //line-vertex, cell-vertex connectivity => vertex-cell, vertex-line, cell-line, line-cell connectivity
    void oneCellLine(Organ* p_g){
        //add line index to cell
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            for(int vi=0;vi<(int)p_g->p_c[ci]->vi.size();vi++){
                for(int li=0; li<(int)p_g->p_l.size();li++){
                    if(p_g->p_l[li]->vi[0]==p_g->p_c[ci]->vi[vi]){
                        for(int vii=0;vii<(int)p_g->p_c[ci]->vi.size();vii++){
                            if(p_g->p_l[li]->vi[1]==p_g->p_c[ci]->vi[vii]){
                                p_g->p_c[ci]->li.push_back(li);
                            }
                        }
                    }
                }
            }
        }

        //debug: cell-line connnectivity information
        /*
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            std::cout<<ci;
            for(int li=0;li<(int)p_g->p_c[ci]->li.size();li++){
                std::cout<<" "<<p_g->p_c[ci]->li[li];
            }
            std::cout<<std::endl;
        }
        */

        //add cell index to line and vertex
        for(int i=0; i<(int)p_g->p_c.size(); i++){
            for(int j=0; j<p_g->p_c[i]->vi.size();j++){
                p_g->p_l[p_g->p_c[i]->li[j]]->ci.push_back(i);
                p_g->p_v[p_g->p_c[i]->vi[j]]->ci.push_back(i);
            }
        }

        //add line index to vertex 
        for(int i=0; i<(int)p_g->p_l.size();i++){
            //a line can only connect two vertices
            p_g->p_v[p_g->p_l[i]->vi[0]]->li.push_back(i);
            p_g->p_v[p_g->p_l[i]->vi[1]]->li.push_back(i);
        }
    }

    //all information of an Organ in one time step
    void oneVTK(Organ* p_g,string CellVTK,string LineVTK){
        oneCell(p_g,CellVTK);
        oneLine(p_g,LineVTK);
        oneCellLine(p_g);
    }

    //read the final VTK inside a filefolder
    void final_VTK(Organ* p_g,string VTK_file_path){
        int length_tmp = VTK_file_path.length();
        char directory_tmp[length_tmp+1];
        strcpy(directory_tmp, VTK_file_path.c_str());
        
        chdir(directory_tmp);
        int file_index = file_process::getFileNum(VTK_file_path);
        cout<<"file_index: "<<file_index<<endl;
        if(file_index==0){
            cout<<"Error: "<<directory_tmp <<"does not contain any files"<<endl;
            exit(1);
        }
        //cout<<"file_index "<<file_index<<endl;
        char str_cell[100], str_line[100];
        sprintf(str_cell,"2dv_cell%.5d00000.vtk",file_index/2-2);
        sprintf(str_line,"2dv_line%.5d00000.vtk",file_index/2-2);
        cout<<"Reading "<<str_cell<<" and "<<str_line<<endl;
        string CellVTK = VTK_file_path+"/"+str_cell;
        string LineVTK = VTK_file_path+"/"+str_line;

        readV::oneVTK(p_g, CellVTK, LineVTK);
    }

    //read the outline information from imageJ extracted xy txt file
    vector<Vertex*> xy_txt_to_vertexXY(string FileName){
        vector<Vertex*> vv;
        ifstream fin(FileName, ios::in);
        if(!fin.is_open()){
            cout<<"Error: missing "<<FileName<<endl;
            exit(1);
        }

        while(!fin.eof()){
            Vertex * v_tmp = new Vertex;
            int tmp;
            fin>>tmp;
            fin>>v_tmp->loc.x;
            fin>>v_tmp->loc.y;
            v_tmp->loc.z=0;
            vv.push_back(v_tmp);
        }
        return vv;
    }

    vector<double> read_vd(string FileName){
        vector<double> vd;
        string nothing;
        
        ifstream fin(FileName, ios::in);
        if(!fin.is_open()){            
            cout<<"Error: missing "<<FileName<<endl;
            exit(1);
        }
        fin>>nothing;
        fin>>nothing;
        while(!fin.eof()){
            double tmp;
            fin>>tmp;
            fin>>tmp;
            vd.push_back(tmp);
        }
        return vd;
    }

    vector<Vertex*> read_vv(string FileName){
        vector<Vertex*> vv;
        ifstream fin(FileName, ios::in);
        if(!fin.is_open()){            
            cout<<"Error: missing "<<FileName<<endl;
            exit(1);
        }
        string s_tmp;
        fin>>s_tmp;
        fin>>s_tmp;
        fin>>s_tmp;
        fin>>s_tmp;
        while(!fin.eof()){
            double tmp;
            Vertex* v_tmp = new Vertex;
            fin>>tmp;
            fin>>tmp;
            v_tmp->loc.x = tmp;
            fin>>tmp;
            v_tmp->loc.y = tmp;
            fin>>tmp;
            v_tmp->loc.z = tmp;
            vv.push_back(v_tmp);
        }
        if(vv[vv.size()-1]->loc.x==0&&vv[vv.size()-1]->loc.y==0)
            vv.pop_back();

        return vv;
    }
/*
    Geometrics_analysis_class read_geo_output_txt(string FileName){
        Geometrics_analysis_class ga;
        ifstream fin(FileName, ios::in);
        if(!fin.is_open()){
            cout<<"Error: missing "<<FileName<<endl;
            exit(1);
        }
        
    }
*/
}

namespace output{

void VTK(Organ* p_g){
    {
        char fname[100];
        sprintf(fname,"2dv_line%010u.vtk", p_g->step);
        std::ofstream fout(fname);

        //vtkファイルのヘッダ
        fout << "# vtk DataFile Version 3.0" << std::endl;
        fout << "2D-vertex" << std::endl;
        fout << "ASCII" << std::endl;
        fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
        fout<<std::endl;
        //点の位置を書き込む
        fout << "POINTS " << p_g->p_v.size() << " float" << std::endl;
        for (int i = 0; i < (int)p_g->p_v.size(); i++) {
        Vertex *vp = p_g->p_v[i];
        fout << vp->loc.x << " ";
        fout << vp->loc.y << " ";
        fout << vp->loc.z << std::endl;
        }
        fout<<std::endl;
        //線の要素を書き込む
        fout << "CELLS " << p_g->p_l.size() << " " << p_g->p_l.size() * 3 << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); i++) {
        Line *lp = p_g->p_l[i];
        fout << "2 " << lp->vi[0] << " " << lp->vi[1] << std::endl;
        }
        fout<<std::endl;
        //CELL_TYPES
        fout << "CELL_TYPES " << p_g->p_l.size() << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); i++) {
        fout << "3" << std::endl;
        }
        fout<<std::endl;

        fout << "CELL_DATA " <<  p_g->p_l.size() << std::endl;
        fout << "SCALARS current_Length float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); ++i) {
            fout <<p_g->p_l[i]->length<<std::endl;
        }
        fout<<std::endl;
        fout << "SCALARS edge_force float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); ++i) {
            fout <<p_g->p_l[i]->edgeForce<<std::endl;
        }
        fout.close();
    }

    {
        char fname[100];
        sprintf(fname,"2dv_cell%010u.vtk", p_g->step);
        std::ofstream fout(fname);

        //vtkファイルのヘッダ
        fout << "# vtk DataFile Version 2.0" << std::endl;
        fout << "2D-vertex" << std::endl;
        fout << "ASCII" << std::endl;
        fout << "DATASET UNSTRUCTURED_GRID" << std::endl;

        //点の位置を書き込む
        fout << "POINTS " << p_g->p_v.size() << " float" << std::endl;
        for (int i = 0; i < (int)p_g->p_v.size(); i++) {
        Vertex *vp = p_g->p_v[i];
        fout << vp->loc.x << " ";
        fout << vp->loc.y << " ";
        fout << vp->loc.z << std::endl;
        }
        int cells_size = 0;

        for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
        Cell *cp = p_g->p_c[i];
        cells_size += cp->vi.size() + 1;
        }

        fout << "CELLS " << p_g->p_c.size() << " " << cells_size << std::endl;
        for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
        Cell *cp = p_g->p_c[i];
        fout << cp->vi.size();
        for (int j = 0; j < (int)cp->vi.size(); j++) {
            fout << " " << cp->vi[j];
        }
        fout << std::endl;
        }
        //CELL_TYPES
        fout << "CELL_TYPES " << p_g->p_c.size() << std::endl;
        for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
        fout << "7" << std::endl;
        }
        
        //Field 
        fout << "CELL_DATA " <<  p_g->p_c.size() << std::endl;
        //output the cell division count
        fout << "SCALARS cell_division_count float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->cellDivisionCount << std::endl;
        }
        fout << "SCALARS cell_time float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->cellTime << std::endl;
        }
        
        fout << "SCALARS cell_area float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->area << std::endl;
        }
        fout << "SCALARS cell_perimeter float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->perimeter << std::endl;
        }
        fout << "SCALARS cell_frequency_modifier float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->frequency_modifier<< std::endl;
        }
        if(division_control=="Gaussian_area"||division_control=="area"){
            fout << "SCALARS cell_area_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->area_modifier<< std::endl;
            }
        }
        if(division_control=="Gaussian_area"){
            fout << "SCALARS cell_Gaussian_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->Gaussian_modifier<< std::endl;
            }
        }
        fout.close();
    }
}

void geo_initial(void){
    string fname = "geometrics_record.txt";
    std::ofstream fout(fname);
    fout<<"Step inner_cell_number epidermal_cell_number center.x center.y cell_layer_number area epiArea inArea epiArea_averaged inArea_averaged area_averaged perimeter perimeter_averaged circularity similarity_index"<<endl;
    fout.close();
}

void geo(Organ* p_g){
    string fname = "geometrics_record.txt";
    ofstream outfile(fname, std::ios::app);
    if(!outfile){
        std::cerr<<"Failed to open "<<fname<<std::endl;
        exit(1);
    }
    outfile<<p_g->step<<" "<<p_g->N_inner_cell<<" "<<p_g->N_epi_cell<<" ";
    outfile<<p_g->center.x<<" "<<p_g->center.y<<" ";
    outfile<<p_g->cell_layer_number<<" ";
    outfile<<p_g->area<<" "<<p_g->epiArea<<" "<<p_g->inArea<<" "<<p_g->epiArea_averaged<<" "<<p_g->inArea_averaged<<" "<<p_g->area_averaged<<" ";
    outfile<<p_g->perimeter<<" "<<p_g->perimeter_averaged<<" ";
    outfile<<p_g->circularity<<" ";
    outfile<<p_g->similarity_index<<" ";
    outfile<<endl;
    outfile.close();
}


void division(Organ *p_g){
    //char fname[100] ="divisionRecord.txt";
    char fname[100];
    sprintf(fname,"divisionRecord.txt");
    std::ofstream fout(fname);
    fout<<"division_index "<<"time_step "<<"divided_cell_index "<<"divided_cell_IsEpidermal "<<"division_axisTheta "<<"divided_cell_x "<<"divided_cell_y "<<"cell_division_count "<<"av_in_frequency_modifier "<<"av_epi_frequency_modifier "<<std::endl;
    for(int di=0;di<(int)p_g->d_r.size();di++){
        fout<<di<<" "<<p_g->d_r[di]->time<<" "<<p_g->d_r[di]->cidx<<" "<<p_g->d_r[di]->IsEpidermal<<" "<<p_g->d_r[di]->axisTheta<<" "<<p_g->d_r[di]->center_x<<" "<<p_g->d_r[di]->center_y<<" "<<p_g->d_r[di]->division_count<<" "<<p_g->d_r[di]->av_in_frequency_modifier<<" "<<p_g->d_r[di]->av_epi_frequency_modifier<<std::endl;
    }
    fout.close();
}

double averaged_vec_double(vector<double> vec){
    double sum = 0;
    for(int i=0;i<vec.size();i++){
        sum+=vec[i];
    }
    return sum/vec.size();
};

double averaged_vec_int(vector<int> vec){
    double sum = 0;
    for(int i=0;i<vec.size();i++){
        sum+=vec[i];
    }
    return sum/vec.size();
};

void batch(Batch *ba){
    char fname[100];
    sprintf(fname,"batch_information.txt");
    std::ofstream fout(fname);

    fout<<"simulation_index sigma_O sigma_L sigma_O/sigma_L N_in N_epi organ_area organ_perimeter circularity regularity_in regularity_epi overlap_index"<<endl;
    
    for(int i=0;i<(int)ba->N_in.size();i++){
        fout<<i<<" "<<sigma_O<<" "<<sigma_L<<" "<<sigma_O/sigma_L<<" "<<ba->N_in[i]<<" "<<ba->N_epi[i]<<" "<<ba->organ_area[i]<<" "<<ba->organ_perimeter[i]<<" "<<4*3.1415926*ba->organ_area[i]/(ba->organ_perimeter[i]*ba->organ_perimeter[i])<<" "<<ba->averaged_regularity_in[i]<<" "<<ba->averaged_regularity_epi[i]<<" "<<ba->overlap_index[i]<<endl;
    }

    /*
    fout<<"simulation_index N_in N_epi N_all organ_area organ_perimeter averaged_in_area averaged_epi_area averaged_perimeter";
    fout<<"regularity_in_av regularity_epi_av real_area overlap_area overlap_index"<<std::endl;
    for(int i=0;i<(int)ba->N_in.size();i++){
        fout<<i<<" "<<ba->N_in[i]<<" "<<ba->N_epi[i]<<" "<<ba->N_in[i]+ba->N_epi[i]<<" "<<ba->organ_area[i]<<" "<<ba->organ_perimeter[i]<<" "<<ba->averaged_inner_area[i]<<" "<<ba->averaged_epi_area[i]<<" "<<ba->averaged_perimeter[i];
        fout<<" "<<ba->averaged_regularity_in[i]<<" "<<ba->averaged_regularity_epi[i]<<" "<<ba->real_area[i]<<" "<<ba->overlap_area[i]<<" "<<ba->overlap_index[i] <<std::endl;
    }
    */
    //output the averaged value and variance for each parameter
    fout<<" "<<endl;
    fout<<"averaged value"<<endl;
    fout<<"N_in N_epi N_all organ_area organ_perimeter averaged_in_area averaged_epi_area averaged_perimeter ";
    fout<<"regularity_in_av regularity_epi_av real_area overlap_area overlap_index";
    fout<<"sigma_O sigma_L sigma_O/sigma_L"<<std::endl;
    fout<<averaged_vec_int(ba->N_in)<<" "<<averaged_vec_int(ba->N_epi)<<" "<<averaged_vec_int(ba->N_in)+averaged_vec_int(ba->N_epi)<<" "<<averaged_vec_double(ba->organ_area)<<" "<<averaged_vec_double(ba->organ_perimeter)<<" "<<averaged_vec_double(ba->averaged_inner_area)<<" "<<averaged_vec_double(ba->averaged_epi_area)<<" "<<averaged_vec_double(ba->averaged_perimeter);
    fout<<" "<<averaged_vec_double(ba->averaged_regularity_in)<<" "<<averaged_vec_double(ba->averaged_regularity_epi)<<" "<<averaged_vec_double(ba->real_area)<<" "<<averaged_vec_double(ba->overlap_area)<<" "<<averaged_vec_double(ba->overlap_index);
    fout<<" "<<sigma_O<<" "<<sigma_L<<" "<<sigma_O/sigma_L<<std::endl;
    fout.close();

}

void batch_final_analysis(vector<double> parameter_1, vector<double> parameter_2,vector<double> vec_result,string output_filename){
    std::ofstream fout(output_filename);
    for(int i=0;i<parameter_1.size();i++){
        for(int j=0;j<parameter_2.size();j++){
            fout<<i*parameter_2.size()+j<<" "<<parameter_1[i]<<" "<<parameter_2[j]<<" "<<vec_result[i*parameter_2.size()+j]<<endl;
        }
        fout<<endl;
    }
    fout.close();
}

void simulation_log(vector<time_t> initial_time, vector<time_t> terminal_time){
    char fname[100];
    sprintf(fname,"simulationLog.txt");
    std::ofstream fout(fname);

    tm * initial_tm;
    char initial_time_string [100];
    initial_tm = localtime(&initial_time[0]);
    strftime(initial_time_string, 100, "The whole simulation started at %T, %Y %B %d", initial_tm);

    tm * terminal_tm;
    char terminal_time_string [100];
    terminal_tm = localtime(&terminal_time[terminal_time.size()-1]);
    strftime(terminal_time_string, 100, "The whole simulation terminated at %T, %Y %B %d", terminal_tm);
    

    vector<time_t> simulation_time;
    for(int i=0; i<initial_time.size();i++){
        if(i==0){
            time_t simulation_time_tmp = difftime(terminal_time[terminal_time.size()-1],initial_time[0]);
            simulation_time.push_back(simulation_time_tmp);
        }
        else{
            time_t simulation_time_tmp = difftime(terminal_time[i-1],initial_time[i]);
            simulation_time.push_back(simulation_time_tmp);
        }
    }

    vector<int> simulation_time_h;
    vector<int> simulation_time_m;
    vector<int> simulation_time_s;
    for(int i=0;i<simulation_time.size();i++){
        simulation_time_h.push_back(simulation_time[i]/3600);
        simulation_time_m.push_back((simulation_time[i]%3600)/60);
        simulation_time_s.push_back(simulation_time[i]%60);
    }

    cout<<"The whole simulation (or analysis) calculation took "<<simulation_time_h[0]<<"h "<<simulation_time_m[0]<<"m "<<simulation_time_s[0]<<"s."<<endl;
    fout<<initial_time_string<<endl;
    fout<<terminal_time_string<<endl;
    fout<<"The whole simulation (or analysis) calculation took "<<simulation_time_h[0]<<"h "<<simulation_time_m[0]<<"m "<<simulation_time_s[0]<<"s."<<endl;

    if(simulation_time.size()>1){
        for(int i=1;i<simulation_time.size();i++){
            fout<<"For the "<<i<<"th batch, it took "<<simulation_time_h[i]<<"h "<<simulation_time_m[i]<<"m "<<simulation_time_s[i]<<"s."<<endl;
        }
    }

    fout.close();
}

}

namespace file_process{

int getFileNum(const string &path) {   //需要用到<dirent.h>头文件
    int fileNum=0;
    DIR *pDir;
    struct dirent* ptr;
    if(!(pDir = opendir(path.c_str())))
        return fileNum;
    while((ptr=readdir(pDir))!=0){
        if(strcmp(ptr->d_name,".")!=0&&strcmp(ptr->d_name,"..")!=0 )
           fileNum++;
    }
    closedir(pDir);
    return fileNum;
}


}

namespace cout_fout_debug{
    void cout_vector_vertex(vector<Vertex*> vv){
        cout<<"index x y"<<endl;
        for(int vi=0; vi<vv.size();vi++){
            cout<<vi<<" "<<vv[vi]->loc.x<<" "<<vv[vi]->loc.y<<endl;
        }
    }

    void fout_vector_vertex(vector<Vertex*> vv, string FileName){
        std::ofstream fout(FileName);
        fout<<"index x y z"<<endl;
        for(int vi=0;vi<vv.size();vi++){
            fout<<vi<<" "<<vv[vi]->loc.x<<" "<<vv[vi]->loc.y<<" "<<vv[vi]->loc.z<<endl;
        }
        fout.close();
    }

    void cout_vector_double(vector<double> vd){
        cout<<"index value"<<endl;
        for(int di=0;di<vd.size();di++){
            cout<<di<<" "<<vd[di]<<endl;
        }
    }

    void cout_vector_int(vector<int> vi){
        cout<<"index value"<<endl;
        for(int i=0; i<vi.size();i++){
            cout<<i<<" "<<vi[i]<<endl;
        }
    }

    void cout_vector_string(vector<string> vs){
        cout<<"index string"<<endl;
        for(int i=0; i<vs.size();i++){
            cout<<i<<" "<<vs[i]<<endl;
        }
    }

    void cout_vector_line(vector<Line*> vl){
        cout<<"index slope intercept"<<endl;
        for(int li=0; li<vl.size();li++){
            cout<<li<<" "<<vl[li]->slope<<" "<<vl[li]->intercept<<endl;
        }
    }

    void fout_vector_double(vector<double> vd, string FileName){
        std::ofstream fout(FileName);
        fout<<"index value"<<endl;
        for(int di=0;di<vd.size();di++){
            fout<<di<<" "<<vd[di]<<endl;
        }
        fout.close();
    }

    void fout_vector_int(vector<int> vd, string FileName){
        std::ofstream fout(FileName);
        fout<<"index value"<<endl;
        for(int di=0;di<vd.size();di++){
            fout<<di<<" "<<vd[di]<<endl;
        }
        fout.close();
    }

    void fout_pair_double(vector<pair<double,double>> pdd,string FileName)
    {
        ofstream fout(FileName);
        fout<<"index x y"<<endl;
        for(int pi=0;pi<pdd.size();pi++){
            fout<<pi<<" "<<pdd[pi].first<<" "<<pdd[pi].second<<endl;
        }
        fout.close();
    }

    void fout_tuple_double(vector<tuple<double,double,double>> tddd,string FileName)
    {
        ofstream fout(FileName);
        fout<<"index x y z"<<endl;
        int index=0;
        for(const auto& ti : tddd){
            fout<<index<<" "<<get<0>(ti)<<" "<<get<1>(ti)<<" "<<get<2>(ti)<<endl;
            index++;
        }
        fout.close();
    }

}
