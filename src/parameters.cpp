#include "parameters.h"
#include <fstream>
#include <iostream>

using namespace std;


int Parameters::str_dist(string s, string t)
{

    //The function then fills in the first row and the first column of the matrix. 
    //The value at the i-th row and j-th column of the matrix represents the Levenshtein distance between the first i characters of s and the first j characters of t.
    ulong len_s = s.length();
    ulong len_t = t.length();

    /* base case: empty strings */
    if (len_s == 0) return int(len_t);
    if (len_t == 0) return int(len_s);

    if(len_s == 1 && len_t ==1)
        return s[0] != t[0];

    Eigen::MatrixXd costos(len_s,len_t);

    for(unsigned i = 0 ; i < s.size(); i++){
        for (unsigned j = 0 ; j < t.size(); j++){
            costos(i,j) = 0;
            costos(0,j) = j;
        }
        costos(i,0) = i;
    }

    int cost;

    for(unsigned i = 1 ; i < s.size(); i++){
        for (unsigned j = 1 ; j < t.size(); j++){
            /* test if last characters of the strings match */
            if (s[i] == t[j])
                cost = 0;
            else
                cost = 1;

            /* return minimum of delete char from s, delete char from t, and delete char from both */
            costos(i,j) =  min(min( costos(i-1,j) + 1,
                                    costos(i,j-1) + 1),
                               costos(i-1,j-1) + cost);
        }
    }

    return costos(s.length()-1,t.length()-1);
}

void Parameters::readConfFile(std::string conf_file_path)
{

    ifstream in(conf_file_path);

    if(!in){
        cout << "[ERROR] Configuration file" << endl;
        return;
    }

    string tmp="";
    while((in >> tmp) && (str_dist(tmp,"<END>") > 2) ){

        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"data_directory") <= 2){
            in >> data_directory;
        }
        if (str_dist(tmp,"filename") <= 2){
            in >> filename;
        }
        else if (str_dist(tmp,"repetitions") <= 2){
            in >> repetitions;
        }
        else if (str_dist(tmp,"beading_variation") <= 2){
            in >> beading_variation;
        }
        else if (str_dist(tmp,"nbr_threads") <= 2){
            in >> nbr_threads;
        }
        else if (str_dist(tmp,"nbr_axons_populations") <= 2){
            in >> nbr_axons_populations;
        }
        else if (str_dist(tmp,"crossing_fibers_type") <= 2){
            in >> crossing_fibers_type;
        }
        else if (str_dist(tmp,"tortuous") <= 2){
            in >> tmp;
            if (std::stod(tmp) > 0){
                tortuous = true;
            }
            else{
                tortuous = false;
            }
        }
        else if (str_dist(tmp,"can_shrink") <= 2){
            in >> tmp;
            if (std::stod(tmp) > 0){
                can_shrink = true;
            }
            else{
                can_shrink = false;
            }
        }

        else if (str_dist(tmp,"alpha") <= 2){
            in >> alpha;
        }
        else if (str_dist(tmp,"beta") <= 2){
            in >> beta;
        }
        else if (str_dist(tmp,"std_dev") <= 2){
            in >> std_dev;
        }
        else if (str_dist(tmp,"ondulation_factor") <= 2){
            in >> ondulation_factor;
        }
        else if (str_dist(tmp,"beading_period") <= 2){
            in >> beading_period;
        }
        else if (str_dist(tmp,"regrow_thr") <= 2){
            in >> regrow_thr;
        }
        else if (str_dist(tmp,"min_rad") <= 2){
            in >> min_rad;
        }
        else if (str_dist(tmp,"mean_glial_process_length") <= 2){
            in >> mean_glial_process_length;
        }
        else if (str_dist(tmp,"std_glial_process_length") <= 2){
            in >> std_glial_process_length;
        }
        else if ((str_dist(tmp,"<vox_sizes>") == 0)){
            vox_sizes.clear();

            in >> tmp;
            while((str_dist(tmp, "</vox_sizes>")) >= 2){

                vox_sizes.push_back(std::stod(tmp));
                in >> tmp;
            }
        }
        
        else if ((str_dist(tmp,"glial_pop1_icvf_soma") <= 2)){
            in >> glial_pop1_icvf_soma;
        }
        else if ((str_dist(tmp,"glial_pop1_icvf_branches") <= 2)){
            in >> glial_pop1_icvf_processes;
        }
        else if ((str_dist(tmp,"glial_pop2_icvf_soma") <= 2)){
            in >> glial_pop2_icvf_soma;
        }
        else if ((str_dist(tmp,"glial_pop2_icvf_branches") <= 2)){
            in >> glial_pop2_icvf_processes;
        }
        else if ((str_dist(tmp,"axons_without_myelin_icvf") <= 2)){
            in >> axons_wo_myelin_icvf;
        }
        else if ((str_dist(tmp,"axons_with_myelin_icvf") <= 2)){
            in >> axons_w_myelin_icvf;
        }
        else if ((str_dist(tmp,"spheres_overlap_factor") <= 2)){

            in >> spheres_overlap_factor;
        }
        else if ((str_dist(tmp,"c2") <= 2)){

            in >> cosPhiSquared;

        }

    }
    in.close();

    return;
}

