#include "parameters.h"
#include <fstream>
#include <iostream>

using namespace std;


int Parameters::str_dist(string s, string t)
{
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
    while((in >> tmp) && (str_dist(tmp,"<END>") >= 2) ){


        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"data_directory") <= 2){
            in >> data_directory;
        }
        else if (str_dist(tmp,"repetitions") <= 2){
            in >> repetitions;
        }
        else if (str_dist(tmp,"beading_variation") <= 2){
            in >> beading_variation;
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
        else if (str_dist(tmp,"draw") <= 2){
            in >> tmp;
            if (std::stod(tmp) > 0){
                draw = true;
            }
            else{
                draw = false;
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
        else if (str_dist(tmp,"beading_frequency") <= 2){
            in >> beading_frequency;
        }
        else if (str_dist(tmp,"regrow_thr") <= 2){
            in >> regrow_thr;
        }
        else if (str_dist(tmp,"min_rad") <= 2){
            in >> min_rad;
        }
        else if ((str_dist(tmp,"<vox_sizes>") == 0)){

            in >> tmp;
            while((str_dist(tmp, "</vox_sizes>")) >= 2){

                vox_sizes.push_back(std::stod(tmp));
                in >> tmp;
            }
        }
        
        else if ((str_dist(tmp,"<capacities>") == 0)){

            in >> tmp;
            while((str_dist(tmp, "</capacities>")) >= 2){

                capacities.push_back(std::stod(tmp));
                in >> tmp;
            }
        }
        else if ((str_dist(tmp,"<icvf>") == 0)){

            in >> tmp;
            while((str_dist(tmp, "</icvf>")) >= 2){

                icvf.push_back(std::stod(tmp));
                in >> tmp;
            }
    
        }
        else if ((str_dist(tmp,"<spheres_overlap_factors>") == 0)){

            in >> tmp;
            while((str_dist(tmp, "</spheres_overlap_factors>")) >= 2){

                spheres_overlap_factors.push_back(std::stod(tmp));
                in >> tmp;
            }
        }
        
        
    }
    in.close();

    cout << "*********** Finished Reading Configuration File ***********" << endl;
    cout << "icvf :" << icvf[0] << endl;
    cout << "beading_variation :" << beading_variation << endl;
    cout << "beading_frequency :" << beading_frequency << endl;
    return;
}

