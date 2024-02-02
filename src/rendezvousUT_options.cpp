#include "../include/rendezvousUT_options.h"


Options_rendezvousUT_t::Options_rendezvousUT_t(const std::string &filename){
    // Read YAML file
    options_yaml = YAML::LoadFile(filename);

    nLeg =  options_yaml["solver"]["nLeg"].as<int>();
    sf =  options_yaml["solver"]["sf"].as<int>();
    Max_iter = options_yaml["solver"]["Max_iter"].as<int>();

    E_DV_cstr =  options_yaml["solver"]["E_DV_cstr"].as<bool>();
    E_Pf_constraint =  options_yaml["solver"]["E_Pf_constraint"].as<bool>();
    
    // E_fakeUT = options_yaml["solver"]["E_fakeUT"].as<bool>();
    E_nav_std = options_yaml["solver"]["E_nav_std"].as<bool>();
    v_RV_free = options_yaml["solver"]["v_RV_free"].as<bool>();
    DV_RV_double = options_yaml["solver"]["DV_RV_double"].as<bool>();
    Fixed_ToF_Leg = options_yaml["solver"]["Fixed_ToF_Leg"].as<bool>();
    Limited_ToF = options_yaml["solver"]["Limited_ToF"].as<bool>();
    
    /* Mission */
    output_folder = options_yaml["mission"]["output_folder"].as<string>();
    firstguess_folder = options_yaml["mission"]["firstguess_folder"].as<string>();
    mu_primary= options_yaml["mission"]["muPrimary"].as<double>();
    mu_FB = options_yaml["mission"]["mu_FB"].as<double>();

    sigma2_r0 =  options_yaml["mission"]["sigma2_r0"].as<double>();
    sigma2_v0 =  options_yaml["mission"]["sigma2_v0"].as<double>();
    sigma2_rf =  options_yaml["mission"]["sigma2_rf"].as<double>();
    sigma2_vf =  options_yaml["mission"]["sigma2_vf"].as<double>();
    sigma2_rRV =  options_yaml["mission"]["sigma2_rRV"].as<double>();
    sigma2_vRV =  options_yaml["mission"]["sigma2_vRV"].as<double>();

    Qd_level =  options_yaml["mission"]["Qd_level"].as<int>();
    
    // int n_r0_dim = options_yaml["mission"]["r0_dim"].size();
    for(int i = 0; i < 3; i++){  
        r0(i) = options_yaml["mission"]["r0"][i].as<double>(); 
        v0(i) = options_yaml["mission"]["v0"][i].as<double>();
        rf(i) = options_yaml["mission"]["rf"][i].as<double>(); 
        vf(i) = options_yaml["mission"]["vf"][i].as<double>(); 
    }

    FB_Legs = options_yaml["mission"]["FB_Legs"].as<vector<int>>();
    
    // Nondimensionalization
    rconv = r0.norm();
    vconv = sqrt(mu_primary/rconv);
    aconv = vconv*vconv/rconv;
    tconv = rconv/vconv;

    int nSeg;
    double tfin;
    string namenSeg, nametfin, namer, namev, namerRV, namevRV;
    for(int i = 0; i < nLeg; i++){
        namenSeg = "nSeg" + std::to_string(i + 1);
        nametfin = "tfin" + std::to_string(i + 1);

        nSeg =  options_yaml["solver"][namenSeg].as<int>();
        tfin =  options_yaml["mission"][nametfin].as<double>()/tconv;

        namerRV = "r0_RV" + std::to_string(i + 1);
        namevRV = "v0_RV" + std::to_string(i + 1);

        for(int i = 0; i < 3; i++){    
            r0_RV(i) = options_yaml["mission"][namerRV][i].as<double>()/rconv;
            v0_RV(i) = options_yaml["mission"][namevRV][i].as<double>()/vconv;
        }
        nSeg_vector.push_back(nSeg);

        rRV_vector.push_back(rRV); 
        vRV_vector.push_back(vRV);
        r0_RV_vector.push_back(r0_RV);
        v0_RV_vector.push_back(v0_RV);

        tfin_vector.push_back(tfin);
    }
    Max_ToF = options_yaml["mission"]["Max_ToF"].as<double>()/tconv;
    r_min = options_yaml["mission"]["r_min"].as<double>();

    // double DVmax;   
    DVtot_single_max = options_yaml["mission"]["DVtot_single_max"].as<double>();   // [km/s]

    r0 = r0/rconv;
    v0 = v0/vconv;
    rf = rf/rconv;
    vf = vf/vconv;

    //DVtot_single_max /= vconv;
}

void Options_rendezvousUT_t::emit(const std::string &filename){

    YAML::Emitter out;
    string namenSeg, nametfin;

    out << YAML::BeginMap;
        out << YAML::Key << "solver";
        out << YAML::BeginMap;
        for(int i = 0; i < nLeg; i++){
            namenSeg = "nSeg" + std::to_string(i + 1);
            nametfin = "tfin" + std::to_string(i + 1);

            out << YAML::Key << namenSeg << YAML::Value << nSeg_vector[i];
            out << YAML::Key << nametfin << YAML::Value << tfin_vector[i];
        }

            out << YAML::Key << "nLeg" << YAML::Value << nLeg;
            out << YAML::Key << "E_Pf_constraint" << YAML::Value << E_Pf_constraint;

            out << YAML::Key << "E_DV_cstr" << YAML::Value << E_DV_cstr;
            out << YAML::Key << "E_navigation_std" << YAML::Value << E_nav_std;
            out << YAML::Key << "v_RV_free" << YAML::Value << v_RV_free;
            out << YAML::Key << "DV_RV_double" << YAML::Value << DV_RV_double;
            out << YAML::Key << "Fixed_DT" << YAML::Value << Fixed_ToF_Leg;
        out << YAML::EndMap;
        out << YAML::Newline << YAML::Newline;

        out << YAML::Key << "mission";
        out << YAML::BeginMap;
            out << YAML::Key << "output_folder" << YAML::Value << output_folder;
            out << YAML::Key << "firstguess_folder" << YAML::Value << firstguess_folder;
            out << YAML::Key << "muPrimary" << YAML::Value << mu_primary;
            out << YAML::Newline << YAML::Newline;

            vector<double> r0_dim = {r0[0]*rconv, r0(1)*rconv, r0(2)*rconv};
            vector<double> v0_dim = {v0[0]*vconv, v0(1)*vconv, v0(2)*vconv};
            vector<double> rf_dim = {rf[0]*rconv, rf(1)*rconv, rf(2)*rconv};
            vector<double> vf_dim = {vf[0]*vconv, vf(1)*vconv, vf(2)*vconv};

            out << YAML::Key << "r0_dim" << YAML::Value << r0_dim;
            out << YAML::Key << "v0_dim" << YAML::Value << v0_dim;
            out << YAML::Key << "rf_dim" << YAML::Value << rf_dim;
            out << YAML::Key << "vf_dim" << YAML::Value << vf_dim;
            out << YAML::Key << "Max_ToF" << YAML::Value << Max_ToF;
            out << YAML::Key << "r_min" << YAML::Value << r_min;
            out << YAML::Key << "DVtot_single_max" << YAML::Value << DVtot_single_max;
            
            out << YAML::Newline << YAML::Newline;
            out << YAML::Key << "sigma2_r0" << YAML::Value << sigma2_r0;
            out << YAML::Key << "sigma2_v0" << YAML::Value << sigma2_v0;
            out << YAML::Key << "sigma2_rf" << YAML::Value << sigma2_rf;
            out << YAML::Key << "sigma2_vf" << YAML::Value << sigma2_vf;
            out << YAML::Key << "sigma2_rRV" << YAML::Value << sigma2_rRV;
            out << YAML::Key << "sigma2_vRV" << YAML::Value << sigma2_vRV;

            out << YAML::Key << "Qd_level" << YAML::Value << Qd_level;
            out << YAML::Key << "sf" << YAML::Value << sf;
            
        out << YAML::EndMap;
    out << YAML::EndMap;

    //-----------------------------------------------//
    ofstream file_out(filename);
    file_out << out.c_str() << endl;
    file_out.close();

}

 