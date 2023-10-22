#include "../include/rendezvousUT_options.h"


Options_rendezvousUT_t::Options_rendezvousUT_t(const std::string &filename)
{
    // Read YAML file
    options_yaml = YAML::LoadFile(filename);

    //--------------------------------------------------//
    // Parsing
    //

    // Solver                     // Population size
    Nsub =  options_yaml["solver"]["Nsub"].as<int>();
    MaxIter =  options_yaml["solver"]["MaxIter"].as<int>();
    beta = options_yaml["solver"]["beta"].as<double>();
    
    nLeg =  options_yaml["solver"]["nLeg"].as<int>();
    sf =  options_yaml["solver"]["sf"].as<int>();
    planar =  options_yaml["solver"]["planar"].as<bool>();
    E_Pf_constraint =  options_yaml["solver"]["E_Pf_constraint"].as<bool>();
    E_DV_cstr =  options_yaml["solver"]["E_DV_cstr"].as<bool>();
    E_uniform_time =  options_yaml["solver"]["E_uniform_time"].as<bool>();
    E_rendezvous = options_yaml["solver"]["E_rendezvous"].as<bool>();
    DVstd_model = options_yaml["solver"]["DVstd_model"].as<int>();  
    sigma_constraint = options_yaml["solver"]["sigma_constraint"].as<int>();
    obj_func_switch = options_yaml["solver"]["objective"].as<int>();
    
    // E_fakeUT = options_yaml["solver"]["E_fakeUT"].as<bool>();
    E_nav_std = options_yaml["solver"]["E_nav_std"].as<bool>();
    v_RV_free = options_yaml["solver"]["v_RV_free"].as<bool>();
    DV_RV_double = options_yaml["solver"]["DV_RV_double"].as<bool>();
    
    objective_std= options_yaml["solver"]["objective_std"].as<bool>();
    cov_collocation_mode= options_yaml["solver"]["cov_collocation_mode"].as<int>();
    cov_propagation_mode= options_yaml["solver"]["cov_propagation"]["mode"].as<int>();
    cov_propagation_nSub= options_yaml["solver"]["cov_propagation"]["nSub"].as<int>();
    
    /* Mission */
    output_folder = options_yaml["mission"]["output_folder"].as<string>();
    firstguess_folder = options_yaml["mission"]["firstguess_folder"].as<string>();
    mu = options_yaml["mission"]["mu"].as<double>();
    // tfin1 =  options_yaml["mission"]["tfin1"].as<double>();
    // tfin2 =  options_yaml["mission"]["tfin2"].as<double>();

    sigma_r0 =  options_yaml["mission"]["sigma_r0"].as<double>();
    sigma_v0 =  options_yaml["mission"]["sigma_v0"].as<double>();
    sigma_rf_des =  options_yaml["mission"]["sigma_rf_des"].as<double>();
    sigma_vf_des =  options_yaml["mission"]["sigma_vf_des"].as<double>();
    sigma_rRV =  options_yaml["mission"]["sigma_rRV"].as<double>();
    sigma_vRV =  options_yaml["mission"]["sigma_vRV"].as<double>();
    Qd_r =  options_yaml["mission"]["Qd_r"].as<double>();
    Qd_v =  options_yaml["mission"]["Qd_v"].as<double>();
    Qd_level =  options_yaml["mission"]["Qd_level"].as<int>();
    nav_sigma_r =  options_yaml["mission"]["nav_sigma_r"].as<double>();
    nav_sigma_v =  options_yaml["mission"]["nav_sigma_v"].as<double>();
    
    // int n_r0_dim = options_yaml["mission"]["r0_dim"].size();
    string name;
    for (int i=0; i<3; i++){
        r0(i) = options_yaml["mission"]["r0"][i].as<double>(); 
        v0(i) = options_yaml["mission"]["v0"][i].as<double>(); 

        rf(i) = options_yaml["mission"]["rf"][i].as<double>(); 
        vf(i) = options_yaml["mission"]["vf"][i].as<double>(); 

        rRV1(i) = options_yaml["mission"]["rRV1"][i].as<double>(); 
        vRV1(i) = options_yaml["mission"]["vRV1"][i].as<double>(); 
    }

    int nSeg;
    Vector3d rRV, vRV;
    double tfin;
    string namenSeg, nametfin, namer, namev;
    for(int i = 0; i < nLeg; i++){
        namenSeg = "nSeg" + std::to_string(i + 1);
        nametfin = "tfin" + std::to_string(i + 1);

        nSeg =  options_yaml["solver"][namenSeg].as<int>();
        tfin =  options_yaml["mission"][nametfin].as<double>();

        if(i < nLeg - 1){
            namer = "rRV" + std::to_string(i + 1);
            namev = "vRV" + std::to_string(i + 1);

            for (int i=0; i<3; i++){    
                rRV(i) = options_yaml["mission"][namer][i].as<double>(); 
                vRV(i) = options_yaml["mission"][namev][i].as<double>();
            }
        }

        nSeg_vector.push_back(nSeg);

        rRV_vector.push_back(rRV); 
        vRV_vector.push_back(vRV);

        tfin_vector.push_back(tfin);
    }

    // double DVmax;   
    DVtot_max = options_yaml["mission"]["DVtot_max"].as<double>(); 
    DVtot_single_max = options_yaml["mission"]["DVtot_single_max"].as<double>();   // [km/s]
    amrif = 1000; //kg

    // Nondimensionalization
    // rconv = r0.norm();
    // vconv = sqrt(amu_dim/rconv);
    // aconv = vconv*vconv/rconv;
    // tconv = rconv/vconv;

    rconv = 384400;
    vconv = 1.024546847245897;
    tconv = rconv/vconv;
    aconv = vconv*vconv/rconv;

    kstd = sqrt(2*log(1./beta) + sqrt(3.));

    r0 = r0;
    v0 = v0;
    rf = rf;
    vf = vf;
    rRV1 = rRV1;
    vRV1 = vRV1;
    tfin1 = tfin1;
    tfin2 = tfin2;

    DVtot_max /= vconv;
    //DVtot_single_max /= vconv;
}
 
void Options_rendezvousUT_t::emit(const std::string &filename){

    YAML::Emitter out;
    string namenSeg, nametfin, namer, namev;

    out << YAML::BeginMap;
        out << YAML::Key << "solver";
        out << YAML::BeginMap;
            //out << YAML::Key << "nSeg" << YAML::Value << nSeg;
            out << YAML::Key << "nLeg" << YAML::Value << nLeg;

            for(int i = 0; i < nLeg; i++){
                namenSeg = "nSeg" + std::to_string(i + 1);

                out << YAML::Key << namenSeg << YAML::Value << nSeg_vector[i];
            }
            out << YAML::Key << "sf" << YAML::Value << sf;
            out << YAML::Key << "Nsub" << YAML::Value << Nsub;
            out << YAML::Key << "beta" << YAML::Value << beta;
            out << YAML::Key << "kstd" << YAML::Value << kstd;
            out << YAML::Newline << YAML::Newline;

            out << YAML::Key << "planar" << YAML::Value << planar;
            out << YAML::Key << "E_Pf_constraint" << YAML::Value << E_Pf_constraint;
            out << YAML::Key << "E_uniform_time" << YAML::Value << E_uniform_time;
            out << YAML::Key << "E_rendezvous" << YAML::Value << E_rendezvous;
            out << YAML::Key << "DVstd_model" << YAML::Value << DVstd_model;
            out << YAML::Key << "sigma_constraint" << YAML::Value << sigma_constraint;
            out << YAML::Newline << YAML::Newline;

            out << YAML::Key << "objective" << YAML::Value << obj_func_switch;
            out << YAML::Key << "E_DV_cstr" << YAML::Value << E_DV_cstr;
            out << YAML::Key << "E_Pf_constraint" << YAML::Value << E_Pf_constraint;
            out << YAML::Key << "objective_std" << YAML::Value << objective_std;
            out << YAML::Key << "E_navigation_std" << YAML::Value << E_nav_std;
            out << YAML::Key << "v_RV_free" << YAML::Value << v_RV_free;
            out << YAML::Key << "DV_RV_double" << YAML::Value << DV_RV_double;
            out << YAML::Newline << YAML::Newline;

            out << YAML::Key << "cov_collocation_mode" << YAML::Value << cov_collocation_mode;
            out << YAML::Key << "cov_propagation";
            out << YAML::BeginMap;
                out << YAML::Key << "mode" << YAML::Value << cov_propagation_mode;
                out << YAML::Key << "nSub" << YAML::Value << cov_propagation_nSub;
            out << YAML::EndMap;
        out << YAML::EndMap;
        out << YAML::Newline << YAML::Newline;

        out << YAML::Key << "mission";
        out << YAML::BeginMap;
            out << YAML::Key << "output_folder" << YAML::Value << output_folder;
            out << YAML::Key << "firstguess_folder" << YAML::Value << firstguess_folder;
            out << YAML::Key << "mu" << YAML::Value << mu;
            out << YAML::Newline << YAML::Newline;

            vector<double> r0_dim={r0[0]*rconv, r0[1]*rconv, r0[2]*rconv};
            vector<double> v0_dim={v0[0]*vconv, v0[1]*vconv, v0[2]*vconv};
            vector<double> rf_dim={rf[0]*rconv, rf[1]*rconv, rf[2]*rconv};
            vector<double> vf_dim={vf[0]*vconv, vf[1]*vconv, vf[2]*vconv};
            vector<double> rRV_dim={rRV1[0]*rconv, rRV1[1]*rconv, rRV1[2]*rconv};
            vector<double> vRV_dim={vRV1[0]*vconv, vRV1[1]*vconv, vRV1[2]*vconv};
            //vector<double> tf_vector_dim(nLeg);
            // for (int i=0; i<nLeg; i++){
            //     tf_vector_dim(i) = tf_vector[0]*tconv; 
            // }

            out << YAML::Key << "r0_dim" << YAML::Value << r0_dim;
            out << YAML::Key << "v0_dim" << YAML::Value << v0_dim;
            out << YAML::Newline << YAML::Newline;
            for(int i = 0; i < nLeg; i++){
                nametfin = "tfin" + std::to_string(i + 1);
                namer = "rRV" + std::to_string(i + 1);
                namev = "vRV" + std::to_string(i + 1);

                out << YAML::Key << nametfin << YAML::Value << tfin_vector[i]*tconv;
                
                if(i < nLeg - 1){
                    namer = "rRV" + std::to_string(i + 1);
                    namev = "vRV" + std::to_string(i + 1);

                    vector<double> rRV_dim = {rRV_vector[i][0]*rconv, rRV_vector[i][1]*rconv, rRV_vector[i][2]*rconv};
                    vector<double> vRV_dim = {vRV_vector[i][0]*vconv, vRV_vector[i][1]*vconv, vRV_vector[i][2]*vconv};
                    out << YAML::Key << namer << YAML::Value << rRV_dim;
                    out << YAML::Key << namev << YAML::Value << vRV_dim;
                }
                out << YAML::Newline << YAML::Newline;
            }
            out << YAML::Key << "rf_dim" << YAML::Value << rf_dim;
            out << YAML::Key << "vf_dim" << YAML::Value << vf_dim;
            // out << YAML::Key << "DVtot_max" << YAML::Value << DVtot_max;
            out << YAML::Key << "DVtot_single_max" << YAML::Value << DVtot_single_max;
            
            out << YAML::Newline << YAML::Newline;
            out << YAML::Key << "sigma_r0" << YAML::Value << sigma_r0;
            out << YAML::Key << "sigma_v0" << YAML::Value << sigma_v0;
            out << YAML::Key << "sigma_rf_des" << YAML::Value << sigma_rf_des;
            out << YAML::Key << "sigma_vf_des" << YAML::Value << sigma_vf_des;
            out << YAML::Key << "sigma_rRV" << YAML::Value << sigma_rRV;
            out << YAML::Key << "sigma_vRV" << YAML::Value << sigma_vRV;
            out << YAML::Key << "Qd_r" << YAML::Value << Qd_r;
            out << YAML::Key << "Qd_v" << YAML::Value << Qd_v;
            out << YAML::Key << "Qd_level" << YAML::Value << Qd_level;
            out << YAML::Key << "nav_sigma_r" << YAML::Value << nav_sigma_r;
            out << YAML::Key << "nav_sigma_v" << YAML::Value << nav_sigma_v;
            out << YAML::Key << "sf" << YAML::Value << sf;
            out << YAML::Newline << YAML::Newline;

            out << YAML::Key << "rconv" << YAML::Value << rconv;
            out << YAML::Key << "vconv" << YAML::Value << vconv;
            out << YAML::Key << "tconv" << YAML::Value << tconv;
            out << YAML::Key << "aconv" << YAML::Value << aconv;
            out << YAML::Newline;

        out << YAML::EndMap;
    out << YAML::EndMap;

    //-----------------------------------------------//
    ofstream file_out(filename);
    file_out << out.c_str() << endl;
    file_out.close();

}




 