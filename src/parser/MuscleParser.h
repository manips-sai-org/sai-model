#ifndef MUSCLE_PARSER_H_
#define MUSCLE_PARSER_H_

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <functional>
#include <Eigen/Dense>
#include "tinyxml2.h"

namespace {
    inline std::vector<double> parseCSV(const char* text) {
        std::vector<double> vec;
        if (!text) return vec;
        std::stringstream ss(text);
        std::string item;
        while (std::getline(ss, item, ',')) {
            try {
                vec.push_back(std::stod(item));
            } catch (...) {}
        }
        return vec;
    }

    inline Eigen::Vector3d parseVector3d(const char* text) {
        Eigen::Vector3d v = Eigen::Vector3d::Zero();
        if (!text) return v;
        std::stringstream ss(text);
        std::string val;
        try {
            if (std::getline(ss, val, ',')) v.x() = std::stod(val);
            if (std::getline(ss, val, ',')) v.y() = std::stod(val);
            if (std::getline(ss, val, ',')) v.z() = std::stod(val);
        } catch (...) {}
        return v;
    }
}

namespace SaiModel {

struct Waypoint {
    std::string link_name;
    Eigen::Vector3d point;
};

struct ActivatorNode {
    double act_time_constant = 0.0;
    double deact_time_constant = 0.0;
};

struct ContractorNode {
    std::vector<Waypoint> muscle_tendon_path;
    double opt_fiber_length = 0.0;
    double slack_tendon_length = 0.0;
    double max_fiber_velocity = 0.0;
    double peak_iso_force = 0.0;
    double pennation_angle = 0.0;
    double force_vel_curvature = 0.0;
    double passive_damping = 0.0;
};

struct MuscleNode {
    std::string muscle_name;
    std::string comments;
    int id = 0;
    ActivatorNode activator;
    ContractorNode contractor;
};

struct MuscleSystemNode {
    std::string muscle_system_name;
    std::string comments;
    std::string robot_name;
    std::vector<double> fiber_length_norm;
    std::vector<double> iso_force_norm;
    std::vector<MuscleNode> muscles;
};

inline MuscleSystemNode parseMuscleXML(
    const std::string& filename,
    const std::function<bool(const std::string&)>& link_exists = nullptr) {
    MuscleSystemNode system;
    tinyxml2::XMLDocument doc;
    
    if (doc.LoadFile(filename.c_str()) != tinyxml2::XML_SUCCESS) {
        std::cerr << "Error: Could not load XML file " << filename << std::endl;
        return system;
    }

    tinyxml2::XMLElement* root = doc.FirstChildElement("muscleSystemNode");
    if (!root) return system;

    // Root metadata
    if (auto e = root->FirstChildElement("muscleSystemName")) 
        system.muscle_system_name = e->GetText() ? e->GetText() : "";
    if (auto e = root->FirstChildElement("comments")) 
        system.comments = e->GetText() ? e->GetText() : "";
    if (auto e = root->FirstChildElement("robotName")) 
        system.robot_name = e->GetText() ? e->GetText() : "";
    if (auto e = root->FirstChildElement("fiberLengthNorm")) 
        system.fiber_length_norm = parseCSV(e->GetText());
    if (auto e = root->FirstChildElement("isoForceNorm")) 
        system.iso_force_norm = parseCSV(e->GetText());

    // Iterate through all MuscleNodes
    for (tinyxml2::XMLElement* mElem = root->FirstChildElement("muscleNode"); mElem; mElem = mElem->NextSiblingElement("muscleNode")) {
        MuscleNode m;
        
        if (auto e = mElem->FirstChildElement("muscleName")) m.muscle_name = e->GetText() ? e->GetText() : "";
        if (auto e = mElem->FirstChildElement("comments")) m.comments = e->GetText() ? e->GetText() : "";
        if (auto e = mElem->FirstChildElement("ID")) m.id = std::stoi(e->GetText());

        // Activator Node
        if (tinyxml2::XMLElement* aElem = mElem->FirstChildElement("activatorNode")) {
            if (auto e = aElem->FirstChildElement("actTimeConstant")) 
                m.activator.act_time_constant = std::stod(e->GetText());
            if (auto e = aElem->FirstChildElement("deactTimeConstant")) 
                m.activator.deact_time_constant = std::stod(e->GetText());
        }

        // Contractor Node
        if (tinyxml2::XMLElement* cElem = mElem->FirstChildElement("contractorNode")) {
            // Path Parsing
            if (tinyxml2::XMLElement* pathElem = cElem->FirstChildElement("muscleTendonPath")) {
                std::vector<Waypoint> parsed_path;
                bool invalid_path = false;
                std::string invalid_link_name = "";

                for (tinyxml2::XMLElement* child = pathElem->FirstChildElement(); child != nullptr; ) {
                    if (std::string(child->Name()) == "linkName") {
                        Waypoint wp;
                        wp.link_name = child->GetText() ? child->GetText() : "";

                        if (link_exists && !link_exists(wp.link_name)) {
                            invalid_path = true;
                            if (invalid_link_name.empty()) {
                                invalid_link_name = wp.link_name;
                            }
                        }
                        
                        tinyxml2::XMLElement* next = child->NextSiblingElement();
                        if (next && std::string(next->Name()) == "point") {
                            wp.point = parseVector3d(next->GetText());
                            parsed_path.push_back(wp);
                            child = next->NextSiblingElement();
                        } else {
                            child = next;
                        }
                    } else {
                        child = child->NextSiblingElement();
                    }
                }

                if (invalid_path) {
                    parsed_path.clear();
                    m.contractor.muscle_tendon_path.clear();
                    std::cerr
                        << "Warning: Ignoring muscle tendon path for muscle '"
                        << m.muscle_name
                        << "' because link '" << invalid_link_name
                        << "' was not found in the robot model." << std::endl;
                } else {
                    m.contractor.muscle_tendon_path = std::move(parsed_path);
                }
            }

            // Contractor constants
            if (auto e = cElem->FirstChildElement("optFiberLength")) m.contractor.opt_fiber_length = std::stod(e->GetText());
            if (auto e = cElem->FirstChildElement("slackTendonLength")) m.contractor.slack_tendon_length = std::stod(e->GetText());
            if (auto e = cElem->FirstChildElement("maxFiberVelocity")) m.contractor.max_fiber_velocity = std::stod(e->GetText());
            if (auto e = cElem->FirstChildElement("peakIsoForce")) m.contractor.peak_iso_force = std::stod(e->GetText());
            if (auto e = cElem->FirstChildElement("pennationAngle")) m.contractor.pennation_angle = std::stod(e->GetText());
            if (auto e = cElem->FirstChildElement("forceVelCurvature")) m.contractor.force_vel_curvature = std::stod(e->GetText());
            if (auto e = cElem->FirstChildElement("passiveDamping")) m.contractor.passive_damping = std::stod(e->GetText());
        }

        system.muscles.push_back(m);
    }

    return system;
}

} // namespace

#endif

// // ==========================================
// // 4. Main Execution
// // ==========================================

// int main() {
//     MuscleSystemNode mySystem = parseMuscleXML("full_body_1_foot_muscle.xml");

//     std::cout << "System: " << mySystem.muscle_system_name << " | Robot: " << mySystem.robot_name << std::endl;
    
//     if (!mySystem.muscles.empty()) {
//         const auto& m = mySystem.muscles[0];
//         std::cout << "Parsed Muscle: " << m.muscle_name << " with " 
//                   << m.contractor.muscle_tendon_path.size() << " waypoints." << std::endl;
//     }

//     return 0;
// }
