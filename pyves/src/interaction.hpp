#pragma once

#include "definitions.hpp"
#include "box.hpp"
#include "particle.hpp"
#include <cmath>
#include <unordered_map>
#include <Eigen/Core>

#include <pybind11/pybind11.h>
namespace py = pybind11;


namespace _pyves
{


    namespace __detail
    {
        template<class A, class B>
        using map_t = std::map<A,B>;
        
    }
    typedef __detail::map_t<std::string, __detail::map_t<std::string, __detail::map_t<std::string, REAL>>> LookupTable_t;



    template<PBC P>
    REAL interaction(const Particle& p1, const Particle& p2, const Box<P>& box, REAL cutoff)
    {
        
        CARTESIAN distance_vec = box.distanceVectorParticle(p1, p2);
        // std::cout << "\n\ndistance_vec " << distance_vec.format(ROWFORMAT) << std::endl;

        const REAL mean_sigma = (p1.sigma+p2.sigma)/2;
        // std::cout << "mean_sigma " << mean_sigma << std::endl;

        REAL r2 = 1.f/(distance_vec.squaredNorm());
        // std::cout << "distance_vec.squaredNorm() " << distance_vec.squaredNorm() << std::endl;
        // std::cout << "1.f/(distance_vec.squaredNorm()) " << 1.f/(distance_vec.squaredNorm()) << std::endl;

        if(r2 < 1.f/(cutoff*cutoff)) return 0.f;
        r2 *= mean_sigma*mean_sigma;
        // std::cout << "r2 " << r2 << std::endl;

        distance_vec.normalize();
        const CARTESIAN p1_orien_kappa = p1.getOrientation().normalized()*p1.kappa/2;
        // const CARTESIAN p1_orien_kappa = p1.orientation.normalized()*p1.kappa/2;
        // std::cout << "p1_orien_kappa " << p1_orien_kappa.format(ROWFORMAT) << std::endl;
        const CARTESIAN p2_orien_kappa = p2.getOrientation().normalized()*p2.kappa/2;
        // const CARTESIAN p2_orien_kappa = p2.orientation.normalized()*p2.kappa/2;
        // std::cout << "p2_orien_kappa " << p2_orien_kappa.format(ROWFORMAT) << std::endl;

        // const REAL mean_kappa = (p1.kappa+p2.kappa)/2;

        const Eigen::AngleAxis<REAL> rotation1 (p1.gamma, CARTESIAN::UnitZ());
        const Eigen::AngleAxis<REAL> rotation2 (p2.gamma, -CARTESIAN::UnitZ());
        const REAL a  = (-(rotation1*(CARTESIAN::UnitY()*p1.kappa/2)) + CARTESIAN::UnitX() + (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
        const REAL b  = (  rotation1*(CARTESIAN::UnitY()*p1.kappa/2)  + CARTESIAN::UnitX() - (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
        const REAL c1 = (-(rotation1*(CARTESIAN::UnitY()*p1.kappa/2)) + CARTESIAN::UnitX() - (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
        const REAL c2 = (  rotation1*(CARTESIAN::UnitY()*p1.kappa/2)  + CARTESIAN::UnitX() + (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
        // std::cout << "a  " << a << std::endl;
        // std::cout << "b  " << b << std::endl;
        // std::cout << "c1 " << c1 << std::endl;
        // std::cout << "c2 " << c2<< std::endl;
        const REAL chi = static_cast<REAL>(
              std::pow(CARTESIAN( -p1_orien_kappa + distance_vec + p2_orien_kappa ).norm() - a,  2)
            + std::pow(CARTESIAN(  p1_orien_kappa + distance_vec - p2_orien_kappa ).norm() - b,  2)
            + std::pow(CARTESIAN( -p1_orien_kappa + distance_vec - p2_orien_kappa ).norm() - c1, 2)
            + std::pow(CARTESIAN(  p1_orien_kappa + distance_vec + p2_orien_kappa ).norm() - c2, 2));
            
        // std::cout << "chi " << chi << std::endl;
        const REAL mean_epsilon = (p1.epsilon+p2.epsilon)/2;
        const REAL r6 = r2*r2*r2;
        // std::cout << "r6 " << r6 << std::endl;
        
        return 4.f*mean_epsilon*(r6*r6 - (1.f-chi)*r6);
    }



    template<PBC P>
    REAL interactionWithLookup(const Particle& p1, const Particle& p2, const Box<P>& box, REAL cutoff, const LookupTable_t& table)
    {
        CARTESIAN distance_vec = box.distanceVectorParticle(p1, p2);

        const REAL mean_sigma = (p1.sigma+p2.sigma)/2;

        REAL r2 = 1.f/(distance_vec.squaredNorm());

        if(r2 < 1.f/(cutoff*cutoff)) return 0.f;
        r2 *= mean_sigma*mean_sigma;

        distance_vec.normalize();
        const CARTESIAN p1_orien_kappa = p1.getOrientation().normalized()*p1.kappa/2;
        const CARTESIAN p2_orien_kappa = p2.getOrientation().normalized()*p2.kappa/2;
        const auto these_values = table.at(p1.name).at(p2.name);
        
        const REAL chi = static_cast<REAL>(
              std::pow(CARTESIAN( -p1_orien_kappa + distance_vec + p2_orien_kappa ).norm() - these_values.at("a"),  2)
            + std::pow(CARTESIAN(  p1_orien_kappa + distance_vec - p2_orien_kappa ).norm() - these_values.at("b"),  2)
            + std::pow(CARTESIAN( -p1_orien_kappa + distance_vec - p2_orien_kappa ).norm() - these_values.at("c1"), 2)
            + std::pow(CARTESIAN(  p1_orien_kappa + distance_vec + p2_orien_kappa ).norm() - these_values.at("c2"), 2));
            
        const REAL mean_epsilon = (p1.epsilon+p2.epsilon)/2;
        const REAL r6 = r2*r2*r2;
        
        return 4.f*mean_epsilon*(r6*r6 - (1.f-chi)*r6);
    }


    
    inline REAL calculateInteractionOptimumA(const Particle& p1, const Particle& p2)
    {
        const Eigen::AngleAxis<REAL> rotation1 (p1.gamma, CARTESIAN::UnitZ());
        const Eigen::AngleAxis<REAL> rotation2 (p2.gamma, -CARTESIAN::UnitZ());
        return (-(rotation1*(CARTESIAN::UnitY()*p1.kappa/2)) + CARTESIAN::UnitX() + (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
    }


    
    inline REAL calculateInteractionOptimumB(const Particle& p1, const Particle& p2)
    {
        const Eigen::AngleAxis<REAL> rotation1 (p1.gamma, CARTESIAN::UnitZ());
        const Eigen::AngleAxis<REAL> rotation2 (p2.gamma, -CARTESIAN::UnitZ());
        return (  rotation1*(CARTESIAN::UnitY()*p1.kappa/2)  + CARTESIAN::UnitX() - (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
    }


    
    inline REAL calculateInteractionOptimumC1(const Particle& p1, const Particle& p2)
    {
        const Eigen::AngleAxis<REAL> rotation1 (p1.gamma, CARTESIAN::UnitZ());
        const Eigen::AngleAxis<REAL> rotation2 (p2.gamma, -CARTESIAN::UnitZ());
        return (-(rotation1*(CARTESIAN::UnitY()*p1.kappa/2)) + CARTESIAN::UnitX() - (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
    }


    
    inline REAL calculateInteractionOptimumC2(const Particle& p1, const Particle& p2)
    {
        const Eigen::AngleAxis<REAL> rotation1 (p1.gamma, CARTESIAN::UnitZ());
        const Eigen::AngleAxis<REAL> rotation2 (p2.gamma, -CARTESIAN::UnitZ());
        return (  rotation1*(CARTESIAN::UnitY()*p1.kappa/2)  + CARTESIAN::UnitX() + (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
    }



    inline auto calculateInteractionLookupTable(std::vector<Particle> particles)
    {
        LookupTable_t table;

        for(const Particle& p1 : particles)
        {
            for(const Particle& p2 : particles)
            {
                REAL a = calculateInteractionOptimumA(p1, p2);
                REAL b = calculateInteractionOptimumB(p1, p2);
                REAL c1 = calculateInteractionOptimumC1(p1, p2);
                REAL c2 = calculateInteractionOptimumC2(p1, p2);
                table[p1.name][p2.name]["a"] = a;
                table[p1.name][p2.name]["b"] = b;
                table[p1.name][p2.name]["c1"] = c1;
                table[p1.name][p2.name]["c2"] = c2;
            }
        }
        return table;
    }



    inline void bind_interaction(py::module& m) 
    {
        m.def("interaction", &interaction<PBC::ON>);
        m.def("interaction", &interaction<PBC::OFF>);
    }
}