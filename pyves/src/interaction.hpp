#pragma once

#include "definitions.hpp"
#include "box.hpp"
#include "particle.hpp"
#include <cmath>
#include <Eigen/Core>

#include <pybind11/pybind11.h>
namespace py = pybind11;


namespace _pyves
{
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



    inline void bind_interaction(py::module& m) 
    {
        m.def("interaction", &interaction<PBC::ON>);
        m.def("interaction", &interaction<PBC::OFF>);
        // m.def("interaction",  [](const Particle& p1, const Particle& p2, const Box<PBC::ON>& box) { return interaction(p1, p2, box);});
    }
}