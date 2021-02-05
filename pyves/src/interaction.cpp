#include "interaction.hpp"


namespace _pyves::interaction
{
    REAL chi(const Particle& p1, const Particle& p2, const Box<PBC::ON>& box, REAL cutoff)
    {
        CARTESIAN distance_vec = box.distanceVector(p1.position, p2.position);
        REAL r2 = 1.f/(distance_vec.squaredNorm());

        if(r2 < 1.f/(cutoff*cutoff)) return 0.f;

        distance_vec.normalize();
        const CARTESIAN p1_orien_kappa = p1.getOrientation().normalized()*p1.kappa/2;
        const CARTESIAN p2_orien_kappa = p2.getOrientation().normalized()*p2.kappa/2;

        const Eigen::AngleAxis<REAL> rotation1 (p1.gamma, CARTESIAN::UnitZ());
        const Eigen::AngleAxis<REAL> rotation2 (p2.gamma, -CARTESIAN::UnitZ());
        const REAL a  = (-(rotation1*(CARTESIAN::UnitY()*p1.kappa/2)) + CARTESIAN::UnitX() + (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
        const REAL b  = (  rotation1*(CARTESIAN::UnitY()*p1.kappa/2)  + CARTESIAN::UnitX() - (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
        const REAL c1 = (-(rotation1*(CARTESIAN::UnitY()*p1.kappa/2)) + CARTESIAN::UnitX() - (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
        const REAL c2 = (  rotation1*(CARTESIAN::UnitY()*p1.kappa/2)  + CARTESIAN::UnitX() + (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
        
        return static_cast<REAL>(
              std::pow(CARTESIAN( -p1_orien_kappa + distance_vec + p2_orien_kappa ).norm() - a,  2)
            + std::pow(CARTESIAN(  p1_orien_kappa + distance_vec - p2_orien_kappa ).norm() - b,  2)
            + std::pow(CARTESIAN( -p1_orien_kappa + distance_vec - p2_orien_kappa ).norm() - c1, 2)
            + std::pow(CARTESIAN(  p1_orien_kappa + distance_vec + p2_orien_kappa ).norm() - c2, 2));
    }



    REAL potentialEnergy(const Particle& p1, const Particle& p2, const Box<PBC::ON>& box, const REAL cutoff)
    {
        CARTESIAN distance_vec = box.distanceVector(p1.position, p2.position);
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
        const REAL mean_epsilon = std::sqrt(p1.epsilon*p2.epsilon);
        const REAL r6 = r2*r2*r2;
        // std::cout << "r6 " << r6 << std::endl;

        const REAL scaling = (p1.name == p2.name) ? (p1.self_affinity) : ((p1.other_affinity + p2.other_affinity)/2);
        
        return scaling*4*mean_epsilon*(r6*r6 - (1.f-chi)*r6);
    }


    
    REAL calculateInteractionOptimumA(const Particle& p1, const Particle& p2)
    {
        const Eigen::AngleAxis<REAL> rotation1 (p1.gamma, CARTESIAN::UnitZ());
        const Eigen::AngleAxis<REAL> rotation2 (p2.gamma, -CARTESIAN::UnitZ());
        return (-(rotation1*(CARTESIAN::UnitY()*p1.kappa/2)) + CARTESIAN::UnitX() + (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
    }


    
    REAL calculateInteractionOptimumB(const Particle& p1, const Particle& p2)
    {
        const Eigen::AngleAxis<REAL> rotation1 (p1.gamma, CARTESIAN::UnitZ());
        const Eigen::AngleAxis<REAL> rotation2 (p2.gamma, -CARTESIAN::UnitZ());
        return (  rotation1*(CARTESIAN::UnitY()*p1.kappa/2)  + CARTESIAN::UnitX() - (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
    }


    
    REAL calculateInteractionOptimumC1(const Particle& p1, const Particle& p2)
    {
        const Eigen::AngleAxis<REAL> rotation1 (p1.gamma, CARTESIAN::UnitZ());
        const Eigen::AngleAxis<REAL> rotation2 (p2.gamma, -CARTESIAN::UnitZ());
        return (-(rotation1*(CARTESIAN::UnitY()*p1.kappa/2)) + CARTESIAN::UnitX() - (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
    }


    
    REAL calculateInteractionOptimumC2(const Particle& p1, const Particle& p2)
    {
        const Eigen::AngleAxis<REAL> rotation1 (p1.gamma, CARTESIAN::UnitZ());
        const Eigen::AngleAxis<REAL> rotation2 (p2.gamma, -CARTESIAN::UnitZ());
        return (  rotation1*(CARTESIAN::UnitY()*p1.kappa/2)  + CARTESIAN::UnitX() + (rotation2*(CARTESIAN::UnitY()*p2.kappa/2))).norm();
    }



    LookupTable_t calculateInteractionLookupTable(std::vector<Particle> particles)
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



    void bind_interaction(py::module& m) 
    {
        m.def("interaction", &potentialEnergy);
        m.def("chi", &chi);
        // m.def("interactionWithLookup", &interactionWithLookup);
        // m.def("interaction", &interaction<PBC::OFF>);
    }
}