#include "external_potential.hpp"


namespace _pyves::interaction
{
    REAL _z_direction_energy(REAL z, const Box<PBC::ON>& box, REAL surface_width, REAL cutoff)
    {
        if(z <= cutoff)
        {
            return 999'999;
        }
        else if(z > box.getLengthZ()-surface_width)
        {
            return -std::abs(z-(box.getLengthZ()-surface_width)) / surface_width;
        }
        else
        {
            return 0;
        }
    }



    REAL _angle_pow2_penalty(const Particle& p)
    {
        const REAL angle = ::_pyves::directed_angle(p.getOrientation(), CARTESIAN::UnitZ());
        return angle*angle;
    }



    REAL surface_potential(const Particle& p, const Box<PBC::ON>& box, REAL surface_width, REAL cutoff)
    {
        auto z = p.getz();
        while(z > box.getLengthZ())
        {
            z -= box.getLengthZ();
        }

        if((z > cutoff) and (z <= (box.getLengthZ()-surface_width))) // z between |---_________v|
        {
            // std::cout << "z " << z 
            //           << "  return 0 from (z > cutoff) " << std::boolalpha << (z > cutoff) 
            //           << "  (z <= (box.getLengthZ()-surface_width)) " << (z <= (box.getLengthZ()-surface_width)) 
            //           << "  ((z > cutoff) and (z <= (box.getLengthZ()-surface_width))) " << ((z > cutoff) and (z <= (box.getLengthZ()-surface_width)))
            //           << "  all " << (((z > cutoff) and (z <= (box.getLengthZ()-surface_width))) or (surface_width < 1e-3)) << "\n";
            return 0;
        }
        else if (surface_width < 1e-3)
        {
            // std::cout << "surface width " << surface_width << "  " << (surface_width < 1e-3) << "\n";
            return 0;
        }
        else if(z <= cutoff)
        {
            return 999'999;
        }

        const auto _z_pot = p.surface_affinity_translation * _z_direction_energy(z, box, surface_width, cutoff);
        const auto _a_pot = p.surface_affinity_rotation * _angle_pow2_penalty(p);
        
        return _z_pot + _a_pot * (std::abs(z-(box.getLengthZ()-surface_width)) / surface_width);
    }



    REAL external_potential(const Particle& p, const Box<PBC::ON>& box, REAL surface_width, REAL cutoff)
    {
        return surface_potential(p, box, surface_width, cutoff);
    }
    


    void bind_external_potential(py::module& m)
    {        
        m.def("surface_potential", &surface_potential);
        m.def("external_potential", &external_potential);
        m.def("internal_angle_pow2_penalty", &_angle_pow2_penalty);
        m.def("internal_z_direction_energy", &_z_direction_energy);
    }
}