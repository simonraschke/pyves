#include "external_potential.hpp"


namespace _pyves::interaction
{
    REAL surface_potential(const Particle& p, const Box<PBC::ON>& box, REAL surface_width, REAL cutoff)
    {
        auto z = p.getz();
        while(z > box.getLengthZ())
        {
            z -= box.getLengthZ();
        }
        
        REAL translation_potential = 0;
        if(z < cutoff)
        {
            translation_potential += 1e7;
        }
        else if(z > surface_width)
        {
            translation_potential += (box.getLengthZ()-z) / surface_width;
        }
        
        const REAL angle = ::_pyves::absolute_angle(p.getOrientation(), CARTESIAN::UnitZ());
        const REAL rotation_potential = angle*angle;

        return p.surface_affinity_translation * translation_potential + p.surface_affinity_rotation * rotation_potential;
    }



    REAL external_potential(const Particle& p, const Box<PBC::ON>& box, REAL surface_width, REAL cutoff)
    {
        return surface_potential(p, box, surface_width, cutoff);
    }
    

    void bind_external_potential(py::module& m)
    {        
        m.def("surface_potential", &surface_potential);
        m.def("external_potential", &external_potential);
    }
}