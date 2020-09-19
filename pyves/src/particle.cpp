#include "particle.hpp"


namespace _pyves{
CARTESIAN::Scalar Particle::getx() const 
{ 
    return position(0); 
}



CARTESIAN::Scalar Particle::gety() const 
{ 
    return position(1); 
}



CARTESIAN::Scalar Particle::getz() const 
{ 
    return position(2); 
}



CARTESIAN Particle::getPosition() const 
{ 
    return CARTESIAN(position); 
}



CARTESIAN Particle::getInitialPosition() const
{
    return *initial_position;
}



void Particle::setPosition(CARTESIAN_CREF v)
{
    position = v;
}



bool Particle::trySetPosition(CARTESIAN_CREF v)
{
    if(!initial_position)
    {
        setInitialPosition(v);
    }

    if((v - getInitialPosition()).squaredNorm() < position_bound_radius_squared)
    {
        setPosition(v);
        return true;
    }
    else
    {
        return false;
    }
}



void Particle::setInitialPosition(CARTESIAN p)
{
    initial_position.reset(new CARTESIAN(p));
}



CARTESIAN::Scalar Particle::getux() const 
{ 
    return orientation(0); 
}



CARTESIAN::Scalar Particle::getuy() const 
{ 
    return orientation(1); 
}



CARTESIAN::Scalar Particle::getuz() const 
{ 
    return orientation(2); 
}


CARTESIAN Particle::getOrientation() const
{
    return CARTESIAN(orientation); 
}



void Particle::setOrientation(CARTESIAN_CREF v)
{
    orientation = v; 
    orientation.normalize(); 
}



bool Particle::trySetOrientation(CARTESIAN_CREF v)
{
    if(!initial_orientation)
    {
        setInitialOrientation(v);
    }
    if(std::acos(v.normalized().dot(getInitialOrientation())) < orientation_bound_radiant)
    {
        setOrientation(v);
        return true;
    }
    else
    {
        return false;
    }
}



CARTESIAN Particle::getInitialOrientation() const
{
    return *initial_orientation;
}



void Particle::setInitialOrientation(CARTESIAN o)
{
    initial_orientation = std::make_unique<CARTESIAN>(o.normalized());
}



void Particle::updateNeighborList(const ParticleRefContainer& particles, const Box<PBC::ON>& box, REAL cutoff)
{
    neighbors.clear();
    const REAL cutoff_sq = cutoff*cutoff;
    std::copy_if(std::begin(particles), std::end(particles), std::back_inserter(neighbors), [&](const Particle& p){
        return p != *this && box.squaredDistance(position, p.position) < cutoff_sq;
    });
    std::shuffle(std::begin(neighbors), std::end(neighbors), RandomEngine.pseudo_engine);
    // std::cout << __func__ << " got "<< neighbors.size() << " neighbors from " << particles.size() << " particles" <<"\n";
}



REAL Particle::potentialEnergy(const Box<PBC::ON>& box, REAL cutoff) const
{
    return std::accumulate(std::begin(neighbors), std::end(neighbors), REAL(0), [&](REAL val, const Particle& other){
        return val + interaction(*this, other, box, cutoff);
    });
}




bool Particle::operator==(const Particle& other) const 
{ 
    return this == &other; 
}



bool Particle::operator!=(const Particle& other) const 
{
    return !(this == &other); 
}



bool Particle::assertIntegrity() const
{
    return all(
        all(std::isfinite(position[0]), std::isfinite(position[1]), std::isfinite(position[2])),
        all(std::isfinite(orientation[0]), std::isfinite(orientation[1]), std::isfinite(orientation[2])),
        std::isfinite(sigma),
        std::isfinite(epsilon),
        std::isfinite(kappa),
        std::isfinite(gamma),
        name != "UNDEFINED"
    );
}



std::string Particle::repr() const
{
    std::stringstream ss;
    ss  << "<Particle " << name << " (REAL=" << type_name<REAL>() << ") at " << position.format(VECTORFORMAT)
        << " in " << orientation.format(VECTORFORMAT) << " direction>";
    return ss.str();
}  



std::string Particle::detailed_repr() const
{
    return repr() + 
    #if __cplusplus <= 201703L
        string_format(") direction with sigma=%f eps=%f kappa=%f gamma=%f >", sigma, epsilon, kappa, gamma);
    #else
        std::format(") direction with sigma={} eps={} kappa={} gamma={} >", sigma, epsilon, kappa, gamma);
    #endif
}  



Particle::Particle(const Particle& other)
    : position(other.position)
    , orientation(other.orientation)
    , sigma(other.sigma)
    , epsilon(other.epsilon)
    , kappa(other.kappa)
    , gamma(other.gamma)
    , name(other.name)
    , position_bound_radius_squared(other.position_bound_radius_squared)
    , orientation_bound_radiant(other.orientation_bound_radiant)
{
    if (other.initial_position)
    {
        CARTESIAN copy = CARTESIAN(*other.initial_position.get());
        initial_position = std::make_unique<CARTESIAN>(CARTESIAN(copy));
    }
    if (other.initial_orientation)
    {
        CARTESIAN copy = CARTESIAN(*other.initial_orientation.get());
        initial_orientation = std::make_unique<CARTESIAN>(CARTESIAN(copy));
    }
}



Particle::Particle(CARTESIAN_CREF pos, CARTESIAN_CREF orien)
    : position(pos)
    , orientation(orien.normalized())
{
    setInitialPosition(pos);
    setInitialOrientation(orien);
}



Particle::Particle(CARTESIAN_CREF pos, CARTESIAN_CREF orien, REAL sigma, REAL eps, REAL kappa, REAL gamma, std::string name)
    : position(pos)
    , orientation(orien.normalized())
    , sigma(sigma)
    , epsilon(eps)
    , kappa(kappa)
    , gamma(gamma)
    , name(name)
{
    setInitialPosition(pos);
    setInitialOrientation(orien);
    // if(!assertIntegrity())
    // {
    //     throw std::logic_error(std::string("Particle should be completely initialized, but is not: ")+detailed_repr());
    // }
}


Particle& Particle::operator=(const Particle& other)
{
    position = other.position;
    orientation = other.orientation;
    sigma = other.sigma;
    epsilon = other.epsilon;
    kappa = other.kappa;
    gamma = other.gamma;
    name = other.name;
    position_bound_radius_squared = other.position_bound_radius_squared;
    orientation_bound_radiant = other.orientation_bound_radiant;
    
    if(other.initial_position)
    {
        CARTESIAN copy = CARTESIAN(*other.initial_position.get());
        initial_position = std::make_unique<CARTESIAN>(CARTESIAN(copy));
    }
    if(other.initial_orientation)
    {
        CARTESIAN copy = CARTESIAN(*other.initial_orientation.get());
        initial_orientation = std::make_unique<CARTESIAN>(CARTESIAN(copy));
    }
    return *this;
}



void bind_particle(py::module& m) 
{
    py::bind_vector<ParticleRefContainer>( m, "ParticleRefContainer" )
        .def(py::init<>())
        .def("__len__", [](const ParticleRefContainer& v) { return v.size(); })
        .def("__iter__", [](ParticleRefContainer& v) 
        {
            return py::make_iterator(std::begin(v), std::end(v));
        }, py::keep_alive<0, 1>())
        .def("size", &ParticleRefContainer::size)
        .def("__repr__", [](const ParticleRefContainer& v) {
            return "ParticleRefContainer\n[\n"
                + std::accumulate(std::begin(v), std::end(v), std::string(""), [](std::string s, const Particle& p) 
                { 
                    return s+"\t"+p.repr()+",\n";
                })
                + "]";
        })
    ;
        
    py::class_<Particle>(m, "Particle", py::dynamic_attr())
        .def(py::init<>())
        .def(py::init<const Particle&>())
        .def(py::init<CARTESIAN_CREF, CARTESIAN_CREF>())
        .def(py::init<CARTESIAN_CREF, CARTESIAN_CREF, REAL, REAL, REAL, REAL, std::string>(), 
                py::arg("pos"), py::arg("orien"), py::arg("sigma"), py::arg("eps"), py::arg("kappa"), py::arg("gamma"), py::arg("name"))
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("assertIntegrity", &Particle::assertIntegrity)
        .def("forceSetPosition", &Particle::setPosition)
        .def_property("position", &Particle::getPosition, &Particle::trySetPosition)
        .def_property("orientation", &Particle::getOrientation, &Particle::trySetOrientation)
        .def("forceSetOrientation", &Particle::setOrientation)
        .def_property("initial_position", &Particle::getInitialPosition, &Particle::setInitialPosition)
        .def_property("initial_orientation", &Particle::getInitialOrientation, &Particle::setInitialOrientation)
        .def_readwrite("translation_bound_sq", &Particle::position_bound_radius_squared)
        .def_readwrite("rotation_bound", &Particle::orientation_bound_radiant)
        .def_readwrite("sigma", &Particle::sigma)
        .def_readwrite("epsilon", &Particle::epsilon)
        .def_readwrite("kappa", &Particle::kappa)
        .def_readwrite("gamma", &Particle::gamma)
        .def_readwrite("name", &Particle::name)
        .def_readonly("neighbors", &Particle::neighbors)
        .def_property_readonly("x", &Particle::getx)//, &Particle::setx)
        .def_property_readonly("y", &Particle::gety)//, &Particle::sety)
        .def_property_readonly("z", &Particle::getz)//, &Particle::setz)
        .def_property_readonly("ux", &Particle::getux)
        .def_property_readonly("uy", &Particle::getuy)
        .def_property_readonly("uz", &Particle::getuz)
        .def("__repr__", &Particle::repr)
    ;
}

} // namespace _pyves