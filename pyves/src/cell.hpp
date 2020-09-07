#pragma once

#include "definitions.hpp"
#include "box.hpp"
#include "particle.hpp"
#include "interaction.hpp"
#include <atomic>
#include <memory>
#include <Eigen/Geometry>

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>



namespace _pyves 
{ 
    enum class CellState{ UNDEFINED, IDLE, BLOCKED, FINISHED };

    static const CARTESIAN CellBoundOffset
    {
        static_cast<REAL>(1e-3),
        static_cast<REAL>(1e-3),
        static_cast<REAL>(1e-3)
    };

    struct Cell; 
    typedef std::vector<std::reference_wrapper<Cell>> CellRefContainer;
}



PYBIND11_MAKE_OPAQUE(_pyves::CellRefContainer)



struct _pyves::Cell
{
    std::atomic<CellState> state {CellState::UNDEFINED};
    Eigen::AlignedBox<REAL,3> bounding_box;

    CellRefContainer proximity;
    CellRefContainer region;
    ParticleRefContainer particles;


    Cell(CARTESIAN_CREF min, CARTESIAN_CREF max);
    Cell(const Cell &other);

    Cell& operator=(const Cell &other);
    bool operator==(const Cell& other) const;

    bool isNeighbourOf(const Cell& other, const Box<PBC::ON>& b) const;
    bool contains(const Particle& p) const;
    bool assertIntegrity() const;
    void shuffle();
    REAL potentialEnergy(const Particle& p, const Box<PBC::ON>& b, REAL cutoff) const;
    
    template<CellState S> bool proximityAllInState() const;
    template<CellState S> bool proximityNoneInState() const;

    template<CellState S> bool regionAllInState() const;
    template<CellState S> bool regionNoneInState() const;

    std::string repr() const;
};



namespace _pyves
{
    template<CellState S>
    bool Cell::proximityAllInState() const
    {
        return std::all_of(std::begin(proximity),std::end(proximity), [](const Cell& cell){ return cell.state == S; } );
    }



    template<CellState S>
    bool Cell::proximityNoneInState() const
    {
        return std::none_of(std::begin(proximity),std::end(proximity), [](const Cell& cell){ return cell.state == S; } );
    }



    template<CellState S>
    bool Cell::regionAllInState() const
    {
        return std::all_of(std::begin(region),std::end(region), [](const Cell& cell){ return cell.state == S; } );
    }



    template<CellState S>
    bool Cell::regionNoneInState() const
    {
        return std::none_of(std::begin(region),std::end(region), [](const Cell& cell){ return cell.state == S; } );
    }



    void bind_cell(py::module& m);
} // namespace _pyves