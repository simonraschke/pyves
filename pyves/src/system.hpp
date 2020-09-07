#pragma once
#include "box.hpp"
#include "stepwidth_alignment_unit.hpp"
#include "particle.hpp"
#include "cell.hpp"
#include "interaction.hpp"
#include "metropolis.hpp"
#include "utility.hpp"

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>
#include <random>
#include <numeric>
#include <tbb/task_arena.h>
#include <tbb/parallel_for_each.h>
#include <tbb/task_group.h>



PYBIND11_MAKE_OPAQUE(std::vector<_pyves::Particle>)
PYBIND11_MAKE_OPAQUE(std::vector<_pyves::Cell>)



namespace _pyves
{
    typedef std::vector<_pyves::Particle> ParticleContainer;
    typedef std::vector<_pyves::Cell> CellContainer;
    



    struct System
    {
        Box<PBC::ON> box;
        ParticleContainer particles;
        CellContainer cells;
        StepwidthAlignmentUnit translation_alignment;
        StepwidthAlignmentUnit rotation_alignment;
    
        std::size_t cores = make_nan<std::size_t>();
        std::size_t threads = make_nan<std::size_t>();
        
        std::mt19937_64 pseudo_engine{std::random_device{}()};

        REAL temperature = make_nan<REAL>();
        std::size_t time_max = make_nan<std::size_t>();

        void setThreads(std::size_t);
        bool particleIsFree(const Particle&) const;
        bool particleIsFree(const Particle&, REAL cutoff) const;
        // void prepare_simulation();
        bool assertIntegrity() const;
        void cellStep(const Cell&);
        void shuffle();
        void singleSimulationStep();

        template<typename FUNCTOR> void cellBasedApplyFunctor(FUNCTOR&& func);

        template<CellState S> bool allCellsInState() const;
        template<CellState S> bool noCellsInState() const;

    private:
        tbb::task_arena task_arena;
    };



    template<CellState S>
    bool System::allCellsInState() const
    {
        return std::all_of(std::begin(cells), std::end(cells), [](const Cell& cell){ return cell.state == S; } );
    }



    template<CellState S>
    bool System::noCellsInState() const
    {
        return std::none_of(std::begin(cells), std::end(cells), [](const Cell& cell){ return cell.state == S; } );
    }



    template<typename FUNCTOR>
    void System::cellBasedApplyFunctor(FUNCTOR&& func)
    {
        std::shuffle(std::begin(cells), std::end(cells), pseudo_engine);
        
        std::cout << __PRETTY_FUNCTION__ << "\n";
        tbb::task_group tg;
        task_arena.execute([&]()
        {
            std::cout << "task_arena execute" << "\n";
            while(! (allCellsInState<CellState::FINISHED>()) )
            {
                std::cout << "while" << "\n";
                for(Cell& cell: cells)
                {
                    if( cell.regionNoneInState<CellState::BLOCKED>() && 
                        cell.state == CellState::IDLE
                    )
                    {
                        cell.state = CellState::BLOCKED;
                        std::cout << cell.repr() << "\n";
                        
                        assert( cell.state == CellState::BLOCKED );
                        assert( cell.proximityNoneInState<CellState::BLOCKED>() );
                        
                        tg.run( [&]
                        {
                            assert( cell.state == CellState::BLOCKED );
                            func( cell ); 
                            cell.state = CellState::FINISHED;
                            assert( cell.state == CellState::FINISHED );
                        } );
                    }
                }
            }
        });
        
        assert( allInState<CellState::FINISHED>() );
    }

    

    void bind_system(py::module& m);
}