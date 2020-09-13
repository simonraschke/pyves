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
#include <chrono>
#include <thread>
#include <mutex>
#include <tbb/task_arena.h>
#include <tbb/parallel_for_each.h>
// #include <tbb/task_group.h>
#include <tbb/global_control.h>
#undef __TBB_show_deprecated_header_message
#include "parallel.hpp"
#ifndef __TBB_show_deprecated_header_message
#define __TBB_show_deprecated_header_message
#endif



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
        
        // std::mt19937_64 pseudo_engine{std::random_device{}()};

        REAL temperature = make_nan<REAL>();
        // std::size_t time_max = make_nan<std::size_t>();
        // std::mutex mutex;

        void setThreads(std::size_t);
        bool particleIsFree(const Particle&) const;
        bool particleIsFree(const Particle&, REAL cutoff) const;
        void prepareSimulationStep();
        bool assertIntegrity();
        void cellStep(const Cell&);
        void shuffle();
        void singleSimulationStep();
        void multipleSimulationSteps(const unsigned long);
        std::size_t numParticlesInCells() const;

        template<typename FUNCTOR> void applyToCells(FUNCTOR&& func);
        template<typename FUNCTOR> void applyToCellsSlowAndSafe(FUNCTOR&& func);

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



    template<typename FUNC>
    void System::applyToCells(FUNC&& func)
    {
        // std::shuffle(std::begin(cells), std::end(cells), pseudo_engine);
        
        std::vector<CellContainer::iterator> iterators(cells.size());
        std::iota(std::begin(iterators), std::end(iterators), std::begin(cells));
        std::shuffle(std::begin(iterators), std::end(iterators), RandomEngine.pseudo_engine);
        
        task_arena.execute([&]
        {
            enhance::scoped_root_dummy ROOT; 
                
            // std::atomic<std::size_t> i = 0;
            while(! (allCellsInState<CellState::FINISHED>()) )
            {
                // i++;
                for(auto& cell_it: iterators)
                {
                    if( Cell& cell = *cell_it;
                        cell.state == CellState::IDLE &&
                        cell.regionNoneInState<CellState::BLOCKED>()
                    )
                    {
                        cell.state = CellState::BLOCKED;
                        
                        assert( cell.state == CellState::BLOCKED );
                        assert( cell.proximityNoneInState<CellState::BLOCKED>() );
                        
                        ROOT.enqueue_child( [&]
                        {
                            assert( cell.state == CellState::BLOCKED );
                            func( cell ); 
                            cell.state = CellState::FINISHED;
                            assert( cell.state == CellState::FINISHED );
                        } );
                    }
                }
            }
            // std::cout << "ran " << i << "\n";
        });
        
        if( !allCellsInState<CellState::FINISHED>() )
        {
            throw std::runtime_error("not all cells finished");
        }
    }

    

    void bind_system(py::module& m);
}