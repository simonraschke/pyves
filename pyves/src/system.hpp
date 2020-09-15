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
#include <taskflow/taskflow.hpp>



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
    
        REAL interaction_cutoff = make_nan<REAL>();
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

        template<CellState S> std::size_t numCellsInState() const;
        template<CellState S> bool allCellsInState() const;
        template<CellState S> bool noCellsInState() const;

    private:
        // tbb::task_arena task_arena;
        std::unique_ptr<tf::Executor> executor;
    };



    template<CellState S>
    std::size_t System::numCellsInState() const
    {
        return std::accumulate(std::begin(cells), std::end(cells), static_cast<std::size_t>(0), [](std::size_t i, const Cell& cell){ return (cell.state == S) ? ++i : i; } );
    }



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
        CellRefContainer cellrefs;
                
        std::copy_if(std::begin(cells), std::end(cells), std::back_inserter(cellrefs), [](Cell& cell){
            return cell.particles.size() > 0;
        });

        // for(auto cell_it = std::begin(cells); cell_it != std::end(cells); cell_it++)
        // {
        //     if(cell_it->particles.size())
        //     {
        //         iterators.push_back(cell_it);
        //     }
        // }
        // std::iota(std::begin(iterators), std::end(iterators), std::begin(cells));

        // auto starttask = taskflow.placeholder();

        // auto B = taskflow.emplace([&](tf::Subflow& subflow){ 
        //     std::cout << "TaskB\n";
        //     auto B1 = subflow.emplace([&](){ std::cout << "TaskB1\n"; }).name("B1");
        //     auto B2 = subflow.emplace([&](){ std::cout << "TaskB2\n"; }).name("B2");
        //     auto B3 = subflow.emplace([&](){ std::cout << "TaskB3\n"; }).name("B3");
        //     B1.precede(B3); 
        //     B2.precede(B3);
        // }).name("B");

        // executor.run_until(taskflow, [] { return allCellsInState<CellState::FINISHED>(); });

        tf::Taskflow taskflow;
        tf::Task join_task = taskflow.placeholder().name("join_task");
        
        for(std::size_t i = 0; i < executor->num_workers(); ++i)
        {
            tf::Task worker_task = taskflow.emplace([&, cellrefs] () mutable
            {   
                std::shuffle(std::begin(cellrefs), std::end(cellrefs), RandomEngine.pseudo_engine);
                
                while(! (allCellsInState<CellState::FINISHED>()) )
                // while(! (numCellsInState<CellState::FINISHED>() <= 1) )
                {
                    for(Cell& cell: cellrefs)
                    {
                        if(// Cell& cell = cell_it.get();
                            cell.state == CellState::IDLE &&
                            cell.regionNoneInState<CellState::BLOCKED>() &&
                            cell.algorithm_mutex.try_lock()
                        )
                        {
                            cell.state = CellState::BLOCKED;
                            assert( cell.state == CellState::BLOCKED );
                            assert( cell.proximityNoneInState<CellState::BLOCKED>() );
                            assert( cell.state == CellState::BLOCKED );
                            func( cell ); 
                            cell.state = CellState::FINISHED;
                            assert( cell.state == CellState::FINISHED );
                            cell.algorithm_mutex.unlock();
                        }
                    }
                }
            }).name(std::string("worker")+std::to_string(i));
            worker_task.precede(join_task);
        }
        // taskflow.dump(std::cout);
        executor->run(taskflow).wait();

        // executor->wait_for_all();
        if( !allCellsInState<CellState::FINISHED>() )
        {
            throw std::runtime_error("not all cells finished");
        }
    }



    // template<typename FUNC>
    // void System::applyToCells(FUNC&& func)
    // { 
    //     std::vector<CellContainer::iterator> iterators(cells.size());
    //     std::iota(std::begin(iterators), std::end(iterators), std::begin(cells));
    //     std::shuffle(std::begin(iterators), std::end(iterators), RandomEngine.pseudo_engine);
        
    //     // auto find_cell = [&iterators] -> CellContainer::iterator { 
    //     //     return std::find_if(std::begin(iterators), std::end(iterators), [](const auto cell_it){
    //     //         return cell.state == CellState::IDLE && cell.regionNoneInState<CellState::BLOCKED>();
    //     //     });
    //     // }
        
    //     while(! (allCellsInState<CellState::FINISHED>()) )
    //     {
    //         for(auto& cell_it: iterators)
    //         {
    //             if( Cell& cell = *cell_it;
    //                 cell.state == CellState::IDLE &&
    //                 cell.regionNoneInState<CellState::BLOCKED>()
    //             )
    //             {
    //                 cell.state = CellState::BLOCKED;
                    
    //                 assert( cell.state == CellState::BLOCKED );
    //                 assert( cell.proximityNoneInState<CellState::BLOCKED>() );  
                    
    //                 executor->async( [&]
    //                 {
    //                     assert( cell.state == CellState::BLOCKED );
    //                     func( cell ); 
    //                     cell.state = CellState::FINISHED;
    //                     assert( cell.state == CellState::FINISHED );
    //                 } );
    //             }
    //         }
    //     }
    //     executor->wait_for_all();
        
    //     if( !allCellsInState<CellState::FINISHED>() )
    //     {
    //         throw std::runtime_error("not all cells finished");
    //     }
    // }

    

    void bind_system(py::module& m);
}