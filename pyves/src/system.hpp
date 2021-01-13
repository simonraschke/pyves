#pragma once
#include "box.hpp"
#include "stepwidth_alignment_unit.hpp"
#include "particle.hpp"
#include "cell.hpp"
#include "interaction.hpp"
#include "external_potential.hpp"
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

#ifdef PYVES_USE_TBB
#include "parallel.hpp"
#else
#include <taskflow/taskflow.hpp>
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
    
        REAL interaction_cutoff = make_nan<REAL>();
        bool interaction_surface = false;
        REAL interaction_surface_width = 0;
        REAL temperature = make_nan<REAL>();
        REAL neighbor_cutoff = make_nan<REAL>();
        REAL global_exchange_ratio = make_nan<REAL>();
        REAL global_exchange_epot_theshold = make_nan<REAL>();
        std::size_t threads = make_nan<std::size_t>();
        std::size_t _cell_update_step_count = 0;
        std::size_t cell_update_interval = make_nan<std::size_t>();
        std::size_t _neighbor_update_step_count = 0;
        std::size_t neighbor_update_interval = make_nan<std::size_t>();

        LookupTable_t lookup_table;

        void setThreads(std::size_t);
        bool particleIsFree(const Particle&) const;
        bool particleIsFree(const Particle&, REAL cutoff) const;
        void prepareSimulationStep();
        bool assertIntegrity();
        void cellStep(const Cell&);
        Cell& cellOfParticle(const Particle&);
        void exchangeParticles(Particle&, Particle&);
        ParticleRefContainer randomParticles(std::size_t num);
        void globalExchange();
        void shuffle();
        void reorderCells();
        void makeNeighborLists();
        void singleSimulationStep();
        void multipleSimulationSteps(const unsigned long);
        void makeInteractionLookupTable(ParticleContainer);
        std::size_t numParticlesInCells() const;
        REAL totalEnergy(const Particle&) const;
        REAL totalEnergy();
        REAL potentialEnergy() const;
        REAL potentialEnergyConcurrent();
        REAL potentialEnergyConcurrentBruteForce();
        Eigen::Matrix<REAL, Eigen::Dynamic, 1> particleEnergies() const;
        Eigen::Matrix<REAL, Eigen::Dynamic, 1> particleChiValues() const;
        Eigen::Matrix<REAL, Eigen::Dynamic, 1> particleSurfacePotentialValues() const;
        Eigen::Matrix<REAL, Eigen::Dynamic, 1> particleExternalPotentialValues() const;

        void benchmark(std::size_t);

        template<typename FUNCTOR> void applyToCells(FUNCTOR&& func);
        template<typename FUNCTOR> void applyToCellsSlowAndSafe(FUNCTOR&& func);

        template<CellState S> std::size_t numCellsInState() const;
        template<CellState S> bool allCellsInState() const;
        template<CellState S> bool noCellsInState() const;
        

    private:

        #ifdef PYVES_USE_TBB
            tbb::task_arena task_arena;
        #else
            std::unique_ptr<tf::Executor> executor;
        #endif
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
            cell.state = (cell.particles.size() > 0) ? CellState::IDLE : CellState::FINISHED;
            return cell.particles.size() > 0;
        });

        // std::cout << "calculating " << cellrefs.size() <<  " from " << cells.size() << "\n";

#ifdef PYVES_USE_TBB
        {   
            pyves::scoped_root_dummy ROOT; 
            for(int i = 0; i < task_arena.max_concurrency(); ++i)
            {
                ROOT.enqueue_child( [&, cellrefs] () mutable
                {   
                    std::shuffle(std::begin(cellrefs), std::end(cellrefs), RandomEngine.pseudo_engine);
                    
                    while(! (allCellsInState<CellState::FINISHED>()) )
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
                });
            }
        }



        // std::shuffle(std::begin(cellrefs), std::end(cellrefs), RandomEngine.pseudo_engine);
        // tbb::parallel_for_each(std::begin(cellrefs), std::end(cellrefs), [&] (Cell& cell)
        // {
        //     while(cell.state != CellState::FINISHED)
        //     {
        //         if(// Cell& cell = cell_it.get();
        //             cell.state == CellState::IDLE &&
        //             cell.regionNoneInState<CellState::BLOCKED>() &&
        //             cell.algorithm_mutex.try_lock()
        //         )
        //         {
        //             cell.state = CellState::BLOCKED;
        //             assert( cell.state == CellState::BLOCKED );
        //             assert( cell.proximityNoneInState<CellState::BLOCKED>() );
        //             assert( cell.state == CellState::BLOCKED );
        //             func( cell ); 
        //             cell.state = CellState::FINISHED;
        //             assert( cell.state == CellState::FINISHED );
        //             cell.algorithm_mutex.unlock();
        //         }
        //     }
        // });

        // // OLD VERSION, slightly slower at 
        // {
        //     pyves::scoped_root_dummy ROOT; 
        //     while(! (allCellsInState<CellState::FINISHED>()) )
        //     {
        //         for(Cell& cell: cellrefs)
        //         {
        //             if(// Cell& cell = cell_it.get();
        //                 cell.state == CellState::IDLE &&
        //                 cell.regionNoneInState<CellState::BLOCKED>()
        //             )
        //             {
        //                 cell.state = CellState::BLOCKED;
        //                 assert( cell.state == CellState::BLOCKED );
        //                 assert( cell.proximityNoneInState<CellState::BLOCKED>() );
        //                 assert( cell.state == CellState::BLOCKED );
        //                 ROOT.enqueue_child( [&]
        //                 {
        //                     assert( cell.state == CellState::BLOCKED );
        //                     func( cell ); 
        //                     cell.state = CellState::FINISHED;
        //                     assert( cell.state == CellState::FINISHED );
        //                 } );
        //             }
        //         }
        //     }
        // }
#else
        tf::Taskflow taskflow;
        tf::Task join_task = taskflow.placeholder().name("join_task");
        
        for(std::size_t i = 0; i < executor->num_workers(); ++i)
        {
            tf::Task worker_task = taskflow.emplace([&, cellrefs] () mutable
            {   
                std::shuffle(std::begin(cellrefs), std::end(cellrefs), RandomEngine.pseudo_engine);
                
                while(! (allCellsInState<CellState::FINISHED>()) )
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
        executor->run(taskflow).wait();
#endif
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