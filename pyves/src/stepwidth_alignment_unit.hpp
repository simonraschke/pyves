#pragma once

#include "definitions.hpp"
#include <cstdlib>
#include <exception>
#include <iostream>
#include <atomic>

#include <pybind11/pybind11.h>
namespace py = pybind11;



namespace _pyves { class StepwidthAlignmentUnit; }



class _pyves::StepwidthAlignmentUnit
{
public:
    StepwidthAlignmentUnit& operator=(const StepwidthAlignmentUnit&);

    void setup(std::size_t, REAL, REAL, REAL);
    void setup(std::size_t, REAL, REAL, REAL, REAL);

    REAL operator()() const;
    void accepted();
    void rejected();
    REAL getTarget() const;
    std::size_t getAccepted() const;
    std::size_t getRejected() const;
    REAL getRatio() const;
    void setAlignmentEvery(std::size_t);

protected:
    void alignment_check();
    void do_aligment();

    REAL stepwidth;
    REAL stepwidth_min;
    REAL stepwidth_max;
    REAL ratio_target;
    REAL ratio_old {0};
    std::atomic<std::size_t> accepted_count {0};
    std::atomic<std::size_t> rejected_count {0};
    std::size_t alignment_every {0};
    bool setup_flag {false};
};



namespace _pyves
{
    void bind_sw_alignment_unit(py::module& m);
}