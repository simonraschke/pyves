#include "stepwidth_alignment_unit.hpp"


namespace _pyves
{
    StepwidthAlignmentUnit& StepwidthAlignmentUnit::operator=(const StepwidthAlignmentUnit& other)
    {
        alignment_every = other.alignment_every;   
        stepwidth = other.stepwidth;
        stepwidth_min = other.stepwidth_min;   
        stepwidth_max = other.stepwidth_max;   
        ratio_target = other.ratio_target;   
        setup_flag = other.setup_flag;   

        accepted_count = other.accepted_count.load();
        rejected_count = other.rejected_count.load();
        return *this; 
    }



    void StepwidthAlignmentUnit::setup(std::size_t interval, REAL min, REAL max, REAL target)
    {
        // tbb::spin_mutex::scoped_lock lock(mutex);
        alignment_every = interval;
        stepwidth = 0.5f*(min+max);
        stepwidth_min = min;
        stepwidth_max = max;
        ratio_target = target;
        setup_flag = true;
    }



    void StepwidthAlignmentUnit::setup(std::size_t interval, REAL min, REAL max, REAL actual_stepwidth, REAL target)
    {
        // tbb::spin_mutex::scoped_lock lock(mutex);
        alignment_every = interval;
        stepwidth = actual_stepwidth;
        stepwidth_min = min;
        stepwidth_max = max;
        ratio_target = target;
        setup_flag = true;
    }

    
    
    REAL StepwidthAlignmentUnit::getMin() const
    {
        return stepwidth_min;
    }



    REAL StepwidthAlignmentUnit::getMax() const
    {
        return stepwidth_max;
    }



    REAL StepwidthAlignmentUnit::operator()() const
    {
        return stepwidth;
    }



    void StepwidthAlignmentUnit::accepted()
    {
        ++accepted_count;
        alignment_check();
    }



    void StepwidthAlignmentUnit::rejected()
    {
        ++rejected_count;
        alignment_check();
    }



    REAL StepwidthAlignmentUnit::getTarget() const
    {
        if( !setup_flag )
            throw std::logic_error("StepwidthAlignmentUnit::setup was not executed");
        else
            return ratio_target;
    }



    REAL StepwidthAlignmentUnit::getRatio() const
    {   
        if( accepted_count == 0 && rejected_count == 0)
            return ratio_old;
        else
            return static_cast<REAL>(accepted_count) / (accepted_count+rejected_count);
    }



    std::size_t StepwidthAlignmentUnit::getAccepted() const
    {
        return accepted_count;
    }



    std::size_t StepwidthAlignmentUnit::getRejected() const
    {
        return rejected_count;
    }



    void StepwidthAlignmentUnit::setAlignmentEvery(std::size_t interval)
    {
        alignment_every = interval;
    }



    void StepwidthAlignmentUnit::alignment_check()
    {
        if( !setup_flag )
            throw std::logic_error("StepwidthAlignmentUnit::setup was not executed");
        else if( accepted_count + rejected_count >= alignment_every ) 
        {
            // this will only be executed once per circle
            // if(!mutex.try_lock())
            //     return;
            ratio_old = getRatio();
            do_aligment();
            accepted_count = rejected_count = 0;
            // mutex.unlock();
        }
    }


    void StepwidthAlignmentUnit::do_aligment()
    {
        REAL deviation = getRatio() - ratio_target;
        if( deviation < 0.0 )
        {
            if( stepwidth > stepwidth_min ) 
            {
                stepwidth *= (1.f+deviation);
            }
        }
        else
        {
            if( stepwidth >= stepwidth_max ) 
            {
                stepwidth = stepwidth_max;
            }
            else if( stepwidth < stepwidth_max )
            {
                stepwidth *= (1.f+deviation);
                if( stepwidth >= stepwidth_max ) 
                {
                    stepwidth = stepwidth_max;
                }
            }
        }
    }


    void bind_sw_alignment_unit(py::module& m)
    {
        py::class_<StepwidthAlignmentUnit>(m, "StepwidthAlignmentUnit", py::dynamic_attr())
            .def(py::init<>()) 
            .def("setup", static_cast<void (StepwidthAlignmentUnit::*)(std::size_t, REAL, REAL, REAL)>(&StepwidthAlignmentUnit::setup))
            .def("setup", static_cast<void (StepwidthAlignmentUnit::*)(std::size_t, REAL, REAL, REAL, REAL)>(&StepwidthAlignmentUnit::setup))
            // .def("__call__", &StepwidthAlignmentUnit::operator())
            .def("ratio", &StepwidthAlignmentUnit::getRatio)
            ;
    }
        
} // namespace _pyves
