#pragma once



namespace _pyves
{
    // singleton implementaton
    // to derive class from
    // USE ONLY if DERIVED must not be constructed more than once
    // constructing DERIVED twice is impossible!
    template<typename DERIVED>
    struct Singleton
    {
        // get the singleton instance
        // DERIVED will be constructed on first call
        static DERIVED& getInstance()
        {
            if( INSTANCE == nullptr)
                INSTANCE = new DERIVED();
            return *INSTANCE;
        }
        
        // manually destroy DERIVED
        static void destroyInstance()
        {
            delete INSTANCE;
            INSTANCE = nullptr;
        }
        
    protected:
        // only to derive from
        explicit Singleton() { }
        
        // destroy if derived is destroyed
        virtual ~Singleton()
        {
            INSTANCE = nullptr;
        }
        
    private:
        static DERIVED* INSTANCE;
        

        // do not copy
        Singleton(const Singleton&) = delete;
        Singleton& operator= (const Singleton&) = delete;
    };



    // define Singleton<DERIVED> type outside of class
    template<typename DERIVED> DERIVED* Singleton<DERIVED>::INSTANCE = nullptr;
}