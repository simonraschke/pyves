#pragma once

#include <memory>
#include <string>
#include <stdexcept>
#include <limits>



namespace _pyves
{
    template<typename... Args>
    constexpr bool all(Args... args)
    {
        return (... && args);
    }

    template<typename... Args>
    constexpr bool none(Args... args)
    {
        return !(... && args);
    }

    template<typename... Args>
    constexpr bool any(Args... args)
    {
        return (... || args);
    }

    template<typename ... Args>
    std::string string_format( const std::string& format, Args ... args )
    {
        std::size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
        if( size <= 0 )
        { 
            throw std::runtime_error( "Error during formatting." ); 
        }
        std::unique_ptr<char[]> buf( new char[ size ] ); 
        snprintf( buf.get(), size, format.c_str(), args ... );
        return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
    }
    
    
    
    template<typename T>
    constexpr auto make_nan() -> T 
    { 
        return std::numeric_limits<T>::signaling_NaN(); 
    }



#if __cplusplus < 201703L

    #include <type_traits>
    #include <typeinfo>
    #ifndef _MSC_VER
    #   include <cxxabi.h>
    #endif
    #include <memory>
    #include <string>
    #include <cstdlib>

    template <class T>
    std::string
    type_name()
    {
        typedef typename std::remove_reference<T>::type TR;
        std::unique_ptr<char, void(*)(void*)> own
            (
    #ifndef _MSC_VER
                    abi::__cxa_demangle(typeid(TR).name(), nullptr,
                                            nullptr, nullptr),
    #else
                    nullptr,
    #endif
                    std::free
            );
        std::string r = own != nullptr ? own.get() : typeid(TR).name();
        if (std::is_const<TR>::value)
            r += " const";
        if (std::is_volatile<TR>::value)
            r += " volatile";
        if (std::is_lvalue_reference<T>::value)
            r += "&";
        else if (std::is_rvalue_reference<T>::value)
            r += "&&";
        return r;
    }

#else

    #include <string_view>
    #include <string>

    template <typename T>
    constexpr auto type_name() noexcept 
    {
        std::string_view name, prefix, suffix;
    #ifdef __clang__
        name = __PRETTY_FUNCTION__;
        prefix = "auto _pyves::type_name() [T = ";
        suffix = "]";
    #elif defined(__GNUC__)
        name = __PRETTY_FUNCTION__;
        prefix = "constexpr auto _pyves::type_name() [with T = ";
        suffix = "]";
    #elif defined(_MSC_VER)
        name = __FUNCSIG__;
        prefix = "auto __cdecl type_name<";
        suffix = ">(void) noexcept";
    #endif
        name.remove_prefix(prefix.size());
        name.remove_suffix(suffix.size());
        return static_cast<std::string>(name);
    }
#endif
}