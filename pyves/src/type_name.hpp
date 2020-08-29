#pragma once


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
        prefix = "auto type_name() [T = ";
        suffix = "]";
    #elif defined(__GNUC__)
        name = __PRETTY_FUNCTION__;
        prefix = "constexpr auto type_name() [with T = ";
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