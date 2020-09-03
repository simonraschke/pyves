#pragma once

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
}