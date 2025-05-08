#pragma once

#include <Eigen/Dense>
#include <format>
#include <ranges>

static constexpr double Infinity = std::numeric_limits<double>::infinity();

using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;

auto operator<<(std::ostream& ostream, Vec2 const& vec) -> std::ostream&;
auto operator<<(std::ostream& ostream, Vec3 const& vec) -> std::ostream&;

// Handle std::format_string for Eigen vector types
template <>
struct std::formatter<Vec2> : std::formatter<std::string> {
    auto format(Vec2 p, std::format_context& ctx) const
    {
        return std::formatter<string>::format(std::format("{}, {}", p.x(), p.y()), ctx);
    }
};

template <>
struct std::formatter<Vec3> : std::formatter<std::string> {
    auto format(Vec3 p, std::format_context& ctx) const
    {
        return std::formatter<string>::format(std::format("{}, {}, {}", p.x(), p.y(), p.z()), ctx);
    }
};

#if __cpp_lib_ranges_enumerate != 202302L
/**
    std::views::enumerate is a C++23 feature (P2164R4).
    However, this feature has yet to be implemented into upstream clang (PR 73617).
 */
#    define enumerate(v) std::views::zip(std::views::iota(0), v)
#else
#    define enumerate(v) std::views::enumerate(v)
#endif

#ifdef DEBUG_FPE_TRAP
// See: https://stackoverflow.com/questions/69059981
#    include <fenv.h>
#    include <math.h>
#    include <signal.h>
#    include <stdio.h>
#    include <stdlib.h>

static void fpe_signal_handler(int sig, siginfo_t* sip, void* scp)
{
    int fe_code = sip->si_code;

    printf("In signal handler : ");

    if (fe_code == ILL_ILLTRP)
        printf("Illegal trap detected\n");
    else
        printf("Code detected : %d\n", fe_code);

    abort();
}

void enable_floating_point_exceptions()
{
    fenv_t env;
    fegetenv(&env);

    env.__fpcr = env.__fpcr | __fpcr_trap_invalid;
    fesetenv(&env);

    struct sigaction act;
    act.sa_sigaction = fpe_signal_handler;
    sigemptyset(&act.sa_mask);
    act.sa_flags = SA_SIGINFO;
    sigaction(SIGILL, &act, NULL);
}
#endif
