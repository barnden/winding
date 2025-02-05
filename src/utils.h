#pragma once
#ifndef UTILS_H
#    define UTILS_H

#    include <Eigen/Dense>
#    include <ranges>

using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;

constexpr double PI = 3.141592653589793;

#define __enumerate(_1, _2, NAME, ...) NAME
#define __enumerate_start_index(v, start) std::views::zip(std::views::iota(start), v)
#define __enumerate_zero(v) __enumerate_start_index(v, 0)
#define enumerate(...) __enumerate(__VA_ARGS__, __enumerate_start_index, __enumerate_zero)(__VA_ARGS__)

struct Options {
    std::string data_path = "./reference";
    std::string file_stem = "vase";
    std::string out_path = "./output";
    std::string experiment = "reference";

    // If pull_in and push_out then alternate pull_in - push_out
    bool push_out = false;
    bool pull_in = false;

    bool use_bvh = false;
};

#    ifdef DEBUG_FPE_TRAP
#        include <fenv.h>
#        include <math.h>
#        include <signal.h>
#        include <stdio.h>
#        include <stdlib.h>

static void
fpe_signal_handler(int sig, siginfo_t* sip, void* scp)
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
#    endif
#endif
