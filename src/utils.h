#pragma once
#ifndef UTILS_H
#    define UTILS_H

#    include <Eigen/Dense>
#    include <ranges>

using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;

constexpr float PI = 3.14159265359;

#    define enumerate(v) std::views::zip(std::views::iota(0), v)
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
