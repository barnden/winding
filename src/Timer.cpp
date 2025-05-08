#include <print>

#include "Timer.h"

using namespace std::chrono;

Timer::Timer(std::string_view name, bool print_on_destruct)
    : m_name(name)
    , m_print_on_destruct(print_on_destruct)
{
    m_start = steady_clock::now();
}

auto Timer::elapsed() const -> double
{
    auto const now = steady_clock::now();
    auto const duration = duration_cast<microseconds>(now - m_start);

    return duration.count();
}

void Timer::print() const
{
    std::println("[Timer] {}: {} us", m_name, elapsed());
}

Timer::~Timer()
{
    if (m_print_on_destruct)
        print();
}