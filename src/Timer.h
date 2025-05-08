/*
 * Copyright (c) 2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include <chrono>
#include <string_view>

class Timer {
    decltype(std::chrono::steady_clock::now()) m_start;
    std::string m_name;
    bool m_print_on_destruct;

public:
    Timer(std::string_view name, bool print_on_destruct = true);

    auto elapsed() const;
    void print() const;
    ~Timer();
};