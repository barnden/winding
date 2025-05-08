#include <chrono>
#include <string_view>

class Timer {
    decltype(std::chrono::steady_clock::now()) m_start;
    std::string m_name;
    bool m_print_on_destruct;

public:
    Timer(std::string_view name, bool print_on_destruct = true);

    auto elapsed() const -> double;
    void print() const;
    ~Timer();
};