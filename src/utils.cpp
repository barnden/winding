#include "utils.h"

auto operator<<(std::ostream& ostream, Vec2 const& vec) -> std::ostream&
{
    return ostream << std::format("{}, {}", vec.x(), vec.y());
}

auto operator<<(std::ostream& ostream, Vec3 const& vec) -> std::ostream&
{
    return ostream << std::format("{}, {}, {}", vec.x(), vec.y(), vec.z());
}