#ifndef DeviceUtils_hpp
#define DeviceUtils_hpp

namespace device {
enum class Dim
{
    X,
    Y,
    Z
};

template <Dim d>
__device__ inline auto
globalThreadIndex() -> unsigned int
{
    return 0;
}

template <Dim d>
__device__ inline auto
globalThreadCount() -> unsigned int
{
    return 0;
}
}  // namespace device

#endif  // DeviceUtils_hpp
