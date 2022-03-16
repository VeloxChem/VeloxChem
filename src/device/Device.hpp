#ifndef Device_hpp
#define Device_hpp

#include <cstddef>
#include <string>
#include <tuple>
#include <type_traits>

/** @file Compatibility layer for device acceleration. */

#if defined VLX_HAS_HIP
#include "hip/HIPDevice.hpp"
#elif defined VLX_HAS_CUDA
#include "cuda/CUDADevice.hpp"
#else
#include "DeviceCompat.hpp"
#endif

#include "DeviceUtils.hpp"
#include "range.hpp"

namespace mem {
/**
 * Tag struct for device backend.
 */
struct Device
{
    inline const static std::string name{"Device"};
};

template <typename B>
inline constexpr bool is_on_device_v = std::is_same_v<B, Device>;
}  // namespace mem

// namespace device {
// template <typename T>
// using step_range = util::lang::step_range_proxy<T>;
//
///**
// * adapted from: https://github.com/harrism/hemi/blob/master/hemi/grid_stride_range.h
// */
// template <typename T, Dim dir = Dim::X, typename U = T, typename V = std::common_type_t<T, U>>
//__host__ __device__ inline auto
// grid_stride_range(T begin, U end) -> step_range<V>
//{
//    // begin at the global thread index
//    begin += globalThreadIndex<dir>();
//    // step of global thread count
//    return step_range<V>(static_cast<V>(begin), static_cast<V>(end), globalThreadCount<dir>());
//}
//
// template <typename T, typename U = T>
//__host__ __device__ inline auto
// grid_stride_range_x(T begin, U end) -> step_range<std::common_type_t<T, U>>
//{
//    return grid_stride_range<T, Dim::X>(begin, end);
//}
//
// template <typename T, typename U = T>
//__host__ __device__ inline auto
// grid_stride_range_y(T begin, U end) -> step_range<std::common_type_t<T, U>>
//{
//    return grid_stride_range<T, Dim::Y>(begin, end);
//}
//
// template <typename T, typename U = T>
//__host__ __device__ inline auto
// grid_stride_range_z(T begin, U end) -> step_range<std::common_type_t<T, U>>
//{
//    return grid_stride_range<T, Dim::Z>(begin, end);
//}
//
// template <typename T>
//__global__ auto
// full1D(T *devPtr, T fill_value, size_t sz) -> void
//{
//    for (auto i : grid_stride_range(0, sz))
//    {
//        devPtr[i] = fill_value;
//    }
//}
//
// template <typename T>
//__global__ auto
// full2D(T *devPtr, T fill_value, size_t nRows, size_t nCols, size_t nPad) ->
//{
//    for (auto i : grid_stride_range_x(0, nRows))
//    {
//        for (auto j : grid_stride_range_y(0, nCols))
//        {
//            devPtr[i * nPad + j] = fill_value;
//        }
//    }
//}
//}  // namespace device

#endif  // Device_hpp
