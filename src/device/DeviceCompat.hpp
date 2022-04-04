#ifndef DeviceCompat_hpp
#define DeviceCompat_hpp

#include <cstddef>
#include <tuple>

/** Empty definitions of device functions. */

#define __device__
#define __host__
#define __global__
#define DEVICE_CHECK(cmd)

/* @{ Memory management API */
#define deviceMemcpy
#define deviceMemcpy2D
#define H2D
#define D2H
#define D2D
/* @} */

/* @{ Event API */
#define deviceEvent
#define deviceEvenCreate
#define deviceEventRecord
#define deviceEventSynchronize
#define deviceEventElapsedTime
#define deviceEventDestroy
/* @} */

/* @{ Stream API */
#define deviceStream
#define deviceStreamSynchronize
/* @} */

/* @{ Device API */
#define _Device
#define deviceSynchronize
/* @} */

#define deviceLaunch

namespace mem {
namespace detail {
/** 1-dimensional allocation of given type on the device.
 *
 * @tparam T scalar type of the allocation.
 * @param[in] count number of element in allocation.
 * @return pointer to the allocation.
 */
template <typename T>
auto
device_allocate(size_t /* count */) -> T *
{
    return nullptr;
}

/** 2-dimensional, aligned allocation of given type on the device.
 *
 * @tparam T scalar type of the allocation.
 * @param[in] height number of rows in allocation.
 * @param[in] width number of columns in allocation.
 * @return pitch of allocation and pointer to it.
 *
 * We use `hipMallocPitch` for allocations of 2-dimensional memory chunks on the GPU.
 * This function aligns the allocation and returns the pitch, _i.e._ the number
 * of columns plus padding (if any), *and* the pointer to the allocated chunk.
 * This "drawing" should exemplify the terms used:
 *
 *              width      padding
 *          /           \ /        \
 *        / [o o o o o o | x x x x ]
 * height   [o o o o o o | x x x x ]
 *        \ [o o o o o o | x x x x ]
 *          \                      /
 *                    pitch
 */
template <typename T>
auto
device_allocate(size_t /* height */, size_t /* width */) -> std::tuple<size_t, T *>
{
    return {0, nullptr};
}

/** Deallocate memory chunk of given type on device.
 *
 * @tparam T type of memory block.
 * @param p the pointer to memory block.
 */
template <typename T>
auto
device_deallocate(T * /* p */) noexcept -> void
{
    // do nothing
}
}  // namespace detail
}  // namespace mem

#endif  // DeviceCompat_hpp
