//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "DeviceFunc.hpp"

#include "ErrorHandler.hpp"

namespace gpu {  // gpu namespace

void
set_device(const int32_t iDevice)
{
#ifdef ENABLE_GPU
    auto cerr = cudaSetDevice(iDevice);

    errors::assertMsgCritical(cerr == cudaSuccess, {"set_device"});
#endif
}

void
synchronize_device()
{
#ifdef ENABLE_GPU
    auto cerr = cudaDeviceSynchronize();

    errors::assertMsgCritical(cerr == cudaSuccess, {"synchronize_device"});
#endif
}

void
allocate_device_memory(void**  pointer,
                       size_t* dataPitch,
                       size_t  dataWidth,
                       size_t  dataHeight)
{
#ifdef ENABLE_GPU
    auto cerr = cudaMallocPitch(pointer, dataPitch, dataWidth, dataHeight);

    errors::assertMsgCritical(cerr == cudaSuccess, {"allocate_device_memory"});
#endif
}

void
free_device_memory(void* pointer)
{
#ifdef ENABLE_GPU
    auto cerr = cudaFree(pointer);

    errors::assertMsgCritical(cerr == cudaSuccess, {"free_device_memory"});
#endif
}

void
copy_to_device_memory(      void* destination,
                          size_t  destinationPitch,
                      const void* source,
                          size_t  sourcePitch,
                          size_t  dataWidth,
                          size_t  dataHeight)
{
#ifdef ENABLE_GPU
    auto cerr = cudaMemcpy2D(destination, destinationPitch, source, sourcePitch, dataWidth, dataHeight,
                             cudaMemcpyHostToDevice);

    errors::assertMsgCritical(cerr == cudaSuccess, {"copy_to_device_memory"});
#endif
}

void
copy_from_device_memory(      void* destination,
                            size_t  destinationPitch,
                        const void* source,
                            size_t  sourcePitch,
                            size_t  dataWidth,
                            size_t  dataHeight)
{
#ifdef ENABLE_GPU
    auto cerr = cudaMemcpy2D(destination, destinationPitch, source, sourcePitch, dataWidth, dataHeight,
                             cudaMemcpyDeviceToHost);

    errors::assertMsgCritical(cerr == cudaSuccess, {"copy_from_device_memory"});
#endif
}

}  // namespace gpu
