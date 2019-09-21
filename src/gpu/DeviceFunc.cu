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
setDevice(const int32_t iDevice)
{
#ifdef ENABLE_GPU
    auto cerr = cudaSetDevice(iDevice);

    errors::assertMsgCritical(cerr == cudaSuccess, {"setCudaDevice"});
#endif
}

void
allocateDeviceMemory(void**  pointer,
                     size_t* dataPitch,
                     size_t  dataWidth,
                     size_t  dataHeight)
{
#ifdef ENABLE_GPU
    auto cerr = cudaMallocPitch(pointer, dataPitch, dataWidth, dataHeight);

    errors::assertMsgCritical(cerr == cudaSuccess, {"allocateDeviceMemory"});
#endif
}

void
freeDeviceMemory(void* pointerToMemory)
{
#ifdef ENABLE_GPU
    auto cerr = cudaFree(pointerToMemory);

    errors::assertMsgCritical(cerr == cudaSuccess, {"freeDeviceMemory"});
#endif
}

void
copyToDeviceMemory(      void*  destination,
                         size_t destinationPitch,
                   const void*  source,
                         size_t sourcePitch,
                         size_t dataWidth,
                         size_t dataHeight)
{
#ifdef ENABLE_GPU
    auto cerr = cudaMemcpy2D(destination, destinationPitch, source, sourcePitch, dataWidth, dataHeight,
                             cudaMemcpyHostToDevice);

    errors::assertMsgCritical(cerr == cudaSuccess, {"copyToDeviceMemory"});
#endif
}

void
copyFromDeviceMemory(     void*  destination,
                          size_t destinationPitch,
                    const void*  source,
                          size_t sourcePitch,
                          size_t dataWidth,
                          size_t dataHeight)
{
#ifdef ENABLE_GPU
    auto cerr = cudaMemcpy2D(destination, destinationPitch, source, sourcePitch, dataWidth, dataHeight,
                             cudaMemcpyDeviceToHost);

    errors::assertMsgCritical(cerr == cudaSuccess, {"copyFromDeviceMemory"});
#endif
}




}  // namespace gpu
