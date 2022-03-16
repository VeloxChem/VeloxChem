#ifndef BufferMacros_hpp
#define BufferMacros_hpp

#define _X_XY_BUFFERS_AND_MDSPANS_WITH_TYPE(SizeTag, Type, TypeTag) \
    template <typename B>                                           \
    using Buffer##SizeTag##TypeTag = Buffer##SizeTag<Type, B>;      \
    template <typename B>                                           \
    using MDSpan##SizeTag##TypeTag = MDSpan##SizeTag<Type, B>;      \
    template <typename B>                                           \
    using ConstMDSpan##SizeTag##TypeTag = ConstMDSpan##SizeTag<Type, B>;

#define _X_XY_GENERATORS_WITH_TYPE(SizeTag, Type, TypeTag)                                   \
    template <typename B>                                                                    \
    auto constantBuffer(size_t nElements, Type fill_value)->Buffer##SizeTag##TypeTag<B>      \
    {                                                                                        \
        return Buffer##SizeTag##TypeTag<B>(nElements).setConstant(fill_value);               \
    }                                                                                        \
    template <typename B>                                                                    \
    auto zeroBuffer(size_t nElements)->Buffer##SizeTag##TypeTag<B>                           \
    {                                                                                        \
        return Buffer##SizeTag##TypeTag<B>(nElements).setConstant(Type{0});                  \
    }                                                                                        \
    template <typename B>                                                                    \
    auto randomBuffer(size_t nElements, Type lower, Type upper)->Buffer##SizeTag##TypeTag<B> \
    {                                                                                        \
        return Buffer##SizeTag##TypeTag<B>(nElements).setRandom(lower, upper);               \
    }

#define _X_XY_BUFFERS_AND_MDSPANS_WITH_BACKEND(SizeTag, Backend, BackendTag) \
    template <typename T>                                                    \
    using Buffer##BackendTag##SizeTag = Buffer##SizeTag<T, Backend>;         \
    template <typename T>                                                    \
    using MDSpan##BackendTag##SizeTag = MDSpan##SizeTag<T, Backend>;         \
    template <typename T>                                                    \
    using ConstMDSpan##BackendTag##SizeTag = ConstMDSpan##SizeTag<T, Backend>;

#define _X_XY_GENERATORS_WITH_BACKEND(SizeTag, Backend, BackendTag)                       \
    template <typename T>                                                                 \
    auto constantBuffer(size_t nElements, T fill_value)->Buffer##BackendTag##SizeTag<T>   \
    {                                                                                     \
        return Buffer##BackendTag##SizeTag<T>(nElements).setConstant(fill_value);         \
    }                                                                                     \
    template <typename T>                                                                 \
    auto zeroBuffer(size_t nElements)->Buffer##BackendTag##SizeTag<T>                     \
    {                                                                                     \
        return Buffer##BackendTag##SizeTag<T>(nElements).setConstant(T{0});               \
    }                                                                                     \
    template <typename T>                                                                 \
    auto randomBuffer(size_t nElements, T lower, T upper)->Buffer##BackendTag##SizeTag<T> \
    {                                                                                     \
        return Buffer##BackendTag##SizeTag<T>(nElements).setRandom(lower, upper);         \
    }

#define _X_XY_BUFFERS_AND_MDSPANS_WITH_BACKEND_AND_TYPE(SizeTag, Type, TypeTag, Backend, BackendTag) \
    using Buffer##BackendTag##SizeTag##TypeTag      = Buffer##SizeTag<Type, Backend>;                \
    using MDSpan##BackendTag##SizeTag##TypeTag      = MDSpan##SizeTag<Type, Backend>;                \
    using ConstMDSpan##BackendTag##SizeTag##TypeTag = ConstMDSpan##SizeTag<Type, Backend>;

#define _X_XY_GENERATORS_WITH_BACKEND_AND_TYPE(SizeTag, Type, TypeTag, Backend, BackendTag)           \
    auto constantBuffer(size_t nElements, Type fill_value)->Buffer##BackendTag##SizeTag##TypeTag      \
    {                                                                                                 \
        return Buffer##BackendTag##SizeTag##TypeTag(nElements).setConstant(fill_value);               \
    }                                                                                                 \
    auto zeroBuffer(size_t nElements)->Buffer##BackendTag##SizeTag##TypeTag                           \
    {                                                                                                 \
        return Buffer##BackendTag##SizeTag##TypeTag(nElements).setConstant(Type{0});                  \
    }                                                                                                 \
    auto randomBuffer(size_t nElements, Type lower, Type upper)->Buffer##BackendTag##SizeTag##TypeTag \
    {                                                                                                 \
        return Buffer##BackendTag##SizeTag##TypeTag(nElements).setRandom(lower, upper);               \
    }

#define VLX_X_XY_DEFINE_BUFFERS_MDSPANS(SizeTag)                                         \
    template <typename T, typename B>                                                    \
    using MDSpan##SizeTag = typename Buffer##SizeTag<T, B>::mdspan_type;                 \
    template <typename T, typename B>                                                    \
    using ConstMDSpan##SizeTag = typename Buffer##SizeTag<T, B>::const_mdspan_type;      \
    _X_XY_BUFFERS_AND_MDSPANS_WITH_TYPE(SizeTag, double, d)                              \
    _X_XY_BUFFERS_AND_MDSPANS_WITH_BACKEND(SizeTag, mem::Host, Host)                     \
    _X_XY_BUFFERS_AND_MDSPANS_WITH_BACKEND(SizeTag, mem::Device, Device)                 \
    _X_XY_BUFFERS_AND_MDSPANS_WITH_BACKEND_AND_TYPE(SizeTag, double, d, mem::Host, Host) \
    _X_XY_BUFFERS_AND_MDSPANS_WITH_BACKEND_AND_TYPE(SizeTag, double, d, mem::Device, Device)

#define _N_MY_XN_BUFFERS_AND_MDSPANS_WITH_TYPE(Size, SizeTag, Type, TypeTag) \
    template <typename B, size_t Size>                                       \
    using Buffer##SizeTag##TypeTag = Buffer##SizeTag<Type, B, Size>;         \
    template <typename B, size_t Size>                                       \
    using MDSpan##SizeTag##TypeTag = MDSpan##SizeTag<Type, B, Size>;         \
    template <typename B, size_t Size>                                       \
    using ConstMDSpan##SizeTag##TypeTag = ConstMDSpan##SizeTag<Type, B, Size>;

#define _N_MY_XN_BUFFERS_AND_MDSPANS_WITH_BACKEND(Size, SizeTag, Backend, BackendTag) \
    template <typename T, size_t Size>                                                \
    using Buffer##BackendTag##SizeTag = Buffer##SizeTag<T, Backend, Size>;            \
    template <typename T, size_t Size>                                                \
    using MDSpan##BackendTag##SizeTag = MDSpan##SizeTag<T, Backend, Size>;            \
    template <typename T, size_t Size>                                                \
    using ConstMDSpan##BackendTag##SizeTag = ConstMDSpan##SizeTag<T, Backend, Size>;

#define _N_MY_XN_BUFFERS_AND_MDSPANS_WITH_BACKEND_AND_TYPE(Size, SizeTag, Type, TypeTag, Backend, BackendTag) \
    template <size_t Size>                                                                                    \
    using Buffer##BackendTag##SizeTag##TypeTag = Buffer##SizeTag<Type, Backend, Size>;                        \
    template <size_t Size>                                                                                    \
    using MDSpan##BackendTag##SizeTag##TypeTag = MDSpan##SizeTag<Type, Backend, Size>;                        \
    template <size_t Size>                                                                                    \
    using ConstMDSpan##BackendTag##SizeTag##TypeTag = ConstMDSpan##SizeTag<Type, Backend, Size>;

#define VLX_N_MY_XN_DEFINE_BUFFERS_AND_MDSPANS(Size, SizeTag)                                     \
    template <typename T, typename B, size_t Size>                                                \
    using MDSpan##SizeTag = typename Buffer##SizeTag<T, B, Size>::mdspan_type;                    \
    template <typename T, typename B, size_t Size>                                                \
    using ConstMDSpan##SizeTag = typename Buffer##SizeTag<T, B, Size>::const_mdspan_type;         \
    _N_MY_XN_BUFFERS_AND_MDSPANS_WITH_TYPE(Size, SizeTag, double, d)                              \
    _N_MY_XN_BUFFERS_AND_MDSPANS_WITH_BACKEND(Size, SizeTag, mem::Host, Host)                     \
    _N_MY_XN_BUFFERS_AND_MDSPANS_WITH_BACKEND(Size, SizeTag, mem::Device, Device)                 \
    _N_MY_XN_BUFFERS_AND_MDSPANS_WITH_BACKEND_AND_TYPE(Size, SizeTag, double, d, mem::Host, Host) \
    _N_MY_XN_BUFFERS_AND_MDSPANS_WITH_BACKEND_AND_TYPE(Size, SizeTag, double, d, mem::Device, Device)

#define _MN_BUFFERS_AND_MDSPANS_WITH_TYPE(Type, TypeTag)       \
    template <typename B, size_t MRows, size_t NCols>          \
    using BufferMN##TypeTag = BufferMN<Type, B, MRows, NCols>; \
    template <typename B, size_t MRows, size_t NCols>          \
    using MDSpanMN##TypeTag = MDSpanMN<Type, B, MRows, NCols>; \
    template <typename B, size_t MRows, size_t NCols>          \
    using ConstMDSpanMN##TypeTag = ConstMDSpanMN<Type, B, MRows, NCols>;

#define _MN_BUFFERS_AND_MDSPANS_WITH_BACKEND(Backend, BackendTag)      \
    template <typename T, size_t MRows, size_t NCols>                  \
    using Buffer##BackendTag##MN = BufferMN<T, Backend, MRows, NCols>; \
    template <typename T, size_t MRows, size_t NCols>                  \
    using MDSpan##BackendTag##MN = MDSpanMN<T, Backend, MRows, NCols>; \
    template <typename T, size_t MRows, size_t NCols>                  \
    using ConstMDSpan##BackendTag##MN = ConstMDSpanMN<T, Backend, MRows, NCols>;

#define _MN_BUFFERS_AND_MDSPANS_WITH_BACKEND_AND_TYPE(Type, TypeTag, Backend, BackendTag) \
    template <size_t MRows, size_t NCols>                                                 \
    using Buffer##BackendTag##MN##TypeTag = BufferMN<Type, Backend, MRows, NCols>;        \
    template <size_t MRows, size_t NCols>                                                 \
    using MDSpan##BackendTag##MN##TypeTag = MDSpanMN<Type, Backend, MRows, NCols>;        \
    template <size_t MRows, size_t NCols>                                                 \
    using ConstMDSpan##BackendTag##MN##TypeTag = ConstMDSpanMN<Type, Backend, MRows, NCols>;

#define VLX_MN_DEFINE_BUFFERS_AND_MDSPANS()                                         \
    template <typename T, typename B, size_t MRows, size_t NCols>                   \
    using MDSpanMN = typename BufferMN<T, B, MRows, NCols>::mdspan_type;            \
    template <typename T, typename B, size_t MRows, size_t NCols>                   \
    using ConstMDSpanMN = typename BufferMN<T, B, MRows, NCols>::const_mdspan_type; \
    _MN_BUFFERS_AND_MDSPANS_WITH_TYPE(double, d)                                    \
    _MN_BUFFERS_AND_MDSPANS_WITH_BACKEND(mem::Host, Host)                           \
    _MN_BUFFERS_AND_MDSPANS_WITH_BACKEND(mem::Device, Device)                       \
    _MN_BUFFERS_AND_MDSPANS_WITH_BACKEND_AND_TYPE(double, d, mem::Host, Host)       \
    _MN_BUFFERS_AND_MDSPANS_WITH_BACKEND_AND_TYPE(double, d, mem::Device, Device)

#endif  // BufferMacros_hpp
