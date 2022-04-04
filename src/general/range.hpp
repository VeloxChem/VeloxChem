/* Copyright 2012–2020 Konrad Rudolph
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef UTIL_LANG_RANGE_HPP
#define UTIL_LANG_RANGE_HPP

#include <cmath>
#include <iterator>
#include <type_traits>

namespace util {
namespace lang {

namespace detail {

template <typename T>
struct range_iter_base : std::iterator<std::input_iterator_tag, T>
{
    __host__ __device__
    range_iter_base(T current)
        : current(current)
    {
    }

    __host__ __device__ T
    operator*() const
    {
        return current;
    }

    __host__ __device__ T const*
    operator->() const
    {
        return &current;
    }

    __host__ __device__ range_iter_base&
    operator++()
    {
        ++current;
        return *this;
    }

    __host__ __device__ range_iter_base
    operator++(int)
    {
        auto copy = *this;
        ++*this;
        return copy;
    }

    __host__ __device__ bool
    operator==(range_iter_base const& other) const
    {
        return current == other.current;
    }

    __host__ __device__ bool
    operator!=(range_iter_base const& other) const
    {
        return not(*this == other);
    }

   protected:
    T current;
};

}  // namespace detail

template <typename T>
struct step_range_proxy
{
    struct iterator : detail::range_iter_base<T>
    {
        __host__ __device__
        iterator(T current, T step)
            : detail::range_iter_base<T>(current), step_(step)
        {
        }

        using detail::range_iter_base<T>::current;

        __host__ __device__ iterator&
        operator++()
        {
            current += step_;
            return *this;
        }

        __host__ __device__ iterator
        operator++(int)
        {
            auto copy = *this;
            ++*this;
            return copy;
        }

        // Loses commutativity. Iterator-based ranges are simply broken. :-(
        __host__ __device__ bool
        operator==(iterator const& other) const
        {
            return step_ > 0 ? current >= other.current : current < other.current;
        }

        __host__ __device__ bool
        operator!=(iterator const& other) const
        {
            return not(*this == other);
        }

        T step_;
    };

    __host__ __device__
    step_range_proxy(T begin, T end, T step)
        : begin_(begin, step), end_(end, step)
    {
    }

    __host__ __device__ iterator
    begin() const
    {
        return begin_;
    }

    __host__ __device__ iterator
    end() const
    {
        return end_;
    }

    __host__ __device__ std::size_t
                        size() const
    {
        if (*end_ >= *begin_)
        {
            // Increasing and empty range
            if (begin_.step_ < T{0}) return 0;
        }
        else
        {
            // Decreasing range
            if (begin_.step_ > T{0}) return 0;
        }
        return std::ceil(std::abs(static_cast<double>(*end_ - *begin_) / begin_.step_));
    }

   private:
    iterator begin_;
    iterator end_;
};

template <typename T>
struct range_proxy
{
    struct iterator : detail::range_iter_base<T>
    {
        __host__ __device__
        iterator(T current)
            : detail::range_iter_base<T>(current)
        {
        }
    };

    __host__ __device__
    range_proxy(T begin, T end)
        : begin_(begin), end_(end)
    {
    }

    __host__ __device__ step_range_proxy<T>
                        step(T step)
    {
        return {*begin_, *end_, step};
    }

    __host__ __device__ iterator
    begin() const
    {
        return begin_;
    }

    __host__ __device__ iterator
    end() const
    {
        return end_;
    }

    __host__ __device__ std::size_t
                        size() const
    {
        return *end_ - *begin_;
    }

   private:
    iterator begin_;
    iterator end_;
};

template <typename T>
struct step_inf_range_proxy
{
    struct iterator : detail::range_iter_base<T>
    {
        __host__ __device__
        iterator(T current = T(), T step = T())
            : detail::range_iter_base<T>(current), step(step)
        {
        }

        using detail::range_iter_base<T>::current;

        __host__ __device__ iterator&
        operator++()
        {
            current += step;
            return *this;
        }

        __host__ __device__ iterator
        operator++(int)
        {
            auto copy = *this;
            ++*this;
            return copy;
        }

        __host__ __device__ bool
        operator==(iterator const&) const
        {
            return false;
        }

        __host__ __device__ bool
        operator!=(iterator const&) const
        {
            return true;
        }

       private:
        T step;
    };

    __host__ __device__
    step_inf_range_proxy(T begin, T step)
        : begin_(begin, step)
    {
    }

    __host__ __device__ iterator
    begin() const
    {
        return begin_;
    }

    __host__ __device__ iterator
    end() const
    {
        return iterator();
    }

   private:
    iterator begin_;
};

template <typename T>
struct infinite_range_proxy
{
    struct iterator : detail::range_iter_base<T>
    {
        __host__ __device__
        iterator(T current = T())
            : detail::range_iter_base<T>(current)
        {
        }

        __host__ __device__ bool
        operator==(iterator const&) const
        {
            return false;
        }

        __host__ __device__ bool
        operator!=(iterator const&) const
        {
            return true;
        }
    };

    __host__ __device__
    infinite_range_proxy(T begin)
        : begin_(begin)
    {
    }

    __host__ __device__ step_inf_range_proxy<T>
                        step(T step)
    {
        return {*begin_, step};
    }

    __host__ __device__ iterator
    begin() const
    {
        return begin_;
    }

    __host__ __device__ iterator
    end() const
    {
        return iterator();
    }

   private:
    iterator begin_;
};

template <typename T, typename U>
__host__ __device__ auto
range(T begin, U end) -> range_proxy<typename std::common_type<T, U>::type>
{
    using C = typename std::common_type<T, U>::type;
    return {static_cast<C>(begin), static_cast<C>(end)};
}

template <typename T>
__host__ __device__ infinite_range_proxy<T>
                    range(T begin)
{
    return {begin};
}

namespace traits {

template <typename C>
struct has_size
{
    template <typename T>
    static auto check(T*) -> typename std::is_integral<decltype(std::declval<T const>().size())>::type;

    template <typename>
    static auto check(...) -> std::false_type;

    using type                  = decltype(check<C>(0));
    static constexpr bool value = type::value;
};

}  // namespace traits

template <typename C, typename = typename std::enable_if<traits::has_size<C>::value>>
__host__ __device__ auto
indices(C const& cont) -> range_proxy<decltype(cont.size())>
{
    return {0, cont.size()};
}

template <typename T, std::size_t N>
__host__ __device__ range_proxy<std::size_t> indices(T (&)[N])
{
    return {0, N};
}

template <typename T>
__host__ __device__ range_proxy<typename std::initializer_list<T>::size_type>
                    indices(std::initializer_list<T>&& cont)
{
    return {0, cont.size()};
}

}  // namespace lang
}  // namespace util

#endif  // ndef UTIL_LANG_RANGE_HPP
