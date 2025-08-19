//
// Created by Cem Akarsubasi on 8/18/25.
//

#ifndef GEL_RANDOM_H
#define GEL_RANDOM_H

#include <cassert>
#include <type_traits>
#include <cstdint>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

/// Random number generation
namespace CGLA::Random
{
namespace detail
{
    /// PCG Random Number Generation. Adapted from code by Melissa O'Neill <oneill@pcg-random.org>
    /// under the Apache 2.0 License.
    ///
    /// For additional information about the PCG random number generation scheme,
    /// including its license and other licensing options, visit
    ///
    ///     http://www.pcg-random.org
    struct Pcg32Random {
    private:
        std::uint64_t state = 0x853c49e6748fea9bULL;
        std::uint64_t inc = 0xda3e39cb94b95bdbULL;

    public:
        constexpr Pcg32Random() = default;
        constexpr Pcg32Random(const std::uint64_t seed, const std::uint64_t seq) { pcg32_seed(seed, seq); }

        /// Seed the rng.  Specified in two parts, state initializer and a
        /// sequence selection constant (a.k.a. stream id)
        constexpr void pcg32_seed(std::uint64_t initstate, std::uint64_t initseq)
        {
            this->state = 0U;
            this->inc = (initseq << 1u) | 1u;
            pcg32_rand();
            this->state += initstate;
            pcg32_rand();
        }

        /// Generate a uniformly distributed 32-bit random number
        constexpr std::uint32_t pcg32_rand()
        {
            const std::uint64_t oldstate = this->state;
            this->state = oldstate * 6364136223846793005ULL + this->inc;
            const std::uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
            const std::uint32_t rot = oldstate >> 59u;
            return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        }

        /// Generate a uniformly distributed number, r, where 0 <= r < bound
        constexpr std::uint32_t pcg32_rand_bounded(std::uint32_t bound)
        {
            // To avoid bias, we need to make the range of the RNG a multiple of
            // bound, which we do by dropping output less than a threshold.
            // A naive scheme to calculate the threshold would be to do
            //
            //     uint32_t threshold = 0x100000000ull % bound;
            //
            // but 64-bit div/mod is slower than 32-bit div/mod (especially on
            // 32-bit platforms).  In essence, we do
            //
            //     uint32_t threshold = (0x100000000ull-bound) % bound;
            //
            // because this version will calculate the same modulus, but the LHS
            // value is less than 2^32.

            std::uint32_t threshold = -bound % bound;

            // Uniformity guarantees that this loop will terminate.  In practice, it
            // should usually terminate quickly; on average (assuming all bounds are
            // equally likely), 82.25% of the time, we can expect it to require just
            // one iteration.  In the worst case, someone passes a bound of 2^31 + 1
            // (i.e., 2147483649), which invalidates almost 50% of the range.  In
            // practice, bounds are typically small and only a tiny amount of the range
            // is eliminated.
            for (;;) {
                std::uint32_t r = pcg32_rand();
                if (r >= threshold)
                    return r % bound;
            }
        }

        constexpr std::uint32_t pcg32_rand_bounded_two_sided(std::uint32_t lower, std::uint32_t upper)
        {
            assert(lower <= upper);
            return lower + pcg32_rand_bounded(upper - lower);
        }
    };
}

/// Fast and high-quality PRNG that implements std::uniform_random_bit_generator
template <std::uint64_t lower = 0ULL, std::uint64_t upper = ~0ULL>
struct GelPrngBase {
    using result_type = std::uint64_t;
    static_assert(lower <= upper);

private:
    detail::Pcg32Random rng[2];

public:
    /// Initialize a random generator with the default seed
    constexpr GelPrngBase()
    {
        pcg32x2_seed(42, 42, 42, 42);
    }

    /// Initialize a random generator with the given seed
    explicit constexpr GelPrngBase(const std::array<std::uint64_t, 4>& seed)
    {
        pcg32x2_seed(seed[0], seed[1], seed[2], seed[3]);
    }

    [[nodiscard]]
    static constexpr result_type min() { return lower; }

    [[nodiscard]]
    static constexpr result_type max() { return upper; }

    /// Generate the next random number
    [[nodiscard]]
    constexpr result_type operator()()
    {
        return pcg32x2_rand_bounded_two_sided(lower, upper);
    }

private:
    constexpr void pcg32x2_seed(std::uint64_t seed1, std::uint64_t seed2, std::uint64_t seq1, std::uint64_t seq2)
    {
        uint64_t mask = ~0ull >> 1;
        // The stream for each of the two generators *must* be distinct
        if ((seq1 & mask) == (seq2 & mask))
            seq2 = ~seq2;
        rng[0].pcg32_seed(seed1, seq1);
        rng[1].pcg32_seed(seed2, seq2);
    }

    constexpr std::uint64_t pcg32x2_rand()
    {
        return (static_cast<std::uint64_t>(rng[0].pcg32_rand()) << 32)
            | rng[1].pcg32_rand();
    }

    constexpr std::uint64_t pcg32x2_rand_bounded(std::uint64_t bound)
    {
        std::uint64_t threshold = -bound % bound;
        for (;;) {
            std::uint64_t r = pcg32x2_rand();
            if (r >= threshold)
                return r % bound;
        }
    }

    constexpr std::uint64_t pcg32x2_rand_bounded_two_sided(std::uint64_t lower_bound, std::uint64_t upper_bound)
    {
        assert(lower_bound <= upper_bound);
        return lower_bound + pcg32x2_rand_bounded(upper_bound - lower_bound);
    }
};

namespace Distributions
{
    /// The distribution concept
    template <typename T>
    concept Distribution =
        std::is_default_constructible_v<T> && std::is_move_constructible_v<T> && requires(GelPrngBase<> rng, T t)
        {
            typename T::result_type;
            { t.min() } -> std::same_as<typename T::result_type>;
            { t.max() } -> std::same_as<typename T::result_type>;
            { t(rng) } -> std::same_as<typename T::result_type>;
        };

    struct UniformUint {
        using result_type = std::uint64_t;

    private:
        result_type m_min = 0ULL;
        result_type m_max = ~0ULL;

    public:
        constexpr UniformUint() = default;

        constexpr UniformUint(const result_type _lower, const result_type _upper) : m_min(_lower), m_max(_upper)
        {
            if (_lower > _upper)
                throw std::invalid_argument("lower bound must be less than upper bound");
        }

        [[nodiscard]]
        constexpr result_type min() const { return m_min; }

        [[nodiscard]]
        constexpr result_type max() const { return m_max; }

        [[nodiscard]]
        constexpr result_type operator()(GelPrngBase<>& rng) const
        {
            return rng();
        }
    };

    static_assert(Distribution<UniformUint>);

    struct UniformDouble {
        using result_type = double;

        [[nodiscard]]
        static constexpr result_type min() { return 0.0; }

        [[nodiscard]]
        static constexpr result_type max() { return 1.0; }

        [[nodiscard]]
        result_type operator()(GelPrngBase<>& rng) { return std::ldexp(rng(), -64); }
    };

    static_assert(Distribution<UniformDouble>);

    enum struct NormalGenerator {
        BoxMueller
    };

    template <NormalGenerator gen = NormalGenerator::BoxMueller>
    struct Gaussian {
        using result_type = double;
        constexpr Gaussian() = default;
        static constexpr result_type min() { return -std::numeric_limits<double>::infinity(); }

        static constexpr result_type max()
        {
            return std::numeric_limits<double>::infinity();;
        }

        result_type operator()(GelPrngBase<>& rng) const
        {
            UniformDouble uinf;
            const auto x = uinf(rng);
            const auto y = uinf(rng);
            const auto z = std::sqrt(-2 * std::log(x)) * std::cos(2 * std::numbers::pi * y);
            return z;
        }
    };

    static_assert(Distribution<Gaussian<>>);

    template <NormalGenerator gen = NormalGenerator::BoxMueller>
    struct Normal {
        using result_type = double;
        double mean = 0.0;
        double stddev = 1.0;
        constexpr Normal() = default;
        constexpr Normal(const double mean, const double stddev) : mean(mean), stddev(stddev) {}
        static constexpr result_type min() { return -std::numeric_limits<double>::infinity(); }

        static constexpr result_type max()
        {
            return std::numeric_limits<double>::infinity();;
        }

        result_type operator()(GelPrngBase<>& rng) const
        {
            Gaussian<gen> g;
            return mean + stddev * g(rng);
        }
    };

    static_assert(Distribution<Normal<>>);
}

/// Core template for all GEL pseudo-random generators
template <Distributions::Distribution Distribution>
// We inherit from Distribution to make use of empty-base optimization.
struct GelPrngWithDistribution : Distribution {
    using result_type = Distribution::result_type;

private:
    GelPrngBase<> rng;

public:
    explicit GelPrngWithDistribution() = default;
    explicit GelPrngWithDistribution(Distribution&& distribution) : Distribution(std::move(distribution)) {}
    explicit GelPrngWithDistribution(const std::array<std::uint64_t, 4>& seed) : rng(seed) {}

    explicit GelPrngWithDistribution(const std::array<std::uint64_t, 4>& seed,
                                     Distribution&& distribution) : Distribution(std::move(distribution)), rng(seed) {}

    [[nodiscard]]
    result_type min() const { return Distribution::min(); }
    [[nodiscard]]
    result_type max() const { return Distribution::max(); }
    /// Generate the next random number
    result_type operator()() { return Distribution::operator()(rng); }
};

/// Uniform double number generator [0.0, 1.0)
using GelPrngDouble = GelPrngWithDistribution<Distributions::UniformDouble>;

/// Uniform uint generator
using GelPrng = GelPrngWithDistribution<Distributions::UniformUint>;

/// Gaussian double generator
using GelPrngGaussian = GelPrngWithDistribution<Distributions::Gaussian<>>;

/// Normal double generator
using GelPrngNormal = GelPrngWithDistribution<Distributions::Normal<>>;

/// Bit generator for std algorithms
using GelBitGenerator = GelPrngBase<>;

/// Function that seeds the GEL pseudo-random number generator
//[[deprecated]]
void gel_srand(unsigned int seed);

//// GEL provides a PCG pseudo-random number
/// generator which is optimized for speed. This version allows
/// an integer argument which is useful for grid-based noise
/// functions.
unsigned int gel_rand(unsigned int k);

/// GEL provides a PCG pseudo-random number
/// generator which is optimized for speed. This means
/// that GEL_RAND_MAX==UINT_MAX.
//[[deprecated]]
unsigned int gel_rand();
}

namespace CGLA
{
// Reexports

using UniformIntegerDistribution = Random::Distributions::UniformUint;
using UniformDoubleDistribution = Random::Distributions::UniformDouble;
using GaussianDistribution = Random::Distributions::Gaussian<>;
using NormalDistribution = Random::Distributions::Normal<>;
using Random::GelPrng;
using Random::GelPrngDouble;
using Random::GelPrngGaussian;
using Random::GelPrngNormal;
using Random::GelPrngWithDistribution;
using Random::GelBitGenerator;
using Random::gel_rand;
using Random::gel_srand;
}

#endif //GEL_RANDOM_H
