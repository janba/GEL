/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/CGLA/CGLA-util.h>
#include <cstdint>
#include <array>
#include <random>

namespace CGLA
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

struct GelPrng {
    using result_type = std::uint64_t;

private:
    Pcg32Random rng[2];
    result_type m_min = 0;
    result_type m_max = std::numeric_limits<result_type>::max();;

public:
    /// Initialize a random generator with the default seed
    constexpr GelPrng()
    {
        pcg32x2_seed(42, 42, 42, 42);
    }

    /// Initialize a random generator with the given seed
    explicit constexpr GelPrng(const std::array<std::uint64_t, 4>& seed)
    {
        pcg32x2_seed(seed[0], seed[1], seed[2], seed[3]);
    }

    /// Initialize a random number generator with the given bounds and seed
    constexpr GelPrng(const std::array<std::uint64_t, 2> bounds,
                      const std::array<std::uint64_t, 4>& seed) : m_min(bounds[0]), m_max(bounds[1])
    {
        pcg32x2_seed(seed[0], seed[1], seed[2], seed[3]);
    }

    [[nodiscard]]
    constexpr result_type min() const { return m_min; }

    [[nodiscard]]
    constexpr result_type max() const { return m_max; }

    [[nodiscard]]
    constexpr result_type operator()()
    {
        return pcg32x2_rand_bounded_two_sided(m_min, m_max);
    }


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

// TODO: uniform_random_bit_generator requires min() and max() to be constant expressions, so we need a templated version
// of this still
// static_assert(std::uniform_random_bit_generator<GelPrng>);

template <std::uint64_t lower_bound, std::uint64_t upper_bound>
struct GelPrngBits {
private:
    GelPrng rng;

public:
    constexpr GelPrngBits() : rng({lower_bound, upper_bound}, {42, 42, 42, 42}) {}
    explicit constexpr GelPrngBits(const std::array<std::uint64_t, 4>& seed) : rng({lower_bound, upper_bound}, seed) {}
    using result_type = std::uint64_t;

    [[nodiscard]]
    static constexpr result_type min() { return lower_bound; }

    [[nodiscard]]
    static constexpr result_type max() { return upper_bound; }

    constexpr result_type operator()()
    {
        return rng();
    }
};

static_assert(std::uniform_random_bit_generator<GelPrngBits<0, 1>>);

/// Uniform double number generator [0.0, 1.0)
struct GelPrngUniformDouble {
    using result_type = double;
    GelPrng rng;
    constexpr GelPrngUniformDouble() = default;
    explicit constexpr GelPrngUniformDouble(const std::array<std::uint64_t, 4>& seed) : rng(seed) {}
    static constexpr result_type min() { return 0.0; }
    static constexpr result_type max() { return 1.0; }
    result_type operator()() { return std::ldexp(rng(), -64); }
};

thread_local Pcg32Random pcg32_local = {0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL};

void gel_srand(unsigned int s)
{
    pcg32_local.pcg32_seed(s, s);
}

unsigned int gel_rand(unsigned int k)
{
    return gel_rand();
}

unsigned int gel_rand()
{
    return pcg32_local.pcg32_rand();
}
}
