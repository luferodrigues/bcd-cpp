// Random numbers for the simulation.
//
// The simulation code never talks to a generator directly; it asks an `Rng` for numbers. `Rng`
// only describes WHICH functions exist (like an abstract base class in Python, abc.ABC); a
// concrete class supplies them:
//   - MtRng (below): the real thing, a seeded Mersenne Twister. Used by patchy_bcd.
//   - ReplayRng (tests/check_tool.cpp): hands out numbers recorded from a run of the Python code,
//     so C++ and Python can be compared step by step with identical random numbers.
#pragma once

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

class Rng {
public:
    virtual ~Rng() = default;
    // A number in [low, high), like np.random.uniform(low, high).
    virtual double uniform(double low, double high) = 0;
    // A normally distributed number, like np.random.normal(mean, sigma).
    virtual double normal(double mean, double sigma) = 0;
    // Random order of v's elements, like np.random.shuffle(v).
    virtual void shuffle(std::vector<int>& v) = 0;
};

class MtRng : public Rng {
public:
    explicit MtRng(std::uint64_t seed) : gen_(seed) {}

    double uniform(double low, double high) override {
        return std::uniform_real_distribution<double>(low, high)(gen_);
    }
    double normal(double mean, double sigma) override {
        if (!(sigma > 0)) return mean; // np.random.normal(m, 0) returns m
        return mean + sigma * std_normal_(gen_);
    }
    void shuffle(std::vector<int>& v) override { std::shuffle(v.begin(), v.end(), gen_); }

private:
    std::mt19937_64 gen_;
    std::normal_distribution<double> std_normal_{0.0, 1.0};
};
