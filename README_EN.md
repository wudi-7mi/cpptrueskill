# C++ TrueSkill Through Time

En / [中文](README.md)

This C++17 header-only library provides both classic TrueSkill and TrueSkill Through Time. The TTT implementation follows [TrueSkillThroughTime.py](https://github.com/glandfried/TrueSkillThroughTime.py) 2.0, including temporal inference, Ordinal/Continuous/Discrete observations, online updates, and model evidence.

## Usage Example

Please refer to the [example.cpp](example.cpp) file.

### 1v1 Match

```cpp
// Initialize
trueskill::Rating p1(80, 2.5);
trueskill::Rating p2(75, 2.5);

// Calculate match quality
double quality1v1 = trueskill::quality_1vs1(p1, p2);
std::cout << "quality1v1: " << quality1v1 << std::endl;

// Calculate match result
auto result1v1 = trueskill::rate_1vs1(p1, p2);
for (const auto& team : result1v1) {
    for (const auto& rating : team) {
        std::cout << rating << " ";
    }
    std::cout << std::endl;
}

// If the match is a draw
auto result1v1_drawn = trueskill::rate_1vs1(p1, p2, true);
for (const auto& team : result1v1_drawn) {
    for (const auto& rating : team) {
        std::cout << rating << " ";
    }
    std::cout << std::endl;
}
```

### Multi-team Match

```cpp
// Initialize
trueskill::TrueSkill ts(25.0, 8.333333, 4.2, 0.3, 0.1);
std::cout << ts << std::endl;

trueskill::Rating alice(34.23, 7.22);
trueskill::Rating bob(24.23, 3.22);   
trueskill::Rating charlie(12.23, 3.9);
trueskill::Rating david(26.23, 6.11);
trueskill::Rating eve(41.2, 4.22);
trueskill::Rating frank(24.23, 5.08);
trueskill::Rating gloria(4.23, 3.08);

// Assign weights and teams
std::vector<std::vector<double>> weights = {{1.0, 1.0}, {0.8, 1.0}, {1.0, 1.0, 1.0}};
std::vector<std::vector<trueskill::Rating>> teams = {{alice, bob}, {charlie, david}, {eve, frank, gloria}};

// Calculate match quality
double quality = ts.quality(teams, weights);
std::cout << "quality: " << quality << std::endl;

// Calculate match result
auto newTeams = ts.rate(teams, {2, 1, 3}, weights);
for (const auto& team : newTeams) {
    for (const auto& rating : team) {
        std::cout << rating << " ";
    }
    std::cout << std::endl;
}
```

### TTT Game Model

`Game` exposes each player's likelihood and posterior together with the evidence of the observed result. These values form the basis of temporal message propagation.

```cpp
trueskill::Player alice(trueskill::math::Gaussian(25, 25.0 / 3), 25.0 / 6, 0.03);
trueskill::Player bob(trueskill::math::Gaussian(25, 25.0 / 3), 25.0 / 6, 0.03);

// As in Python TrueSkillThroughTime, a larger result means a better placement.
trueskill::Game game({{alice}, {bob}}, {1, 0});
auto posteriors = game.posteriors();
double evidence = game.evidence;
```

Note that classic `TrueSkill::rate()` accepts ranks where a smaller value is better, while TTT `Game` results use larger values for better placements.

### Historical Inference

```cpp
trueskill::History::Composition composition = {
    {{"alice"}, {"bob"}},
    {{"bob"}, {"charlie"}},
    {{"charlie"}, {"alice"}}
};
trueskill::History::Results results = {{1, 0}, {1, 0}, {1, 0}};

trueskill::History history(composition, results, {1, 2, 3});
history.convergence();
auto curves = history.learning_curves();
double log_evidence = history.log_evidence();
```

`History` groups games by time, propagates uncertainty using `gamma * sqrt(elapsed_time)`, and iterates forward and backward messages. Custom priors, per-player `beta/gamma`, team weights, and simultaneous games are supported.

### Mixed Observations and Incremental Updates

```cpp
std::vector<trueskill::Observation> observations = {
    trueskill::Observation::Ordinal,
    trueskill::Observation::Continuous,
    trueskill::Observation::Discrete
};

trueskill::History online_history(
    composition, results, {1, 2, 3}, {},
    0, 3, 1, 0.15, 0, true, {}, observations);
online_history.convergence();
online_history.add_history({{{"alice"}, {"charlie"}}}, {{2, 0}});
online_history.convergence();
```

- `Ordinal`: wins, rankings, and draws.
- `Continuous`: continuous differences such as time or distance.
- `Discrete`: Poisson counts such as goals or points, using a fixed-point Gaussian approximation.
- `online=true`: freezes online estimates per revealed time batch while retaining smoothed posteriors.

### State Persistence

`History::save(std::ostream&)` and `History::load(std::istream&)` preserve the complete converged state, including forward/backward messages, likelihoods, online estimates, and evidence. The current versioned native binary format is intended for restoring state within the same build ecosystem.

## Compile and Run

Use xmake to compile and run:

```shell
# Test
xmake build test_main
xmake run test_main

# Example
xmake build example
xmake run example
```
