# cpp TrueSkill

En / [中文](README.md)

This is a header-only C++ implementation of TrueSkill, fully based on a [Python implementation](https://github.com/sublee/trueskill).

[中文版](README.md)

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