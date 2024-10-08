# cpp TrueSkill

[En](README_EN.md) / 中文

这是一个 C++ 实现的 TrueSkill 纯头文件库，完全参考了一个 [python 实现](https://github.com/sublee/trueskill)。

## 使用示例

参考 [example.cpp](example.cpp) 文件。

### 1对1对战

```cpp
// 初始化
trueskill::Rating p1(80, 2.5);
trueskill::Rating p2(75, 2.5);

// 计算比赛质量
double quality1v1 = trueskill::quality_1vs1(p1, p2);
std::cout << "quality1v1: " << quality1v1 << std::endl;

// 计算比赛结果
auto result1v1 = trueskill::rate_1vs1(p1, p2);
for (const auto& team : result1v1) {
    for (const auto& rating : team) {
        std::cout << rating << " ";
    }
    std::cout << std::endl;
}

// 如果比赛是平局
auto result1v1_drawn = trueskill::rate_1vs1(p1, p2, true);
for (const auto& team : result1v1_drawn) {
    for (const auto& rating : team) {
        std::cout << rating << " ";
    }
    std::cout << std::endl;
}
```

### 多队伍对战

```cpp
// 初始化
trueskill::TrueSkill ts(25.0, 8.333333, 4.2, 0.3, 0.1);
std::cout << ts << std::endl;

trueskill::Rating alice(34.23, 7.22);
trueskill::Rating bob(24.23, 3.22);   
trueskill::Rating charlie(12.23, 3.9);
trueskill::Rating david(26.23, 6.11);
trueskill::Rating eve(41.2, 4.22);
trueskill::Rating frank(24.23, 5.08);
trueskill::Rating gloria(4.23, 3.08);

// 分配权重和队伍
std::vector<std::vector<double>> weights = {{1.0, 1.0}, {0.8, 1.0}, {1.0, 1.0, 1.0}};
std::vector<std::vector<trueskill::Rating>> teams = {{alice, bob}, {charlie, david}, {eve, frank, gloria}};

// 计算比赛质量
double quality = ts.quality(teams, weights);
std::cout << "quality: " << quality << std::endl;

// 计算比赛结果
auto newTeams = ts.rate(teams, {2, 1, 3}, weights);
for (const auto& team : newTeams) {
    for (const auto& rating : team) {
        std::cout << rating << " ";
    }
    std::cout << std::endl;
}
```

## 编译和运行

使用xmake进行编译和运行：

```shell
# 测试
xmake build test_main
xmake run test_main

# 示例
xmake build example
xmake run example
```