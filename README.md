# C++ TrueSkill Through Time

[English](README_EN.md) | 中文

一个 C++17 纯头文件库，同时提供经典 TrueSkill 和 TrueSkill Through Time（TTT）。TTT 实现对齐 [TrueSkillThroughTime.py](https://github.com/glandfried/TrueSkillThroughTime.py) 2.0，支持历史消息传播、时间动态、混合观测、在线更新和模型证据。

## 功能

- 经典 TrueSkill：1v1、多队伍、平局、玩家权重和比赛质量。
- TrueSkill Through Time：前向/后向传播、收敛迭代和完整学习曲线。
- Ordinal、Continuous 和 Discrete 三种观测模型。
- 自定义玩家先验、`beta`、`gamma`、比赛时间和权重。
- 历史增量追加、在线估计和二进制状态恢复。
- C++17 header-only，无运行时第三方依赖。

## 选择 API

| 使用场景 | 主要类型 | 更新方式 |
| --- | --- | --- |
| 比赛结束后立即更新当前分数 | `Rating`、`TrueSkill` | 对每场比赛调用 `rate()` |
| 利用完整历史回推每个时刻的技能 | `Player`、`Game`、`History` | 构建历史后调用 `convergence()` |

两套 API 的结果方向不同，不能直接复用同一组名次值：

| API | 结果含义 | 两队示例 |
| --- | --- | --- |
| `TrueSkill::rate()` | rank 越小越好 | `{0, 1}` 表示第一队获胜 |
| `Game` / `History` | result 越大越好 | `{1, 0}` 表示第一队获胜 |

从 1.x 升级时，请先阅读 [迁移指南](MIGRATION.md)。

## 引入

将 `src` 加入 include path，然后包含统一入口头文件：

```cpp
#include "trueskill.hpp"
```

## 经典 TrueSkill

### 1v1

```cpp
#include "trueskill.hpp"

trueskill::Rating alice;
trueskill::Rating bob;

double quality = trueskill::quality_1vs1(alice, bob);
auto updated = trueskill::rate_1vs1(alice, bob);       // alice 获胜
auto drawn = trueskill::rate_1vs1(alice, bob, true);   // 平局

alice = updated[0][0];
bob = updated[1][0];
```

### 多队伍和权重

```cpp
trueskill::TrueSkill environment;

std::vector<std::vector<trueskill::Rating>> teams = {
    {trueskill::Rating(30, 7), trueskill::Rating(27, 6)},
    {trueskill::Rating(25, 8)},
    {trueskill::Rating(32, 5)}
};
std::vector<std::vector<double>> weights = {
    {1.0, 0.5},
    {1.0},
    {1.0}
};

// rank 越小越好：第二队第一，第一队第二，第三队第三。
auto updated = environment.rate(teams, {1, 0, 2}, weights);
double quality = environment.quality(teams, weights);
```

相同 rank 表示平局。经典接口中的零权重会被限制为内部最小值；如果需要让玩家完全不参与某场比赛，请使用 TTT 权重语义。

## TrueSkill Through Time

### 单局推断

`Player` 保存技能先验以及玩家自己的表现噪声 `beta` 和时间动态 `gamma`。`Game` 计算单局 likelihood、posterior 和 evidence。

```cpp
trueskill::Player alice(
    trueskill::math::Gaussian(25, 25.0 / 3),
    25.0 / 6,
    0.03);
trueskill::Player bob(
    trueskill::math::Gaussian(25, 25.0 / 3),
    25.0 / 6,
    0.03);

// result 越大越好：alice 获胜。
trueskill::Game game({{alice}, {bob}}, {1, 0});
auto posteriors = game.posteriors();
double evidence = game.evidence;
```

### 历史推断和学习曲线

`History` 使用字符串玩家 ID 描述比赛。每个 composition 元素是一场比赛，每场比赛包含多个队伍，每个队伍包含多个玩家。

```cpp
#include <iostream>
#include "trueskill.hpp"

trueskill::History::Composition composition = {
    {{"alice"}, {"bob"}},
    {{"bob"}, {"charlie"}},
    {{"charlie"}, {"alice"}}
};
trueskill::History::Results results = {
    {1, 0},
    {1, 0},
    {1, 0}
};
trueskill::History::Times times = {1, 2, 3};

trueskill::History history(composition, results, times);
auto convergence = history.convergence();

for (const auto& [name, curve] : history.learning_curves()) {
    for (const auto& [time, skill] : curve) {
        std::cout << name << " @ " << time
                  << ": mu=" << skill.mu()
                  << ", sigma=" << skill.sigma() << '\n';
    }
}

double log_evidence = history.log_evidence();
```

时间相同的比赛属于同一个 batch。跨 batch 传播时，模型会叠加标准差为 `gamma * sqrt(elapsed_time)` 的动态噪声。未传 `times` 时，比赛依次使用递增时间。

### 自定义先验和权重

```cpp
trueskill::History::Priors priors = {
    {"alice", trueskill::Player(
        trueskill::math::Gaussian(25, 8), 4.0, 0.05)}
};
trueskill::History::HistoryWeights weights = {
    {{1.0}, {0.0}}
};

trueskill::History history(
    {{{"alice"}, {"bob"}}},
    {{1, 0}},
    {1},
    priors,
    0, 3, 1, 0.15, 0, false,
    weights);
```

没有自定义先验的玩家使用 `History` 的全局参数。TTT 中权重为 `0` 表示该玩家对该局的队伍表现没有贡献。

## 观测模型

| 模型 | 用途 | result 示例 |
| --- | --- | --- |
| `Observation::Ordinal` | 胜负、排名和平局 | `{1, 0}`、`{0, 0}` |
| `Observation::Continuous` | 时间、距离等连续结果 | `{12.4, 11.8}` |
| `Observation::Discrete` | 进球、击杀等计数结果 | `{3, 1}` |

三种模型都遵循 result 越大表现越好的约定。耗时等原始值越小越好的指标，需要先取负数或做其他单调反向变换。

可为每场比赛指定模型；只传一个元素时，该模型应用于全部比赛：

```cpp
trueskill::History::Observations observations = {
    trueskill::Observation::Ordinal,
    trueskill::Observation::Continuous,
    trueskill::Observation::Discrete
};

trueskill::History history(
    composition, results, times, {},
    0, 3, 1, 0.15, 0, false, {}, observations);
```

## 在线和增量更新

`online=true` 会按时间 batch 依次揭示历史，保存进入每个 batch 前、仅依赖更早 batch 的在线先验快照，同时仍可读取包含当前比赛和前后消息的平滑 posterior。

```cpp
trueskill::History history(
    composition, results, times, {},
    0, 3, 1, 0.15, 0, true);

history.convergence();
auto smoothed = history.learning_curves();
auto online = history.learning_curves({}, true);

history.add_history(
    {{{"alice"}, {"charlie"}}},
    {{2, 0}},
    {4});
history.convergence();
```

在线模式只能追加晚于最后一个已揭示 batch 的时间。离线模式允许追加历史并重新平滑。

## 模型证据

- `Game::evidence`：单场比赛在当前先验下的证据。
- `History::log_evidence()`：历史中所有已处理比赛的 log evidence 之和。
- `History::geometric_mean()`：比赛 evidence 的几何平均值。

这些值适合比较参数或观测模型；它们不是经典 API 的 `quality()`，含义不能互换。

## 状态保存

```cpp
#include <fstream>

std::ofstream output("history.bin", std::ios::binary);
history.save(output);
output.close();

std::ifstream input("history.bin", std::ios::binary);
auto restored = trueskill::History::load(input);
```

保存内容包括前向/后向消息、likelihood、在线估计和 evidence。当前格式是带版本号的本机二进制格式，用于同一构建体系内恢复；不要将其作为跨平台或长期稳定的交换格式。

## 默认参数

| API | `mu` | `sigma` | `beta` | 动态参数 | 平局概率 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 经典 `TrueSkill` | 25 | 25/3 | 25/6 | `tau = 25/300` | 0.1 |
| TTT `Player` | 0 | 6 | 1 | `gamma = 0.03` | - |
| TTT `History` | 0 | 3 | 1 | `gamma = 0.15` | 0 |

## 构建和测试

```sh
xmake build example
xmake run example

xmake build test_main
xmake run test_main
```

更多完整用法见 [example.cpp](example.cpp)，算法说明见 [TS_TSTT_GUIDE.md](TS_TSTT_GUIDE.md)，从旧版接口升级见 [MIGRATION.md](MIGRATION.md)。
