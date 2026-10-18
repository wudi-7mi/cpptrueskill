# 从 cpptrueskill 1.x 迁移到 2.x

2.x 在保留经典 TrueSkill 的基础上新增了 TrueSkill Through Time（TTT）。如果项目只需要比赛结束后的即时评分，可以继续使用经典接口；只有需要历史平滑、时间动态、学习曲线或混合观测时，才需要迁移到 `History`。

## 需要注意的兼容性变化

| 1.x / 经典用法 | 2.x 对应项 | 影响 |
| --- | --- | --- |
| `Player(id, Rating)` | `MatchPlayer(id, Rating)` | `Player` 名称现在属于 TTT 玩家模型 |
| `Rating` | `Rating` | 保持不变 |
| `TrueSkill::rate()` | `TrueSkill::rate()` | 保持不变，rank 越小越好 |
| 逐场重复调用 `rate()` | `History` | 可回推历史中每个时刻的技能 |
| `quality()` | `Game::evidence` / `History::log_evidence()` | 含义不同，不能直接替换 |
| 全局 `tau` | 玩家或历史的 `gamma` | TTT 按实际时间间隔传播不确定性 |

此外，底层 `math::Gaussian` 现在支持独立高斯的加减、标量乘法以及 `sigma == 0` 的退化高斯。这主要服务于 TTT 的消息传播；旧代码不应依赖此前对零方差输入的拒绝行为。

## 路径一：继续使用经典 TrueSkill

`Rating`、`TrueSkill::rate()`、`TrueSkill::quality()`、`rate_1vs1()` 和 `quality_1vs1()` 的调用方式不变。

```cpp
#include "trueskill.hpp"

trueskill::Rating alice;
trueskill::Rating bob;
trueskill::TrueSkill environment;

// 经典接口：rank 越小越好。
auto updated = environment.rate({{alice}, {bob}}, {0, 1});
alice = updated[0][0];
bob = updated[1][0];
```

如果旧代码使用了用于保存 ID 和评分的 `Player`，只需改名：

```cpp
#include "match.hpp"

// 1.x
// trueskill::Player player(42, trueskill::Rating());

// 2.x
trueskill::MatchPlayer player(42, trueskill::Rating());
```

`MatchPlayer` 没有通过 `trueskill.hpp` 导出，需要显式包含 `match.hpp`。新建代码通常也可以直接在业务对象中保存玩家 ID 和 `Rating`，不必依赖这个包装类型。

## 路径二：迁移到历史推断

### 1. 将评分转换为 TTT 玩家先验

经典 `Rating` 和 TTT `Player` 不是同一种类型。迁移已有评分时，显式复制均值和标准差，并选择 `beta` 与 `gamma`：

```cpp
trueskill::Rating rating(25, 25.0 / 3);

trueskill::Player player(
    trueskill::math::Gaussian(rating.mu(), rating.sigma()),
    25.0 / 6,  // beta：单局表现噪声
    0.03);     // gamma：单位时间的技能动态
```

不要无意中依赖 TTT 的默认值。经典 `Rating` 默认是 `mu=25, sigma=25/3`，而独立构造的 TTT `Player` 默认是 `mu=0, sigma=6, beta=1, gamma=0.03`。

### 2. 反转结果方向

这是迁移时最重要的行为差异：

```cpp
// 经典 API：第一队获胜，rank 越小越好。
auto updated = environment.rate(teams, {0, 1});

// TTT API：第一队获胜，result 越大越好。
trueskill::Game game(ttt_teams, {1, 0});
```

多队比赛同理。例如经典 rank `{1, 0, 2}` 可以转换成保持同一排序的 TTT result `{1, 2, 0}`。相同值在两套 API 中都表示平局。

### 3. 用玩家 ID 描述比赛历史

迁移前的逐场更新通常类似：

```cpp
for (const Match& match : matches) {
    auto updated = environment.rate(match.teams, match.ranks);
    // 将 updated 写回当前 Rating。
}
```

迁移后先构建完整历史，再运行消息传播：

```cpp
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
trueskill::History::Times times = {10, 20, 30};

trueskill::History history(composition, results, times);
history.convergence();

auto curves = history.learning_curves();
const auto& alice_curve = curves.at("alice");
const auto& [time, posterior] = alice_curve.back();
```

`learning_curves()` 返回每名玩家在其参赛 batch 上的 `(time, Gaussian)` 序列，而不只是一个当前评分。

### 4. 保留已有玩家参数

通过 `custom_priors` 为部分或全部玩家设置已有评分和个体参数：

```cpp
trueskill::History::Priors priors = {
    {"alice", trueskill::Player(
        trueskill::math::Gaussian(25, 25.0 / 3),
        25.0 / 6,
        0.03)},
    {"bob", trueskill::Player(
        trueskill::math::Gaussian(30, 6),
        25.0 / 6,
        0.03)}
};

trueskill::History history(
    composition, results, times, priors,
    25, 25.0 / 3, 25.0 / 6, 0.03);
```

构造函数中的全局参数用于没有出现在 `priors` 中的玩家。

## 时间语义

经典 `tau` 在每次 `rate()` 时注入固定动态噪声。TTT 的 `gamma` 使用时间差：每次传播会叠加标准差为 `gamma * sqrt(elapsed_time)` 的动态噪声。

- 应使用单调、含义稳定的时间单位；更换单位时同步调整 `gamma`。
- 时间相同的比赛属于同一个 batch，会联合处理。
- 不传 `times` 时，每场比赛自动使用递增时间。
- 在线模式追加的时间必须晚于最后一个已处理 batch。

因此，把数据库时间戳直接换成秒通常会显著放大动态噪声。可以改用天、周等业务单位，或相应缩放 `gamma`。

## 权重语义

两套接口对零权重的处理不同：

- 经典 `TrueSkill` 会把零权重限制为内部最小值，玩家仍有极小贡献。
- TTT 的零权重表示玩家对该场队伍表现完全没有贡献。

非零权重仍表示玩家表现对队伍表现的相对贡献。迁移后应特别复核替补、掉线和未出场玩家的处理。

```cpp
trueskill::History::HistoryWeights weights = {
    {{1.0, 0.0}, {1.0}}
};

trueskill::History history(
    {{{"alice", "substitute"}, {"bob"}}},
    {{1, 0}},
    {1},
    {}, 0, 3, 1, 0.15, 0, false,
    weights);
```

## 在线处理和追加历史

需要按比赛到达顺序保留实时估计时，启用在线模式：

```cpp
trueskill::History history(
    composition, results, times, {},
    0, 3, 1, 0.15, 0, true);

history.convergence();

auto online_curves = history.learning_curves({}, true);
auto smoothed_curves = history.learning_curves();

history.add_history(
    {{{"alice"}, {"charlie"}}},
    {{1, 0}},
    {40});
history.convergence();
```

`learning_curves({}, true)` 读取进入每个 batch 前的在线先验快照，只使用更早 batch 的数据；默认调用读取包含当前比赛和后续信息的平滑 posterior。

## 迁移观测模型

经典 TrueSkill 只表达顺序结果。TTT 可按比赛选择观测模型：

| 模型 | 适用数据 | 结果解释 |
| --- | --- | --- |
| `Ordinal` | 胜负、排名、平局 | 只使用顺序关系 |
| `Continuous` | 时间、距离、连续分差 | 使用连续结果值 |
| `Discrete` | 比分、进球、击杀 | 使用 Poisson 计数近似 |

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

如果所有比赛使用同一种模型，可只传一个 observation。无论哪种模型，`result` 都是值越大表现越好；耗时等原始值越小越好的指标需要先取负数或做其他单调反向变换。

## `quality` 与 evidence

经典 `quality()` 衡量双方先验下产生接近结果的概率，常用于匹配质量。TTT 的 evidence 衡量实际观测结果在模型下的概率，常用于参数和模型比较。

迁移时不要将以下代码机械替换：

```cpp
double match_quality = environment.quality(teams);
double game_evidence = game.evidence;
double history_log_evidence = history.log_evidence();
```

它们解决的问题不同，数值范围和业务阈值也不能沿用。

## 状态持久化

`History::save()` 和 `History::load()` 保存完整推断状态，恢复后可以继续读取曲线或追加比赛：

```cpp
#include <fstream>

std::ofstream output("history.bin", std::ios::binary);
history.save(output);
output.close();

std::ifstream input("history.bin", std::ios::binary);
auto restored = trueskill::History::load(input);
```

这是带版本号的本机二进制格式，适合同一构建体系内恢复。它不是跨架构、跨编译器或长期存档的稳定交换格式；长期数据仍应保存原始比赛和模型参数，并在需要时重建 `History`。

## 验证清单

- 所有经典 `Player(id, rating)` 已改为 `MatchPlayer`，或替换为业务自己的玩家类型。
- 所有传给 `Game` / `History` 的结果都已改成“值越大越好”。
- 已显式确认 TTT 的 `mu`、`sigma`、`beta`、`gamma` 和平局概率，而非误用默认值。
- 时间单位与 `gamma` 匹配，相同时间确实应该作为同一 batch。
- 已复核零权重玩家的语义。
- 业务没有把 `quality()` 的阈值用于 evidence。
- 在线与平滑学习曲线读取方式符合产品需求。
- 使用固定数据集比较了迁移前后的排序、均值和不确定性；TTT 会利用未来比赛回推历史，因此数值不应期待逐场完全一致。
