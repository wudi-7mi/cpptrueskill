# TrueSkill 与 TrueSkill Through Time：从直觉到公式

本文介绍经典 TrueSkill（TS）和 TrueSkill Through Time（TSTT，也常简称 TTT）。全文分为两个相对独立的部分：

- 第一部分使用具体比赛和 C++ 代码解释“它解决什么问题、结果怎么看、API 怎么用”。阅读这一部分不要求概率论基础。
- 第二部分使用统一符号、概率模型和数学公式说明“为什么这样更新、消息如何传播、当前实现具体计算了什么”。

第一次接触 TrueSkill 时，可以只读第一部分。需要调参、验证实现或研究模型时，再继续阅读第二部分。

## 阅读路径

- 想快速接入：阅读第一部分的第 1、3、5、6、9 节。
- 想理解 TSTT 为什么能修正历史：阅读第一部分第 4 节和第二部分第 17 至 19 节。
- 想校验 TS 更新公式：阅读第二部分第 11 至 16 节。
- 想理解 Continuous、Discrete 和 evidence：阅读第二部分第 20 至 22 节。
- 想调参：先读第一部分第 3.5、6.3 节，再读第二部分第 23、24 节。

---

# 第一部分：用例子理解算法和用法

## 1. 从 Alice 和 Bob 的比赛说起

假设 Alice 和 Bob 都是新玩家。系统不可能立刻知道他们的真实水平，因此不会只给每人保存一个固定分数，而是保存两个数：

```text
Alice: mu = 25, sigma = 8.33
Bob:   mu = 25, sigma = 8.33
```

可以先把它们理解为：

- `mu`：目前认为玩家大概有多强。
- `sigma`：对这个判断有多没把握。

两人的 `mu` 相同，表示当前看起来水平相近；`sigma` 很大，表示系统其实还不太了解他们。

如果 Alice 赢了，系统通常会：

- 提高 Alice 的 `mu`。
- 降低 Bob 的 `mu`。
- 降低两人的 `sigma`，因为这场比赛提供了新信息。

但一次获胜不会证明 Alice 永远比 Bob 强。比赛包含临场状态、地图、队友和随机事件，所以算法只会进行有限度的调整。

这就是 TrueSkill 的基本想法：同时估计“水平”和“对水平估计的可信程度”。

## 2. `sigma` 和发挥稳定性不是一回事

这是理解 TrueSkill 最重要的区别之一。

假设有两位玩家：

```text
Carol: mu = 30, sigma = 1
David: mu = 30, sigma = 8
```

系统认为他们的平均水平都在 30 左右，但对 Carol 很有把握，对 David 还缺少信息。

玩家每场发挥的随机程度由另一个参数 `beta` 表示：

- `sigma` 大：系统不了解这个玩家。
- `beta` 大：即使了解这个玩家，他每场发挥仍可能有较大波动。

所以一个参赛很多但发挥不稳定的老玩家，可能是“小 `sigma`、大 `beta`”；一个只打过一场的新玩家，可能是“大 `sigma`、普通 `beta`”。

## 3. TS：每场结束后更新当前评分

经典 TS 适合实时天梯。每场比赛结束后，把当前评分和比赛名次交给 `rate()`，再保存返回的新评分。

### 3.1 1v1

```cpp
#include "trueskill.hpp"

trueskill::Rating alice;
trueskill::Rating bob;

// 当前阵容的匹配质量。
double quality = trueskill::quality_1vs1(alice, bob);

// Alice 获胜。
auto updated = trueskill::rate_1vs1(alice, bob);
alice = updated[0][0];
bob = updated[1][0];
```

使用默认参数时，这个例子的近似结果是：

```text
比赛前：Alice = (mu=25.000, sigma=8.333)
        Bob   = (mu=25.000, sigma=8.333)

比赛后：Alice = (mu=29.396, sigma=7.171)
        Bob   = (mu=20.604, sigma=7.171)
```

Alice 和 Bob 原本看起来一样强，所以 Alice 获胜提供了明显信息：双方均值向相反方向移动，同时 `sigma` 从约 8.33 降到 7.17。系统变得更有把握，但一次比赛后仍保留了相当大的不确定性。

平局时传入第三个参数：

```cpp
auto updated = trueskill::rate_1vs1(alice, bob, true);
```

### 3.2 多队比赛

TS 使用 `rank` 表示名次，并且数值越小越好：

```cpp
trueskill::TrueSkill environment;

std::vector<std::vector<trueskill::Rating>> teams = {
    {trueskill::Rating(30, 7), trueskill::Rating(27, 6)},
    {trueskill::Rating(25, 8)},
    {trueskill::Rating(32, 5)}
};

// 第二队第一、第一队第二、第三队第三。
auto updated = environment.rate(teams, {1, 0, 2});
```

相同 rank 表示平局：

```cpp
// 第一队和第二队并列第一，第三队第三。
auto updated = environment.rate(teams, {0, 0, 2});
```

返回值保持输入的队伍顺序和玩家顺序，所以通常可以直接按原位置写回业务对象。

### 3.3 玩家权重

假设第一队的第二名玩家只参与了一半时间，可以使用权重：

```cpp
std::vector<std::vector<double>> weights = {
    {1.0, 0.5},
    {1.0},
    {1.0}
};

auto updated = environment.rate(teams, {1, 0, 2}, weights);
double quality = environment.quality(teams, weights);
```

权重会同时影响队伍表现和结果对玩家评分的反馈。经典实现会把零权重限制为一个内部最小值，因此零权重玩家仍有极小贡献。

### 3.4 `quality()` 是什么

`quality()` 回答的是：“按当前评分，这组阵容产生接近结果的可能性有多高？”

它适合用于辅助匹配，但不能简单理解为某一队的胜率，也不是比赛结束后的模型证据。

两个玩家的胜率可以单独计算：

```cpp
double p_alice_wins = environment.calculateWinProbability(alice, bob);
```

### 3.5 TS 参数的直觉

```cpp
trueskill::TrueSkill environment(
    25.0,        // mu：新玩家评分中心
    25.0 / 3,    // sigma：新玩家初始不确定性
    25.0 / 6,    // beta：单场发挥噪声
    25.0 / 300,  // tau：每次更新允许的技能变化
    0.1);        // draw_probability：平局概率
```

直观上：

- 增大 `sigma`：新玩家需要更多比赛才能稳定。
- 增大 `beta`：单场结果更可能是偶然，评分变化通常更保守。
- 增大 `tau`：系统更愿意相信玩家技能会在相邻更新间发生变化。
- 增大 `draw_probability`：更宽的表现差会被当作平局。

## 4. 为什么还需要 TSTT

现在假设 Alice 有三场比赛：

```text
第 1 天：Alice 击败 Bob
第 30 天：Alice 击败 Carol
第 60 天：Alice 击败 David
```

逐场 TS 可以不断更新 Alice 的当前评分，但它存在两个限制。

第一，它只向前看。处理完第 1 天的比赛后，就形成了当时的评分；第 60 天的新信息不会回头修正第 1 天的判断。

第二，它不知道真实时间间隔。连续两局和相隔一年，在经典 TS 中都只是“又调用了一次 `rate()`”。

TSTT 会为 Alice 在不同时间建立不同的技能节点：

```text
Alice@day1 <-> Alice@day30 <-> Alice@day60
```

相邻节点之间允许技能随时间变化。推断时，信息既从过去向未来传播，也从未来向过去传播。因此 Alice 后来的持续获胜可以提高模型对她早期技能的判断。

TSTT 最终返回的不是一个当前分数，而是一条学习曲线：

```text
Alice: [(day1, Gaussian), (day30, Gaussian), (day60, Gaussian)]
```

## 5. 用 `Game` 理解 TSTT 单局

TSTT 的 `Player` 与经典 TS 的 `Rating` 不同。它包含技能先验、单场发挥噪声 `beta` 和时间动态参数 `gamma`：

```cpp
trueskill::Player alice(
    trueskill::math::Gaussian(25, 25.0 / 3),
    25.0 / 6,
    0.03);

trueskill::Player bob(
    trueskill::math::Gaussian(25, 25.0 / 3),
    25.0 / 6,
    0.03);
```

可以用 `Game` 推断一场比赛：

```cpp
// TSTT 的 result 越大越好，因此 {1, 0} 表示 Alice 获胜。
trueskill::Game game({{alice}, {bob}}, {1, 0});

auto likelihoods = game.likelihoods;
auto posteriors = game.posteriors();
double evidence = game.evidence;
```

`Game` 只处理单场比赛，不进行跨时间传播。完整历史应使用 `History`。

## 6. 用 `History` 推断完整时间线

### 6.1 描述比赛历史

`History` 使用字符串玩家 ID 描述比赛：

```cpp
#include <iostream>
#include "trueskill.hpp"

trueskill::History::Composition composition = {
    {{"alice"}, {"bob"}},
    {{"bob"}, {"charlie"}},
    {{"charlie"}, {"alice"}}
};

// 每个 result 对应一场比赛；数值越大，表现越好。
trueskill::History::Results results = {
    {1, 0},
    {1, 0},
    {1, 0}
};

trueskill::History::Times times = {1, 2, 3};

trueskill::History history(composition, results, times);
history.convergence();
```

`Composition` 的四层嵌套关系是：

```text
历史 -> 比赛 -> 队伍 -> 玩家 ID
```

例如：

```cpp
{{{"alice", "bob"}, {"charlie"}}}
```

表示历史中有一场比赛，第一队是 Alice 和 Bob，第二队是 Charlie。

### 6.2 读取学习曲线

```cpp
auto curves = history.learning_curves();

for (const auto& [name, curve] : curves) {
    for (const auto& [time, skill] : curve) {
        std::cout << name << " @ " << time
                  << ": mu=" << skill.mu()
                  << ", sigma=" << skill.sigma() << '\n';
    }
}
```

只读取部分玩家：

```cpp
auto selected = history.learning_curves({"alice", "bob"});
```

### 6.3 时间和 batch

时间相同的比赛会被放进同一个 batch：

```cpp
trueskill::History::Times times = {1, 1, 2};
trueskill::History history(composition, results, times);

std::size_t games = history.size();          // 3 场比赛
std::size_t batches = history.batch_count(); // 2 个时间点
```

同一个玩家在同一 batch 内只有一个技能状态。这个 batch 中的多场比赛会共同更新它。

跨时间时，模型根据 `gamma * sqrt(elapsed_time)` 加入动态噪声。时间单位必须稳定：如果 `gamma` 是按“天”校准的，就不要直接传 Unix 秒时间戳。

### 6.4 自定义玩家先验

```cpp
trueskill::History::Priors priors = {
    {"alice", trueskill::Player(
        trueskill::math::Gaussian(25, 8),
        4.0,
        0.05)}
};

trueskill::History history(
    composition,
    results,
    times,
    priors,
    0,     // 其他玩家的 mu
    3,     // 其他玩家的 sigma
    1,     // 其他玩家的 beta
    0.15,  // 其他玩家的 gamma
    0);    // draw_probability
```

出现在 `priors` 中的玩家使用自己的参数，其他玩家使用 `History` 的全局参数。

## 7. TSTT 的三种比赛观测

### 7.1 Ordinal：只关心胜负和排名

```cpp
trueskill::Game game(
    teams, {1, 0}, 0.0, {},
    trueskill::Observation::Ordinal);
```

Ordinal 只关心大小关系，所以 `{1, 0}` 和 `{100, 0}` 表示相同结果。相同值表示平局。

适合：胜负、排名、晋级结果。

### 7.2 Continuous：使用实际差值

```cpp
trueskill::Game game(
    teams, {12.4, 11.8}, 0.0, {},
    trueskill::Observation::Continuous);
```

Continuous 会使用结果差 `12.4 - 11.8`，而不只是判断谁更大。

适合：连续得分、距离或经过变换的完成时间。

### 7.3 Discrete：使用离散计数差

```cpp
trueskill::Game game(
    teams, {3, 1}, 0.0, {},
    trueskill::Observation::Discrete);
```

Discrete 使用 Poisson 计数模型处理排序后相邻队伍的非负整数结果差。

适合：进球、击杀等计数。

三种模型都遵循“result 越大表现越好”。耗时等原始值越小越好的指标，需要先取负数或进行其他单调反向变换。

历史中可以混合观测模型：

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

只传一个 observation 时，它会应用于所有比赛；否则数量必须与比赛数量相同。

## 8. TSTT 权重、在线模式和增量追加

### 8.1 权重

```cpp
trueskill::History::HistoryWeights weights = {
    {{1.0, 0.5}, {1.0}}
};

trueskill::History history(
    {{{"alice", "substitute"}, {"bob"}}},
    {{1, 0}},
    {1},
    {}, 0, 3, 1, 0.15, 0, false,
    weights);
```

TSTT 中零权重表示玩家完全不贡献该场队伍表现，这一点与经典 TS 不同。

### 8.2 在线先验和平滑后验

```cpp
trueskill::History history(
    composition, results, times, {},
    0, 3, 1, 0.15, 0, true);

history.convergence();

auto smoothed = history.learning_curves();
auto online_priors = history.learning_curves({}, true);
```

当前实现中：

- `smoothed` 包含当前 batch 的比赛信息以及过去、未来的消息。
- `online_priors` 是进入当前 batch 前保存的先验快照，只包含更早 batch 的信息。

因此 `online_priors` 适合回答“比赛开始前系统如何评价玩家”，它不是包含当前比赛结果的 filtering posterior。

### 8.3 追加比赛

```cpp
history.add_history(
    {{{"alice"}, {"charlie"}}},
    {{2, 0}},
    {4});

history.convergence();
```

在线模式只能追加晚于最后一个已揭示 batch 的时间。离线模式可以插入其他时间，之后重新运行收敛和平滑。

## 9. evidence、保存状态和选择建议

### 9.1 evidence

evidence 衡量实际比赛结果在当前模型和先验下有多合理：

```cpp
double one_game = game.evidence;
double total_log = history.log_evidence();
double geometric_mean = history.geometric_mean();
```

它可以用于比较参数或观测模型，但必须在相同数据和相同结果单位下比较。它不是 TS 的 `quality()`。

### 9.2 保存完整历史状态

```cpp
#include <fstream>

std::ofstream output("history.bin", std::ios::binary);
history.save(output);
output.close();

std::ifstream input("history.bin", std::ios::binary);
auto restored = trueskill::History::load(input);
```

这是带版本号的本机二进制格式，适合同一构建体系内恢复，不应当作跨平台长期交换格式。

### 9.3 应该选 TS 还是 TSTT

选择 TS，当：

- 只需要比赛结束后更新当前评分。
- 主要目标是实时天梯或阵容匹配。
- 不需要利用未来比赛修正历史评分。
- 固定的 `tau` 已足够描述技能变化。

选择 TSTT，当：

- 需要每名玩家的完整学习曲线。
- 需要历史平滑，让未来比赛修正过去判断。
- 比赛间隔不均匀，需要显式建模时间。
- 需要 Continuous、Discrete 观测或 evidence。
- 需要比较赛前在线先验和事后平滑后验。

### 9.4 最容易犯的错误

两套 API 的结果方向相反：

```text
TS rate():         rank 越小越好，第一队胜使用 {0, 1}
TSTT Game/History: result 越大越好，第一队胜使用 {1, 0}
```

此外还应注意：

- `sigma` 是估计不确定性，不是单场发挥波动。
- TS 的零权重仍有极小贡献；TSTT 的零权重表示完全不贡献。
- `quality()` 和 evidence 含义不同。
- TSTT 会使用未来比赛回推历史，不应期待它与逐场 TS 数值一致。
- `Player` 和 `History` 的 TSTT 默认 prior 不同。

读到这里，已经可以正确选择和使用本仓库的两套 API。下面开始形式化描述算法。

---

# 第二部分：严谨的概率模型和数学推导

## 10. 记号约定

本文使用以下记号：

| 符号 | 含义 |
| --- | --- |
| `i` | 玩家索引 |
| `k` | 队伍索引 |
| `e` | 比赛事件索引 |
| `t` | 时间或时间 batch |
| `s_i`、`s_i(t)` | 玩家潜在技能 |
| `p_i`、`p_i(t)` | 玩家单场表现 |
| `T_k` | 队伍表现 |
| `d_k = T_k - T_(k+1)` | 排序后相邻队伍的表现差 |
| `mu_i` | 技能分布均值 |
| `sigma_i` | 技能分布标准差 |
| `beta_i` | 单场表现噪声标准差 |
| `tau` | TS 每次更新的动态标准差 |
| `gamma_i` | TSTT 单位时间动态参数 |
| `w_i` | 玩家在队伍中的权重 |
| `phi(x)` | 标准正态概率密度函数 |
| `Phi(x)` | 标准正态累积分布函数 |

下面的 `tau` 指 TS 动态参数。高斯自然参数中的精度均值将写作 `eta`，以免与代码中 `Gaussian::tau()` 的命名混淆。

## 11. 高斯分布和消息运算

### 11.1 均值和方差形式

一维高斯分布写作：

```text
N(x; mu, sigma^2)
  = 1 / sqrt(2*pi*sigma^2)
    * exp(-(x-mu)^2 / (2*sigma^2))
```

独立高斯随机变量的线性组合仍是高斯。若：

```text
X ~ N(mu_x, sigma_x^2)
Y ~ N(mu_y, sigma_y^2)
```

则：

```text
X + Y ~ N(mu_x + mu_y, sigma_x^2 + sigma_y^2)
X - Y ~ N(mu_x - mu_y, sigma_x^2 + sigma_y^2)
aX    ~ N(a*mu_x, a^2*sigma_x^2)
```

这三条规则解释了为什么个人表现可以聚合成队伍表现，以及为什么可以直接得到队伍表现差的分布。

### 11.2 自然参数形式

高斯分布也可以写成自然参数形式：

```text
pi  = 1 / sigma^2
eta = mu / sigma^2
```

忽略归一化常数后：

```text
N(x) proportional to exp(-0.5*pi*x^2 + eta*x)
```

两个高斯消息相乘时，自然参数直接相加：

```text
pi_product  = pi_1 + pi_2
eta_product = eta_1 + eta_2
```

消息相除时自然参数相减：

```text
pi_quotient  = pi_1 - pi_2
eta_quotient = eta_1 - eta_2
```

本仓库的 `math::Gaussian` 使用 `pi` 和名为 `tau` 的精度均值保存概率密度乘除。这里的 `Gaussian::tau()` 对应本文的 `eta`，与经典 TS 的动态参数 `TrueSkill::_tau` 不是同一个概念。

## 12. TS 的生成模型

### 12.1 技能动态

比赛更新前，TS 给已有技能加入固定动态噪声：

```text
s_i' | s_i ~ N(s_i, tau^2)
```

如果更新前：

```text
s_i ~ N(mu_i, sigma_i^2)
```

则动态传播后：

```text
s_i' ~ N(mu_i, sigma_i^2 + tau^2)
```

记传播后的标准差为：

```text
sigma_tilde_i = sqrt(sigma_i^2 + tau^2)
```

### 12.2 个人表现

给定技能后，个人单场表现为：

```text
p_i | s_i' ~ N(s_i', beta^2)
```

对技能先验积分后：

```text
p_i ~ N(mu_i, sigma_tilde_i^2 + beta^2)
```

### 12.3 队伍表现

队伍 `k` 的表现定义为：

```text
T_k = sum(i in team k, w_i * p_i)
```

假设各玩家表现条件独立，则：

```text
E[T_k] = sum(i in team k, w_i * mu_i)

Var[T_k]
  = sum(i in team k,
        w_i^2 * (sigma_tilde_i^2 + beta^2))
```

相邻两队表现差：

```text
d_k = T_k - T_(k+1)
```

仍为高斯：

```text
E[d_k]   = E[T_k] - E[T_(k+1)]
Var[d_k] = Var[T_k] + Var[T_(k+1)]
```

## 13. TS 的结果因子和截断高斯

### 13.1 获胜和平局约束

设平局边界为 `epsilon >= 0`。排序后前一队击败后一队时：

```text
d_k > epsilon
```

两队平局时：

```text
-epsilon <= d_k <= epsilon
```

`draw_probability` 会结合双方表现噪声换算成 `epsilon`。两队共有 `n` 名玩家且使用相同 `beta` 时，本仓库使用：

```text
epsilon
  = beta * sqrt(n)
    * Phi_inverse((draw_probability + 1) / 2)
```

### 13.2 单侧截断

若：

```text
d ~ N(mu_d, sigma_d^2)
```

并观测到 `d > epsilon`，定义：

```text
z = (mu_d - epsilon) / sigma_d
v(z) = phi(z) / Phi(z)
w(z) = v(z) * (v(z) + z)
```

则截断后的矩为：

```text
E[d | d > epsilon]
  = mu_d + sigma_d * v(z)

Var[d | d > epsilon]
  = sigma_d^2 * (1 - w(z))
```

`v` 控制均值移动，`w` 控制方差收缩。结果越出乎先验预料，`v` 通常越大，评分调整也越明显。

### 13.3 双侧截断和平局

平局条件是 `-epsilon <= d <= epsilon`。定义：

```text
a = (-epsilon - mu_d) / sigma_d
b = ( epsilon - mu_d) / sigma_d
Z = Phi(b) - Phi(a)
```

则：

```text
E[d | draw]
  = mu_d
    + sigma_d * (phi(a) - phi(b)) / Z

Var[d | draw]
  = sigma_d^2 * [
      1
      + (a*phi(a) - b*phi(b)) / Z
      - ((phi(a) - phi(b)) / Z)^2
    ]
```

实际实现通过等价的 `v_draw` 和 `w_draw` 计算这些矩。

## 14. TS 1v1 无平局更新公式

考虑 Alice 击败 Bob，权重均为 1，忽略平局边界。动态传播后的方差记作：

```text
v_A = sigma_A^2 + tau^2
v_B = sigma_B^2 + tau^2
```

表现差的总方差为：

```text
c^2 = v_A + v_B + 2*beta^2
```

标准化均值差：

```text
z = (mu_A - mu_B) / c
v = phi(z) / Phi(z)
w = v * (v + z)
```

更新后的均值为：

```text
mu_A' = mu_A + (v_A / c) * v
mu_B' = mu_B - (v_B / c) * v
```

更新后的方差为：

```text
sigma_A'^2 = v_A * (1 - (v_A / c^2) * w)
sigma_B'^2 = v_B * (1 - (v_B / c^2) * w)
```

这些公式展示了三个重要性质：

1. 胜者均值上升，负者均值下降。
2. 不确定性更大的玩家通常调整更多。
3. 爆冷时 `z` 较小甚至为负，`v` 更大，因此调整更明显。

多玩家、多队和平局不能仅靠这一组 1v1 公式处理。本仓库为技能、个人表现、队伍表现、相邻队伍差和结果约束建立因子图，并迭代传递高斯消息。

## 15. TS 因子图和消息传播

一场多队比赛的因子图可以概括为：

```text
skill priors
    |
    v
individual performances
    |
    v
weighted team performances
    |
    v
adjacent team differences
    |
    v
win/draw truncation factors
```

向下传播计算各中间变量的先验分布；结果因子进行截断矩匹配；再向上传播 likelihood，更新每个技能变量。

多队伍排序后，相邻差变量会共享中间队伍，因此彼此耦合。本仓库在差值层进行前向和反向迭代，直到消息变化不超过 `min_delta` 或达到内部迭代上限。

## 16. TS 的匹配质量

设共有 `n` 名玩家和 `K` 支队伍，所有玩家均值组成列向量 `mu`，技能协方差为对角矩阵 `Sigma`。令设计矩阵 `A` 的形状为 `n x (K-1)`，使 `A^T * x` 将玩家表现列向量 `x` 映射为相邻队伍表现差；玩家权重已包含在 `A` 中。令：

```text
B = beta^2 * A^T * A
C = A^T * Sigma * A
```

本实现的质量计算可写成：

```text
quality
  = sqrt(det(B) / det(B + C))
    * exp(-0.5 * mu^T * A * (B + C)^(-1) * A^T * mu)
```

第一项惩罚不确定性，第二项惩罚队伍间预期表现差。它衡量阵容产生接近比赛的程度，不依赖实际比赛结果。

## 17. TSTT 的状态空间模型

### 17.1 时间转移

TSTT 为玩家 `i` 在其参与的每个时间 batch 建立技能变量 `s_i(t)`。相邻参赛时间之间使用高斯随机游走：

```text
s_i(t) | s_i(t_prev)
  ~ N(s_i(t_prev), gamma_i^2 * delta_t)

delta_t = |t - t_prev|
```

因此，如果上一时刻的消息为：

```text
s_i(t_prev) ~ N(mu_prev, sigma_prev^2)
```

传播到当前时刻后：

```text
s_i(t)
  ~ N(mu_prev,
      sigma_prev^2 + gamma_i^2 * delta_t)
```

动态噪声标准差是：

```text
drift_i = gamma_i * sqrt(delta_t)
```

当前实现使用：

```text
drift_i
  = min(gamma_i * sqrt(delta_t),
        1.67 * History::sigma)
```

这个上限避免超长间隔让单次传播的不确定性无限增大。上限使用 `History` 的全局 `sigma`，即使该玩家有自定义 prior。

### 17.2 比赛模型

给定当前技能：

```text
p_i(t) | s_i(t) ~ N(s_i(t), beta_i^2)
```

队伍表现为：

```text
T_k(t) = sum(i in team k, w_i * p_i(t))
```

因此 TSTT 与 TS 使用相同的单局表现结构；区别在于技能 prior 来自时间链上的前向和后向消息，而不是一个孤立的当前 `Rating`。

## 18. TSTT 的 batch 后验

对于玩家 `i` 在时间 `t` 的技能节点，记：

- `f_i,t(s)`：从更早时间传播来的 forward 消息。
- `b_i,t(s)`：从更晚时间传播来的 backward 消息。
- `l_i,t,e(s)`：同一 batch 中比赛 `e` 对该玩家产生的 likelihood。

平滑后验满足：

```text
p(s_i(t) | all observations)
  proportional to
  f_i,t(s_i(t))
  * b_i,t(s_i(t))
  * product(e in events_i,t, l_i,t,e(s_i(t)))
```

由于这些消息都由高斯表示，乘法可以通过自然参数相加完成。

向未来发送的消息不能重复包含来自未来的 backward 信息，因此：

```text
forward_posterior
  proportional to
  forward * product(likelihoods)
```

向过去发送时同理：

```text
backward_posterior
  proportional to
  backward * product(likelihoods)
```

它们分别经过时间转移模型传播到相邻技能节点。

## 19. 前向过滤、后向平滑和在线快照

### 19.1 前向传播

按时间从早到晚处理 batch：

1. 将上一参赛时间的 `forward_posterior` 加入动态噪声。
2. 得到当前技能节点的 forward prior。
3. 用当前 batch 的比赛更新各 likelihood。
4. 将 `forward * likelihoods` 发送到下一个参赛时间。

### 19.2 后向传播

按时间从晚到早执行对称过程：

1. 将下一参赛时间的 `backward_posterior` 反向传播。
2. 更新当前技能节点的 backward 消息。
3. 重新计算与该节点相连的比赛 likelihood。
4. 将 `backward * likelihoods` 发送到上一个参赛时间。

### 19.3 当前实现的 online 值

`History` 的 `online` 字段在 batch 首次接收 forward 消息时保存：

```text
online_i,t = forward prior before games at time t
```

所以：

```cpp
history.learning_curves({}, true)
```

返回赛前先验快照，不包含当前 batch 的 likelihood。默认的 `learning_curves()` 返回完整平滑后验。

## 20. 三种 TSTT likelihood 的形式化定义

队伍按 result 从大到小排序。对排序后相邻队伍定义：

```text
d_k = T_k - T_(k+1)
r_k = result_k - result_(k+1)
```

### 20.1 Ordinal likelihood

非平局时：

```text
L(d_k) = I(d_k > epsilon_k)
```

平局时：

```text
L(d_k) = I(-epsilon_k <= d_k <= epsilon_k)
```

将该 likelihood 与高斯 prior 相乘会产生截断高斯。本实现用截断分布的一阶、二阶矩构造高斯近似。

如果 `d_k ~ N(mu_d, sigma_d^2)`，非平局 evidence 为：

```text
P(d_k > epsilon_k)
  = 1 - Phi((epsilon_k - mu_d) / sigma_d)
```

平局 evidence 为：

```text
P(-epsilon_k <= d_k <= epsilon_k)
  = Phi((epsilon_k - mu_d) / sigma_d)
    - Phi((-epsilon_k - mu_d) / sigma_d)
```

### 20.2 Continuous likelihood

当前实现把观测结果差当作精确的队伍表现差：

```text
L(d_k | r_k) = delta(d_k - r_k)
```

其中 `delta` 是 Dirac delta。对应代码中的退化高斯：

```text
Gaussian(r_k, 0)
```

若 prior 为 `N(mu_d, sigma_d^2)`，则 evidence 是该点的高斯密度：

```text
p(r_k)
  = N(r_k; mu_d, sigma_d^2)
```

因此 Continuous result 的单位会直接影响 evidence。当前 API 没有额外的连续观测噪声参数；需要测量噪声时，应在业务层调整结果尺度或扩展 likelihood。

### 20.3 Discrete likelihood

Discrete 模型把非负整数差 `r_k` 看作 Poisson 观测，并使用表现差作为 log rate：

```text
r_k | d_k ~ Poisson(exp(d_k))
```

其 likelihood 为：

```text
L(d_k | r_k)
  = exp(r_k*d_k - exp(d_k)) / r_k!
```

由于高斯 prior 与该 likelihood 不共轭，后验不再是精确高斯。本实现通过固定点迭代匹配近似后验的高斯参数。

设差值 prior 为 `N(mu_d, sigma_d^2)`，记 `q = sigma_d^2`、观测计数为 `r`。当前实现从 `kappa_0 = 1` 开始迭代：

```text
u_n = kappa_n - mu_d - r*q - 1

kappa_(n+1)
  = log(
      [mu_d + r*q - 1 - kappa_n
       + sqrt(u_n^2 + 2*q)]
      / (2*q)
    )
```

当相邻两次 `kappa` 的变化不超过容差，或达到 16 次迭代后，令：

```text
lambda = exp(kappa)
```

近似 posterior 的参数为：

```text
mu_post
  = mu_d + q * (r - lambda)

sigma_post^2
  = q / (1 + q*lambda)
```

evidence 为高斯与 Poisson likelihood 的积分：

```text
p(r_k)
  = integral[
      Poisson(r_k; exp(d_k))
      * N(d_k; mu_d, sigma_d^2)
    ] dd_k
```

当前实现使用有限区间上的 Simpson 数值积分近似该值，并将极小结果限制为 `1e-12`。

## 21. TSTT 的迭代与收敛

比赛 likelihood 依赖技能 prior，而技能 prior 又会因相邻时间和其他比赛消息改变，所以完整历史通常不能通过一次前向、后向传播达到固定点。

`History::convergence()` 重复迭代，追踪所有 likelihood 更新中的最大变化：

```text
step = (max absolute mu change,
        max absolute sigma change)
```

当：

```text
max(step.mu, step.sigma) <= epsilon
```

或达到最大迭代次数时停止。

```cpp
auto [step, iterations] = history.convergence(
    8,     // 最大迭代次数
    1e-5); // epsilon
```

多队比赛内部还会对相邻队伍差的消息进行固定次数的前后向更新。外层 `History` 收敛处理的是比赛之间和时间之间的耦合。

## 22. evidence 与 log evidence

对比赛 `e`，设其观测为 `y_e`，比赛前技能信息为 `D_<e`。单场 evidence 的一般定义是：

```text
Z_e = p(y_e | D_<e, model)
```

它通过对潜在表现变量积分得到：

```text
Z_e
  = integral[
      p(y_e | latent performances)
      * p(latent performances | priors)
    ] d(latent performances)
```

完整历史的 log evidence 为：

```text
log Z_history = sum(e, log Z_e)
```

使用 log 可以避免大量小概率相乘导致浮点下溢。几何平均 evidence 为：

```text
geometric_mean
  = exp(log Z_history / number_of_games)
```

比较不同参数时，应使用相同数据、时间单位和观测模型。Continuous evidence 是密度而非离散概率，结果单位变化会改变其数值。

## 23. 参数尺度与辨识

### 23.1 平移和缩放

Ordinal 模型主要由技能差决定。若同时对 `mu`、`sigma`、`beta` 和结果边界进行一致的尺度变换，相对排序可以保持不变。因此评分绝对数值通常是展示约定，而非可直接解释的物理量。

### 23.2 `beta` 与结果信息量

增大 `beta` 会增加队伍表现差的先验方差。同一个胜负结果相对于更大的表现噪声提供的信息更少，因此技能 posterior 通常变化更保守。

### 23.3 `gamma` 与时间单位

动态方差为：

```text
gamma^2 * delta_t
```

如果时间单位从天改为秒：

```text
delta_seconds = 86400 * delta_days
```

为了保持动态方差不变，应使用：

```text
gamma_per_second
  = gamma_per_day / sqrt(86400)
```

所以时间单位和 `gamma` 必须成对保存、成对校准。

### 23.4 `sigma` 与冷启动

初始 `sigma` 越大，新玩家 prior 越宽，早期比赛越容易推动均值大幅变化。随着 likelihood 累积，posterior 方差通常下降；时间动态又会在长期间隔后增加不确定性。

## 24. 本仓库默认参数

| API | `mu` | `sigma` | `beta` | 动态参数 | 平局概率 |
| --- | ---: | ---: | ---: | ---: | ---: |
| TS `TrueSkill` | 25 | 25/3 | 25/6 | `tau = 25/300` | 0.1 |
| TSTT `Player` | 0 | 6 | 1 | `gamma = 0.03` | - |
| TSTT `History` | 0 | 3 | 1 | `gamma = 0.15` | 0 |

`Player` 的默认 prior 与 `History` 自动创建玩家时使用的默认 prior 不同。构造历史时应明确以 `History` 参数或 `custom_priors` 为准。

## 25. 从模型到 API 的对应关系

| 数学对象 | C++ 类型或成员 |
| --- | --- |
| 技能高斯分布 | `trueskill::math::Gaussian` |
| TS 当前评分 | `trueskill::Rating` |
| TS 环境参数和更新 | `trueskill::TrueSkill` |
| TSTT 技能先验、`beta`、`gamma` | `trueskill::Player` |
| 单场 likelihood 和 posterior | `trueskill::Game` |
| 时间 batch 和消息传播 | `trueskill::History` |
| 单个时间节点的消息容器 | `trueskill::Skill` |
| 平滑后验 | `History::learning_curves()` |
| 赛前 forward 先验快照 | `History::learning_curves({}, true)` |
| 单场模型证据 | `Game::evidence` |
| 历史 log evidence | `History::log_evidence()` |

## 26. 相关文档

- [README.md](README.md)：API 速查、默认参数和构建方式。
- [MIGRATION.md](MIGRATION.md)：从 1.x 经典实现迁移到 2.x。
- [example.cpp](example.cpp)：可直接构建运行的示例。
