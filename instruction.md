你是一个专业的 python 和 c++ 程序员。现在有一个 TrueSkill python 的实现，需要你用 c++ 实现一个相同功能的库。
现在你已经在这个项目中了，项目包含有：

- `py` 目录，里面是 python 的实现
- `src` 目录，里面放了 C++ 实现的一些文件，`mathematics.hpp` 中定义了高斯分布类、矩阵类和一些数学函数，`factors.hpp` 中定义了因子图的一些基本操作，`trueskill.hpp` 中定义了 TrueSkill 算法的一些基本操作。
- `example.cpp` 是一个使用这个库的例子

目前我们已经完成了 `src/mathematics.hpp` 和 `src/factors.hpp` 的实现，接下来需要你根据 `py/__init__.py` 中的实现，完善 `src/trueskill.hpp` 中的实现，注意涉及全局环境的那些都不需要。
