-- 设置项目名称和版本
set_project("TrueSkill")
set_version("1.0.0")

-- 设置最小xmake版本
set_xmakever("2.5.0")

-- 设置C++标准
set_languages("c++17")

-- 添加编译选项
add_cxflags("-Wall")

-- 设置构建模式
add_rules("mode.debug", "mode.release")

-- 添加 doctest 依赖
add_requires("doctest")

-- 定义为 header-only 库
target("trueskill")
    set_kind("headeronly")  -- 标记为 header-only 库
    add_headerfiles("src/*.hpp")  -- 添加所有头文件
    add_includedirs("src", {public = true})  -- 设置 include 路径，并且是 public 的，以便在其他项目中使用

-- 定义示例程序
target("example")
    set_kind("binary")
    add_files("example.cpp")
    add_deps("trueskill")  -- 依赖 trueskill
    add_includedirs("src")

-- 定义测试目标
target("test_main")
    set_kind("binary")
    add_files("tests/test_main.cpp")
    add_deps("trueskill")  -- 依赖 trueskill
    add_packages("doctest")  -- 使用 doctest 包
    add_includedirs("src", "tests")
