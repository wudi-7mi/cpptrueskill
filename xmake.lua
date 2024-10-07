-- 设置项目名称和版本
set_project("TrueSkill")
set_version("1.0.0")

-- 设置最小xmake版本
set_xmakever("2.5.0")

-- 设置C++标准
set_languages("c++17")

-- 添加编译选项
add_cxflags("-Wall", "-Wextra", "-pedantic")

-- 设置构建模式
add_rules("mode.debug", "mode.release")

-- 添加 doctest 依赖
add_requires("doctest")

-- 定义目标
target("trueskill")
    set_kind("static")
    add_headerfiles("src/*.hpp")
    add_includedirs("src")

-- 定义示例程序
target("example")
    set_kind("binary")
    add_files("example.cpp")
    add_deps("trueskill")
    add_includedirs("src")

-- 定义测试目标
target("test_main")
    set_kind("binary")
    add_files("tests/test_main.cpp")
    add_deps("trueskill")
    add_packages("doctest")
    add_includedirs("src", "tests")