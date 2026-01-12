#!/usr/bin/env python3
"""测试 examples 模块的功能。"""

import sys
from pathlib import Path

# 添加 src 目录到路径
sys.path.insert(0, str(Path(__file__).parent / "src"))

from cpptree.examples import list_examples, get_example, get_example_path

def test_list_examples():
    """测试 list_examples() 函数。"""
    print("=" * 80)
    print("测试 list_examples()")
    print("=" * 80)
    
    examples = list_examples()
    print(f"找到 {len(examples)} 个示例文件:")
    for i, name in enumerate(examples, 1):
        print(f"  {i}. {name}")
    
    return examples


def test_get_example():
    """测试 get_example() 函数。"""
    print("\n" + "=" * 80)
    print("测试 get_example()")
    print("=" * 80)
    
    # 测试获取第一个示例
    examples = list_examples()
    if examples:
        example_name = examples[0]
        print(f"\n获取示例: {example_name}")
        
        try:
            source_code, tree = get_example(example_name)
            print(f"✓ 成功读取文件")
            print(f"  源代码长度: {len(source_code)} 字符")
            print(f"  AST 根节点类型: {tree.root_node.type}")
            print(f"  AST 根节点子节点数: {len(tree.root_node.children)}")
            
            # 显示源代码的前几行
            print(f"\n源代码前 10 行:")
            print("-" * 80)
            lines = source_code.split('\n')[:10]
            for i, line in enumerate(lines, 1):
                print(f"{i:3}: {line}")
            if len(source_code.split('\n')) > 10:
                print("  ...")
            
        except Exception as e:
            print(f"✗ 错误: {e}")
            return False
    
    return True


def test_get_example_path():
    """测试 get_example_path() 函数。"""
    print("\n" + "=" * 80)
    print("测试 get_example_path()")
    print("=" * 80)
    
    examples = list_examples()
    if examples:
        example_name = examples[0]
        try:
            path = get_example_path(example_name)
            print(f"示例文件路径: {path}")
            print(f"文件存在: {path.exists()}")
            print(f"文件大小: {path.stat().st_size} 字节")
        except Exception as e:
            print(f"✗ 错误: {e}")
            return False
    
    return True


def test_get_example_without_extension():
    """测试不带扩展名的文件名。"""
    print("\n" + "=" * 80)
    print("测试 get_example() - 不带扩展名")
    print("=" * 80)
    
    try:
        source_code, tree = get_example("01_simple_defines")  # 不带 .f90
        print("✓ 成功：不带扩展名也能工作")
        print(f"  AST 根节点类型: {tree.root_node.type}")
    except Exception as e:
        print(f"✗ 错误: {e}")
        return False
    
    return True


if __name__ == "__main__":
    print("cpptree.examples 模块测试\n")
    
    try:
        test_list_examples()
        test_get_example()
        test_get_example_path()
        test_get_example_without_extension()
        
        print("\n" + "=" * 80)
        print("所有测试完成！")
        print("=" * 80)
        
    except Exception as e:
        print(f"\n✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

