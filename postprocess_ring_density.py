#!/usr/bin/env python3
"""
向后兼容包装脚本，用于调用 postprocessing 模块中的 postprocess_ring_density.py
"""

import sys
import os

# 添加当前目录到路径，确保可以导入postprocessing模块
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from postprocessing.postprocess_ring_density import main
except ImportError as e:
    print(f"错误: 无法导入postprocessing模块: {e}")
    print("请确保postprocessing目录存在且包含postprocess_ring_density.py")
    sys.exit(1)

if __name__ == '__main__':
    main()