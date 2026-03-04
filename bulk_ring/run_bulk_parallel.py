#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
并行运行多个化学势的体相环状聚合物密度模拟。
基于 test_bulk_density_final.cpp 编译的可执行文件。

使用说明：
1. 修改参数：在 main() 函数开头修改硬编码参数
2. 运行脚本：conda run -n base --no-capture-output python run_bulk_parallel.py
3. 查看结果：输出文件在当前目录生成

注意：
- 配置文件路径在C++代码中硬编码，如需修改需要重新编译
- 默认使用Release版本可执行文件，确保已编译
- 根据系统资源调整并行进程数
"""

import os
import sys
import subprocess
import multiprocessing
import time
from typing import List, Tuple, Optional
import logging

# 设置日志格式
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def run_single_mu(executable: str, mu: float, config_file: Optional[str] = None) -> Tuple[bool, str]:
    """
    运行单个化学势的模拟。

    Args:
        executable: 可执行文件路径
        mu: 化学势值
        config_file: 可选，配置文件路径，如果为None则使用C++默认路径

    Returns:
        (成功标志, 输出消息)
    """
    try:
        # 构建命令：可执行文件 + 化学势（作为字符串）
        cmd = [executable, str(mu)]

        # 如果需要覆盖配置文件，可以通过环境变量传递？
        # 由于C++代码中配置文件是硬编码的，这里无法直接传递。
        # 如果需要修改配置文件，可以考虑修改C++程序或使用符号链接。
        # 目前先保持原样。

        logger.info(f"开始运行化学势 mu = {mu}")
        start_time = time.time()

        # 运行子进程，捕获输出
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='replace'
        )

        elapsed = time.time() - start_time

        if result.returncode == 0:
            logger.info(f"化学势 mu = {mu} 运行完成，耗时 {elapsed:.1f} 秒")
            # 检查输出文件是否生成
            output_file = f"bulk_density_results_mu_{mu:.3f}.txt"
            if os.path.exists(output_file):
                logger.debug(f"结果文件已生成: {output_file}")
            else:
                logger.warning(f"未找到结果文件: {output_file}")

            # 返回成功和最后几行输出
            output_lines = result.stdout.strip().split('\n')[-5:] if result.stdout else []
            return True, f"成功: {mu}\n输出摘要:\n" + "\n".join(output_lines)
        else:
            error_msg = f"化学势 mu = {mu} 运行失败，返回码: {result.returncode}"
            logger.error(error_msg)
            if result.stderr:
                logger.error(f"错误输出:\n{result.stderr[:500]}")
            return False, f"失败: {mu}\n错误:\n{result.stderr[:500] if result.stderr else '无错误信息'}"

    except Exception as e:
        error_msg = f"运行化学势 mu = {mu} 时发生异常: {e}"
        logger.exception(error_msg)
        return False, error_msg


def worker(args: Tuple[str, float, Optional[str]]) -> Tuple[float, bool, str]:
    """用于进程池的工作函数"""
    executable, mu, config_file = args
    success, message = run_single_mu(executable, mu, config_file)
    return mu, success, message


def main():
    # ==================== 硬编码参数设置 ====================
    # 使用前请根据需求修改以下参数：

    # 1. 化学势列表（直接在此处指定）
    #    示例：[-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5]
    #mu_values = [3.5, 4.0, 4.5]

    # 或者使用范围模式（注释掉上面的mu_values，取消下面的注释）：
    start, stop, step = -1,1.0, 0.5  # 起始值，结束值，步长
    mu_values = []
    current = start
    while current <= stop + 1e-10:
        mu_values.append(round(current, 6))
        current += step

    mu_values = [5.5,5.75]
    # 2. 可执行文件路径
    #    默认使用Release版本，如需Debug版本请修改
    executable_path = './Release/test_bulk_density.exe'  # 默认Release版本
    # executable_path = './Debug/test_bulk_density.exe'  # Debug版本

    # 3. 并行进程数
    #    默认使用CPU核心数，可根据系统资源调整
    jobs = multiprocessing.cpu_count()  # 默认使用CPU核心数
    jobs = 4  # 或指定固定值，如4个并行进程

    # 4. 配置文件路径（目前C++代码中硬编码，此参数暂未使用）
    #    C++程序中配置文件路径为硬编码，如需修改需要重新编译C++代码
    config_path = None

    # 5. 日志文件（可选）
    #    不保存日志文件到磁盘（仅控制台输出）
    log_file = None  # 不保存日志文件
    # log_file = 'bulk_parallel.log'  # 保存到文件

    # 6. Dry-run模式（只打印命令不执行）
    #    用于预览将要执行的命令，不实际运行模拟
    dry_run = False
    # dry_run = True  # 启用Dry-run模式，只打印命令
    # =====================================================

    # 配置日志文件
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
        logger.addHandler(file_handler)

    logger.info(f"化学势列表: {mu_values}")
    logger.info(f"共 {len(mu_values)} 个任务")

    # 解析可执行文件路径
    if not os.path.isabs(executable_path):
        # 相对于脚本所在目录
        script_dir = os.path.dirname(os.path.abspath(__file__))
        executable_path = os.path.join(script_dir, executable_path)

    if not os.path.exists(executable_path):
        logger.error(f"可执行文件不存在: {executable_path}")
        sys.exit(1)

    logger.info(f"使用可执行文件: {executable_path}")
    logger.info(f"并行进程数: {jobs}")

    # 配置文件路径（目前仅作记录）
    if config_path and not os.path.exists(config_path):
        logger.warning(f"配置文件不存在: {config_path}")
        config_path = None

    # 准备参数列表
    task_args = [(executable_path, mu, config_path) for mu in mu_values]

    if dry_run:
        logger.info("Dry-run 模式，将执行以下命令:")
        for executable, mu, _ in task_args:
            cmd = [executable, str(mu)]
            logger.info(f"mu={mu}: {' '.join(cmd)}")
        logger.info(f"总共 {len(task_args)} 个任务，最大并发数: {jobs}")
        return

    # 开始并行执行
    logger.info(f"开始并行执行，最大并发数: {jobs}")
    start_time = time.time()

    results = []
    success_count = 0
    failure_count = 0

    # 使用进程池
    with multiprocessing.Pool(processes=jobs) as pool:
        # 使用imap_unordered以便实时获取结果
        for i, (mu, success, message) in enumerate(pool.imap_unordered(worker, task_args), 1):
            results.append((mu, success, message))
            if success:
                success_count += 1
                logger.info(f"任务进度: {i}/{len(mu_values)} [mu={mu}] 成功")
            else:
                failure_count += 1
                logger.error(f"任务进度: {i}/{len(mu_values)} [mu={mu}] 失败")

    total_time = time.time() - start_time

    # 汇总结果
    logger.info("=" * 60)
    logger.info("任务完成汇总:")
    logger.info(f"总任务数: {len(mu_values)}")
    logger.info(f"成功: {success_count}")
    logger.info(f"失败: {failure_count}")
    logger.info(f"总耗时: {total_time:.1f} 秒")
    logger.info(f"平均每个任务: {total_time/len(mu_values):.1f} 秒")

    # 输出失败任务详情
    if failure_count > 0:
        logger.warning("失败任务详情:")
        for mu, success, message in results:
            if not success:
                logger.warning(f"化学势 mu={mu}: {message[:200]}...")

    # 生成结果摘要文件
    summary_file = "bulk_parallel_summary.txt"
    try:
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(f"# 体相环状聚合物并行模拟摘要\n")
            f.write(f"# 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# 总任务数: {len(mu_values)}\n")
            f.write(f"# 成功: {success_count}\n")
            f.write(f"# 失败: {failure_count}\n")
            f.write(f"# 总耗时: {total_time:.1f} 秒\n")
            f.write(f"# 平均每个任务: {total_time/len(mu_values):.1f} 秒\n\n")

            f.write("化学势, 状态, 备注\n")
            for mu, success, message in results:
                status = "成功" if success else "失败"
                # 提取消息的第一行
                first_line = message.split('\n')[0] if message else ""
                f.write(f"{mu:.6f}, {status}, {first_line}\n")

        logger.info(f"结果摘要已保存到: {summary_file}")
    except Exception as e:
        logger.error(f"保存摘要文件失败: {e}")

    # 如果有失败任务，返回非零退出码
    if failure_count > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()