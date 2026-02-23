#!/usr/bin/env python3

import os
import sys
import glob
import json
import argparse
import subprocess
import time
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


# 假设 run_simulation.py 在同级目录
SCRIPT_DIR = Path(__file__).resolve().parent
SIMULATION_SCRIPT = SCRIPT_DIR / 'run_simulation.py'

def check_config_conflicts(config_files):
    """检测输出目录冲突"""
    output_map = {}
    conflicts = {}

    for cfg_path in config_files:
        try:
            with open(cfg_path, 'r') as f:
                data = json.load(f)
            # 兼容性处理：防止key不存在报错
            out_dir = data.get('output_params', {}).get('output_dir')
            if not out_dir:
                print(f"Warning: No output_dir in {cfg_path}, skipping conflict check.")
                continue
                
            # 标准化路径以准确比对
            abs_out_dir = str(Path(out_dir).resolve())

            if abs_out_dir in output_map:
                if abs_out_dir not in conflicts:
                    conflicts[abs_out_dir] = [output_map[abs_out_dir]]
                conflicts[abs_out_dir].append(cfg_path)
            else:
                output_map[abs_out_dir] = cfg_path
        except Exception as e:
            print(f"Error reading {cfg_path}: {e}")

    if conflicts:
        print(f"\n❌ 发现输出目录冲突 ({len(conflicts)} 处):")
        for d, files in conflicts.items():
            print(f"  目录: {d}")
            for f in files:
                print(f"    - {f}")
        return False
    return True

def run_single_task(config_path):
    """
    单个任务执行函数，将在线程/进程池中运行
    """
    config_path = Path(config_path)
    try:
        # 读取配置获取输出目录
        with open(config_path, 'r') as f:
            config = json.load(f)
        output_dir = Path(config['output_params']['output_dir'])
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 复制配置文件到输出目录
        src_config = config_path
        dest_config = output_dir / config_path.name
        shutil.copy2(src_config, dest_config)
        print(f"Copied config file to: {dest_config}")
        
        log_file = output_dir / 'simulation.log'
        
        start_time = time.time()
        
        # 使用 with open 确保文件流正确关闭
        with open(log_file, 'w') as f_log:
            # 启动子进程
            # 注意：这里 cmd 用列表形式更安全
            cmd = [sys.executable, str(SIMULATION_SCRIPT), '--config', str(config_path)]
            
            # 使用 subprocess.run 替代 Popen，它是同步阻塞的，
            # 但因为我们在 ThreadPoolExecutor 里运行，所以实际上是并行的
            result = subprocess.run(
                cmd,
                stdout=f_log,
                stderr=subprocess.STDOUT, # 将 stderr 合并到 stdout
                text=True,
                check=False # 不自动抛出异常，由我们检查 returncode
            )
            
        duration = time.time() - start_time
        return {
            'config': config_path.name,
            'status': 'SUCCESS' if result.returncode == 0 else 'FAILED',
            'code': result.returncode,
            'duration': duration,
            'output': str(output_dir)
        }

    except Exception as e:
        return {
            'config': config_path.name,
            'status': 'ERROR',
            'error': str(e),
            'code': -1
        }

def main():
    parser = argparse.ArgumentParser(description="Parallel Simulation Runner")
    parser.add_argument('--input-dir', '-i', type=str, default=str(SCRIPT_DIR/ 'input/Linear_configs'))# input/Ring_configs
    parser.add_argument('--max-processes', '-m', type=int, default=3)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"Error: Input directory {input_dir} not found.")
        sys.exit(1)

    config_files = sorted(list(input_dir.glob('*.json')))
    if not config_files:
        print("No config files found.")
        sys.exit(1)

    # 1. 预检查
    print(f"检查 {len(config_files)} 个配置文件的冲突...")
    if not check_config_conflicts(config_files):
        sys.exit(1)

    print(f"\n🚀 开始执行任务 (并发数: {args.max_processes})...")
    
    results = []
    # 使用 ThreadPoolExecutor 管理子进程调用
    # 为什么用 ThreadPool 而不是 ProcessPool？
    # 因为我们的任务主要是 IO 等待（等待 subprocess 完成），Python 的 GIL 不会影响 subprocess 的运行。
    # ThreadPool 开销比 ProcessPool 小。
    with ThreadPoolExecutor(max_workers=args.max_processes) as executor:
        # 提交所有任务
        future_to_config = {executor.submit(run_single_task, f): f for f in config_files}
        
        for i, future in enumerate(as_completed(future_to_config)):
            res = future.result()
            progress = f"[{i+1}/{len(config_files)}]"
            
            if res['status'] == 'SUCCESS':
                print(f"{progress} ✅ {res['config']} (耗时: {res['duration']:.1f}s)")
            elif res['status'] == 'FAILED':
                print(f"{progress} ❌ {res['config']} (Return Code: {res['code']}) -> 查看日志: {res['output']}")
            else:
                print(f"{progress} ⚠️ {res['config']} Exception: {res.get('error')}")
            
            results.append(res)

    # 统计结果
    success_count = sum(1 for r in results if r['status'] == 'SUCCESS')
    print(f"\n=== 🏁 执行完成 ===")
    print(f"总计: {len(config_files)} | 成功: {success_count} | 失败: {len(config_files) - success_count}")

if __name__ == "__main__":
    main()