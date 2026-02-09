#!/usr/bin/env python3
"""
优化版批量配置文件生成脚本
特性：Pathlib路径处理、自动字典递归更新、逻辑解耦
"""

import json
import numpy as np
from pathlib import Path
from itertools import product
from copy import deepcopy

# Import external potential functions
from external_potential import generate_vext_params

# from batch_config_generator_backup import INDEPENDENT_VARS

# ==============================
# 1. 核心配置与参数
# ==============================

CURRENT_SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = CURRENT_SCRIPT_DIR / "input" / "Ring_configs"



BASE_CONFIG = {
    "input_params": 
    {
        "polymer_type": "Linear", 
        "knot_type": "Trivial", 
        "configuration": "input/test_1.dat",
        "mu_b": 0.870821, 
        "M": 40, 
        "init_N": 64, 
        "rho_b": 0.4,
        "box_xy": 17.8886, 
        "H": 20.0, 
        "rcut": 1.0, 
        "max_N": 1280,
        "external_potential": "custom"
    },
    "simulation_params": 
    {
        "EPS_TRANS": 0.02, 
        "ROT_RATIO": 0.3, 
        "K_MAX": 10,
        "sample_interval": 5, 
        "sample_time": 10000, 
        "sample_block": 4, 
        "dz": 0.05
    },
    "Vext_params": {},
    "output_params": 
    { 
        "output_dir": "output/1.0", 
        "output_prefix": "mc_simulation" 
    }
}

# 耦合变量组 (H, rho_b, M, EPS_TRANS)

#必须置入的变量

MUST_INPUT = {"H":8.0,"polymer_type":"Ring","knot_type":"Trivial","init_N":64,"external_potential":"None"}
# MUST_INPUT = {"H":6.0,"knot_type":"Linear","init_N":64,"external_potential":"custom"} 

# 耦合变量
COUPLED_GROUPS = [{"EPS_TRANS":0.05,"K_MAX":10}]
for item in COUPLED_GROUPS:
     item.update(MUST_INPUT)


INDEPENDENT_VARS = {
    "M":[8],
    "rho_b":[0.1],
    "mu_b":[1.225130559]
}


"""
COUPLED_GROUPS = [
    {"H": 20.0, "rho_b": 0.1, "M": 40, "EPS_TRANS": 0.1},
    {"H": 20.0, "rho_b": 0.3, "M": 40, "EPS_TRANS": 0.05},
    {"H": 20.0, "rho_b": 0.5, "M": 40, "EPS_TRANS": 0.02}
]
INDEPENDENT_VARS = {
    "knot_type": ["Trefoil", "4_1Knot", "5_1Knot"],
}

"""

# 独立变量

# "external_potential": ["None", "custom"] # 可以把外势类型也加在这里

# ==============================
# 2. 工具函数 (核心逻辑)
# ==============================

def recursive_update(d, u):
    """递归更新字典，消除大量 if-else 判断"""
    for k, v in u.items():
        if isinstance(v, dict):
            d[k] = recursive_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def flatten_config(config):
    """将嵌套字典扁平化，用于查找参数在哪一层 (辅助函数)"""
    flat = {}
    for section, params in config.items():
        if isinstance(params, dict):
            for k in params:
                flat[k] = section
    return flat

# 预计算参数所在的层级结构，避免在循环中重复查找
PARAM_MAP = flatten_config(BASE_CONFIG)


def calculate_derived_params(params):
    """计算依赖参数 (物理公式)"""
    p = params.copy()
    
    # 提取基础变量
    N, M = p["init_N"], p["M"]
    rho, H = p["rho_b"], p["H"]

    knot = p["knot_type"]
    mu_b = p["mu_b"]

    # 物理计算
    p["box_xy"] = np.sqrt(N * M / (rho * H))
    p["max_N"] = int(p["box_xy"]**2 * H / M)
    
    # 路径更新
    config_fname = f"{knot}N{N}M{M}H{H:.2f}rho{rho:.2f}.dat"
    p["configuration"] = f"input/init_config_Z/{config_fname}"

    #p["output_prefix"] = f"mc_sim_{knot}_H{H:.1f}_rho{rho:.2f}_M{M}"
    p["output_prefix"] = f"mc_sim_{knot}_H{H:.1f}_mub{mu_b:.2f}_M{M}"
    
    # 外势处理
    if p.get("external_potential") == "custom":
        # 传入 seed 确保可复现，这里简单用 Hash 模拟
        # seed = hash(config_fname) % (2**32)
        # 将生成的 Vext 参数直接放入一个特殊的 key，稍后处理
        p["_Vext_payload"] = generate_vext_params(H)
        
    return p

# ==============================
# 3. 主流程
# ==============================

def main():
    print(f"=== 批量生成开始 | 目标: {OUTPUT_ROOT} ===")
    
    # 1. 生成所有参数组合 (笛卡尔积)
    # 将独立变量的 key 和 value 分离
    indep_keys = list(INDEPENDENT_VARS.keys())
    indep_vals = list(INDEPENDENT_VARS.values())
    
    
    count = 0
    # 遍历：耦合组 x 独立变量组合 x 外势类型
    for group_params in COUPLED_GROUPS:
        for indep_values in product(*indep_vals):
            # 合并独立变量
            current_params = group_params.copy()
            current_params.update(dict(zip(indep_keys, indep_values)))
            

                
                # --- 生成单个配置 ---
                # 1. 计算衍生参数
            full_params = calculate_derived_params(current_params)
                
                # 2. 准备配置对象 (深拷贝)
            config = deepcopy(BASE_CONFIG)
                
                # 3. 智能更新 (核心简化步骤)
                # 自动将 full_params 分发到 input/simulation/output 等部分
            for key, val in full_params.items():
                if key == "_Vext_payload":
                        config["Vext_params"] = val
                elif key in PARAM_MAP:
                        section = PARAM_MAP[key]
                        config[section][key] = val
                elif key == "output_dir":
                        config["output_params"]["output_dir"] = val
                
                # 4. 生成文件名和路径
                # 格式化文件名：去掉多余的0，保留必要精度
            fname_parts = [
                    f"config_{count:04d}",
                    f"M{full_params['M']:d}",
                    f"{full_params['knot_type']}",
                    f"H{full_params['H']:.1f}",
                    f"mu{full_params['mu_b']:.2f}"
            ]
            filename = "_".join(fname_parts) + ".json"
                
            # 确保输出目录存在
            # 动态生成子目录名
            sub_dir_name = f"OUT_{filename.replace('.json', '')}"
            config["output_params"]["output_dir"] = f"output/{sub_dir_name}" # 写入配置文件的路径
                
            file_path = OUTPUT_ROOT / filename
            file_path.parent.mkdir(parents=True, exist_ok=True)
                
            # 5. 写入文件
            with open(file_path, 'w') as f:
                json.dump(config, f, indent=4)
                
            print(f"[{count+1}] 生成: {file_path.name}")
            count += 1

    print(f"\n=== 完成: 共生成 {count} 个文件 ===")

if __name__ == "__main__":
    main()