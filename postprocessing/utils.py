"""
后处理共享工具函数

包含多个后处理脚本共用的函数，避免代码重复。
"""

import os
import json
import re
import numpy as np
from pathlib import Path
try:
    from scipy.interpolate import CubicSpline
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def natural_sort_key(s):
    """
    自然排序键函数，用于按数字顺序排序文件。

    参数:
        s: 字符串

    返回:
        list: 用于排序的键列表
    """
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', s)]


def find_block_files(output_dir, pattern='block_*_rho_profile.dat'):
    """
    在输出目录中查找block文件并按自然顺序排序。

    参数:
        output_dir: 输出目录路径
        pattern: 文件匹配模式

    返回:
        list: 排序后的文件路径列表
    """
    import glob
    files = glob.glob(os.path.join(output_dir, pattern))
    files.sort(key=natural_sort_key)
    return files


def load_total_density_profile(file_path):
    """
    加载密度分布文件。

    参数:
        file_path: 密度文件路径

    返回:
        np.ndarray: 密度分布数组
    """
    return np.loadtxt(file_path)


def extract_config_name(filename):
    """
    从文件名提取配置名。
    支持两种格式：
    1. 新格式: avg_density_config_XXXX.npz
    2. 旧格式: config_XXXX_avg_density.npz

    参数:
        filename: 文件名或完整路径

    返回:
        str: 配置名
    """
    base_name = os.path.basename(filename)

    # 尝试匹配新格式: avg_density_config_XXXX.npz
    if base_name.startswith('avg_density_') and base_name.endswith('.npz'):
        config_name = base_name[len('avg_density_'):-len('.npz')]
    # 尝试匹配旧格式: config_XXXX_avg_density.npz
    elif base_name.endswith('_avg_density.npz'):
        config_name = base_name[:-len('_avg_density.npz')]
    else:
        # 否则使用不带扩展名的文件名
        config_name = os.path.splitext(base_name)[0]

    return config_name


def resolve_output_dir(json_config, json_path):
    """
    解析JSON配置中的output_dir为绝对路径。

    参数:
        json_config: JSON配置字典
        json_path: JSON文件路径

    返回:
        str: 绝对路径的输出目录
    """
    output_dir = json_config['output_params']['output_dir']

    # 确保output_dir是绝对路径（相对于项目根目录）
    if not os.path.isabs(output_dir):
        # 假设JSON文件路径相对于项目根目录，计算绝对路径
        # 项目根目录是JSON文件父目录的父目录的父目录
        project_root = Path(json_path).parent.parent.parent
        output_dir = os.path.join(project_root, output_dir)

    return output_dir


def get_npz_path_from_json(json_path):
    """
    根据JSON文件路径获取对应的NPZ文件路径。

    参数:
        json_path: JSON配置文件路径

    返回:
        str: NPZ文件路径
    """
    # 读取JSON文件获取output_dir
    with open(json_path, 'r') as f:
        config = json.load(f)

    # 获取output_dir
    output_dir = resolve_output_dir(config, json_path)

    # 从JSON文件名提取配置名
    config_name = os.path.splitext(os.path.basename(json_path))[0]

    # 构建两种可能的npz文件路径（新格式和旧格式）
    npz_file_new = os.path.join(output_dir, f"avg_density_{config_name}.npz")
    npz_file_old = os.path.join(output_dir, f"{config_name}_avg_density.npz")

    # 检查文件是否存在（优先新格式）
    if os.path.exists(npz_file_new):
        return npz_file_new
    elif os.path.exists(npz_file_old):
        return npz_file_old
    else:
        # 如果都不存在，返回新格式路径（用于创建新文件）
        return npz_file_new


def ensure_directory_exists(dir_path):
    """
    确保目录存在，如果不存在则创建。

    参数:
        dir_path: 目录路径
    """
    os.makedirs(dir_path, exist_ok=True)


def load_json_file(json_path):
    """
    加载JSON配置文件。

    参数:
        json_path: JSON文件路径

    返回:
        dict: JSON配置字典
    """
    with open(json_path, 'r') as f:
        return json.load(f)


def save_npz_data(file_path, data_dict):
    """
    保存数据到NPZ文件。

    参数:
        file_path: 输出文件路径
        data_dict: 要保存的数据字典
    """
    ensure_directory_exists(os.path.dirname(file_path))
    np.savez(file_path, **data_dict)


def load_npz_data(file_path):
    """
    从NPZ文件加载数据。

    参数:
        file_path: NPZ文件路径

    返回:
        dict: 加载的数据字典
    """
    return np.load(file_path)


def load_mu_ex_interpolator(csv_path=None):
    """
    加载mu_b到mu_ex的映射数据并创建三次样条插值器。

    参数:
        csv_path: CSV文件路径，默认为"input/Ring_mu_ex.csv"

    返回:
        callable: 插值函数，输入mu_b返回mu_ex_b，如果CSV文件不存在或scipy不可用则返回None
    """
    if csv_path is None:
        csv_path = os.path.join("input", "Ring_mu_ex.csv")

    if not os.path.exists(csv_path):
        print(f"警告: CSV文件不存在: {csv_path}")
        return None

    try:
        # 加载CSV数据
        data = np.loadtxt(csv_path, delimiter=',', skiprows=1)
        mu_b_values = data[:, 0]
        mu_ex_values = data[:, 1]

        # 检查数据有效性
        if len(mu_b_values) < 2:
            print(f"警告: CSV文件中数据点不足: {csv_path}")
            return None

        # 使用三次样条插值（如果scipy可用）
        if SCIPY_AVAILABLE:
            interpolator = CubicSpline(mu_b_values, mu_ex_values)
            # print(f"已创建三次样条插值器，数据点: {len(mu_b_values)}")
            return interpolator
        else:
            print("警告: scipy不可用，使用线性插值代替")
            # 使用numpy的线性插值作为备选
            from numpy import interp
            # 返回一个包装函数
            def linear_interp(mu_b):
                return np.interp(mu_b, mu_b_values, mu_ex_values)
            return linear_interp
    except Exception as e:
        print(f"加载CSV文件或创建插值器时出错: {e}")
        import traceback
        traceback.print_exc()
        return None


def calculate_mu_ex_b(mu_b, interpolator):
    """
    使用插值器计算mu_ex_b。

    参数:
        mu_b: 体相化学势
        interpolator: 插值函数

    返回:
        float: 插值得到的mu_ex_b，如果插值器为None则返回None
    """
    if interpolator is None:
        return None

    try:
        return float(interpolator(mu_b))
    except Exception as e:
        print(f"插值计算mu_ex_b时出错 (mu_b={mu_b}): {e}")
        return None