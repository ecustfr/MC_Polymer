# 温度调度和Step External Potential配置说明

## 概述

本文档提供了关于MC_Polymer模拟系统中温度调度（Temperature Schedule）和阶梯势能（Step External Potential）的配置说明。系统已经移除了向后兼容的退火系统，现在只支持新的温度调度系统。

## 文件备份

本文件夹包含以下备份文件：
1. `run_simulation_backup.py` - 修改后的主模拟脚本备份
2. `batch_config_generator_backup.py` - 批量配置文件生成脚本备份
3. 示例JSON配置文件
4. 本说明文档

## 温度调度配置

### 基本概念

温度调度用于控制模拟过程中外势能项的温度（beta值）变化。beta = 1/(kT)，其中：
- 低beta值对应高温（高动能）
- 高beta值对应低温（低动能）

### 配置参数

温度调度配置位于JSON文件的`simulation_params`部分：

```json
"temperature_schedule": {
    "beta_values": [0.5, 1.0, 2.0, 1.0],
    "step_counts": [1000, 2000, 3000],
    "schedule_type": "linear",
    "schedule_param": 1.0
}
```

#### 参数说明

1. **beta_values** (必需)
   - 类型：浮点数数组
   - 说明：定义温度调度各阶段的beta值
   - 要求：数组长度必须比`step_counts`长度多1
   - 示例：`[0.5, 1.0, 2.0, 1.0]` 表示有4个温度点

2. **step_counts** (必需)
   - 类型：整数数组
   - 说明：定义每个温度段的模拟步数
   - 要求：数组长度必须比`beta_values`长度少1
   - 示例：`[1000, 2000, 3000]` 表示：
     - 第1段：从beta=0.5到beta=1.0，持续1000步
     - 第2段：从beta=1.0到beta=2.0，持续2000步
     - 第3段：从beta=2.0到beta=1.0，持续3000步

3. **schedule_type** (可选，默认："linear")
   - 类型：字符串
   - 可选值：
     - `"linear"`：线性温度变化
     - `"exponential"`：指数温度变化
   - 说明：定义温度在每段内的变化方式

4. **schedule_param** (可选，默认：1.0)
   - 类型：浮点数
   - 说明：调度参数，对于指数调度控制衰减速率
   - 对于指数调度：`decay = exp(-schedule_param * progress)`

### 调度类型详解

#### 线性调度 (linear)
- 公式：`current_beta = start_beta + (end_beta - start_beta) * progress`
- progress：在当前段内的进度（0.0到1.0）
- 温度均匀变化

#### 指数调度 (exponential)
- 公式：`decay = exp(-schedule_param * progress)`
- `current_beta = end_beta + (start_beta - end_beta) * decay`
- 温度变化速率逐渐减慢
- `schedule_param`越大，变化越快

### 使用示例

#### 示例1：简单升温再降温
```json
"temperature_schedule": {
    "beta_values": [0.5, 2.0, 0.5],
    "step_counts": [5000, 5000],
    "schedule_type": "linear"
}
```
- 从高温（beta=0.5）降温到低温（beta=2.0），5000步
- 再从低温升温回高温，5000步

#### 示例2：复杂温度循环
```json
"temperature_schedule": {
    "beta_values": [0.1, 0.5, 1.0, 2.0, 1.0, 0.5],
    "step_counts": [1000, 2000, 3000, 2000, 1000],
    "schedule_type": "exponential",
    "schedule_param": 2.0
}
```

## Step External Potential配置

### 基本概念

阶梯势能（Step Potential）是一种分段常数的外部势能，在z方向的不同区间内具有不同的势能值。

### 配置参数

Step Potential配置位于JSON文件的`Vext_params`部分，当`external_potential`为`"step"`时必需：

```json
"Vext_params": {
    "potential_type": "step",
    "boundaries": [0.0, 2.0, 4.0, 6.0],
    "potentials": [1.5, -0.8, 2.3],
    "C": 0.5
}
```

#### 参数说明

1. **potential_type** (必需)
   - 类型：字符串
   - 值：必须为`"step"`
   - 说明：标识势能类型

2. **boundaries** (必需)
   - 类型：浮点数数组
   - 说明：定义阶梯区间的边界点
   - 要求：
     - 必须包含0.0和盒子高度H
     - 必须严格递增
     - 长度必须比`potentials`长度多1
   - 示例：`[0.0, 2.0, 4.0, 6.0]` 表示3个区间：
     - [0.0, 2.0)
     - [2.0, 4.0)
     - [4.0, 6.0]

3. **potentials** (必需)
   - 类型：浮点数数组
   - 说明：每个区间的势能值
   - 要求：长度必须比`boundaries`长度少1
   - 示例：`[1.5, -0.8, 2.3]` 表示：
     - 区间[0.0, 2.0)：势能=1.5
     - 区间[2.0, 4.0)：势能=-0.8
     - 区间[4.0, 6.0)：势能=2.3

4. **C** (可选，默认：1.0)
   - 类型：浮点数
   - 说明：势能幅度缩放因子
   - 公式：最终势能 = C × 原始势能

### 物理意义

- 正势能值：排斥区域，聚合物倾向于避开
- 负势能值：吸引区域，聚合物倾向于聚集
- 势能为0：中性区域

### 使用示例

#### 示例1：简单三阶梯势能
```json
"Vext_params": {
    "potential_type": "step",
    "boundaries": [0.0, 3.0, 6.0],
    "potentials": [2.0, -1.0],
    "C": 0.5
}
```
- 区间[0.0, 3.0)：排斥势能 1.0 (2.0 × 0.5)
- 区间[3.0, 6.0)：吸引势能 -0.5 (-1.0 × 0.5)

#### 示例2：复杂多阶梯势能
```json
"Vext_params": {
    "potential_type": "step",
    "boundaries": [0.0, 1.0, 2.5, 4.0, 5.5, 8.0],
    "potentials": [3.0, -2.0, 1.5, -1.0, 0.5],
    "C": 0.3
}
```

## 批量配置文件生成

### 使用batch_config_generator.py

`batch_config_generator.py`脚本可以自动生成包含step potential和温度调度的配置文件。

#### 关键配置参数

在`batch_config_generator.py`中：

1. **设置external_potential为"step"**：
```python
MUST_INPUT = {"H":6.0, "polymer_type":"Linear", "knot_type":"Linear",
              "init_N":64, "external_potential":"step"}
```

2. **配置step potential参数**：
在`calculate_derived_params`函数中：
```python
if external_potential_type == "step":
    p["_Vext_payload"] = generate_vext_params(
        H,
        potential_type=external_potential_type,
        potential_mean=0.0,   # 平均高度
        potential_std=3.0,    # 高度波动
        C=0.5
    )
```

3. **配置温度调度**：
在`BASE_CONFIG`中：
```python
"temperature_schedule": {
    "beta_values": [0.5, 1.0, 1.0, 0.1, 1.0],
    "step_counts": [1000, 2000, 2000, 2000],
    "schedule_type": "linear"
}
```

## 运行模拟

### 基本命令
```bash
python run_simulation.py --config input/Linear_configs/config_0000_M6_Linear_H6.0_mu0.50.json
```

### 输出文件

模拟会生成以下文件：
1. `{external_potential}_potential_values.txt` - 势能值数据
2. `{external_potential}_potential_plot.png` - 势能分布图
3. `block_{n}_ensemble_data.dat` - 每个block的统计量
4. `block_{n}_rho_two_profile.dat` - 密度分布
5. `config.dat` - 模拟配置摘要

## 注意事项

1. **参数验证**：系统移除了防御性编程，参数错误会直接抛出异常
2. **向后兼容**：旧的`external_annealing_params`参数不再支持
3. **温度调度**：必须正确配置`beta_values`和`step_counts`的长度关系
4. **Step Potential**：`boundaries`必须包含0和H，且严格递增
5. **文件路径**：确保输入配置文件路径正确

## 故障排除

### 常见错误

1. **KeyError: 'temperature_schedule'**
   - 原因：JSON文件中缺少`temperature_schedule`字段
   - 解决：添加温度调度配置或使用默认beta=1.0

2. **ValueError: boundaries must be strictly increasing**
   - 原因：`boundaries`数组不是严格递增
   - 解决：检查并修正boundaries值

3. **ValueError: potentials length must be boundaries length - 1**
   - 原因：`potentials`和`boundaries`长度不匹配
   - 解决：确保`len(potentials) = len(boundaries) - 1`

4. **FileNotFoundError: Configuration file not found**
   - 原因：配置文件路径错误
   - 解决：检查文件路径和名称

## 示例文件

本文件夹包含以下示例配置文件：
1. `example_step_potential_temperature.json` - 包含step potential和温度调度
2. `example_temperature_only.json` - 只包含温度调度
3. `example_step_potential_only.json` - 只包含step potential

使用这些示例作为模板创建自己的配置文件。