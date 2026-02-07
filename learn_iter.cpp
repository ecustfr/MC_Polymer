# include <iostream>
# include <vector>
void grow(int current,int parent ,int direction ,std::vector<int>& order, std::vector<int>& connections,int M);


void grow(int current,int parent ,int direction ,std::vector<int>& order, std::vector<int>& connections,int M)
{
    // 边界检查：如果超出索引范围，停止递归
    if (current < 0 || current >= M) {
        return;
    }

    // 1. 记录当前单体和它的连接父亲
    order.push_back(current);
    connections.push_back(parent);
    for(const auto &monomer:order)
    {
        std::cout << monomer<<" " ;
    }
    std::cout<<"|" ;
    // 2. 向指定方向（-1 或 1）递归生长下一个单体
    // 此时当前的 current 变成了下一个单体的 parent
    grow(current + direction, current, direction, order, connections,M);
}





/**
 * @brief 递归生长分支
 * @param next_idx 当前要插入的单体索引
 * @param parent_idx 作为参照的父单体索引
 * @param step 方向步长（-1 表示向左，1 表示向右）
 */
bool MuVT_MC_LinearPolymer::grow_branch_recursive(
    int next_idx, int parent_idx, int step, 
    double &total_W, 
    std::vector<std::array<double, 3>> &r_new, 
    std::vector<bool> &is_inserted,
    int k_max) 
{
    // 基准情况：超出链边界，说明该方向生长完成
    
}

/**
 * @brief 主入口：从中间点插入并双向生长
 */


bool MuVT_MC_LinearPolymer::grow_branch_recursive(
    int next_idx, int parent_idx, int step, 
    double &total_W, 
    std::vector<std::array<double, 3>> &r_new, 
    std::vector<int> & insert_over_monomer
    int k_max) 
{
    // 基准情况：超出链边界，说明该方向生长完成
    if (next_idx < 0 || next_idx >= this->M) {
        return true;
    }

    std::vector<std::array<double, 3>> candidates(k_max);
    std::vector<double> w_trial(k_max, 0.0);
    double sum_w = 0.0;

    // 1. 尝试 k_max 个候选位置
    for (int k = 0; k < k_max; ++k) {
        candidates[k] = r_new[parent_idx]; // 基于父节点位置
        add_random_unit_vector(candidates[k].data(), this->gen);

        // 2. 碰撞检测：检查其他链 + 检查本链已插入的部分
        if (!this->overlap_all_other_polymer(candidates[k].data()) && 
            !this->overlap_with_inserted(candidates[k], r_new, is_inserted)) {
            
            w_trial[k] = this->Vext_bz(candidates[k].data());
            sum_w += w_trial[k];
        }
    }

    if (sum_w <= 1e-12) return false;

    // 3. 选择并记录位置
    int selected = RouteWheel(w_trial.data(), k_max, this->gen);
    r_new[next_idx] = candidates[selected];
    is_inserted[next_idx] = true; // 关键：标记为已插入
    
    // 4. 更新权重
    total_W *= (sum_w / static_cast<double>(k_max));

    // 5. 继续向同一方向递归
    return grow_branch_recursive(next_idx + step, next_idx, step, total_W, r_new, is_inserted, k_max);
}


int main()
{
Get_G_insert(4 , 10);
}