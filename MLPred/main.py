from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.svm import SVR
from sklearn.ensemble import AdaBoostRegressor
from sklearn.neighbors import KNeighborsRegressor
from keras.src.models import Sequential
from keras.src.layers import Dense, Input
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import networkx as nx
import numpy as np
import time
import sys

def read_graph(graph_file, node_set):
    g = nx.Graph()
    with open(graph_file, 'r') as file:
        for line in file:
            line = line.split()
            src_idx = int(line[0])
            dst_idx = int(line[1])
            if src_idx == dst_idx:
                continue
            node_set.add(src_idx)
            node_set.add(dst_idx)
            #g.add_edge(src_idx, dst_idx)

    g_size = max(node_set) + 1
    return g, g_size, node_set

def train_input2(graph, mode, seedset_num, k):
    train_seedsets = []
    train_inf = []
    k_str = str(k)
    if mode == "IM":
        train_seedsets = [f"../seedset/{graph}/seed_IM_{k}.txt" for k in [k, k, k, k, k]]
        train_inf = [f"../dataset/{graph}/inf_exp_IC_IM_s{k_str}_{i}.txt" for i in range(0, seedset_num, 1)]
    elif mode == "random":
        train_seedsets = [f"../seedset/{graph}/seed_random_{k}.txt" for k in [k, k, k, k, k]]
        train_inf = [f"../dataset/{graph}/inf_exp_IC_random_s{k_str}_{i}.txt" for i in range(0, seedset_num, 1)]

    #print(train_seedsets)
    #print(train_inf)
    return train_seedsets, train_inf

def train_input1(graph, mode, seedset_num, k):
    train_seedsets = []
    train_inf = []
    if graph in ["Email-Core", "Wiki-Vote", "Epinions", "Slashdot", "Email-Euall"]:
        if mode == "IM":
            train_seedsets = [f"../seedset/{graph}/seed_IM_{i}.txt" for i in range(5, 105, 5)]
            train_inf = [f"../dataset/{graph}/inf_exp_IC_IM_s{i}.txt" for i in range(5, 105, 5)]
        elif mode == "random":
            train_seedsets = [f"../seedset/{graph}/seed_random_{i}.txt" for i in range(10, 110, 5)]
            train_inf = [f"../dataset/{graph}/inf_exp_IC_random_s{i}.txt" for i in range(10, 110, 5)]
    elif graph in ["Pokec", "Livejournal", "Wiki-Link", "Twitter"]:
        if mode == "IM":
            train_seedsets = [f"../seedset/{graph}/seed_IM_{i}.txt" for i in range(50, 1050, 50)]
            train_inf = [f"../dataset/{graph}/inf_exp_IC_IM_s{i}.txt" for i in range(50, 1050, 50)]
        elif mode == "random":
            train_seedsets = [f"../seedset/{graph}/seed_random_{i}.txt" for i in range(100, 1100, 50)]
            train_inf = [f"../dataset/{graph}/inf_exp_IC_random_s{i}.txt" for i in range(100, 1100, 50)]

    return train_seedsets, train_inf

# 读取训练种子集文件
def read_seeds(train_seedsets, X):
    i = 0
    for file_name in train_seedsets:
        trainS = set()
        with open(file_name) as f:
            for line in f:
                seed_node = int(line.strip())
                trainS.add(seed_node)

        # 将当前种子集的编码结果存储到 X 中
        for value in trainS:
            X[i, value - 1] = 1

        i = i + 1

    return X

# 读取节点的影响概率文件
def read_inf(train_inf, y):
    sn = 0
    for file_name in train_inf:
        node_inf = {}
        with open(file_name) as f:
            for line in f:
                line = line.split()
                node_idx = int(line[0])
                inf_val = float(line[1])
                node_inf[node_idx] = inf_val

        # 将当前种子集的节点信息存储到数组 y 中
        for j, (node_idx, inf_val) in enumerate(node_inf.items()):
            y[sn, node_idx-1] = inf_val

        sn = sn + 1

    return y

# 读取待预测种子集文件
def read_test(graph, g_size, mode, k):
    S = set()
    k_str = str(k)
    pre_seedfile = "../seedset/" + graph + "/seed_" + mode + "_" + k_str+".txt"
    for line in open(pre_seedfile):
        line = line.split()
        seed_node = int(line[0])
        S.add(seed_node)

    # 创建一个全零的数组，大小为 1xg_size
    S_array = np.zeros((1, g_size), dtype=int)

    # 将集合中的元素位置设为 1
    for element in S:
        S_array[0, element - 1] = 1

    return S_array

# 写入预测结果
def write_result(graph, predictions, node_set, pre_result, _model, mode, k):
    # 打开文件，如果文件不存在则创建
    k_str=str(k)
    output = "result/"+graph+"/preprob_"+_model+"_"+mode+"_"+k_str+".txt"
    with open(output, "w") as f:
        for idx in node_set:
            prob = predictions[0][idx-1]
            f.write(f"{idx} {prob:.3f}\n")
            #print(f"{idx} {prob:.3f}\n")
            pre_result.append(prob)

# 模型评估
def evaluate(graph, pre_result, mode, k, train_time, predict_time):
    # 计算MSE（均方误差（Mean Squared Error，MSE））
    k_str = str(k)
    real_result_file = "../dataset/" + graph + "/inf_exp_IC_"+ mode + "_s" + k_str + ".txt"
    real_result = []
    with open(real_result_file) as f:
        for line in f:
            line = line.split()
            node_idx = int(line[0])
            inf_val = float(line[1])
            real_result.append(inf_val)

    mse = mean_squared_error(real_result, pre_result)
    print("Mean Squared Error (MSE): ", mse)

    # 计算RMSE（均方根误差（Root Mean Squared Error，RMSE））
    rmse = np.sqrt(mse)
    print("Root Mean Squared Error (RMSE): ", rmse)

    # 计算MAE（平均绝对误差（Mean Absolute Error，MAE））
    mae = mean_absolute_error(real_result, pre_result)
    print("Mean Absolute Error (MAE): ", mae)

    # 计算R-Squared（R方（r-squared，决定系数））
    r2 = r2_score(real_result, pre_result)
    print("R-Squared: ", r2)
    print("=" * 80)
    print()

    #outwrite result
    output = "result/" + graph + "/time_preprob_" + _model + "_" + mode + "_" + k_str + ".txt"
    with open(output, "w") as f:
        f.write("Train time: {:.6f} sec".format(train_time) + "\n")
        f.write("Predict time: {:.6f} sec".format(predict_time) + "\n")
        f.write("Mean Squared Error (MSE): " + str(mse) + "\n")
        f.write("Root Mean Squared Error (RMSE): " + str(rmse) + "\n")
        f.write("Mean Absolute Error (MAE): " + str(mae) + "\n")
        f.write("R-Squared: "+ str(r2) + "\n")



def main(graph, seedset_num, _model, mode, k):
    print("read graph...")
    graph_name = "../dataset/"+graph+"/graph.txt"
    node_set = set()
    g, g_size, node_set = read_graph(graph_name, node_set)
    print("graph size:", g_size)

    train_seedsets, train_inf = train_input1(graph, mode, seedset_num, k)

    # 创建一个全零的二维数组
    X = np.zeros((seedset_num, g_size), dtype=int)
    y = np.zeros((seedset_num, g_size), dtype=float)  # 注意数据类型应该为 float

    X = read_seeds(train_seedsets, X)
    print(X)

    y = read_inf(train_inf, y)
    print(y)

    # 从原始数据集中随机选择部分样本
    # sample_indices = np.random.choice(len(X), size=2, replace=False)  # 选择10000个样本，不重复
    # X = X[sample_indices]
    # y = y[sample_indices]
    # # 使用 X_sampled 和 y_sampled 进行模型训练

    print("\n"+"training...")
    start_time_train = time.perf_counter()
    # 使用随机森林模型进行训练
    # (default：'max_depth': None, 'min_samples_leaf': 1, 'min_samples_split': 2, 'n_estimators': 100)
    if _model == 'RF':
        model = RandomForestRegressor(n_estimators=100, random_state=0)  # n_estimators is the number of trees in the forest
        model.fit(X, y)
    # 使用决策树模型进行训练
    elif _model == 'DT':
        model = DecisionTreeRegressor(random_state=0)
        model.fit(X, y)
    # 使用K近邻回归
    elif _model == 'KNN':
        model = KNeighborsRegressor(n_neighbors=2)
        model.fit(X, y)
    # elif _model == 'DNN':
    #     # 定义神经网络模型
    #     model = Sequential()
    #     model.add(Input(shape=(g_size,)))  # 添加输入层，指定输入形状为 g_size
    #     model.add(Dense(g_size, activation='relu'))  # 隐藏层，使用 ReLU 激活函数
    #     model.add(Dense(g_size, activation='linear'))  # 输出层，使用线性激活函数
    #
    #     # 编译模型
    #     model.compile(loss='mean_squared_error', optimizer='adam')
    #
    #     # 训练模型
    #     model.fit(X, y, epochs=100, verbose=0)

    end_time_train = time.perf_counter()
    train_time = end_time_train - start_time_train
    print("Train time: {:.6f} sec".format(train_time) + "\n")

    # 读取种子集文件
    print("input query seed vector...")
    S_array = read_test(graph, g_size, mode, k)
    print(S_array)

    # 预测输入种子集下节点被影响的概率
    print("\n"+"predicting...")
    start_time_predict = time.perf_counter()

    predictions = model.predict(S_array)

    end_time_predict = time.perf_counter()
    predict_time = end_time_predict - start_time_predict
    print("Predicted probabilities: ", predictions[0])

    # 输出预测结果
    #for prob in predictions[0]:
    #    print("{:.3f}".format(prob)+" ")
    print("Predict time: {:.6f} sec".format(predict_time) + "\n")

    # 写入预测结果(只写入随机森林模型的结果)
    pre_result = []
    write_result(graph, predictions, node_set, pre_result, _model, mode, k)
    # if _model == "RF":
    #     write_result(graph, predictions, node_set, pre_result, _model, mode, k)
    # else:
    #     for idx in node_set:
    #         prob = predictions[0][idx]
    #         pre_result.append(prob)

    # 写入训练/预测时间和评估指标

    # 模型评估
    evaluate(graph, pre_result, mode, k, train_time, predict_time)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(1)
    graph = sys.argv[1]
    # 训练集种子集的数量
    seedset_num = int(sys.argv[2])
    # 训练模型：RF(随机森林) DT(决策树) KNN(k近邻)
    _model = sys.argv[3]
    # 训练种子集：IM random
    mode = sys.argv[4]
    # 测试种子集个数
    k = sys.argv[5]
    main(graph, seedset_num, _model, mode, k)
