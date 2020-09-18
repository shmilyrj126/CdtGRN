import os
import sys
import time

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
import tensorflow as tf
from scipy.stats import reciprocal
from sklearn.model_selection import RandomizedSearchCV
#利用sklearn内置函数进行输入层归一化
from sklearn.preprocessing import StandardScaler
from tensorflow import keras

#选定GPU为训练计算器
physical_devices = tf.config.experimental.list_physical_devices('GPU')
assert len(physical_devices) > 0, "Not enough GPU hardware devices available"
tf.config.experimental.set_memory_growth(physical_devices[0], True)

#从csv文件导入数据集
input_dir="datasets"
df_data = pd.read_csv(os.path.join(input_dir, "train"+".csv"))
df_label = pd.read_csv(os.path.join(input_dir, "valid"+".csv"), header=None)
x_train_all, y_train_all = df_data, df_label

#调整训练集，首先将数据集转换为(2行 7 列 的向量，便于CNN识别规律)，接着将验证集与训练集分开，前5000为验证集数据，后面为验证集数据（预处理时已经将阴阳集混匀打散，这里没有进行打散）
x_valid, x_train = x_train_all.to_numpy().reshape(-1,2,7)[:5000], x_train_all.to_numpy().reshape(-1,2,7)[5000:]
y_valid, y_train = y_train_all.to_numpy()[:5000], y_train_all.to_numpy()[5000:]
x_test, y_test = x_train_all.to_numpy().reshape(-1,2,7),y_train_all.to_numpy()

scaler = StandardScaler()
# x_train: [None, 1,14] -> [None, 2,7]
'''
 x = (x - u) / std -> 正态分布
Z-Score:去除均值和方差缩放
transform：公式为：(X-mean)/std  计算时对每个属性/每列分别进行。
将数据按期属性（按列进行）减去其均值，并处以其方差。得到的结果是，对于每个属性/每列来说所有数据都聚集在0附近，方差为1。
'''
x_train_scaled = scaler.fit_transform(
    x_train.astype(np.float32).reshape(-1, 1)).reshape(-1,2,7,1)
x_valid_scaled = scaler.transform(
    x_valid.astype(np.float32).reshape(-1, 1)).reshape(-1,2,7,1)
x_test_scaled = scaler.transform(
    x_test.astype(np.float32).reshape(-1, 1)).reshape(-1,2,7,1)
print(np.max(x_train_scaled), np.min(x_train_scaled))

#构建可以进行超参数优化的模型。

def build_model(hidden_layers = 1,
                layer_size = 30,dropout_rate=0.5,
                learning_rate = 3e-3):
    """
    构建可以进行超参数优化的模型。
    :param hidden_layers: DNN层数
    :param layer_size: DNN单层神经元个数
    :param dropout_rate: dropout的比率
    :param learning_rate: 学习率
    """
    model = keras.models.Sequential()
    model.add(keras.layers.Conv2D(filters=8, kernel_size=3,
                              padding='same',
                              activation='selu',
                              input_shape=(2, 7, 1)))
    model.add(keras.layers.SeparableConv2D(filters=16, kernel_size=3,
                                       padding='same',
                                       activation='selu')) #selu:带有批归一化功能的relu激活函数,padding = “SAME” 输出大小等于输入大小除以步长向上取整
    model.add(keras.layers.MaxPool2D(pool_size=2))
    model.add(keras.layers.Flatten())
    model.add(keras.layers.Dense(128, activation='selu'))
    model.add(keras.layers.Dense(32,activation='selu'))
    for _ in range(hidden_layers - 1):
        model.add(keras.layers.Dense(layer_size,
                                     activation = 'selu'))
        model.add(keras.layers.Dropout(rate=dropout_rate))
    model.add(keras.layers.Dense(1,activation='sigmoid'))
    optimizer = keras.optimizers.SGD(learning_rate)
    model.compile(loss = 'binary_crossentropy', optimizer = optimizer,metrics = ["accuracy",
                                                                                 tf.keras.metrics.AUC()])
    return model

sklearn_model = keras.wrappers.scikit_learn.KerasRegressor(
    build_fn = build_model)

# f(x) = 1/(x*log(b/a)) a <= x <= b
# 定义超参数搜索范围
param_distribution = {
    "hidden_layers":[1, 2, 3, 4],
    "dropout_rate":np.arange(0,1),
    "layer_size": np.arange(1, 100),
    "learning_rate": reciprocal(1e-4, 1e-2),
}

callbacks = [keras.callbacks.EarlyStopping(patience=5, min_delta=1e-2)]
random_search_cv = RandomizedSearchCV(sklearn_model,
                                      param_distribution,
                                      n_iter = 10,
                                      cv = 3,
                                      scoring = 'neg_log_loss',
                                      n_jobs = -1)
#n_iter=300 训练次数为300次，数值越大，获得的参数精度越大，但是搜索时间越长
#n_jobs = -1 使用所有cpu进行训练
#scoring = 'neg_log_loss' 使用negative_log-loss 为评分标准，解决'log_loss'为负数的问题

history = random_search_cv.fit(x_train_scaled, y_train, epochs = 30,
                     validation_data = (x_valid_scaled, y_valid),
                     callbacks = callbacks)
#循环训练30次

#模型保存与数据输出
def plot_learning_curves(history):
    pd.DataFrame(history.history).plot(figsize=(8, 5))
    plt.grid(True)
    plt.gca().set_ylim(0, 1)
    plt.show()

plot_learning_curves(history)
model = random_search_cv.best_estimator_.model
model.evaluate(x_test_scaled, y_test)
print("Saving model to disk \n")
mp = "E://logs/iris_model.h5"
model.save(mp)
