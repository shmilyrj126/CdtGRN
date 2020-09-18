import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sklearn
import pandas as pd
import os
import sys
import time
import tensorflow as tf
from tensorflow import keras

physical_devices = tf.config.experimental.list_physical_devices('GPU')
assert len(physical_devices) > 0, "Not enough GPU hardware devices available"
tf.config.experimental.set_memory_growth(physical_devices[0], True)

input_dir="datasets"
df_data = pd.read_csv(os.path.join(input_dir, "train"+".csv"))
df_label = pd.read_csv(os.path.join(input_dir, "valid"+".csv"), header=None)
x_train_all, y_train_all = df_data, df_label

x_valid, x_train = x_train_all.to_numpy().reshape(-1,2,7)[:5000], x_train_all.to_numpy().reshape(-1,2,7)[5000:]
y_valid, y_train = y_train_all.to_numpy()[:5000], y_train_all.to_numpy()[5000:]
x_test, y_test = x_train_all.to_numpy().reshape(-1,2,7),y_train_all.to_numpy()


from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
# x_train: [None, 1,14] -> [None, 7,7]
x_train_scaled = scaler.fit_transform(
    x_train.astype(np.float32).reshape(-1, 1)).reshape(-1,2,7,1)
x_valid_scaled = scaler.transform(
    x_valid.astype(np.float32).reshape(-1, 1)).reshape(-1,2,7,1)
x_test_scaled = scaler.transform(
    x_test.astype(np.float32).reshape(-1, 1)).reshape(-1,2,7,1)
print(np.max(x_train_scaled), np.min(x_train_scaled))

def build_model(hidden_layers = 1,
                layer_size = 30,
                learning_rate = 3e-3):
    model = keras.models.Sequential()
    model.add(keras.layers.Conv2D(filters=8, kernel_size=3,
                              padding='same',
                              activation='selu',
                              input_shape=(2, 7, 1)))
    model.add(keras.layers.SeparableConv2D(filters=16, kernel_size=3,
                                       padding='same',
                                       activation='selu'))
    model.add(keras.layers.MaxPool2D(pool_size=2))
    model.add(keras.layers.Flatten())
    model.add(keras.layers.Dense(128, activation='selu'))
    model.add(keras.layers.Dense(32,activation='selu'))
    for _ in range(hidden_layers - 1):
        model.add(keras.layers.Dense(layer_size,
                                     activation = 'relu'))
    model.add(keras.layers.Dense(1,activation='sigmoid'))
    optimizer = keras.optimizers.SGD(learning_rate)
    model.compile(loss = 'binary_crossentropy', optimizer = optimizer,metrics = ["accuracy",tf.keras.metrics.AUC()])
    return model

logdir = './callbacks'
if not os.path.exists(logdir):
    os.mkdir(logdir)
output_model_file = os.path.join(logdir,
                                 "model.h5")


sklearn_model = keras.wrappers.scikit_learn.KerasRegressor(
    build_fn = build_model)

from scipy.stats import reciprocal
# f(x) = 1/(x*log(b/a)) a <= x <= b

param_distribution = {
    "hidden_layers":[1, 2, 3, 4],
    "layer_size": np.arange(1, 100),
    "learning_rate": reciprocal(1e-4, 1e-2),
}

from sklearn.model_selection import RandomizedSearchCV
callbacks = [keras.callbacks.ModelCheckpoint(output_model_file,save_best_only = True),keras.callbacks.EarlyStopping(patience=5, min_delta=1e-2)]
random_search_cv = RandomizedSearchCV(sklearn_model,
                                      param_distribution,
                                      n_iter = 10,
                                      cv = 3,
                                      n_jobs = 1)
history = random_search_cv.fit(x_train_scaled, y_train, epochs = 30,
                     validation_data = (x_valid_scaled, y_valid),
                     callbacks = callbacks)

def plot_learning_curves(history):
    pd.DataFrame(history.history).plot(figsize=(8, 5))
    plt.grid(True)
    plt.gca().set_ylim(0, 1)
    plt.show()

plot_learning_curves(history)

model = random_search_cv.best_estimator_.model
model.evaluate(x_test_scaled, y_test)