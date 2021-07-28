import numpy as np
from tensorflow import keras
from tensorflow.keras import layers

# Model / data parameters
num_classes = 10
input_shape = (28, 28, 1)

# the data, split between train and test sets
(x_train, y_train), (x_test, y_test) = keras.datasets.mnist.load_data()

# Scale images to the [0, 1] range
x_train = x_train.astype("float32") / 255
x_test = x_test.astype("float32") / 255
# Make sure images have shape (28, 28, 1)
x_train = np.expand_dims(x_train, -1)
x_test = np.expand_dims(x_test, -1)
print("x_train shape:", x_train.shape)
print(x_train.shape[0], "train samples")
print(x_test.shape[0], "test samples")


# convert class vectors to binary class matrices
y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)

model = keras.Sequential(
    [
        keras.layers.InputLayer(input_shape=input_shape),
        layers.Conv2D(32, kernel_size=(3, 3), activation="relu"),
        layers.MaxPooling2D(pool_size=(2, 2)),
        layers.Conv2D(64, kernel_size=(3, 3), activation="relu"),
        layers.MaxPooling2D(pool_size=(2, 2)),
        layers.Flatten(),
        layers.Dropout(0.5),
        layers.Dense(num_classes, activation="softmax"),
    ]
)

model.summary()

batch_size = 128
epochs = 15

model.compile(loss="categorical_crossentropy", optimizer="adam", metrics=["accuracy"])

model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs, validation_split=0.1)












# from pandas import DataFrame
# from pandas import Series
# from pandas import concat
# from pandas import read_csv
# from pandas import datetime
# from sklearn.metrics import mean_squared_error
# from sklearn.preprocessing import MinMaxScaler
# from keras.models import Sequential
# from keras.layers import Dense
# from keras.layers import LSTM
# from math import sqrt
# from matplotlib import pyplot
# import numpy

# # date-time parsing function for loading the dataset
# def parser(x):
# 	return datetime.strptime('190'+x, '%Y-%m')

# # frame a sequence as a supervised learning problem
# def timeseries_to_supervised(data, lag=1):
# 	df = DataFrame(data)
# 	columns = [df.shift(i) for i in range(1, lag+1)]
# 	columns.append(df)
# 	df = concat(columns, axis=1)
# 	df.fillna(0, inplace=True)
# 	return df

# # create a differenced series
# def difference(dataset, interval=1):
# 	diff = list()
# 	for i in range(interval, len(dataset)):
# 		value = dataset[i] - dataset[i - interval]
# 		diff.append(value)
# 	return Series(diff)

# # invert differenced value
# def inverse_difference(history, yhat, interval=1):
# 	return yhat + history[-interval]

# # scale train and test data to [-1, 1]
# def scale(train, test):
# 	# fit scaler
# 	scaler = MinMaxScaler(feature_range=(-1, 1))
# 	scaler = scaler.fit(train)
# 	# transform train
# 	train = train.reshape(train.shape[0], train.shape[1])
# 	train_scaled = scaler.transform(train)
# 	# transform test
# 	test = test.reshape(test.shape[0], test.shape[1])
# 	test_scaled = scaler.transform(test)
# 	return scaler, train_scaled, test_scaled

# # inverse scaling for a forecasted value
# def invert_scale(scaler, X, value):
# 	new_row = [x for x in X] + [value]
# 	array = numpy.array(new_row)
# 	array = array.reshape(1, len(array))
# 	inverted = scaler.inverse_transform(array)
# 	return inverted[0, -1]

# # fit an LSTM network to training data
# def fit_lstm(train, batch_size, nb_epoch, neurons):
# 	X, y = train[:, 0:-1], train[:, -1]
# 	X = X.reshape(X.shape[0], 1, X.shape[1])
# 	model = Sequential()
# 	model.add(LSTM(neurons, batch_input_shape=(batch_size, X.shape[1], X.shape[2]), stateful=True))
# 	model.add(Dense(1))
# 	model.compile(loss='mean_squared_error', optimizer='adam')
# 	for i in range(nb_epoch):
# 		model.fit(X, y, epochs=1, batch_size=batch_size, verbose=0, shuffle=False)
# 		model.reset_states()
# 	return model

# # make a one-step forecast
# def forecast_lstm(model, batch_size, X):
# 	X = X.reshape(1, 1, len(X))
# 	yhat = model.predict(X, batch_size=batch_size)
# 	return yhat[0,0]

# # load dataset
# series = read_csv('shampoo.csv', header=0, parse_dates=[0], index_col=0, squeeze=True, date_parser=parser)

# # transform data to be stationary
# raw_values = series.values
# diff_values = difference(raw_values, 1)

# # transform data to be supervised learning
# supervised = timeseries_to_supervised(diff_values, 1)
# supervised_values = supervised.values

# # split data into train and test-sets
# train, test = supervised_values[0:-12], supervised_values[-12:]

# # transform the scale of the data
# scaler, train_scaled, test_scaled = scale(train, test)

# # fit the model
# lstm_model = fit_lstm(train_scaled, 1, 3000, 4)
# # forecast the entire training dataset to build up state for forecasting
# train_reshaped = train_scaled[:, 0].reshape(len(train_scaled), 1, 1)
# lstm_model.predict(train_reshaped, batch_size=1)

# # walk-forward validation on the test data
# predictions = list()
# for i in range(len(test_scaled)):
# 	# make one-step forecast
# 	X, y = test_scaled[i, 0:-1], test_scaled[i, -1]
# 	yhat = forecast_lstm(lstm_model, 1, X)
# 	# invert scaling
# 	yhat = invert_scale(scaler, X, yhat)
# 	# invert differencing
# 	yhat = inverse_difference(raw_values, yhat, len(test_scaled)+1-i)
# 	# store forecast
# 	predictions.append(yhat)
# 	expected = raw_values[len(train) + i + 1]
# 	print('Month=%d, Predicted=%f, Expected=%f' % (i+1, yhat, expected))

# # report performance
# rmse = sqrt(mean_squared_error(raw_values[-12:], predictions))
# print('Test RMSE: %.3f' % rmse)
# # line plot of observed vs predicted
# pyplot.plot(raw_values[-12:])
# pyplot.plot(predictions)
# pyplot.show()