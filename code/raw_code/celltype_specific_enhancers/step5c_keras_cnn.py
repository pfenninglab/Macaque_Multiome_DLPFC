from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv2D, MaxPooling2D, Flatten
from tensorflow.keras.metrics import AUC, TruePositives, FalsePositives, TrueNegatives, FalseNegatives, Precision, Recall
from tensorflow.keras.regularizers import l2
from tensorflow.keras.optimizers import SGD
from tensorflow.keras.models import load_model
import tensorflow as tf
import tensorflow.keras.backend as K
from Bio import SeqIO
import numpy as np
import argparse, sys, os
from tqdm import tqdm
os.system("nvidia-smi")
print(tf.__version__)
print(tf.keras.__version__)

# set seed
from numpy.random import seed
seed(1)
from tensorflow.random import set_seed
set_seed(2)

class CyclicLR(tf.keras.callbacks.Callback):
  '''
  Implementation of cyclic learning rate scheduler (https://arxiv.org/abs/1506.01186)
  '''
  def __init__(self,base_lr, max_lr, step_size, base_m, max_m, cyclical_momentum):
    self.base_lr = base_lr
    self.max_lr = max_lr
    self.base_m = base_m
    self.max_m = max_m
    self.cyclical_momentum = cyclical_momentum
    self.step_size = step_size
    self.clr_iterations = 0.
    self.cm_iterations = 0.
    self.trn_iterations = 0.
    self.history = {}

  def clr(self):
    cycle = np.floor(1+self.clr_iterations/(2*self.step_size))
    if cycle == 2:
      x = np.abs(self.clr_iterations/self.step_size - 2*cycle + 1)
      return self.base_lr-(self.base_lr-self.base_lr/100)*np.maximum(0,(1-x))
    else:
      x = np.abs(self.clr_iterations/self.step_size - 2*cycle + 1)
      return self.base_lr + (self.max_lr-self.base_lr)*np.maximum(0,(1-x))

  def cm(self):
    cycle = np.floor(1+self.clr_iterations/(2*self.step_size))
    if cycle == 2:
      x = np.abs(self.clr_iterations/self.step_size - 2*cycle + 1)
      return self.max_m
    else:
      x = np.abs(self.clr_iterations/self.step_size - 2*cycle + 1)
      return self.max_m - (self.max_m-self.base_m)*np.maximum(0,(1-x))

  def on_train_begin(self, logs={}):
    logs = logs or {}
    if self.clr_iterations == 0:
      K.set_value(self.model.optimizer.lr, self.base_lr)
    else:
      K.set_value(self.model.optimizer.lr, self.clr())
    if self.cyclical_momentum == True:
      if self.clr_iterations == 0:
        K.set_value(self.model.optimizer.momentum, self.cm())
      else:
        K.set_value(self.model.optimizer.momentum, self.cm())

  def on_batch_begin(self, batch, logs=None):
    logs = logs or {}
    self.trn_iterations += 1
    self.clr_iterations += 1
    self.history.setdefault('lr', []).append(K.get_value(self.model.optimizer.lr))
    self.history.setdefault('iterations', []).append(self.trn_iterations)
    if self.cyclical_momentum == True:
      self.history.setdefault('momentum', []).append(K.get_value(self.model.optimizer.momentum))
    for k, v in logs.items():
      self.history.setdefault(k, []).append(v)
    K.set_value(self.model.optimizer.lr, self.clr())
    if self.cyclical_momentum == True:
      K.set_value(self.model.optimizer.momentum, self.cm())

def onehot_seq(seq):
  '''
  Arg:
    seq (string) - DNA sequence
  Returns:
    to_return (2D numpy array)- one hot encoding of DNA sequence
  '''
  letter_to_index =  {'A':0, 'a':0,
                      'C':1, 'c':1,
                      'G':2, 'g':2,
                      'T':3, 't':3}
  to_return = np.zeros((len(seq),4), dtype='int8')
  for idx,letter in enumerate(seq):
    if letter in letter_to_index:
      to_return[idx,letter_to_index[letter]] = 1
  return to_return

def encode_sequence(fasta_pos, fasta_neg, shuffleOff = True):
  '''
  Args:
    fasta_pos - path to fasta file of positive class
    fasta_neg - path to fasta file of negative class
    shuffleOff (optional; default = do nothing) - if False, randomly shuffle order of sequences 
  Returns:
    x (3D numpy array) -  positive and negative class one-hot-encoded DNA sequences
    y (2D numpy array) - binary values representing true labels for x
  '''
  x_pos = np.array([onehot_seq(seq) for seq in tqdm(SeqIO.parse(fasta_pos, "fasta"))] +
  [onehot_seq(seq.reverse_complement()) for seq in tqdm(SeqIO.parse(fasta_pos, "fasta"))])
  x_neg = np.array([onehot_seq(seq) for seq in tqdm(SeqIO.parse(fasta_neg, "fasta"))] +
  [onehot_seq(seq.reverse_complement()) for seq in tqdm(SeqIO.parse(fasta_neg, "fasta"))])
  # concatenate positives and negatives
  print(f'There are {x_pos.shape[0]} positives and {x_neg.shape[0]} negatives.')
  x = np.expand_dims(np.concatenate((x_pos, x_neg)), axis=3)
  y = np.concatenate((np.ones(len(x_pos)),np.zeros(len(x_neg))))
  # need to shuffle order of training set for validation splitting last
  if not shuffleOff:
    indices = np.arange(y.shape[0])
    np.random.shuffle(indices)
    x = x[indices,:]
    y = y[indices]
  return x, y

def macro_f1(y, y_hat, thresh=0.5):
  """Compute the macro F1-score on a batch of observations (average F1 across labels)
  Args:
      y (int32 Tensor): labels array of shape (BATCH_SIZE, N_LABELS)
      y_hat (float32 Tensor): probability matrix from forward propagation of shape (BATCH_SIZE, N_LABELS)
      thresh: probability value above which we predict positive

  Returns:
      macro_f1 (scalar Tensor): value of macro F1 for the batch
  """
  y_pred = tf.cast(tf.greater(y_hat, thresh), tf.float32)
  tp = tf.cast(tf.math.count_nonzero(y_pred * y, axis=0), tf.float32)
  fp = tf.cast(tf.math.count_nonzero(y_pred * (1 - y), axis=0), tf.float32)
  fn = tf.cast(tf.math.count_nonzero((1 - y_pred) * y, axis=0), tf.float32)
  f1 = 2*tp / (2*tp + fn + fp + 1e-16)
  macro_f1 = tf.reduce_mean(f1, axis=-1)
  return macro_f1

def get_model(input_shape, options):
  """Define a TensorFlow CNN architecture
  Args:
      input_shape (2D numpy array): dimensionality of any observation in the dataset; typically
  Returns:
      model : TensorFlow model
  """
  model = Sequential()
  # 5 Convolutional Layers with 250 filters
  for i in range(5):
    # first convolutional layer kernel size is (11,4); any following layer kernel is (11, 1)
    model.add(Conv2D(filters = 250, kernel_size = (11,4 if i == 0 else 1),
                    activation='relu', kernel_initializer='he_normal',
                    bias_initializer='he_normal', kernel_regularizer = l2(l=1e-5),
                    input_shape=input_shape))
  # Pool
  model.add(MaxPooling2D(pool_size=(26,1), strides=26))
  model.add(Flatten())
  # Linear Layers
  for _ in range(1):
    model.add(Dense(units = 300, activation = 'relu',
                    kernel_initializer='he_normal', bias_initializer='he_normal',
                    kernel_regularizer = l2(l=1e-5)
                    )
    )
  # output layer
  model.add(Dense(units = 1, activation = 'sigmoid',
                kernel_initializer='he_normal', bias_initializer='he_normal',
                kernel_regularizer = l2(l=1e-5)))
  myoptimizer = SGD(lr=options.baselr, momentum=options.maxmomentum)
  model.compile(optimizer=myoptimizer,
                loss="binary_crossentropy",
                # show number of true positives, false positives, true negatives, false negatives,
                # auroc, and auprc for training and validation sets during training
                metrics=[
                        TruePositives(name='TP'),FalsePositives(name='FP'),
                        TrueNegatives(name='TN'), FalseNegatives(name='FN'),
                        AUC(name='auroc', curve='ROC'), AUC(name='auprc', curve='PR')
                        ]
                )
  model.summary()
  return model

def train_model_clr(x_train, y_train, x_valid, y_valid, options):
  # compute class weights to weight the loss function
  total = y_train.shape[0]
  weight_for_0 = (1 / np.sum(y_train==0))*(total)/2.0
  weight_for_1 = (1 / np.sum(y_train==1))*(total)/2.0
  class_weight = {0: weight_for_0, 1: weight_for_1}
  # An epoch is calculated by dividing the number of training images by the batchsize
  iterPerEpoch = y_train.shape[0] / options.batch 
  # number of training iterations per half cycle.
  # Authors suggest setting step_size = (2-8) x (training iterations in epoch)
  iterations = list(range(0,round(y_train.shape[0]/options.batch*options.epoch)+1))
  step_size = len(iterations)/2.5
  # set cyclic learning rate
  scheduler =  CyclicLR(base_lr=options.baselr,
              max_lr=options.maxlr,
              step_size=step_size,
              max_m=0.99,
              base_m=options.basemomentum,
              cyclical_momentum=True)
  model = get_model(x_train.shape[1:], options)
  # save best model
  model_checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(filepath=options.name,
                                                                 save_weights_only=False,
                                                                 monitor='val_loss',
                                                                 mode='min',
                                                                 save_best_only=True)
  # train model
  hist = model.fit(x_train,
                    y_train,
                    batch_size = options.batch,
                    epochs = options.epoch,
                    verbose = 1,
                    class_weight = class_weight,
                    validation_data=(x_valid, y_valid),
                    callbacks = [scheduler,
                                model_checkpoint_callback])
  return model, scheduler, hist

def load_data(options):
  # encode and shuffle training sequences
  print("loading training data...")
  (x_train, y_train) = encode_sequence(options.postrain, options.negtrain, shuffleOff = False)
  # encode validation sequences
  print("loading validation data...")
  (x_valid, y_valid) = encode_sequence(options.posval, options.negval, shuffleOff = True)
  return x_train, y_train, x_valid, y_valid

def main(options):
  DIR="../../../data/tidy_data/celltype_specific_enhancers/models_cnn/"
  if not os.path.exists(DIR+options.name.split("_")[0]):
    os.mkdir(DIR+options.name.split("_")[0])
  os.chdir(DIR+options.name.split("_")[0])
  print(options.name.split("_")[0])
  # load data
  x_train, y_train, x_valid, y_valid = load_data(options)
  K.clear_session()
  # train model
  model, clr, hist = train_model_clr(x_train, y_train,
                                     x_valid, y_valid,
                                     options)
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-b", "--batch", default=32, type=int, help="Batch size", required=False)
  parser.add_argument("-e", "--epoch", default=50, type=int, help="Number of epochs", required=False)
  parser.add_argument("-n", "--name", default="cnn.h5", type=str, help="Output model name", required=False)
  parser.add_argument("-pt", "--postrain", type=str, help="Path to positive training set input FASTA file", required=True)
  parser.add_argument("-nt", "--negtrain", type=str, help="Path to negative training set input FASTA file", required=True)
  parser.add_argument("-pv", "--posval", type=str, help="Path to positive training set input FASTA file", required=True)
  parser.add_argument("-nv", "--negval", type=str, help="Path to negative training set input FASTA file", required=True)
  parser.add_argument("-bl", "--baselr", default=1e-5, type=float, help="Base learning rate", required=False)
  parser.add_argument("-ml", "--maxlr", default=0.1, type=float, help="Maximum learning rate", required=False)
  parser.add_argument("-bm", "--basemomentum", default=0.875, type=float, help="Base momentum", required=False)
  parser.add_argument("-mm", "--maxmomentum", default=0.99, type=float, help="Maximum momentum", required=False)

  options, args = parser.parse_known_args()

  if (len(sys.argv)==1):
    parser.print_help(sys.stderr)
    sys.exit(1)
  elif (
        options.postrain is None or 
        options.negtrain is None or 
        options.posval is None or 
        options.negval is None
        ):
    parser.print_help(sys.stderr)
    sys.exit(1)
  else:
    main(options)
