import matplotlib.pyplot as plt
import tensorflow as tf
from tf.keras import backend as K

def focal_loss(gamma=2.0, alpha=0.25):
    """
    Focal Loss for binary classification.
    
    Arguments:
        gamma: focusing parameter for modulating factor (1-p)^gamma.
        alpha: balance parameter for positive class.
    """
    def focal_loss_fixed(y_true, y_pred):
        y_true = tf.cast(y_true, tf.float32)
        # Clip predictions to prevent log(0)
        y_pred = tf.clip_by_value(y_pred, K.epsilon(), 1 - K.epsilon())
        # Compute binary crossentropy
        bce = - (y_true * tf.math.log(y_pred) + (1 - y_true) * tf.math.log(1 - y_pred))
        # Compute modulating factor
        modulating_factor = tf.where(tf.equal(y_true, 1), 1 - y_pred, y_pred)
        modulating_factor = tf.pow(modulating_factor, gamma)
        # Compute alpha weighting
        alpha_factor = tf.where(tf.equal(y_true, 1), alpha, 1 - alpha)
        # Compute final focal loss
        loss = alpha_factor * modulating_factor * bce
        return K.mean(loss)
    
    return focal_loss_fixed

def plot_training_history(history, figsize=(12,5)):
    """
    Plots training and validation loss and accuracy curves from a Keras history object.

    Arguments:
    history : the object returned by model.fit()
    figsize : figure size (width, height)
    """
    plt.figure(figsize=figsize)

    # Plot Loss
    plt.subplot(1,2,1)
    plt.plot(history.history['loss'], label='Train Loss', color='blue')
    plt.plot(history.history['val_loss'], label='Validation Loss', color='orange')
    plt.title('Model Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()

    # Plot Accuracy
    plt.subplot(1,2,2)
    plt.plot(history.history['accuracy'], label='Train Accuracy', color='blue')
    plt.plot(history.history['val_accuracy'], label='Validation Accuracy', color='orange')
    plt.title('Model Accuracy')
    plt.xlabel('Epoch')
    plt.ylabel('Accuracy')
    plt.legend()

    plt.tight_layout()
    plt.show()
