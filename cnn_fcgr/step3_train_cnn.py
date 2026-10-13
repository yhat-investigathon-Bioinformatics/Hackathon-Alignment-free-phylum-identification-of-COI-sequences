import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix
import tensorflow as tf
from tensorflow.keras import layers, models, callbacks
from tensorflow.keras.utils import Sequence
import gc
import os

# Constants
BATCH_SIZE = 32  
EPOCHS = 30

class DiskDataGenerator(Sequence):
    def __init__(self, x_path, y_data, batch_size, shuffle=True):
        self.x_path = x_path
        self.y_data = y_data
        self.batch_size = batch_size
        self.shuffle = shuffle
        
        self.x_mmap = np.load(self.x_path, mmap_mode='r')
        self.indices = np.arange(len(self.y_data))
        
        if self.shuffle:
            np.random.shuffle(self.indices)

    def __len__(self):
        return int(np.ceil(len(self.y_data) / self.batch_size))

    def __getitem__(self, index):
        batch_indices = self.indices[index * self.batch_size : (index + 1) * self.batch_size]
        
        X_batch = self.x_mmap[batch_indices]
        y_batch = self.y_data[batch_indices]
        
        return X_batch, y_batch

    def on_epoch_end(self):
        if self.shuffle:
            np.random.shuffle(self.indices)
            gc.collect()

def main():
    print(f"=== TRAINING ===")

    # 1. CHECK INPUT FILES
    path_X_train = "X_train.npy"
    path_X_test = "X_test.npy"
    
    # Abort if input files are missing
    if not os.path.exists(path_X_train) or not os.path.exists(path_X_test):
        print("Error: X_train.npy or X_test.npy not found"); return

    # 2. LOAD LABELS (small, kept in RAM)
    print("\nLoading labels...")
    y_train = np.load("y_train.npy")
    y_test = np.load("y_test.npy")
    print(f"   Train samples: {len(y_train)}, Test samples: {len(y_test)}")

    # 3. ENCODING
    le = LabelEncoder()
    y_train_enc = le.fit_transform(y_train)
    y_test_enc = le.transform(y_test)
    
    class_names = le.classes_
    num_classes = len(class_names)
    print(f"   Classes ({num_classes}): {class_names}")
    
    np.save("classes.npy", class_names)
    print("File 'classes.npy' saved.")

    # 4. GENERATORS
    print("\nInitializing disk data generators...")
    train_gen = DiskDataGenerator(path_X_train, y_train_enc, BATCH_SIZE, shuffle=True)
    # Important: shuffle=False for validation so predictions align with y_test_enc at the end
    val_gen = DiskDataGenerator(path_X_test, y_test_enc, BATCH_SIZE, shuffle=False)

    # 5. Model definition
    model = models.Sequential([
        layers.Input(shape=(64, 64, 1)),
        
        layers.Conv2D(32, (3, 3), padding='same', activation='relu'),
        layers.MaxPooling2D((2, 2)),
        
        layers.Conv2D(64, (3, 3), padding='same', use_bias=False),
        layers.BatchNormalization(),
        layers.Activation('relu'),
        layers.MaxPooling2D((2, 2)),
        
        layers.Conv2D(128, (3, 3), padding='same', use_bias=False),
        layers.BatchNormalization(),
        layers.Activation('relu'),
        layers.MaxPooling2D((2, 2)),
        
        layers.Flatten(),
        layers.Dense(128, activation='relu'),
        layers.Dropout(0.3),
        layers.Dense(num_classes, activation='softmax')
    ])

    model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])

    # Automatically save the best model whenever 'val_loss' improves
    checkpoint = callbacks.ModelCheckpoint(
        'classifier.keras', 
        monitor='val_loss', 
        save_best_only=True, 
        verbose=1
    )
    
    early_stop = callbacks.EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True)
    
    # 7. TRAIN
    history = model.fit(
        train_gen,
        epochs=EPOCHS,
        validation_data=val_gen,
        callbacks=[early_stop, checkpoint]
    )

    # 8. FINAL MANUAL SAVE
    print("\nSaving final model...")
    model.save("classifier.keras")
    print("Model saved as 'classifier.keras'")

    # 9. EVALUATION AND CONFUSION MATRIX
    print("\n" + "="*40)
    print("FINAL EVALUATION")
    print("="*40)
    
    loss, acc = model.evaluate(val_gen, verbose=0)
    print(f"FINAL ACCURACY: {acc*100:.2f}%")
    
    print("Generating confusion matrix...")
    y_pred_probs = model.predict(val_gen)
    y_pred = np.argmax(y_pred_probs, axis=1)
    
    plt.figure(figsize=(10, 8))
    cm = confusion_matrix(y_test_enc, y_pred)
    cm_norm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    
    sns.heatmap(cm_norm, annot=True, fmt='.2f', cmap='Blues', xticklabels=class_names, yticklabels=class_names)
    
    plt.title(f"Confusion Matrix - Acc: {acc:.2%}")
    plt.xlabel('Predicted')  # X axis
    plt.ylabel('True')      # Y axis

    plt.tight_layout()
    plt.savefig("confusion_matrix.png")
    print("Plot saved as 'confusion_matrix.png'")
    plt.show()

if __name__ == "__main__":
    main()