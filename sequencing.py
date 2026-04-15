import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Embedding, Flatten, Conv1D, MaxPooling1D, Dropout, LSTM, Bidirectional
from tensorflow.keras.utils import to_categorical
from Bio.Seq import Seq
from Bio.Data import CodonTable
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# এক্সেল ফাইল থেকে ডেটা সেট লোড করুন
df = pd.read_excel('CRISPR_Mutation_Dataset_Provisional_200.xlsx', usecols=['Reference DNA', 'Patient DNA', 'Disease'])

reference_sequences = df['Reference DNA'].tolist()
patient_sequences = df['Patient DNA'].tolist()
diseases = df['Disease'].tolist()

def one_hot_encode(seq, max_length):
    encoder = {'A': [1, 0, 0, 0], 'T': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'C': [0, 0, 0, 1]}
    encoded = [encoder[base] for base in seq]
    if len(encoded) < max_length:
        encoded += [[0, 0, 0, 0]] * (max_length - len(encoded))
    return encoded[:max_length]

max_length = max(max(len(seq) for seq in reference_sequences), max(len(seq) for seq in patient_sequences))
encoded_refs = np.array([one_hot_encode(seq, max_length) for seq in reference_sequences])
encoded_patients = np.array([one_hot_encode(seq, max_length) for seq in patient_sequences])

X = np.concatenate((encoded_refs, encoded_patients), axis=0)
y = np.array([0] * len(encoded_refs) + [1] * len(encoded_patients))

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

model = Sequential([
    Conv1D(filters=64, kernel_size=3, activation='relu', input_shape=(max_length, 4)),
    MaxPooling1D(pool_size=2),
    Dropout(0.5),
    Bidirectional(LSTM(64, return_sequences=True)),
    Flatten(),
    Dense(64, activation='relu'),
    Dropout(0.5),
    Dense(1, activation='sigmoid')
])

model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

model.fit(X_train, y_train, epochs=20, batch_size=64, validation_data=(X_test, y_test), verbose=1)

y_pred = (model.predict(X_test) > 0.5).astype(int)
print(f"Accuracy: {accuracy_score(y_test, y_pred):.4f}")
print(f"Precision: {precision_score(y_test, y_pred):.4f}")
print(f"Recall: {recall_score(y_test, y_pred):.4f}")
print(f"F1-Score: {f1_score(y_test, y_pred):.4f}")

def find_mutation_position(ref_seq, patient_seq):
    mutations = []
    for i, (r, p) in enumerate(zip(ref_seq, patient_seq)):
        if r != p:
            mutations.append((i, r, p))
    return mutations

def sift_polyphen_prediction(mutation):
    position, ref_base, patient_base = mutation
    return f"Position {position}: {ref_base} -> {patient_base}"

def visualize_mutations_dual(ref_seq, patient_seq, mutations):
    fig, ax = plt.subplots(figsize=(len(ref_seq) // 2, 4))
    
    nucleotide_colors = {
        'A': 'blue',
        'T': 'yellow',
        'G': 'green',
        'C': 'red'
    }

    for i, base in enumerate(ref_seq):
        color = nucleotide_colors[base]
        ax.add_patch(Rectangle((i, 1.5), 1, 1, color=color, alpha=0.7))

    for i, base in enumerate(patient_seq):
        color = nucleotide_colors[base]
        ax.add_patch(Rectangle((i, 0), 1, 1, color=color, alpha=0.7))

    for mutation in mutations:
        position, ref_base, patient_base = mutation
        ax.add_patch(Rectangle((position, 0), 1, 1, fill=False, edgecolor='black', linewidth=1.5))  
        ax.add_patch(Rectangle((position, 1.5), 1, 1, fill=False, edgecolor='black', linewidth=1.5))  

    plt.text(len(ref_seq)/3, 2.6, "Actual DNA Sequence", fontsize=12, color="black", weight="bold")
    plt.text(len(ref_seq)/3, 1.1, "Mutated DNA Sequence", fontsize=12, color="black", weight="bold")

    legend_patches = [Rectangle((0, 0), 1, 1, color=color, alpha=0.7, label=base) 
                      for base, color in nucleotide_colors.items()]
    ax.legend(handles=legend_patches, loc='lower right', title="Nucleotides", fontsize=8)

    ax.set_xlim(0, len(ref_seq))
    ax.set_ylim(-1, 3)
    ax.axis('off')  
    plt.title("DNA Sequence Mutation Visualization")
    plt.show()

print("\nDNA Analysis Results:\n")
for i in range(len(reference_sequences)):
    ref_seq = reference_sequences[i]
    patient_seq = patient_sequences[i]
    mutations = find_mutation_position(ref_seq, patient_seq)
    mutation_descriptions = [sift_polyphen_prediction(mutation) for mutation in mutations]
    print(f"Reference DNA: {ref_seq}")
    print(f"Patient DNA: {patient_seq}")
    print(f"Disease: {diseases[i]}")
    print(f"Mutations: {mutation_descriptions}\n")

    visualize_mutations_dual(ref_seq, patient_seq, mutations)