import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import LabelEncoder
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Embedding, Flatten
from Bio.Seq import Seq
from Bio.Data import CodonTable

# এক্সেল ফাইল থেকে ডেটা সেট লোড করুন
df = pd.read_excel('CRISPR_Mutation_Dataset_Provisional_200.xlsx', usecols=['Reference DNA', 'Patient DNA', 'Disease'])

reference_sequences = df['Reference DNA'].tolist()
patient_sequences = df['Patient DNA'].tolist()
diseases = df['Disease'].tolist()

def find_mutation_position(ref_seq, patient_seq):
    mutations = []
    for i, (r, p) in enumerate(zip(ref_seq, patient_seq)):
        if r != p:
            mutations.append((i, r, p))
    return mutations

def sift_polyphen_prediction(mutation):
    position, ref_base, patient_base = mutation
    return f"Position {position}: {ref_base} -> {patient_base}"

def encode_sequence(seq, max_length):
    encoder = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    encoded = [encoder[base] for base in seq]
    return encoded + [0] * (max_length - len(encoded))

def to_binary_vector(seq, length):
    vector = [0] * (4 * length)
    for i, base in enumerate(seq):
        vector[i * 4 + base] = 1
    return vector

def tanimoto_similarity(vec1, vec2):
    intersection = np.dot(vec1, vec2)
    union = np.sum(vec1) + np.sum(vec2) - intersection
    return intersection / union

max_length = max(max(len(seq) for seq in reference_sequences), max(len(seq) for seq in patient_sequences))
encoded_refs = [encode_sequence(seq, max_length) for seq in reference_sequences]
encoded_patients = [encode_sequence(seq, max_length) for seq in patient_sequences]

binary_refs = [to_binary_vector(seq, max_length) for seq in encoded_refs]
binary_patients = [to_binary_vector(seq, max_length) for seq in encoded_patients]

tanimoto_scores = [tanimoto_similarity(ref, pat) for ref, pat in zip(binary_refs, binary_patients)]

X = np.array(encoded_refs + encoded_patients)
y = np.array([0] * len(encoded_refs) + [1] * len(encoded_patients))

model = Sequential([
    Embedding(input_dim=4, output_dim=8, input_shape=(max_length,)),
    Flatten(),
    Dense(10, activation='relu'),
    Dense(1, activation='sigmoid')
])

model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.fit(X, y, epochs=10, verbose=1)

predictions = [model.predict(np.array([encoded_patient]))[0][0] for encoded_patient in encoded_patients]

print("DNA Analysis Results:\n")
for i in range(len(reference_sequences)):
    mutations = find_mutation_position(reference_sequences[i], patient_sequences[i])
    mutation_descriptions = [sift_polyphen_prediction(mutation) for mutation in mutations]
    print(f"Reference DNA: {reference_sequences[i]}")
    print(f"Patient DNA: {patient_sequences[i]}")
    print(f"Disease: {diseases[i]}")
    print(f"Mutations: {mutation_descriptions}")
    print(f"Tanimoto Similarity: {tanimoto_scores[i]:.2f}")
    print(f"Predicted Mutation Impact: {predictions[i]:.4f}\n")