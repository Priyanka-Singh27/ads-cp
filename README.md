# 🧬 DNA Sequence Matching System (C Project)

## 📌 Overview
This project is a DNA sequence analysis system built purely in C. It supports:
- Fast DNA sequence searching
- Mutation detection
- Sequence alignment
- Similarity scoring
- Multi-species comparison

The system uses advanced data structures:
- Suffix Tree
- Trie
- Hash Table
- Needleman-Wunsch Algorithm
- Skip List

---

# 📂 Project Structure

/src        → Implementation files  
/include    → Header files (interfaces)  
/data       → DNA datasets  

---

# ⚙️ SYSTEM FLOW

1. Load DNA sequences from file  
2. Store sequences using Trie  
3. Build Suffix Tree for fast search  
4. Use Hash Table for quick filtering  
5. Perform alignment using Needleman-Wunsch  
6. Rank results using Skip List  

---

# 👥 TEAM RESPONSIBILITIES (VERY IMPORTANT)

Each member must implement ONE module.

---

## 🌲 MEMBER 1 — SUFFIX TREE

### 📄 File:
`src/suffix_tree.c`

### 🔗 Connected in:
- Called after loading data
- Used in search queries

### 🎯 Responsibilities:
- Build suffix tree for a DNA sequence
- Implement fast pattern search

### 🧠 Functions to implement:
- `build_suffix_tree(const char* text)`
  → Build tree from DNA sequence

- `search_pattern(const char* pattern)`
  → Return 1 if found, else 0

- `free_suffix_tree()`
  → Free memory

### 🚀 Feature:
- Exact DNA matching
- Pattern detection

---

## 🌳 MEMBER 2 — TRIE

### 📄 File:
`src/trie.c`

### 🔗 Connected in:
- Used while loading dataset

### 🎯 Responsibilities:
- Store multiple DNA sequences
- Allow prefix-based search

### 🧠 Functions:
- `create_trie()`
- `insert_sequence(root, seq)`
  → Insert each DNA sequence

- `search_sequence(root, seq)`
  → Check if sequence exists

- `free_trie(root)`

### 🚀 Feature:
- Multi-sequence storage
- Species-level grouping

---

## 🧩 MEMBER 3 — HASH TABLE

### 📄 File:
`src/hashtable.c`

### 🔗 Connected in:
- Used before suffix tree search

### 🎯 Responsibilities:
- Store k-mers (substrings of DNA)
- Provide fast lookup

### 🧠 Functions:
- `insert_kmer(kmer)`
  → Store substring

- `search_kmer(kmer)`
  → Return existence

- `free_table()`

### 🚀 Feature:
- Fast filtering of candidate matches

---

## 🧮 MEMBER 4 — NEEDLEMAN-WUNSCH

### 📄 File:
`src/alignment.c`

### 🔗 Connected in:
- Used after search (when exact match fails)

### 🎯 Responsibilities:
- Compare two DNA sequences
- Return similarity score
- Print alignment

### 🧠 Functions:
- `needleman_wunsch(seq1, seq2)`
  → Return similarity score

- `print_alignment(seq1, seq2)`
  → Show aligned sequences

### 🚀 Feature:
- Mutation detection
- Similarity analysis

---

## 🪜 MEMBER 5 — SKIP LIST

### 📄 File:
`src/skiplist.c`

### 🔗 Connected in:
- Used after alignment

### 🎯 Responsibilities:
- Store results sorted by score
- Return top matches

### 🧠 Functions:
- `insert_skiplist(seq, score)`
- `display_top_matches(k)`
- `free_skiplist()`

### 🚀 Feature:
- Ranking best DNA matches

---

# 🧪 UTILITIES (ALREADY IMPLEMENTED)

## 📄 utils.c

### Functions:
- `validate_sequence(seq)`
  → Checks if sequence contains only A, T, C, G

- `normalize_sequence(seq)`
  → Converts to uppercase

- `load_sample_data(file)`
  → Loads dataset

---

# 🧪 EXPECTED FINAL FEATURES

- Exact DNA search  
- Approximate matching  
- Mutation detection  
- Similarity scoring  
- Multi-species comparison  
- Top-K best matches  

---

# 🚀 HOW TO RUN

```bash
gcc src/*.c -Iinclude -o dna_app
./dna_app