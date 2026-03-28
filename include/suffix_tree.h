#ifndef SUFFIX_TREE_H
#define SUFFIX_TREE_H

typedef struct SuffixTreeNode {
    struct SuffixTreeNode* children[4];
    int start;
    int end;
} SuffixTreeNode;

void build_suffix_tree(const char* text);
int search_pattern(const char* pattern);
void free_suffix_tree();

#endif