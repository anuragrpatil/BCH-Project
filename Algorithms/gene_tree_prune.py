import os
import gzip
from Bio import SeqIO
import re
from newick import read, dump


def clean_node(node_name):
    regex = r"[a-z]_[a-zA-Z]+_[a-zA-Z1-9]+.[a-zA-Z1-9]+.[a-zA-Z1-9]+.[a-zA-Z1-9]+"
    if re.search(regex, str(node_name)):
        node_name.name = re.search(regex, str(node_name)).group(0)

def main():
    filename = ['OG0000002_tree.txt', ]
    for file in filename:
        #Read Dataset
        trees = read('./dataset/{}'.format(file))
        node_needed = trees[0].get_node('n441').get_leaves()
        all_node = trees[0].get_leaves()

        trees[0].prune(node_needed, inverse=True)

        #Rename Nodes
        trees[0].visit(clean_node)

        #Get node names
        leaves = trees[0].get_leaves()

        #Get Node lengths
        unique = {}
        for node in leaves:
            regex = r"[a-z]_[a-zA-Z]+_[a-zA-Z1-9]"
            if re.search(regex, node.name):
                node_new_name = re.search(regex, str(node.name)).group(0)
            if node_new_name in unique:
                unique[node_new_name] = max(unique[node_new_name], node.length)
            else:
                unique[node_new_name] = node.length 

        #prune to just keep the longest unique nodes
        for node in leaves:
            regex = r"[a-z]_[a-zA-Z]+_[a-zA-Z1-9]"
            if re.search(regex, node.name):
                node_new_name = re.search(regex, str(node.name)).group(0)
            if node.length != unique[node_new_name]:
                trees[0].prune_by_names(node.name)

        #remove c elegans ref
        c_elegans_remove = ['c_elegans_ref_protein_PAR-1']
        trees[0].prune_by_names(c_elegans_remove)

        #Remove nodes with no chil
        while trees[0].walk(mode='postorder'):
            atleast_once = True
            for n in trees[0].walk(mode='postorder'):
                regex = r"^n[1-9]+"
                if n.ancestor and len(n.descendants) == 0 and re.search(regex, n.name):
                    trees[0].prune_by_names(n.name)
                    atleast_once = False
            if atleast_once == True:
                break

        #dump the tree
        with open('{}_required_tree.txt'.format(file), 'w') as fobj:
            dump(trees, fobj)


if __name__ == "__main__":
    main()


