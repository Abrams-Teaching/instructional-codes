# Generic tree for enumerating all points in a parameter space
#
# Cameron F Abrams cfa22@drexel.edu
#
class Node:
    def __init__(self,index=0,label='',value=None,fmt=None):
        self.index=index
        self.label=label
        self.value=value
        self.fmt=fmt
        self.children=[]
    def add_child(self,child):
        self.children.append(child)
    def __str__(self):
        # converts Node to a '-flag value' format
        if self.value is not None:
            return '-'+self.label+' '+self.fmt.format(self.value)
        else:
            return ''
    def dfs(self):
        # returns a depth-first search of tree rooted at self
        # as a list of Nodes
        branches=[]
        for c in self.children:
            if c!=None:
                b=c.dfs()
                branches.extend(b)
        return branches+[self]

class OptionsTree:
    def __init__(self):
        self.root=None
        self.nodes=[]
    def from_dict(self,my_dict):
        self.add_node(label='ROOT')
        ncomb=1
        for d in my_dict.values():
            self.grow_tree(dim=d)
            ncomb*=len(d[1])
        return self,ncomb
    def add_node(self,parent=None,label='',value=None,fmt=None):
        if self.root==None:
            self.root=Node(index=0,label=label,value=value,fmt=fmt)
            self.nodes.append(self.root)
        else:
            newNode=Node(len(self.nodes),label=label,value=value,fmt=fmt)
            parent.add_child(newNode)
            self.nodes.append(newNode)
    def grow_tree(self,dim=None):
        # adds children patterned by 'dim' to every leaf
        # dim is three-element list
        # dim[0] is the label for this axis in parameter space
        # dim[1] is a list of values along that axis
        # dim[2] is the format code for outputting values
        for node in self.root.dfs():
            if len(node.children)==0:
                l,v,s=dim
                for n in v:
                    self.add_node(parent=node,label=l,value=n,fmt=s)
    def tree_paths(self, root):
        if root is None: 
            return []
        if (len(root.children)==0):
            return [str(root)]
        # if left/right is None we'll get empty list anyway
        my_list=[]
        for c in root.children:
            my_list.extend(self.tree_paths(c))
        return [str(root) + ' '+ l for l in my_list]
