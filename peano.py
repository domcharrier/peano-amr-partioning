#!/usr/bin/env python3
class TreeNode:
  def __init__(self,dim):
      self.hasChildren = False
      if dim==2:
        self.children = [[[None for x in range(0,3)] for y in range(0,3)] for z in range(0,1)] 
      elif dim==3:
        self.children = [[[None for x in range(0,3)] for y in range(0,3)] for z in range(0,3)] 

  def putChild(self,pos,node):
      self.hasChildren = True
      self.children[pos[2]][pos[1]][pos[0]] = node

  def getChild(self,pos):
      return self.children[pos[2]][pos[1]][pos[0]]

def linearise(pos,dim):
    ix=pos[0]
    iy=pos[1]
    iz=pos[2] if dim==3 else 0
    return ix + 3*iy + 9*iz

def delinearise(ind):
    pos = [0]*3
    pos[0] = ind % 3
    rem = ind - pos[0]
    pos[1] = int ( rem % 9 / 3 )
    rem -= 3*pos[1]
    pos[2] = int ( rem / 9 )
    return pos 

def peanoIndexToPos2D(childPeanoIndex,parentMotif):
    """
    Given a parent index, translate the Peano SFC index of the
    child to a position in the parent's array of children.

    :param childPeanoIndex: Peano index of the child
    :param parentMotif:     Peano SFC motif of the parent (P,Q,R,S).
    
    :see: M. Bader, Space-Filling Curves, vol. 9. Berlin, Heidelberg: Springer Berlin Heidelberg, 2013.

    Details:
    
    P motif:

    The indexing
    
    6 7 8 | 0 1 2 | 6 7 8 
    5 4 3 | 5 4 3 | 5 4 3 
    0 1 2 | 6 7 8 | 0 1 2 
    ---------------------
    8 7 6 | 2 1 0 | 8 7 6
    3 4 5 | 3 4 5 | 3 4 5
    2 1 0 | 8 7 6 | 2 1 0
    ---------------------
    6 7 8 | 0 1 2 | 6 7 8 
    5 4 3 | 5 4 3 | 5 4 3 
    0 1 2 | 6 7 8 | 0 1 2

    is describes the P replacement pattern:

    P ->

    | P | Q | P | 
    -------------
    | S | R | S |
    -------------
    | P | Q | P |

    We can decude the other indices from here.

    """
    indexMap = {}
    indexMap["P"] = [0,1,2,5,4,3,6,7,8] 
    indexMap["Q"] = [6,7,8,5,4,3,0,1,2] 
    indexMap["R"] = [8,7,6,3,4,5,2,1,0] 
    indexMap["S"] = [2,1,0,3,4,5,8,7,6]
    linearisedPos = indexMap[parentMotif].index(childPeanoIndex) 
    return delinearise(linearisedPos);

def posToMotif2D(childPos,parentMotif):
    """
    Refined cells get assigned a motif.
    The motif determines how the Peano
    curve is mapped to a subsquare.

    The motif of a (refined) child cell depends on 
    the motif of its (refined) parent
    plus its position with respect to the parent.

    :param childPos:    position of the child in the parent's array of children.
    :param parentMotif: Peano SFC motif of the child's parent.

    :see: M. Bader, Space-Filling Curves, vol. 9. Berlin, Heidelberg: Springer Berlin Heidelberg, 2013.
    """
    motifs = {}
    motifs["P"] = ["P", "Q", "P", "S", "R", "S", "P", "Q", "P"] 
    motifs["Q"] = ["Q", "P", "Q", "R", "S", "R", "Q", "P", "Q"] 
    motifs["R"] = ["R", "S", "R", "Q", "P", "Q", "R", "S", "R"] 
    motifs["S"] = ["S", "R", "S", "P", "Q", "P", "S", "R", "S"]
    return motifs[parentMotif][linearise(childPos,2)]


def refineTree(node,offset,size,refinementCriterion,l=0):
    if refinementCriterion(offset,size,l):
        dim = len(offset)
        izmax = 3 if dim == 3 else 1
        for iz in range(0,izmax):
            for iy in range(0,3):
                for ix in range(0,3):
                    pos = [ix,iy,iz]
                    child = TreeNode(dim)
                    node.putChild(pos,child)

                    childSize = size.copy()
                    for i,val in enumerate(childSize):
                        childSize[i] /= 3.0
                    childOffset = offset.copy()
                    for i,val in enumerate(childOffset):
                        childOffset[i] += pos[i]*childSize[i];
                    refineTree(child,childOffset,childSize,refinementCriterion,l+1)

def computePeanoSFCCoordinates(centres,offsets,sizes,pos,parentMotif,node,offset,size,l,firstLevel):
    """
    :param centres: list to store the coordinates in
    :param offsets: list to store the coordinates in
    :param sizes:   list to store the coordinates in
    """
    if l >= firstLevel and not node.hasChildren:
        centre = offset.copy()
        for i,val in enumerate(size):
            centre[i] += 0.5 * size[i]
        centres.append(centre)
        offsets.append(offset)
        sizes.append(size)

    for peanoIndex in range(0,3**2):
        childPos = peanoIndexToPos2D(peanoIndex,parentMotif) 
        child    = node.getChild(childPos)
        
        if child != None:
            childSize = size.copy()
            for i,val in enumerate(childSize):
                childSize[i] /= 3.0
            childOffset = offset.copy()
            for i,val in enumerate(childOffset):
                childOffset[i] += childPos[i]*childSize[i]
            childMotif = posToMotif2D(childPos,parentMotif)
            computePeanoSFCCoordinates(centres,offsets,sizes,childPos,childMotif,child,childOffset,childSize,l+1,firstLevel)

if __name__ == "__main__":
    import math

    dim  = 2
    domainOffset = [0.0,0.0] 
    domainSize   = [3.0,3.0]

    def refinementCriterion(offset,size,l):
        """
        Bounding box scaling.
        """
        lmax = 3
        if l>=lmax:
            return False 
        
        centre = [0]*2
        if l > 0:
            for i,dx in enumerate(size):
                centre[i] = offset[i] + 0.5*dx
                if centre[i] < (1.0 - dx) or\
                   centre[i] > (2.0 + dx):
                    return False
            if centre[0] > 1.0 and\
               centre[0] <= 2.0 and\
               centre[1] > 1.0 and\
               centre[1] <= 2.0:
                return False
        return True
    
    def refinementCriterion2(offset,size,l):
        """
        Delayed broadcasts.
        """
        lmax = 3
        if l>=lmax:
            return False 
        elif l < 1:
            return True
        
        centre = [0]*2
        for i,dx in enumerate(size):
           centre[i] = offset[i] + 0.5*dx
        if centre[0] > 2.0 and\
           centre[1] > 2.0:
            return True
        return False
    
    def refinementCriterion3(offset,size,l):
        """
        Enclaves mesh.
        """
        lmax = 3
        if l<lmax-1:
            return True
        elif l==lmax-1:
             centre = [0]*2
             for i,dx in enumerate(size):
                 centre[i] = offset[i] + 0.5*dx
             if centre[0] > 1.0 and centre[0] < 2.0 and\
                centre[1] > 1.0 and centre[1] < 2.0:
                 return True
        return False
    
    def refinementCriterion4(offset,size,l):
        """
        SFC-Cuts Sphere
        """
        lmax = 3
        x0   = [domainSize[0]/2,domainSize[1]/2]
        if l<1:
            return True
        elif l<lmax:
             centre = [0]*2
             for i,dx in enumerate(size):
                 centre[i] = offset[i] + 0.5*dx
             r = math.sqrt( (centre[0]-x0[0])**2 + (centre[1]-x0[1])**2 )
             if r > 0.99-size[0]*2/3 and r < 0.99+size[0]*2/3:
                 return True
        return False

    root = TreeNode(dim)
    refineTree(root,domainOffset,domainSize,refinementCriterion4)

    centres = []
    offsets = []
    sizes   = []
    computePeanoSFCCoordinates(centres,offsets,sizes,[0,0],"P",root,domainOffset,domainSize,0,1)
    
    partitions = 4
    partitionSize = int(len(offsets)/partitions) # floors
    currentPartition=0;
    for i in range(0,len(offsets)):
        currentPartition = min(partitions-1,int(i/partitionSize)) # floors
        lowerLeft  = "%f, %f" % ( offsets[i][0], offsets[i][1] )
        upperRight = "%f, %f" % ( offsets[i][0]+sizes[i][0], offsets[i][1]+sizes[i][1] )
        print("\\draw[fill=c%d,draw=black] (%s) rectangle (%s);" % ( currentPartition, lowerLeft, upperRight) ) 
    # print sfc
    coords = ["(%1.3f,%1.3f)" % (x[0],x[1]) for x in centres]
    print("\\draw[] %s ;" % " -- ".join(coords))
