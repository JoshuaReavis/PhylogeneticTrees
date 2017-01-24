/* 
 * PhyloTree.java
 *
 * Defines a phylogenetic tree, which is a strictly binary tree 
 * that represents inferred hierarchical relationships between species
 * 
 * There are weights along each edge; the weight from parent to left child
 * is the same as parent to right child.
 *
 * Joshua Reavis
 * Spring 2016
 */

import java.lang.*;
import java.util.*;
import java.io.*;
 

public class PhyloTree {
    ArrayList<Species> allSpecies = new ArrayList<Species>(); // used for getAllSpecies()
    public HashMap<String, PhyloTreeNode> forest = new HashMap<String, PhyloTreeNode>(); // initializes our forest of trees
    private PhyloTreeNode overallRoot;    // The actual root of the overall tree
    private int printingDepth;            // How many spaces to indent the deepest 
                                          // node when printing

    // CONSTRUCTOR

    // PhyloTree
    //        - speciesFile contains the path of a valid FASTA input file
    //        - printingDepth is assumed as a positive number

    //        - Creates a linked tree structure representing the inferred hierarchical
    //          species relationship has been created, and overallRoot points to
    //          the root of this tree

    public PhyloTree(String speciesFile, int printingDepth) {
        this.printingDepth = printingDepth;
        buildTree(loadSpeciesFile(speciesFile));        
    }

    // ACCESSORS

    //    - Returns the overall root
    public PhyloTreeNode getOverallRoot() {
        return this.overallRoot;
    }

    //    - Returns a string representation of the tree starting at the overall root node
    public String toString() {       
        return toString(getOverallRoot(), 0, getWeightedHeight());
    }


    // Assumptions:
    //    - node points to the root of a tree you intend to print
    //    - weightedDepth is the sum of the edge weights from the
    //    - overall root to the current root
    //    - maxDepth is the weighted depth of the overall tree

    // Post-conditions:
    //    - Returns a string representation of the subtree starting at the given node
    private String toString(PhyloTreeNode node, double weightedDepth, double maxDepth) {
        StringBuilder stringbuild = new StringBuilder();
        
        if (node.getRightChild() != null){
           String addRightChild = toString(node.getRightChild(), weightedNodeDepth(node.getRightChild()), maxDepth);
           stringbuild.append(addRightChild);
        }
        
        int k = (int)(printingDepth * (weightedDepth / maxDepth));
        while (k != 0){
           stringbuild.append(".");
           k--;
        }
        stringbuild.append(node.toString() + "\n");
        
        if (node.getLeftChild() != null){
           String addLeftChild = toString(node.getLeftChild(), weightedNodeDepth(node.getLeftChild()), maxDepth);
           stringbuild.append(addLeftChild);
        }
        
        return stringbuild.toString();
    }

 
    //    - Returns a string representation in tree format starting at the overall root node
    public String toTreeString() {
        return toTreeString(getOverallRoot());
    }


    //    - Returns a string representation in tree format for a given subtree whose root is the node given
    private String toTreeString(PhyloTreeNode node) {
        StringBuilder stringbuild = new StringBuilder();
        if (node.isLeaf()){
           stringbuild.append(node.getLabel() + ":" + String.format("%.5f", node.getParent().getDistanceToChild())); 
        }else{
           stringbuild.append("(" + toTreeString(node.getRightChild()) + "," + toTreeString(node.getLeftChild()) + ")");
           if (node != getOverallRoot()){
           stringbuild.append(":" + String.format("%.5f", node.getParent().getDistanceToChild()));
           } 
        }
        return stringbuild.toString();
    }


    //    - Returns the tree height as defined in class
    public int getHeight() {
       return nodeHeight(getOverallRoot()); 
    }
    



    //    - Returns the sum of the edge weights along the "longest" (highest weight) path from the overall root to any leaf node.
    public double getWeightedHeight() {      
        return weightedNodeHeight(getOverallRoot());
    }

 
    //    - Returns the number of species in the tree starting from the overall root
    public int countAllSpecies() {
        return countAllSpeciesHelper(getOverallRoot());
    }
    
    //    - Returns the number of species in the tree starting with the given root
    public int countAllSpeciesHelper(PhyloTreeNode Tree){
       return getAllSpecies().size();
             
    }

 
    //    - Returns an ArrayList containing all species in the tree
    public java.util.ArrayList<Species> getAllSpecies() {
       return allSpecies;
    }
 

    //    - label is the label of a tree node you intend to find
    //    - Assumes labels are unique in the tree
    //    - If found: returns the PhyloTreeNode with the specified label
    //    - If not found: returns null
    public PhyloTreeNode findTreeNodeByLabel(String label) {
       return findTreeNodeByLabel(getOverallRoot(), label);
    }



    //    - label1 and label2 are the labels of two species in the tree
    //    - If either node cannot be found: returns null
    //    - If both nodes can be found: returns the PhyloTreeNode of their common ancestor with the largest depth
     public PhyloTreeNode findLeastCommonAncestor(String label1, String label2) {
        ArrayList<PhyloTreeNode> node1Ancestors = new ArrayList<PhyloTreeNode>(); // ArrayList containing ancestors of node 1
        ArrayList<PhyloTreeNode> node2Ancestors = new ArrayList<PhyloTreeNode>(); // ArrayList containing ancestors of node 2
                

        // finds the nodes for our labels
        PhyloTreeNode node1 = findTreeNodeByLabel(getOverallRoot(), label1);
        PhyloTreeNode node2 = findTreeNodeByLabel(getOverallRoot(), label2);
        
        if (node1 == null | node2 == null){
           return null;
        }        
        
        // nodeiterators so we can iterator through the parents of the nodes
        PhyloTreeNode node1Iterator = node1;
        PhyloTreeNode node2Iterator = node2;
        
        // populates node1Ancestor with all of node1's ancestors
        while (node1Iterator.getParent() != null){
           node1Ancestors.add
           (node1Iterator.getParent());
           node1Iterator = node1Iterator.getParent();
        }
        
        node2Ancestors.add(node2Iterator.getParent());
        node2Iterator = node2Iterator.getParent();
        int i = 0;
        
        // check to see if each closest ancestor of node2 is also an ancestor of node1
        while (!node1Ancestors.contains(node2Ancestors.get(0))){
           node2Ancestors.remove(0);
           node2Ancestors.add(node2Iterator.getParent());
           node2Iterator = node2Iterator.getParent();
           i++;   
                
        }
               
        return node2Ancestors.get(0);
    }
    
    // findEvolutionaryDistance
    //    - label1 and label2 are the labels of two species in the tree
    //    - If either node cannot be found: returns POSITIVE_INFINITY
    //    - If both nodes can be found: returns the sum of the weights  along the paths from their least common ancestor to each of the two nodes
     public double findEvolutionaryDistance(String label1, String label2) {
     
        PhyloTreeNode node1 = findTreeNodeByLabel(label1);
        PhyloTreeNode node2 = findTreeNodeByLabel(label2);
        
        if (node1 == null | node2 == null){
           return java.lang.Double.POSITIVE_INFINITY;
        }
        
        PhyloTreeNode ancestor = findLeastCommonAncestor(node1, node2);
        
        double distance1 = 0;
        while (!node1.getLabel().equals(ancestor.getLabel())){
           node1 = node1.getParent();
           distance1 += node1.getDistanceToChild();
        }
        
        double distance2 = 0;
        while (!node2.getLabel().equals(ancestor.getLabel())){
           node2 = node2.getParent();
           distance2 += node2.getDistanceToChild();
        }
        
        return distance1+distance2;
    }

    // MODIFIER

    //    - species contains the set of species for which you want to infer a phylogenetic tree
    //    - Creates a linked tree structure representing the inferred hierarchical species relationship has been created, and overallRoot points to the root of said tree
    private void buildTree(Species[] species) {
       MultiKeyMap<Double> Tdistances = new MultiKeyMap<Double>(); // initializes our tree Tdistances multikeymap
       PhyloTreeNode Tree1 = null;
       PhyloTreeNode Tree2 = null;       
       
       // creates a tree in the forest for each species with null parent
       for (int i = 0; i < species.length; i++){
          forest.put(species[i].getName(), new PhyloTreeNode(null, species[i]));
          this.allSpecies.add(species[i]);

       }
       
       // fills out our Tdistances map with the Tdistances between two trees
       for (String key1 : forest.keySet()){
          for (String key2 : forest.keySet()){
             if (forest.get(key1) != forest.get(key2)){
             // for the third parameter we need to have 2 species, so we get the species associated with the tree in forest hashmap
             double dist = Species.distance(forest.get(key1).getSpecies(), forest.get(key2).getSpecies());
             Tdistances.put(key1, key2, dist);
             }         
          }
       }       
       
       PhyloTreeNode Tnew = null;
       // Consolidates the forest of trees into a single tree containing branched tree nodes
       while (forest.size() > 1){
          // finds the two trees with the smallest distance between them
          double minDistance = java.lang.Double.POSITIVE_INFINITY;
          for (String key1 : forest.keySet()){
             for (String key2 : forest.keySet()){
                if (forest.get(key1) != forest.get(key2) && Tdistances.get(key1, key2) < minDistance){
                   minDistance = Tdistances.get(forest.get(key1).getLabel(), forest.get(key2).getLabel());
                   if (key1.compareTo(key2) <= 0){
                      Tree1 = forest.get(key1);
                      Tree2 = forest.get(key2);
                   }else{
                      Tree1 = forest.get(key2);
                      Tree2 = forest.get(key1);
                   }
                }
             }   
          }
          
          // removes the two trees from the forest
          forest.remove(Tree1.getLabel());
          forest.remove(Tree2.getLabel());
          
          // Puts the new PhyloTreeNode into the forest 
          Tnew = new PhyloTreeNode(Tree1.getLabel() + Tree2.getLabel(), null, Tree1, Tree2, minDistance/2.0);    
          forest.put(Tree1.getLabel() + Tree2.getLabel(), Tnew);
          
          // sets the parents of the copy of our tree nodes
          Tree1.setParent(Tnew);
          Tree2.setParent(Tnew);
          
          // Calculates the distance of the new tree to all of the other trees in the forest
          double count1 = Tree1.getNumLeafs();
          double count2 = Tree2.getNumLeafs();
         
          for (String key : forest.keySet()){
                if (!key.equals(Tnew.getLabel())){
                   //System.out.println(Tdistances.get(key, Tree1.getLabel()));
                   double TotherT1 = Tdistances.get(key, Tree1.getLabel());
                   //Tdistances.remove(key, Tree1.getLabel());
                   double TotherT2 = Tdistances.get(key, Tree2.getLabel());
                   //Tdistances.remove(key, Tree2.getLabel());
                                     
                   double dist1 = ((count1 / (count1 + count2)) * TotherT1);
                   double dist2 = ((count2 / (count2 + count1)) * TotherT2);
                   double Tdist = dist1 + dist2;             
                   Tdistances.put(Tnew.getLabel(), key, Tdist);
                }
          }
          
          this.overallRoot = Tnew;
       
       } 
    }

    // STATIC

    //    - node is null or the root of tree (possibly subtree)
    //    - If null: returns -1
    //    - Else: returns the depth of the node within the overall tree
    public static int nodeDepth(PhyloTreeNode node) {
        if (node == null){
           return -1;   
        }   
        int count = 0;
        
        while (node != null){
           count += 1;
           node = node.getParent();
        }
        
        return count; 
    }
    
    public static double weightedNodeDepth(PhyloTreeNode node){
        if (node == null){
           return java.lang.Double.NEGATIVE_INFINITY;   
        }   
        double weightedDepth = 0.0;
        
        while (node.getParent() != null){
           node = node.getParent();
           weightedDepth += node.getDistanceToChild();
        }
        
        return weightedDepth;        
    }

    //    - node should be null or the root of tree (possibly subtree)
    //    - If null: returns -1
    //    - Else: returns the height subtree rooted at node
    public static int nodeHeight(PhyloTreeNode node) {
       if (node == null){
          return -1;
       }else if (node.isLeaf()){
          return 0;
       }else{
          int countLeft = 1 + nodeHeight(node.getLeftChild());
          int countRight = 1 + nodeHeight(node.getRightChild());
          if (countLeft >= countRight){
             return countLeft;
          }else{
             return countRight;
          }
       }
    }

    //    - node is null or the root of tree (possibly subtree)
    //    - If null: returns NEGATIVE_INFINITY
    //    - Else: returns the weighted height subtree rooted at node
    public static double weightedNodeHeight(PhyloTreeNode node) {
        if (node == null){
           return java.lang.Double.NEGATIVE_INFINITY;
       }else if (node.isLeaf()){
          return 0;
       }else{
          double weightLeft = node.getDistanceToChild() + weightedNodeHeight(node.getLeftChild());
          double weightRight = node.getDistanceToChild() + weightedNodeHeight(node.getRightChild());
          if (weightLeft >= weightRight){
             return weightLeft;
          }else{
             return weightRight;
          }
       }        
   }


    //    - filename contains the path of a valid FASTA input file
    //    - Creates and returns an array of species objects representing all valid species in the input file
    //    - Species without names are skipped
    public static Species[] loadSpeciesFile(String filename) {
        ArrayList<Species> loadIn = new ArrayList<Species>(); // ArrayList which loads in all of our Species
        ArrayList<String> speciesName = new ArrayList<String>(); // ArrayList of species names
        ArrayList<String> sequence = new ArrayList<String>(); // ArrayList for the sequences
        String array[];
        
        try{        
        Scanner input = new Scanner(new File(filename));
        Scanner input2 = new Scanner (new File(filename)); //2 scanners becomes helpful later for checking the line after the one we are on
        String line = input.nextLine();
        String lineNext = input2.nextLine();
        lineNext = input2.nextLine();
        
        // runs through the whole file
        while (input.hasNextLine()){
        
        // if there's a ref in the line, grab the sequence and name
           if (line.contains("ref")){ // if we are on a ref index
              array = line.split("\\|");
              speciesName.add(array[6]);  
              String sequencebuild = "";  
              while (!lineNext.contains(">") && input.hasNextLine()){
                 line = input.nextLine();
                 if (input2.hasNextLine()){
                    lineNext = input2.nextLine();
                 }
                 sequencebuild += line;
                 line = "";
              }
              sequence.add(sequencebuild);
           
           }
           // continue along the lines with the scanner if we did not incur a ref   
           else{
              line = input.nextLine();
              if (input2.hasNextLine()){
                 lineNext = input2.nextLine();
              }
           }
        } 
                                      
           // creates the species and adds it to the master arraylist of species
           for (int i = 0; i < sequence.size(); i++){
              String temp[] = sequence.get(i).split("");
              Species species = new Species(speciesName.get(i), temp);
              loadIn.add(species);               
           }

        
        
        }catch(FileNotFoundException ex){
        }
        
        //transfers the master species arraylist into an array        
        Species result[] = new Species[loadIn.size()];
        for (int i = 0; i < loadIn.size(); i++){
           result[i] = loadIn.get(i);
        }
                
        return result;
    }



    //    - node points to a node in a phylogenetic tree structure
    //    - descendants is a non-null reference variable to an empty arraylist object
    //    - descendants is populated with all species in the subtree rooted at node
    private static void getAllDescendantSpecies(PhyloTreeNode node, java.util.ArrayList<Species> descendants) {
       if (node.getLeftChild() != null && node.getRightChild() !=null){
          descendants.add(node.getLeftChild().getSpecies());
          descendants.add(node.getRightChild().getSpecies());
          getAllDescendantSpecies(node.getLeftChild(), descendants);
          getAllDescendantSpecies(node.getRightChild(), descendants);                
       }      
    }   
    

    //    - node points to a node in a phylogenetic tree structure
    //    - label is the label of a tree node that you intend to locate
    //    - If no node with the label exists in the subtree, return null
    private static PhyloTreeNode findTreeNodeByLabel(PhyloTreeNode node, String label) {
      if (node.getLabel().equals(label)) {
         return node;
      }
      if (node.isLeaf()) {
         return null;
      } 
      else {
         PhyloTreeNode left = findTreeNodeByLabel(node.getLeftChild(), label);
         PhyloTreeNode right = findTreeNodeByLabel(node.getRightChild(), label);
      
         if (left != null) {
            return left;
         } 
         else {
            return right;
         }
      }

    }

    //    - Assume node1 and node2 point to nodes in the phylogenetic tree
    //    - If node1 or node2 are null, return null
    //    - Else: returns the PhyloTreeNode of their common ancestor with the largest depth
     private static PhyloTreeNode findLeastCommonAncestor(PhyloTreeNode node1, PhyloTreeNode node2) {
        ArrayList<PhyloTreeNode> node1Ancestors = new ArrayList<PhyloTreeNode>(); // ArrayList containing ancestors of node 1
        ArrayList<PhyloTreeNode> node2Ancestors = new ArrayList<PhyloTreeNode>(); // ArrayList containing ancestors of node 2
        
        // return null if either node is null       
        if (node1 == null | node2 == null){
           return null;
        }
        if (node1.getLabel().equals(node2.getLabel())){
           return node1;
        }
               
        // populates node1Ancestor with all of node1's ancestors
        while (node1.getParent() != null){
           node1Ancestors.add(node1.getParent());
           node1 = node1.getParent();
        }
        
        node2Ancestors.add(node2.getParent());
        node2 = node2.getParent();                   
        int i = 0;
        
        // check to see if each closest ancestor of node2 is also an ancestor of node1
        while (!node1Ancestors.contains(node2Ancestors.get(i))){
           node2Ancestors.add(node2.getParent());
           node2 = node2.getParent();
           i++;   
                
        }
               
        return node2Ancestors.get(i);
    }
}
