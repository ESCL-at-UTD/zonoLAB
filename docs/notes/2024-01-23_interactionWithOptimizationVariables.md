# zonoLAB - Interacting with Optimization Variables

## MPT3 Notes
**MPT3 Polyhedron does not allow for objects to be dependent on optimization variables.**
YALMIP is used interanlly a lot for opimizations, however the objects are not able to be constructed to depend on optimization variables

### Constructor Implementation:
- Only interaction w/ YALMIP objects is is to convert an LMI into a geometric object describing the feasible region
### mtimes() Implementation:
- MPT3 has multiple individual functions that are then selected between them within the mtimes() function
- If either input is a scaler, then it just scales
- If both are sets (Polyhedrons) then it is the Cartesian product
- If the left is a Matrix then it is an Affine (linear) Map: P.affineMap()
- If the right is a Matrix then it is an inverase affine map: P.invAffineMap()
Note: This ignores the mtimes() overloading issue as it would not allow for non-real matrices to be used as affine maps.

## Potential Implementation
### Interaction with optimization variables
- Proposed: Allow any variables/object-types of the correct size to serve as the individual matrices stored within any of the zonoLAB obetcts
- Alternatives: 
  - a separate class (optimZono?) could exist to allow for this additional functionality but without other operations (plot?, volume?) that wouldn't work if they are optimvars/sdpvars/sym
### mtimes() syntax
- It does not appear to be possible to get MATLAB to call the overloaded mtimes() function from the second object listed...
- Proposed new syntax:
  - define new functions specifying the operations that do not directly overload the operation
  - use the same setup/logic as the MPT3 interaction does that would allow for mtimes to be used for set-operations when called (scaling, affineMap, invAffineMap, cartProd?)
  - If using an optimvar or sdpvar then the direct operation could be called instead and not fail: Z.affineMap(A), Z.invAffineMap(A), etc.
- We can/can do this same thing for the other overloaded operations (and -> generalizedIntersection/labeledIntersection, plus -> minkSum, etc) and 
- Alternatives: 
  - Only do these type of operations with the special object

