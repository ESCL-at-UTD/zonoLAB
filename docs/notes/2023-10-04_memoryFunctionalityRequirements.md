# zonoLAB Memory Functional Requirements/Disscussion
(Meeting Brief: 2023-10-04)

## Memory Functionality Implementation Requirements
### Mandatory Base Changes
- None (did make a few changes that can be removed once past prototyping)

### New Functions/Methods
- Modified Set Operations (account for shared factors)
  - Minkowski Sum 
  - Cartesian Product
  - Generalized Intersection
  - (Union, Compliment, Projection, PontryDiff may have yet to be defined modified versions)
- Fancy/Additional Set Operations
  - Labeled Intersection (rows labeled)
  - AddInfo (different name?... cartProd or intersection)

### Likely Base Changes
- Modifying abstract class to account for memory functionality and/or additional memory-based objects
  - Such as adding the additional switch/case for seperate classes and/or adjusting constructors of each class

### Proposed Usage Syntax
- Overload standard definitions with (...,'labels',{'x','y'})
  - change constructors to use arguments block or similar...
- Define separate memory maintaining set operations 
  - either defaulted to via a setting or as separate methods
- Indexing/Projection
  - define `subsref` and `subsasign`
  - setup projection method (class specific or abstract-level) to accept labels for dims

## Proposed Underlying setup
### Two potential structure approches:
1. Full integration
   - Integrate a set of options at the abstract level (and individual constructors) to enable various (not all yet determined) functionalities
   - Restructure the toolbox to have 
     - either define specific methods that account for memory 
     - or have a property that decides whether to incorporate memory into set operations
2. Seperate/additional classes
   - less scalable/well-integrated but potentially simpler (demonstration is this way)
   - Could either be a separate abstract with sub versions or just have the memory added to the individual -M 

### Data storage
- Book-keeping
  - Currently using built-in table functionality for book keeping of all labels (it's both simpler and more complicated)... may have computational limitations/inefficiencies yet to be compared
  - Direct book-keeping of indicies with associated labels dictionaries/arrays used to store said data

- Modification:
  - Instead of storing the Gc and Gb redundently then hidden, you could define just the G w/ an auxilory vset array defining cont or discrete (or something like that)... Gc and Gb could then just be dependent properties



