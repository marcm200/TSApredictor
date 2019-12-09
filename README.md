# TSApredictor
Fast routine to predict a refinement level at which interior cells in the Figueiredo-algorithm are detectable.

Extension of the articles:
```
"Images of Julia sets that you can trust"
by Luiz-Henrique de Figueiredo, Diego Nehab, Jofge Stolfi, Joao Batista Oliveira-
from 2013

"Rigorous bounds for polynomial Julia sets"
by Luiz-Henrique de Figueiredo, Diego Nehab, Jofge Stolfi, Joao Batista Oliveira-
```

### Disclaimer

Although I tried my best to ensure correctness of the implementation and therefore the results,
the code comes with no warranty.

## README is organized as:

1. Quick start
2. Background information on the algorithm
3. Command-line parameters
4. Limitations
5. Contact and links


## (1) Quick start

Compiling the source code after commenting in and out the desired datatype (#defines _DOUBLE, _LONGDOUBLE, _QUADMATH) with
a suitable C++-compiler, best with optimizations on. Assuming the executable is named `TSApredictor_d` (for double; _ld, _qd equivalently).

The output is written in the text file `TSApredictor.log`.

The fastest example: The basilica
<br>`TSApredictor_d func=z2c c=-1,0`

A quartic polynomial with long double datatype necessary:
<br>`TSApredictor_ld func=z4azc level=10,15 encw=256 c=0,-0.171875 a=1.375,0`
<br>(Lowering the ENCW value to 128 or 64 detects fewer and fewer cycles. See (3)).

## (2) Background

For background on the cell mapping/interval arithmetics approach by Figueiredo et. al,
please see my julia-tsa-core or juliatsa3d projects.

Computing sets at a high refinement level (18 or above) or with float128 is computationally very
expensive, so a quick check whether interior cells of a specific cycle are detectable at a given level 
would be beneficial.

The TSApredictor works in 2 phases: 1) a numerical one, followed by 2) an interval arithmetics one.

Phase 1: The algorithms finds the critical points of the given formula/parameter combination by classic Newton
iterations with starting points on a rectangle circumferencing all roots in 3 times the Lagrangian estimate
(modified from Hubbard, Schleicher, Sutherland. How to find all roots of complex polynomials, 2001).

Afterwards the critical orbits are computed to detect all attracting cycles.

Phase 2: For every cycle a rectangle is placed around every periodic point. The
IA algorithmic phase detects whether a complex interval therein has a path to outside any neighbourhood.
If there is no such path, all orbits of that interval are bounded and part of the cycle. Then this cycle
can be detected in the current resolution. Note, that this phase needs to be done with a datatype providing
sufficient bits of precision to handle the intermediate results (see limitations).

Note also, that the size of the rectangle around the periodic points influences whether a cycle can be detected or not (see limitations).


## (3) Command-line parameters (case-insensitive, order not relevant)

`FUNC=string` (if not provided, standard value is Z2C)
<br>The desired polynomial to use. Implemented are in general `zNazc` for the polynomials z^N+A*z+c, where N is 2..6.
Symbols A and c are complex numbers. For performance reasons a special function Z2C for z^2+c is also implemented which uses a more optimized version of computing the bounding box (time benefit not measured).

'C=double1,double2`
<br>Sets the seed value: double1 as real part, double2 as imaginary part.

`A=double1,double2`
<br>Sets the degree-1 coefficiant accordingly.

<b>Note</b> that real numbers given in the command line as parameters are always read in as C++ double, no matter 
what underlying data type is used in the binary. Values are internally transformed to 
floor(value*2^25)/2^25. 
This implies that c and A values have a fraction bit precision of at most 25.

`PERIODS=N,M` (standard: all cycles)
<br>Only cycles with period length between N and M (both included) are analyzed.

`ÈNCW=N` (standard value: 128, minimum: 32)
<br>The number of pixels the enclosement of each periodic point is enlarged by. A negative value is interpreted as: All cells contained
in the overall enclosement of the periodic points of a cycle (enlarged by |N|) are analyzed (computationally very expensive, i.e. one
tries to compute all immediate basins).

`LEVEL=N,M` (standard 10,24)
<br>Only refinement levels between N and M are checked. 


## (4) Limitations

- The software computes a value R being a power of two so that the filled-in Julia set is definitely contained therein
(Figueiredo 2013, paragraph 7), so refinement levels are always to be read relative to that R. Different R / refinement level combinations
can have the same width a screen pixel is representing.
- A negative answer of the predictor does not rule out interior cells at that level. A different value of ENCW might change that.
- A positive answer means that the julia-tsa-core finds black at latest at that level, but it might do so earlier (a different ENCW value might show that).
- If two cycles share a significant overlap in their rectangles enclosing the union of all the cycle's immediate basins, the predictor will output the same result for both cycles: the faster detectable cycle will dominate. If all cycles were analyzed (PERIOD command line parameter omitted), the code checks for an overlap of such two regions and prints a notification.
- Level 31 is the current maximum supported level for prediction analysis.
- The datatype used for phase 2 needs to provide enough bits to handle all intermediate results. Bounding box calculcations are
done on an expanded form of the real and imaginary component function of the polynomial. 

The following is recommended if the filled-in set is contained
in the 4-square, i.e. |x|<=4, |y|<=4 and |c|<=2, |A|<=2, fractional precision of c, A is 25 (see above):
<br><b>z2c</b> double till level 24
<br><b>z3azc</b> double till 17, long double 18-20
<br><b>z4azc</b> double till 13, long double 14-15
<br><b>z5azc,z6azc</b> float128


## (5) Contact

Please direct any comments to:

marcm200@freenet.de

forum: https://fractalforums.org/fractal-mathematics-and-new-theories/28/julia-sets-true-shape-and-escape-time/2725

Marc Meidlinger, December 2019

