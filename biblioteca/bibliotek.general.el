;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;  All files in this directory are copyright 1997, 1998    ;;;;;;;;
;;;;;;  by Rafael D. Sorkin.     All rights reserved.           ;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


                          ; Time-stamp: < Nov 05 1998 03:05:47 (13889 23643) >

;: This is { bibliotek.general }, containing functions of general utility.

;: Notes 
;
; 1. Most macros are in the separate file { bibliotek.macros }
;

;: Index of functions defined

  ; iff          (the usual predicate of that name)
  ;
  ; firstn       (initial segment of a list)
  ;
  ; omitnth      (a list with the nth element omitted)
  ;
  ; merge-lists  (merges two ordered lists of numbers)
  ;
  ; sort-list   (essentially redundant since `sort' is built in, but this
  ;                  one is non-destructive.)
  ;
  ;     (The menge family remove duplicate elements from a list.
  ;      Partly redundant because `remove-duplicates' is
  ;      built in but latter upsets order and is slow)
  ;
  ; shuffle       (randomly permutes a list) 
  ;
  ; menge   
  ; n-menge   
  ; n-menge%    
  ; menge% (menge%-alt)
  ;
  ; union-of    (forms union of two arbitrary lists-as-sets)
  ;             (cl-19 has `union' as a built-in function, hence this name)
  ;
  ; less/n-less (set difference, relative complement)
  ;
  ; intersect   (returns intersection of two sets, using `eq')
  ;             (CL has `intersection'  built-in)
  ;
  ; symmetric-difference (CL has `set-exclusive-or' built-in)
  ;
  ; equal-as-sets
  ;
  ; union-m      (forms union of a list of sets)
  ;              (the `m' stands for `members' or `menge' or `multiple')
  ;
  ; meetp       (true if two subsets meet each other)
  ; disjointp   (its opposite)
  ;
  ; fibers      (forms equvalence classes)
  ;
  ; memberp     (essentially same as `memq')
  ;
  ; member-member-p
  ;
  ; singleton-p 
  ;
  ; relate
  ;
  ; leaves     (returns the set of "leaves" of the "tree" of a list)
  ;
  ; multiplicity (finds number of times an element occurs in a multiset (list))
  ;              (duplicates CL function `count')
  ;
  ; sum-m == sum           (sums the elements of a list)
  ; product-m == product   (multiplies the elements of a list)
  ;
  ; !       (factorial function for integer argument (for TCL; elisp version is
  ;          in bib.float))
  ; log_2 
  ; ratio*  (just an alias for `/' in TCL, but converts to float in elisp)
  ; square
  ; sgn
  ; reciprocal 
  ;
  ; sup (alias max-m)
  ; inf (alias inf-m)
  ;
  ; median
  ; mean 
  ; sigma (alias std-dev standard-deviation)
  ; std-err-in-mean
  ; stats
  ; stats-xy
  ; correlator
  ; combine-stats (largely useless)
  ;
  ; random-draw
  ;
  ; get-symbolic-name  (finds the ``name'' of its argument)
  ;
  ; sharep, set-of-cells
  ;
  ; defconst-f    (like defconst, but a function not a macro)
  ;
  ; plist-to-alist  
  ; alist-to-plist  
  ;
  ; invert-cons-cell


;: The functions themselves

;=======================================================================
(in-package 'user)                     ; merely to make SCL happy
;=======================================================================

(defun iff (x y) "True when x and y are either both nil or both not"
  (or
   (and x y)
   (not (or x y))))


(defun firstn  (n A) 
  " Returns a (new) list comprising the first n elements of A."
  (reverse (nthcdr (- (length A) n) (reverse A))))  
 ;
 ; This function is almost useless, because you
 ; can get the same sublist using (subseq A 0 n) (which seems to copy too).
 ;
 ; We get a new list because `reverse' is used.
 ;
 ; Localization is unnecessary, because nothing is evaluated inside it other
 ; than the parameters n and A themselves.


(defun omitnth (n L) 
  " Returns a new list omitting element number n from L"
  (append (firstn n L) (nthcdr (1+ n) L) nil))
  ;
  ; Comments
  ;
  ; This has not proved very useful
  ;
  ; The first element is numbered `0'.
  ;
  ; The list returned is "new": it shares no cons-cells with L (which would
  ; not have been the case had `nil' been omitted from the end).
  ;
  ; Alternate plan: copy the list first, then splice out by resetting the
  ; relevant cdr-pointer.  


(deff merge-lists  (A B) 
  "\
 Merges two NUMERICAL lists which are ASSUMED to be in increasing order.
  "    
  (varbind
   x (car A)
   y (car B))   
  (cond 
   ( (null A) B )			; return B if A is empty
   ( (null B) A )			; or vice versa
   ( (< x y)  (cons x (merge-lists (cdr A) B)))
   (    t     (cons y (merge-lists (cdr B) A)))))  
 ;
 ; Comments
 ;
 ; 1. If extend to more general lists than just numerical, remember to
 ;    localize the internal variables.
 ;
 ; 2. In general the merged list will share cons-cells with A and B.
 ;
 ; 3. Has not found much use yet


(deff sort-list (A)     
  "\
 Sorts a list of NUMBERS into increasing order.  
 The built-in function `sort' is much faster (and destructive)."
  (cond
   ((o null cdr A) A)			   ; case of length 0 or 1
   (t					   ; general case
    (varbind N1 (truncate (/ (length A) 2)))
    (merge-lists 
     (sort-list (firstn  N1 A))
     (sort-list (nthcdr  N1 A))))))
   ;
   ; Comments 
   ;
   ;  The built in function `sort' is much faster but destructive.
   ;  To use it instead, just copy arg first, if desire to protect it. 
   ;
   ;  In TCL `truncate' is really needed; in elisp its use is only for style,
   ;  since integer division automatically truncates (unfortunately). 
   ;
   ; In case the list is of length 1 it is not copied.
   ;
   ; Could speed this up by not firstn and nthcdr, which do a lot of copying,
   ; instead could really break the list in the middle.
   ;
   ; Not localized since intended only for numbers


(defun shuffle (L)
  "Not so quick, but simple.  Randomly permutes elts of the list L"
  (if 
      (null L) nil  ; else
    (let
        ((n (Random (length L))))
      (cons
       (nth n L)
       (shuffle (omitnth n L))))))
        ;
  ; IMPROVEMENT TO MAKE: first generate a random perm using the factorial
  ; trick, then just apply it.  This is only one call to random, rather than
  ; many.   Also working with vectors  would be faster, so could make L a
  ; vector temporarily if it isnt already. 



;--------------------------------------------------------------------

;:: comments on menge family
;
; need menge or n-menge for general lists, could keep only latter and redefine
; former in terms of it (the built-in remove-duplicates is too slow in elisp)
;
; For sorted numerical lists, can probably keep only n-menge% (though it is
; maybe 1/3 slower than menge% and menge%-alt when there are few duplicates
;
; Thus can get rid of: menge%-alt menge%
 
(deff menge (L)  
  "\
 Deletes duplicate elements from a list of objects (non-destructively).
 Comparison is done with `eq' NOT with `equal'."
  (cond
   ((null L) ())
   (t
    (varbind tail-menge (o menge cdr L))
    (if 
	(memq (car L) tail-menge)
	tail-menge
      (cons (car L) tail-menge)))))
 
(defun n-menge (L) 
  "\
  Destructive version of `menge', comparison done with `eq'.
  BEWARE: not all equal integers are `eq' "
  (if (null L) ()
    (cons (car L) (n-menge (delq (car L) (cdr L))))))

 ; Notes on `menge' and `n-menge'
 ;
 ; These are general forms of menge, which can operate with any type of object.
 ;
 ; `menge' is non-destructive (since `cons' is), while `n-menge' is
 ; destructive, but intended to be faster.
 ;
 ; `menge' is "stable" in the sense that it preserves order, but from the end,
 ;  NOT the beginning.  In contrast `n-menge' preserves order from the
 ;  beginning, assuming that `delq' does not upset order.
 ;
 ; For lists of numbers, or for poset-elements e0 e1 etc, it would be much
 ; faster to sort them first so that only neighbors would need to be compared.
 ; But even faster is probably to use `remove-duplicates'
 ;  However latter also alters order.
 ;
 ; It appears that `menge(L)' will not share con-cells with L.

;--------------------------------------------------------------------
 
; In following, the name `menge#' is not used, because elisp can't read it
; properly, we use { menge% } instead

; ; ; Alternate version (originally was main version)
;
(deff menge% (S) 
  "\
  A version of `menge' intended for SORTED NUMERICAL lists.  
  (More generally, it will--stably and non-destructively--remove adjacent
  duplicates from any numerical list.)
  Comparison done with `='
  "   
  (cond 
   ((null S) nil)
   (t
    (varbind x (car S))
    (while (and S (= x (car S))) (setq S (cdr S)))
    (cons x (menge% S)))))

(quote
 
(defun menge%-alt (S) 
  "\
  A version of `menge' intended for SORTED NUMERICAL lists.  
  (More generally, it will--stably and non-destructively--remove adjacent
  duplicates from any numerical list.)
  Comparison done with `='
  "
  (cond
   ((o null cdr S)  S)
   ((= (car S) (cadr S)) (o menge%-alt cdr S))
   (t (cons (car S) (menge%-alt (cdr S))))))               )
  ;
  ; This code is most perspicuous and of similar speed
 
 
(deff n-menge% (L) 
  "\
 A DESTRUCTIVE member the `menge%' family.  It  removes duplicate elements
 from a SORTED NUMERICAL list (or adjacent duplicates from any numerical list)
 using only the storage for the list itself.
  "    
  (varbind cons L)
  (while (cdr cons)
    (if (= (car cons) (cadr cons)) 
	(setf (cdr cons) (cddr cons))
      (setq cons (cdr cons))))
  L)
 ;
 ; This "zipper" version uses no recursion.  It is overall much faster than
 ; menge% (except where the list is mostly repetition).
 ; Both versions are much faster than plain `menge' of course.
 ; How does it compare to remove-duplicates?

;----------------------------------------------------------------

(deff union-of (A B)
  "\
  Accepts a pair of ``sets'' (lists bila duplications) and returns their 
  union (as an entirely new set of cons cells).  
  Comparison is done with `eq'. 
  The efect is to delete from A all elements shared with B and then form
  (append A B).
  For collated sets, union% or union%% should be much faster.
  (This was previously called 'union-alt')  \
  "
  (cond
   ((null A)         (copy-list B))
   ((memq (car A) B) (union-of (cdr A) B))
   (t                (cons (car A) (union-of (cdr A) B)))))
   ;
   ; Seems to be about twice as fast as union-alt, as might expect.
   ;
   ; To get a nicer order, do this:
   ;
   ;   (reverse (union-of (reverse B) (reverse A)))

(deff union-m (L) 
  "\
 Accepts a single argument, which should be a list of sets, 
 and returns the union of those sets (as an entirely new set). 
 The ``m'' in the name stands for ``multiple'' or ``menge''.   \
  "
  (varbind 
   Y nil
   length-L (length L))
  (&bind-too j)
  (kwa 
   j from 0 upto length-L
   (setq Y (union-of Y (nth j L))))
  Y)
  ;
  ;  Could also use `reduce' to do this more compactly, but it does more
  ;  rearrangement (see tests)
  ;  (reduce 'union L)
  ;
  ; "Union" might be a nicer name, but you can't use it with case-folding 
  ;  lisps (which unfortunately, are most of them) 


(defun less (A B) 
  "\
 Returns A\\B for sets A and B.  The order of the elements is not affected,
 and the original lists A and B are unchanged.  Comparison done with `eq'.
  " 
  (if (null A) nil
    (if 
        (memq (car A) B)  (less (cdr A) B)
      (cons (car A)  (less (cdr A) B)))))
        
(defun n-less (A B)
  " Destructively returns A\\B for sets A and B"
  (loop for bw in B do (setq A (delq bw A)))
  A)
 ;
 ; This is of same speed, but could be useful since not recursive
 ; aside from that, could be deleted
  

(deff intersect (A B)  
  "\
 Returns the intersection of two ``sets as lists''. 
 Comparison done with `eq'  (not `equal') "
  (&bind-too x)
  (varbind C nil)
  (while A
    (setq 
     x (car A)
     A (cdr A))
    (if 
	(memq x B)
	(setq C (cons x C))))
  (reverse C))
 
(defun intersect-alt (A B)  
  " Returns the intersection of two ``sets as lists''."
  (less A (less A B)))
 ;
 ; This more perspicuous version also works.
 ; It may be slowed down by copying lists and recursion, and in fact it
 ; can sometimes be much slower, though often not.
 ; seems of comparable speed if smaller set comes first (as it clearly should,
 ; so add this feature if use it!)


(defun symmetric-difference (A B)
  " Self-explanatory set function.  Comparison done however `less' does it."
  (append (less A B) (less B A)))

(defun equal-as-sets (A B) 
  "\
 Self-explanatory set function.  
 Comparison of elements done however `symmetric-difference' does it."
  (not (symmetric-difference A B)))

(defun meetp (A B)  
  "Tests whether sets A and B overlap.  Comparison with `eq' "
  (cond 
   ((null A) nil)
   ((memq (car A) B) t)
   (t (meetp (cdr A) B))))

(defun disjointp (A B) "True if A and B are disjoint."  
  (not (meetp A B)))


(deff fibers (S &key ((:test equiv) #'eq))
  "\
 The equivalence classes into which the (nonempty) set S is fibered by the
 equivalence relation specified by the keyword :test (defaults to `eq').  
 Order of elements is preserved.
  "
  (&localize S equiv fibres x f)
  (if (null S) (error "Can't fiber an empty set"))
  (varbind fibres nil)
 ;-------------------------
 ; for each element of S...
 ;-------------------------
  (loop 
   for x in S do
   (cond
   ;----------------------------------------------------------------
   ; if x belongs in an already existing fiber, add it to that fiber
   ;----------------------------------------------------------------
    ((loop for f on fibres
      if (funcall equiv x (caar f))
      do 
      (progn
	(push x (car f))
	(return t))) nil)
   ;----------------------------------
   ; otherwise make a new fiber for x
   ;----------------------------------
    (t (push (list x) fibres))))
 ;------------------------------------
 ; reorder the elements and the fibers
 ;------------------------------------
  (reverse (mapcar #'reverse fibres)))


(defun memberp (x S) 
  "\
 Returns `t' if `x' is an element of `S', otherwise nil.  
 Comparison done with `eq', not `equal'.  
 Exactly same as `memq' except that (in positive case)
 we returns `t' rather than the cons cell where found `x'. 
  "
  (if (memq x S) t nil))


(defun member-member-p (x y) 
  "\
 Returns t if x is an element of an element of y, nil otherwise.
 Membership assessed by `member', which tests equality with `equal'.
  "
  (cond
   ((null y) nil)
   ((atom (car y)) (member-member-p x (cdr y)))
   ((member x (car y)) t)
   ((member-member-p x (cdr y)))))
 ;
 ; Notes: 
 ; This version is not paralellizable, as we might want for big sets.
 ;
 ; No need for localization, it would seem


(defun singleton-p (S) "True if S is a singleton" 
  (and 
   (consp S)
   (null (cdr S))))


(deff relate (R x &key (test 'equal)) 
  "\
 Apply the RELATION R to the element x.  That is, we return the set 
 (NOT the multiset) of all y such that (x . y) is an element of R.
 Arguments:
               R  x  &key test [default: equal] 
  "
  (setq R (Member x R :key (function car) :test test))
  (when R (adjoin (cdar R) (relate (cdr R) x))))


(defun leaves (cc) 
  "\
  The set of leaves of the so-called ``tree'' of the given cons-cell.
  NB: The symbol `nil' can never occurs as a leaf, since it is identified (in
  elisp and TCL) with the empty list `()'.
  " 
  (cond
   ((null cc) nil)
   ((atom cc) (list cc))                ; just a trick for the recursion
   (t 
    (union 
     (leaves (car cc))
     (leaves (cdr cc))))))
  ;
  ;; possisbly using mapcar would be more efficient here


(deff multiplicity (x M) 
  "\
 The multiplicity of x in the list or multiset M (equality defined by `eq'). 
 This essentially duplicates the CL function `count'.
  "
  (varbind tail (memq x M))
  (if (null tail) 0
    (1+ (multiplicity x (cdr tail)))))


(defun sum-m (X) "adds together numbers in a list" (apply '+ X))
(defalias 'sum 'sum-m)


(when *TCL* (unless (fboundp '!)
(defun ! (n) "Factorial function for integers only"
  (cond 
   ((not (integerp n)) (error "This one is for integers only"))
   ((< n 0) Infinity%)
   ((= 0 n) 1)
   (t (* n (! (1- n))))))    ))
  ;
  ; The elisp version is in bib.float rather than here, because it converts
  ; its argument to float. 
  ;
  ; Since `!' is built into some TCL implementations, we don't define it
  ; without checking first that it really is absent.  (It seems not to be
  ; present in CLtL1 or GCL, but it definitely is in Clisp) 


(defun log_2 (x) " Base 2 logarithm" (log x 2))  

(defun square (x) (* x x))

(defun sgn (x) "The signum function, returning -1 0 1"
  (cond
   ((> x 0)  1)
   ((< x 0) -1)
   (t 0)))

;--------------------------------------------------------

;:: Some arithmetic functions that differ between elisp and TCL 
;
; Define some arithmetic fcns differently in elisp, to circumvent "large
; integer nuisance" (i.e. that integer overflow is not detected!), also to
; avoid problem that integer division is truncated in elisp (like fortran).
; To avoid these problems in elisp, we convert to floats where it seems
; appropriate.
; Also note that elisp lacks rationals, unlike TCL.
;-------------------------------------------------------------
; Note concerning the following fcns.
;  The fcn ratio* converts its args to float if we're 
;  in elisp (but if we're in TCL, it is just aliased to { / })
;--------------------------------------------------------------
(cond
 ;
(*elisp*
 ;
(defun ratio* (x y) 
  "\
 Floats its arguments and then takes their ratio.  This is necessary with elisp
 since division of integers truncates the quotient to an integer!
  " 
  (/ (float x) (float y)))
 ;
(defun reciprocal (x) 
  " First convert to float, since elisp truncates integer division" 
  (/ 1 (float x)))
 ;
(defun product-m (X) 
  "\
 The argument should be a list of numbers.  
 We convert them to float (sigh) and return their product." 
  (apply '* (mapcar 'float X)))
 ; 
)
 ;
(*TCL*
 ; 
(defalias 'ratio* '/)
(defun product-m (X) "multiplies together numbers in a list" (apply '* X)) 
(defun reciprocal (x) (/ 1 x)))
 ; 
(t 
 (error "Both `*elisp*' and `*TCL*' are nil!")))
 
(defalias 'product 'product-m)  


;;: Separate definitions of `sup' and `inf' for elisp and GCL
;;; Must add appropriate defs for other types (eg Clisp SCL)
;
(case *lisp-type*
;
((elisp)
 ; 
(defsubst sup (seq) 
  "\
 Maximum of the numbers in a sequence.  This is just like `max' but 
 instead of multiple arguments it takes a single sequence."
  (apply 'max (coerce seq 'list)))
 ;
(defsubst inf (seq) 
  "\
 Minimum of the numbers in a sequence.  This is just like `min' but 
 instead of multiple arguments it takes a single sequence."
  (apply 'min (coerce seq 'list)))  )
 ;
((gcl)
 ; 
(defun sup (seq)
  "\
 Maximum of the numbers in a sequence.  This is just like `max' but 
 instead of multiple arguments it takes a single sequence."
  (loop for x in (coerce seq 'list) maximize x))
 ;
(defun inf (seq)
  " Minimum of the numbers in a sequence. "
  (loop for x in (coerce seq 'list) minimize x))  )  )
 
(defalias 'max-m 'sup)
(defalias 'min-m 'inf)

;---------------------------------------------------

(deff median (L) 
  "\
 Median of a list of numbers: if N is odd we take the middle value, if N is
 even we take the mean of the two middle values.  \
  "
  (varbind N (length L))
  (setq 
   L (sort (copy-list L) '<))
 ;-------------------------------------
 ; handle odd and even cases separately
 ;-------------------------------------
  (if (oddp  N) (nth (/ (1- N) 2) L)
    (mean 
     (list
      (nth (/ N 2)      L)
      (nth (1- (/ N 2)) L)))))


(cond
;-------------------------
(*TCL*
;/case of true common lisp
;-------------------------
(defun mean (L)  " The arithmetic mean of a list of numbers"
  (ratio* 
   (sum-m L)
   (length L))))
;----------------------------------------
(*elisp*
;/case of elisp: first convert to floats
;----------------------------------------
(defun mean (L)  " The arithmetic mean of a list of numbers"
     (setq L (mapcar 'float L))
     (ratio* 
      (sum-m L)
      (length L))))) 

(deff sigma (L) " The ``sample standard deviation'' a list of numbers"
  (varbind 
    N (length L)
    mu (mean L)
    sigma (mean (image of (square (- $ mu)) on L))
    r (ratio* N (1- N)))
  (sqrt (* r sigma)))
  ;
(defalias 'std-dev 'sigma)
(defalias 'standard-deviation 'sigma)

(defun std-err-in-mean (L)   
  "\
 The ``std error in the mean'' of a list of numbers, defined as the
 sample standard deviation divided by sqrt N.  \
  "
  (ratio* (std-dev L) (sqrt (length L))))

(deff stats (x)
  "\
 The arg should be a list of numbers (regarded as independent samplings
 of some random variable x).  We return a plist containing a grab bag of
 statistics derived from the list.  The main plist indicators
 are
           N  mean  variance  sigma  U-mean  U-var  U-sigma ,
 most of
 which estimate parameters of the underlying probability distribution
 of x.  For example `variance' estimates var(x):=Cxx as the mean-square
 variation in the list, weighted by the standard factor of N/N-1.  
 In addition to these statistics, we return in raw form the second third
 and fourth sample moments about the sample mean as `P2' `P3' and `P4'.
 They can be useful for further processing.

 The U-prefix is short for ``uncertainty in'' (or ``standard error in'', to
 use another term). For example  `U-mean' means ``stderr-in-the-mean''.
 
 Some explanation of U-var and U-sigma:
 The statistic `U-var' is an estimate of the amount by which the estimated
 variance deviates from the true variance C(x,x).  This estimate is slightly 
 on the conservative side for finite sample size N, but ``asymptotically just''
 as N--> infty.  More precisely, let var-est be our estimated variance and
 var-true the true variance, C(x,x), of x.  We employ an estimator 
 of U := (var-est - var-true)^2  whose expectation value is asymptotically
 equal to that of U (thus equal to the the mean square fluctuation of var-est
 about var-true), and greater than it for finite N.  This estimator is

                           _________________
               N^2                  ____  2
          --------------   ( x^x^ - x^x^ ) 
          (N-1) (N-2)^2

 where x-hat = x - {bar x} and bar denotes sample mean.  Our value
 for U-var is then the square root of this formula.
 This approach should yield reasonable results in general, though it clearly
 can be ``way off'' for certain cooked up data (eg one can easily concoct 
 data for which U-var = 0).  \
  "
 ;-------------------------------------------
 ;/ count data points and signal error if < 3  
 ;-------------------------------------------
  (varbind N_points (length x))		; number of data points
  (if 
    (< N_points 3) 
      (error "Can't compute all the statistics with fewer than 3 points"))
 ;--------------------------------------------------------------------------
 ;/ define the local function `minus(X y)' to subtract y from every elt of X
 ;--------------------------------------------------------------------------
  (fbind minus (X y) (image on X of (- $ y))) 
 ;-----------------------
 ;/ compute the statistics
 ;-----------------------
  (varbind
   N      (float N_points)	       ; integer division nuisance in elisp!
				       ; (also preempts rational arith in TCL) 
   N/N-1  (atl  N / (N - 1))           ; ubiquitous correction factor 
   fac2   (/
	   (square  N)
	   (atl 
	    (N - 1)*(N - 2)*(N - 2))) ; another correction factor 
   ;
   xbar   (mean x)                    ; sample mean of x, bar(x)
   ;
   x^     (minus x xbar)              ; x-hat = x - xbar
   Pxx    (mean (mapcar #'square x^)) ; 2nd sample moment about mean
   Pxxx   (mean (image 
		   on x^ 
		   of (* $ $ $)))     ; third sample moment about mean
   Pxxxx  (mean 
	   (image 
	    on x^
	    of (o square square $)))  ; fourth sample moment about mean
   ;
   Cxx    (* N/N-1 Pxx)		      ; estimates (variance x)
   sigma  (sqrt Cxx)                  ; estimates (std-deviation x)
   U-var  (sqrt                       ; uncertainty of (variance x) estimate
            (*
              fac2
	      (- Pxxxx (square Pxx))))
   U-sigma (max
             (- 
               (sqrt(+ Cxx U-var))
                sigma)
             (-
               sigma
               (sqrt(max 0 (- Cxx U-var))))))
               ;     ^^^^^ protection for sqrt
 ;----------------------
 ;/ return the statistics
 ;----------------------
  (list
   'N         N_points			; return an integer, not a float
   'mean      xbar
   'U-mean    (sqrt (/ Cxx N))
   'sigma     sigma
   'U-sigma   U-sigma
   'variance  Cxx
   'U-var     U-var
   'P2        Pxx
   'P3        Pxxx
   'P4        Pxxxx))  
   ;
   ; Notes
   ;
   ; It may be that the sqrt above doesn't really need its "protection",
   ; because its arg can never really be negative.  Modulo arithmetic, the
   ; condition is that <y^4> <= <y^2><y^2>(N^2-3N+3)/(N-1) for any vble y of
   ; mean zero and N>2), so for big N it's almost certain.
   ;
   ; We compute the third moment p(xxx) even though it isn't needed for any of
   ; the other statistics, the reason being that we do need it in order to
   ; *combine* two samples.


(deff stats-xy (x y)
  "\
 Input is a pair of numerical lists regarded as corresponding samples 
 of the random variables x and y.  We return a plist beginning with

                    Cxy   <value>  
                    UCxy  <value>  

 where Cxy estimates the correlator

                     C(x,y):=<xy>-<x><y> ,  
 and UCxy 
 estimates the uncertainty in our estimate of C(x,y) (and may be called 
 the ``standard error of Cxy''.   Additionally, the plist contains N and the
 raw sample moments (about the sample means), Pxy, Pxxy, Pxyy, Pxxyy.  
 The moments Pxxxy and Pxyyy are NOT given.
 The moments Pxx, Pxxx, Pxxxx and (x-->y) are not given either, 
 but they can be obtained from `stats'. 

 For the correlator C(x,y) we use the standard unbiased  
 estimate 
                        N      __   _ _
                     ------- ( xy - x y )              (C-est)
                      N - 1
 where bar 
 denotes sample mean.  To estimate the uncertainty UC in C-est we first
 estimate {UC}^2 := (variance C-est).  In order that this latter estimate be
 identically non-negative, we use a formula which is slightly biased in a
 conservative direction: its expectation is slightly greater than the true
 value of (variance C-est).  The bias is only O(1/N) however.  This formula is
 the following expression (where x^ := x - x-bar, and similarly for y^): 

                          _________________
               N^2                 ____  2
          --------------  ( x^y^ - x^y^ )
          (N-1) (N-2)^2

 Finally, to get UC itself, we just take sqrt of this expression.  \
  "
 ;--------------------------------------
 ;/ define the auxiliary function `hat'
 ;--------------------------------------
  (fbind hat (w) "deviation of w from its mean"
    (varbind wbar (mean w))
    (image on w of (- $ wbar)))
 ;-----------------------------------------------
 ;/ count the data points and signal error if < 3
 ;-----------------------------------------------
  (varbind 
   N_points (length x)
   N (float N_points))   ; to avoid integer division in elisp
  (if 
   (< N 3)
   (error 
    (concat
      "Can't compute all the requested correlation "
      "coefficients with < 3 points.")))
 ;------------------------------------
 ;/ check that x and y have same length
 ;------------------------------------
  (unless (= N (length y))
    (error "Can't compute correlator for lists of unequal length."))
 ;------------------------------------------------------------
 ;/ compute x^ y^ etc. and then the desired statistics
 ;------------------------------------------------------------
  (varbind
   N/N-1    (/ N (- N 1))
   fac2     (/
	      (square N)
	      (- N 1)
	      (square (- N 2)))
   x^       (hat x)
   y^       (hat y)
   x^y^     (map 'list (function *) x^ y^)
   Pxy      (mean x^y^)
   Pxxy     (mean (map 'list (function *) x^ x^y^))
   Pxyy     (mean (map 'list (function *) x^y^ y^))
   Pxxyy    (mean (map 'list (function *) x^y^ x^y^))

   Cxy      (*  N/N-1 Pxy)

   UCxy     (sqrt (* fac2 (- Pxxyy (square Pxy)))))

 ;-------------------------------------------------------------
 ;/ return  Cxy and its "std-err", as well as certain moments
 ;-------------------------------------------------------------
  (list
   'Cxy   Cxy
   'UCxy  UCxy
   'Pxy   Pxy
   'Pxxy  Pxxy
   'Pxyy  Pxyy
   'Pxxyy Pxxyy
   'N     N))

(defun correlator (x y) " The correlator of two lists of numbers"
  (getf (stats-xy x y) 'Cxy))

;; following now largely useless
;;
(Defun combine-stats
  (samples &key dispersion-measure   
           &aux N mu spread-type A B D spread)
  "\
 Takes a set of lists of the form (size mean spread) and returns a plist
 containing the corresponding statistics for the combined sample.   
 The entry `spread' is either a sample-std-dev or a std-error, as indicated
 by whether the keyword argument `:dispersion-measure' equals `sigma'or
 `std-err' (several synonyms for each are accepted and whichever is given is
 used in the return plist as well).  

 Arguments:  (samples &key dispersion-measure)
  "  
  (flet
   ((size   (sample) (float (nth 0 sample)))   ; float needed only for elisp
    (mu     (sample) (float (nth 1 sample)))
    (spread (sample) (float (nth 2 sample))))
   ;;
   ;; Find combined size
   (setq N (sum (mapcar 'size samples)))
   ;;
   ;; Compute combined sample-mean
   (setq 
    mu (loop for $ in samples sum (* (size $) (mu $)))
    mu (ratio* mu N))
   ;;
   ;; Find out which of sigma or std-err (= sigma/sqrt N) is required
   (setq 
    spread-type
    (case dispersion-measure
      ((sigma std-dev standard-deviation sample-std-dev) 'sigma)
      ((std-err stderr standard-error)                   'std-err)
      (otherwise (error "Sijapata valid spread-type"))))
   ;;
   ;; Compute sigma or std-err as the case may be
   (setq 
    A (loop for $ in samples sum (* (size $) (square (mu $))))
    A (atl  A - N * (square mu)))
   ;;
   (when (equal spread-type 'sigma)
     (setq 
      B (loop for $ in samples sum (atl ((size $)- 1) * (square (spread $))))
      D (atl N - 1)))
   ;;
   (when (equal spread-type 'std-err)
     (setq 
      B (loop for $ in samples sum 
              (atl (size $) * ((size $)- 1) * (square (spread $))))
      D (atl N * (N - 1))))
   ;;
   (setq spread (sqrt (ratio* (+ A B) D)))
   ;;
   ;; Return results
   (list 'size N 'mean mu dispersion-measure spread)))

 ;; Possible improvements
  ; 1. Make  it return both sigma AND std-err by default
  ; 2. Check for negative spread in input (no problem, but nonsensical)
  ; 3. Check for n < 0 ?


(deff random-draw (seq) "Select at random one element of a sequence"
  (elt seq (o Random length seq)))


(unless *gcl*      
(deff get-symbolic-name (obj &key ((:ignore ignore))) 
  "\
 Returns the ``name'' of its (actual) argument.  The optional second
 argument is a name it will overlook.
 WARNING: this has not been defined for gcl yet.
  "
  (&localize obj majina $ ignore)
  (&bind-too majina)
  (loop initially (setq majina nil)
        for $ being each symbol do
        (if (boundp $)
            (if (eq (symbol-value $) obj) 
                (push $ majina))))
  (setq majina (delq 'obj majina))      ; since this is temporarily a name too!
  (setq majina (delq ignore majina))
  (setq majina (delq 'ignore majina))
  (setq majina (delq 'edebug-previous-result majina))
    (cond
     ((singleton-p majina) (car majina))
     ((null majina) (error "I can't find a name for that object."))
     (t  
      (cond
       (*elisp* (error "That object has more than one name: %s" majina))
       (*TCL*   (error "That object has more than one name: ~s" majina)))))) )
 ;
 ;
 ; Let OBJ be the actual argument to this function.  It looks for a symbol
 ; whose ``value-pointer'' points to OBJ.  If it finds exactly one, it returns
 ; that symbol, else it signals an error.  (It has NOT been made sophisticated
 ; enough to work properly when the correct result would have been one of 
 ; its own internal symbols (whose names are derived from  "obj" "majina" or
 ; "$") or (??) 'edebug-previous-result). 
 ;
 ; It  seems to work, should perhaps extend it to handle its special internal
 ; symbols too, which it could do.  
 ;
 ; The second argument allows you to circumvent the problem
 ; that, when used inside a function it finds the dummy name as well.
 ;
 ;; TCL warning: we haven't distinguished local from global bindings here, one
 ;; could do this with another keyword.  (I think only the global name will be
 ;; found as it is)
 ;
 ; TCL warning: in interactive use, * ** *** are bound to the most recent
 ; printed results (which can confuse things, can easily fix of course).


(defun set-of-cells (L) 
  " Returns a list whose elts are the cells of the (nil-terminated) list L"
  (if (null L) nil (cons L (set-of-cells (cdr L)))))

(defun sharep (A B) " Do these two (undotted) lists share cons-cells?"
  (meetp (set-of-cells A) (set-of-cells B)))


(defun defconst-f (symbol &optional doc-str) 
  " A ``functional'' version of `defconstant' (i.e. args are evaluated)."
  (eval (list 'defconst symbol (list 'quote symbol) doc-str)))
 ;
 ; Here is the back-quoted version of the final form:
 ;
 ;      (eval `(defconst ,symbol ',symbol ,doc-str))
 ;
 ; It is important NOT to write this as a macro, because we want it
 ; to be evaluated at run-time, not compile time.  Doing so entails a clever
 ; use of eval, which compensates for the fact that `defconst' (being itself
 ; a macro) does not evaluate its arguments.
 ;
 ; This function might be more generally useful.


(defun plist-to-alist (P) " Converts a plist to an alist"
  (if (oddp (length P))
      (error 
       "Argument not a valid p-list because its length (= %d) is odd. "
       (length P)))
  (if (null P) nil
    (cons 
     (cons (car P) (cadr P))
     (plist-to-alist (cddr P)))))

(defun alist-to-plist (A) " Converts an alist to a plist"
  (if (null A) nil
    (append
     (list (caar A) (cdar A))
     (alist-to-plist (cdr A)))))


(defun invert-cons-cell (cell) 
  " (X . Y) --> (Y . X)    Apparently not in CLtL1 !"
  (cons (cdr cell) (car cell)))

