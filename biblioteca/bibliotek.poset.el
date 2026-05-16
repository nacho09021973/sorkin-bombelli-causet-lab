;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;  All files in this directory are copyright 1997, 1998,   ;;;;;;;;
;;;;;;  1999 by Rafael D. Sorkin.  All rights reserved.         ;;;;;;;;
;;;;;;                                                          ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;                          Time-stamp: < Dec 16 1998 02:26:50 (13943 24762) >

;        ===========================================
;        Functions Designed Specifically for Causets
;        ===========================================

;: Some nomenclature 
;  =================
;
;     order = poset = partial order = acyclic transitive digraph
;     preorder = preposet = acyclic irreflexive relation
;     pre-pseudo-order = digraph = relation 
;     tonomorphic = isomorphic in the category of orders 
;     link = hasse-link = irreducible relation
;     chain = linearly ordered subset of an order
;     path = saturated chain (all its links are links of the ambient order)
;     past : will be taken wo the *irreflexive* convention
;
;     See Repository for full definitions of these terms and more 

;: Roster of functions put together so far (some still unchecked)
;
;  past  (returns the "absolute" past of an element)
;
;     setf method also provided for `past'
;
;     (The function `past' defined herein assumes that every poset element is
;     a symbol, and that its past is just its symbol-value.   This can be
;     overridden by redefining `past' and its @ setf-method, but you might
;     then have to recompile every function invoking `past'.)
;
;  past-of  (extensions to relative past, past of set, and inclusive past)
;
;  prec    (the fundamental order-relation as a predicate)
;
;  future    (of an element rel a subset)
;  future-of (of an element or a set rel a subset)
;
;  interval  (the exclusive (``open'') interval between two elements)
;
;  jana    (the set of immediate ancestors, "parents")
;
;  kesho   (children) (needed?) (NOT DEFINED HERE, but present in workspace)
;
;  maximal-of-order/maximal  (the maximal layer of an order)
;  n-max+                    (auxiliary for previous)
;  maximal-of-preorder       (does not assume transitivity, useless?)
;
;  minimal             (the set of minimal elements of a SUB-order)
;
;  connected-part/disconnected-part  (of a poset rel a subset)
;
;  compute-futures-for
;
;  compute-links-for
;
;  compute-levels-for
;
;  assign-levels-to-sorted-preorder (could just absorb this into previous)
;
;  level  (includes a defsetf for it)
;
;  count-chains-for/count-chains-of-poset-inayoitwa 
;
;  count-paths-for/count-paths-of-poset-inayoitwa 
;
;  count-spanning-trees-for/count-spanning-trees-of-poset-inayoitwa 
;
;  count-spanning-hasse-trees-for /
;    count-spanning-hasse-trees-of-poset-inayoitwa 
;
;  prepare-substrate   (makes a "blank" substrate)
;
;  create-relation (useful for making small relations by hand)
;  create-order    (useful for making small posets by hand) [being written]
;
;  sortedp == up-sorted-p    ("Is this digraph sorted upward?")
;  down-sorted-p
;
;  posort-upward  (sorts a preorder upward)
;
;  posort    (sorts a preorder upward or downward as specified, using
;             posort-upward) 
;
;  posort-order-upward
;
;  sort-pasts-downward  
;
;  transitive-p           (is this relation transitive?)
;  non-transitive-p       (is this relation NOT transitive?)
;
;  t-close-preorder       (transitively closes a preposet, "destructive")
;
;  t-close-preorder-alt   (alternate version, slower)
;
;  t-close-digraph        (allows for cycles)
;
;  n-delete-from-poset
;
;  copy-subposet/copy-poset
;
;  count-relations
;  count-relations-for  
;
;  count-links-for  
;

;=====================================================================

;: Preparations

(in-package 'user)              ; just to make TCL happy
(provide 'bibliotek-poset)	; GCL dislikes this at end for some reason

;:------------------------------------------------------------------------
;: Define past to be VALUE of symbol, unless past has been defined already
;:------------------------------------------------------------------------   

(unless (fboundp 'past)
 (defalias 'past 'symbol-value)       
 (defsetf past (x) (y) (list 'set x y)))
 ;
 ; Comments
 ;   Using `symbol-value' is not only faster than `eval'(nach CLtL1),
 ;   but it automatically detects when arg is not a symbol. 
 ; 
 ;   Would it be faster to use `defsubst' (but TCL lacks it) or perhaps just
 ;   to directly set function cell of `past' to be that of `symbol-value'?
 ;
 ;   Do we need to localize the defsetf here?
 ;   Seems not since it is a macro, whence by the time it is used the
 ;   temporary bindings of x and y will be gone.


;: And here are the functions themselves

(deff past-of (&key 
	        ((:elt x) "" elt-arg)
		((:set S) "" set-arg)
		((:rel R) "" relativize)
		((:inclusive inclusive) nil))
  "\
 The past of an element or a set, intersected with R when supplied.
 (The past of a set means the union of the pasts of its elements.)
 Keyword args  
              :elt  :set  :rel :inclusive {defaults to nil}
 Usage:
      (pastof :elt x :rel R)
      (pastof :set S :rel R)
      (pastof :elt x :inclusive t), etc.  \
  "
  (&localize x S R inclusive elt-arg set-arg relativize P)
  (&bind-too P)
  (cond 
   (elt-arg 
    (setq P (past x))
    (if inclusive (pushnew x P)))	; see note 3 below
   (set-arg 
    (setq P (union-m (mapcar 'past S)))
    (if inclusive 
        (setq P (union-of S P))))       ; see note 4 below
   (t 
    (error "Must supply either an element or a set to `pastof'")))
  (if relativize (intersect P R) P))
  ;
  ; Comments
  ;
  ; 1. Localization is crucial.
  ;
  ; 2. The seemingly redundant way keywords are put into the "lambda list" is
  ;    necessary, if you want to localize them.  Otherwise the keyword name
  ;    itself would get localized => you would have to refer to it by some
  ;    horrible name.
  ;
  ; 3. The use of `pushnew' rather than `push' in elt-arg case is only a
  ;    precaution, it will be needed if we ever switch over to the inclusive
  ;    convention for pasts. 
  ;
  ; 4. The referenced use of `union-of' rather than `union' is done to
  ;    preserve order, which is not guaranteed in TCL.  If order is not
  ;    important can switch to `union' if it is faster (which it aint in elisp)
  ;
  ; 5. For using pasts as numerically sorted lists, just replace union-of by
  ;    union%% (and same for union-m) see elsewhere for the resulting function,
  ;   `past-of%%'
  ;
(defalias 'pastof 'past-of)


(deff future (x R)
  "\
 Compute the future of x relative to R, ie the set {y in R | x prec y }"
  (&localize x y R)
  (loop for y in R if (prec x y) collect y))
;;;
;;; For now `future' uses the simplest possible method of direct computation.
;;  We do not look whether future has been computed and stored in plist.
;;  Such a more flexible version is being developed in { developing },
;;  but not clear vale la pena so far.


;::: a new future-of is being developed in developing
; 
(deff future-of (&key 
		((:elt x) "" elt-arg)
		((:set S) "" set-arg)
		((:rel R) "" relativized)
		((:inclusive inclusive) nil))
  "\
 The future of an element x or a set S relative to R.
 (The future of a set means the union of the futures of its elements.)
 (x or S should be within R, if not result is undefined.)
 Keyword args are  
                  :elt  :set  :rel :inclusive {defaults to nil}
 Usage:
           (future-of :elt x :rel R)
           (future-of :set S :rel R)
	   (future-of :elt x :inclusive t)  , etc. \
  "
  (unless relativized (error "Must supply R to `future-of'"))
  (&localize x S R elt-arg set-arg relativized inclusive F)
  (&bind-too F)
  (cond 
    (elt-arg 
      (setq F (future x R))
      (if inclusive (pushnew x F) F))
    (set-arg 
     (setq F (union-m (image on S of (future $ R))))
     (if inclusive 
	 (setq F (union S F)) F))
    (t 
     (error "Must supply either an element or a set to `future-of'"))))
     ; 
(defalias 'futureof 'future-of)
 ;
 ; See notes for `past-of'



(deff interval (x y &key ((:inclusive incl) nil)) 
  "\
 The interval with bottom x and top y.
 It is exclusive (``order open'') unless the keyword :inclusive is nonnil"
  (&localize x y incl)
  (if incl
      (append (list x y) (future x (past y)))
    (future x (past y))))
  ;
  ; Note: can't use  `(intersect (future x) (past y))' since we don't define
  ; absolute futures, currently

;-----------------------------------------------------------------

;;; A new version of `jana' is being developed in `developing'
; 
(deff jana (x) 
  "\
 The immediate past of a poset element: its ``parents''.  
 We retrieve it from x's plist if (x is a symbol and) it's there, 
 if not we compute it.  
 NB: We assume that the relation `past' is an *order* (ie we rely on 
     it being transitive). "
  (&localize x W)
  (varbind W "void")
  (cond
   ((and 
     (symbolp x)
     (not (equalp "void" (setf W (Get x 'jana "void")))))
    W)
   (t
    (maximal-of-order (past x)))))
;;
;; Should we make a relativized version?  If so how define (jana x rel Y))?
;; Well, clearly it should mean jana of x regarded as an elt of the subposet Y
;; (it should not just be (intersect Y (jana x)) when these differ)


(deff prec (x y) 
  "\
  Returns `t' if x precedes y with irreflexive convention (or more generally,
  with whatever convention is embodied in `past').
  "
  (&localize x y)
  (if (memq x (past y)) t nil))         
  ;
  ; Would `member' have been better?  Well, not if x and y are guaranteed to be
  ; symbols, as now they are.


(deff prepare-substrate 
        (&key 
	     ((:name name) 'P)
	     ((:elts symbols) 'void elts-supplied)
	     ((:N N) 0 N-supplied) 
	     ((:anon anon) t))
  "\
 The possible args (all keywords) are

  :elts = a list containing the symbols to serve as the poset elements
  :N    = the number of elements of the poset to be created
  :anon = t (default) if symbols should be uninterned (``anonymous'')
  :name = a symbol to serve as the name of the new poset (default=`P')

 Specify either :elts or :N, not both.
 If :elts is given then :N and :anon are ignored.
 If :N is given then the N elements will be 
   anonymous case: uninterned symbols named 0  1  2 ...
   non-anonymous:  the interned symbols    e0 e1 e2 ...
 Further, all pasts will be nil, as will all plists except that 
 each element-symbol will have been given an integer ``key for collation'' as
 its `kfc' property, and `P' itself will have  `poset' as its `type' property.
 We return the symbolic name of the new poset. 
 BEWARE: 
   @ Don't accidentally use name of poset as one of its elements!!
   @ Using a constant list like '(a b c d) as :elts can be dangerous!  \
  "
  (&localize name symbols N anon elts-supplied N-supplied j $)
  (&bind-too $)
 ;----------------------------------------------------------------------
 ;/create the symbols to comprise the poset (unless they were suppplied) 
 ;----------------------------------------------------------------------
  (unless elts-supplied
    (setq 
     symbols
     (cond 
      (anon
       (loop 
	for j from 0 below N collect 
	(make-symbol 
	 (tcl-or-elisp
	  (format nil "~s" j)
	  (format     "%s" j)))))
      (t 
       (loop 
	for j from 0 below N 
	collect 
	(intern 
	 (tcl-or-elisp
	  (format nil "e~s" j)
	  (format     "e%s" j))))))))
 ;-------------------------------------
 ;/set pasts to nil and install labels
 ;-------------------------------------
  (loop 
   for j from 0 below (length symbols) do
   (setq $ (nth j symbols))
   (setf (past $) nil)
   (setf (symbol-plist $) (list 'kfc  j)))
   ;(setf (symbol-plist $) (list 'label (cons name j))))
 ;-------------------------------------
 ;/initialize the plist of name of poset
 ;-------------------------------------
  (setplist name (list 'type 'preorder))
 ;--------------------------
 ;/finish up and return name
 ;--------------------------
  (set name symbols)
  name) 
  ;
  ;  Rationale for offering uninterned symbols as elts:  When the elts are
  ;  uninterned symbols, there can never be clashes with other symbols. 
  ;  Also anonymous symbols are more like unlabelled elts.
  ;
  ; from older version of documentation.:
  ;   Further, all pasts will be nil, as will all plists except that (for
  ;   example) the element-symbol e4 will have as its `label' property the
  ;   ordered pair (P . 4) (assuming the poset is named `P'), and `P' itself
  ;   will carry the information that its `type' is poset.
  ;
(defalias 'prepare-disconnected-poset 'prepare-substrate)

;_____________________________________________________________________

(deff sortedp (P) 
  "\
 Is this digraph (pre-pseudo-order) sorted upward? 
 By definition the test fails if any element precedes itself.  We'd have to
 change this if we wanted to use with reflexive orders. 
 This a completely straightforward - and slow - test of sortedness.
 (A digraph with cycles obviously can't be sorted if it contains cycles, so
 this function is really only useful for preorders.
 It will work for any preorder, but of course it will fail in general for
 SUB-preorders.)  \
  "
  (&localize P j k test)
  (&bind-too j k)
  (varbind test t)
  (kwa j from 0 upto (length P)
  (kwa k from j upto (length P)
  ;     
  (if (prec (nth k P) (nth j P))  (setq test nil))))
  ;
  test)
 
(defalias 'up-sorted-p 'sortedp)
 
(deff down-sorted-p (P) 
  " Is this digraph sorted downward? cf `sortedp' " 
  (&localize P)
  (up-sorted-p (reverse P)))

 ; NOTE: The above is is clearly slow on two counts: it keeps going even after
 ; it has found a pair out of order, and it keeps running thru P to find jth
 ; and kth elts. 
 
;_____________________________________________________________________

;; The sorters here have not been tested all _that_ thoroughly

(deff posort-upward (R)
  "\
 Sorts the preorder with substrate R.
 It is nondestructive, except that it modifies the `sort-label' slots of the
 plists of the poset elts, installing there a natural labeling 0 1 2 ...  of
 the poset.  If desired, these labels can be used to sort the pasts as well.
 The algorithm is just t-close-preorder without the modification of the pasts. 
 BEWARE:
  Depends (via n-max+) on internal operational details of `loop' and `delq'. \
  "
  (&localize R S x counter)
  (varbind counter 0)
 ;------------------------------------------
 ;/ the recursive fcn that does all the work
 ;------------------------------------------
  (fbind inner-sort (S)
    (loop 
     for x in S do
     (unless (get x 'sort-label)
       (when (past x)
	 (inner-sort (past x)))
       (put x 'sort-label counter)
       (incf counter))))
 ;------------------
 ;/ the "outer" part
 ;------------------
  (unless (listp R) 
    (error "Error: argument to posort-upward should be a list"))
 ;----------------------
 ;/erase the sort labels
 ;----------------------
  (loop for x in R do (put x 'sort-label nil))
 ;----------------
 ;/call the sorter
 ;----------------
  (inner-sort R)
 ;-----------------------------------
 ;/sort on the sort labels and return
 ;-----------------------------------
  (Sort (copy-list R) '< :key (lambda (e) (get e 'sort-label))))

 
(deff posort 
            (&key 
	     ((:elements P)  nil  arg-was-substrate)
	     ((:name name)   nil  arg-was-name)
	     ((:direction direction) 'upward))
  "\
 Both arguments are keyword args:   

      {:elements | :name}  --> list of elements or symbolic name of poset
               :direction  --> `upward' (default) or `downward'
 
 The relation to be sorted can be a preorder: it need not be transtive.

 You can supply either a symbol (:name) or the elements themselves (:elements),
 but not both:

  . If substrate is supplied (keyword :elements) then the sorting does not
    affect it but returns a new list of elts in sorted order. 

  . If :name is supplied then, if its plist indicates that P is already sorted
    then nothing is done, except possibly to reverse the substrate if required.
    Otherwise destructive sorting occurs and the resulting sorted substrate is
    made the new value of name.  (Here `destructive' means that the original
    substrate-list is modified, it doesn't refer to the type of rearrangement
    of elements which occurs.)  Finally the plist of name is adjusted to
    reflect the current sorting.   

 In both cases, the sorted substrate is returned.  
 For more, see the documentation of `posort-upward'.  \
  " 
  (&localize P name direction arg-was-substrate arg-was-name currently)
  (if (and arg-was-substrate arg-was-name)
      (error "Can't give both a name and a poset to posort"))
  (if (not (or arg-was-substrate arg-was-name))
      (error "Must supply either a name or a poset to posort"))
  ;
  (when arg-was-substrate                       ; case of no name given
    (return-from posort 
      (case direction
        (upward (posort-upward P))
        (downward (reverse (posort-upward P)))
        (otherwise (error "invalid direction given to posort")))))
  ;
  (when arg-was-name                    ; case where name given
    (let 
        ((currently (get name 'sorted))
         (P (eval name)))
      (set 
       name
       (case direction
         (upward
          (case currently
            (upward P)
            (downward (reverse P))
            (otherwise (posort-upward P))))
         (downward
          (case currently
            (downward P)
            (upward (reverse P))
            (otherwise (reverse (posort-upward P)))))
         (otherwise 
          (error "invalid direction given to posort")))))
    (setf (get name 'sorted) direction)))
   ;
   ; Localization is crucial.
   ;
   ; NB: even if P is already sorted the sorting is not guaranteed to leave it
   ; unchanged, especially for the downward case where the sorting is done by
   ; first sorting upward, then reversing.  


(deff posort-order-upward (P)
  "\
 The argument should be the substrate of an order. 
 It must be TRANSITIVE, a preposet won't do.
 We return a new list of the substrate elements sorted upward.  \
  "
  (&localize P x y)
  (image on P of (put $ 'past-size (o length past $)))
  (sort 
    (copy-list P)
    (lambda (x y) 
      (< 
       (get x 'past-size)
       (get y 'past-size)))))
  ;;
  ;; erase past-size from list?
  ;
  ; method is simply to sort on (o card past x)


(deff sort-pasts-downward (P)
  "\
 First sorts P itself using `posort-upward' and then sorts downward all of its
 pasts, using the resulting sort-labels.
  "
  (&localize P x)
  (posort-upward P)
  (loop 
   for x in P do
   (setf 
     (past x) 
     (Sort (past x) #'> :key (lambda (x) (get x 'sort-label)))))
  "Pasts sorted downward")

;_____________________________________________________________________

(deff maximal-of-preorder (P)
  "\
 We find the maximal layer of a PRE-order in the most direct possible manner.
  "
  (&localize P M x)
  (loop 
   with M = P 
   for x in P
   do (setq M (less M (past x)))
   finally (return M)))
 ;
 ; This will be slow relative to a version which takes advantage of a prior
 ; sorting of the preposet.  On the other hand, the sorting itself has much the
 ; same slowness built in (verdad?)
 ;
 ; However, if we knew that P was an order, and not just a pre-order, then we
 ; could certainly be more efficient.  Check this by comparisons below


(deff n-max+ (S)
  "\
 The argument S should be a preorder (i.e. acyclic).  If it is transitive
 as well (so an order) then we return the set of its maximal elements,
 otherwise we return some superset thereof (hence the `+' in the name).
 BEWARE: We rely heavily on the internal workings of `delq' and `loop'. 
 Advertencia: This function is destructive: it modifies the list S. 
  "
  (&localize S cons y)
  (loop 
    for cons on S do
  (loop 
    for y in (o past car cons)
    do (setq S (delq y S)))) 
  S)
 ; 
 ; Notes
 ;
 ; This relies heavily on the internal workings of `delq', and `loop'.
 ; It needs that `delq' splices out just the right cons when it deletes an
 ; element from S.  It also needs that the loop index-pointer "stay put"
 ; when this happens, eg, it would be bad if `loop' used a numerical index to
 ; find its place, or looped over a *copy* of S. 
 ; We could remove this dependence by expanding the macros and using the
 ; resulting code explicitly!
 ;
 ; The advantage of this algorithm over that of maximal-of-preorder is that
 ; once an element is deleted, we never look at its past.  Obviously, this
 ; saves time.
 ;
 ; FOR CHECKOUT, can put the following just before inner loop
 ;
 ;   (princ (format "S=%s"S))
 ;   (princ (format "[%s], " (car cons)))

(deff maximal-of-order (P) 
  "\
 Returns the maximal layer of the order P.
 The relation P must be acyclic AND transitive.
 This is just the same as n-max+ except it copies P to avoid changing it.
 (So use n-max+ if you want a destructive version of this.)
  "
  (&localize P)
  (n-max+ (copy-list P)))
 
(defalias 'maximal 'maximal-of-order)


(deff minimal (S)
  " Returns the set of minimal elements of the SUB-poset S."
  (&localize S $)
  (loop 
   for $ in S
   unless (meetp S (past $))
   collect $))
 ;
 ;;   if (null (pastof :elt $ :rel S))
 ;
 ; The following alternative is much slower (why? use of disjoint should
 ; probably be faster!)
 ;   (image of  (if (disjointp S (past $)) (setq M (cons $ M))) on S out nil)
;; Try just this alone substituted for if clause in above loop:
 ;   if (disjointp S (past $))

;-------------------------------------------------------------------------

(deff compute-futures-for  (name) 
  "\
 Let *P be the poset named (i.e. the poset pointed to by the symbol to 
 which `name' gets bound) and suppose that that symbol is `P' (so P => *P).  
 We find the futures of the elements of *P, put them in the `future' slots of
 the @ plists, and record in the plist of `P' that its futures have been
 computed. 
  " 
  (&localize name P x y)
  (varbind P (eval name))
 ;------------------------------
 ; initialize all futures to nil
 ;------------------------------
  (image on P of (put $ 'future nil)) 
 ;--------------------------
 ; loop over pairs y < x
 ;--------------------------
  (mapcar                             
   (lambda (x)
     (mapcar                          
      (lambda (y)
	(push x (get y 'future)))
      (past x)))
   P)
  (pushnew 'futures (get name 'already-computed)))
;_________________________________________________________________________

(deff compute-links-for (name)
  "\
 The argument should be a symbol naming a poset (not just a preposet).  
 We compute the ``immediate past'' of each element of the poset and install it
 in the `jana' slot of the element's plist. \
  "
  (&localize name)
  (image on (eval name) of (put $ 'jana (maximal-of-order (past $))))
  (pushnew 'links  (get name 'already-computed)))


(defvar *nontrans* 
  ""
  "\
 Place to return information on failure of transitivity (in lieu of multiple
 return values)")
 ;
 ; This way of proceeding could be improved, see `~/lisp/2:SUGGESTIONS'
 
(deff transitive-p (R &key ((:type type) 'preorder))
  "\
 Usage:  (transitive-p R :type TYPE)  where TYPE can be 
   either 
          relation = digraph = pre-pseudo-order     (for a general relation)
   or 
          preorder = pre-order = preposet = acyclic (for an acyclic relation)

 If :type is not specified, it defaults to `preorder'.
 If R is not transitive, then the first element z which fails the test 
 (meaning that there are x and y such that  { x < y < z but not x < z })
 is placed into the global variable `*nontrans*' (in lieu of having it be a
 second return value, which is impossible in elisp). \
  "
  (&localize R type x)
  (makunbound '*nontrans*)
  (case type
   ;---------------------
   ;/acyclic case
   ;---------------------    
    ((preorder pre-order preposet acyclic) 
     (loop 
      for x in R
      unless (subsetp 
	      (union-m (mapcar #'past (o maximal-of-order past x)))
            ;;(union-m (mapcar #'past (jana x)                   ))
              (past x))
      do 
      (setq *nontrans* x)
      and return nil
      finally (return t)))
   ;---------------------
   ;/general case
   ;---------------------    
    ((digraph pre-pseudo-order relation) 
     (loop 
      for x in R
      unless (subsetp 
              (union-m (mapcar #'past (past x)))
              (past x))
      do 
      (setq *nontrans* x)
      and return nil
      finally (return t)))
   ;-------------
   ;/unknown case
   ;-------------
    (otherwise (error "Unknown type given to `transitive-p'"))))
 ;
 ; In the acyclic case, it works simply by checking for each element whether
 ; its past contains its immediate past's past. (To see that this jana
 ; suffices, just think about a simple chain.)
 ; notice that we use ` maximal-of-order' on a pre-order, which is kosher in
 ; this case, though not good form.
 
(deff non-transitive-p (P)
  "\
 A completely straightforward test for transitivity.  Return nil if P is
 transitive, otherwise return the first non-transitive triple found. 
 In most cases this is much slower than `transitive-p' "
  (&localize P x y z)
  (loop 
   for x in P 
   thereis
   (loop 
    for y in P
    if (prec x y)
    thereis
    (loop 
     for z in P 
     if (and (prec x y) (prec y z) (not(prec x z)))
     return (list x y z)))))

(deff t-close-preorder (R)
  "\
 Transitively closes the preorder with substrate R by adjusting the pasts 
 of its elements.
 Modifies the `done' slots of the plists of the poset elts, installing there
 a natural labeling 0 1 2 ...  of the poset.  If desired, this can be used to 
 sort the poset or its pasts (which in general are *not* sorted automatically).
 Returns only a string saying that preorder has become an order.   
 Calls `n-max+' and `past-of'.
 BEWARE: depends (via n-max+) on internal operational 
         details of `loop' and `delq'.  \      
  "
  (&localize R S x counter)
  (varbind counter 0)
 ;------------------------------------------
 ;/ the recursive fcn that does all the work
 ;------------------------------------------
  (fbind inner-close (S)
    (loop 
     for x in S do
     (unless (get x 'done)
       (when (past x)
	 (inner-close (past x))
	 (setf (past x) (past-of :set (o n-max+ past x) :inclusive t)))
       (put x 'done counter)
       (incf counter))))
 ;------------------
 ;/ the "outer" part
 ;------------------
  (unless (listp R) 
    (error "Error: argument to t-close-preorder should be a list"))
  (loop for x in R do (put x 'done nil))
  (inner-close R)
  "preorder has become order")
  ;
  ; Notes: 
  ;
  ;  Seems likely that current incarnation of this function automatically
  ;  sorts the pasts downward, as a side effect.
  ;
;;; [perhaps add :arg-is name to this?]


(deff t-close-preorder-alt (P)
  "\
 Replaces a(n upwardly) PRESORTED acyclic relation on P by its transitive
 closure.  (If substrate is not sorted then some relations might be left out,
 no spurious ones will be added, in any case.)
 Slower than `t-close-preorder'.
  "
  (&localize P N N-1 i j k x y z)
  (&bind-too i j k x z y)
  (when *carefully*
    (unless (up-sorted-p P)
      (error "Unsorted relation given to `t-close-preorder-alt'")))
  (varbind 
    N (1- (length P))
    N-1 (1- N))
  (kwa i 
     from 0 to (- N 2)  (setq x (nth i P))
  (kwa j 
     from (1+ i) to N-1 (setq y (nth j P))
  (when (prec x y) 
  (kwa k 
     from (1+ j) to N   (setq z (nth k P))
     (if
       (prec y z) (pushnew x (past z)))))))
  (if *carefully*
   "Finished."
   "Finished, but did you presort?"))
  ;
  ; Notes
  ;  It upsets sorting of pasts.
  ;  Seee for explanation of how it works.
  ;
  ;  It may be absurdly slow in some cases, but not in others?  must check. 
  ;  Also coercing the substrate to a vector would speed things up a bit.
  ;
  ;  To avoid sorting altogether can use "?Warschull's alorithm", which is
  ;  simply a triple loop:  {  for y for x for z  if x<y<z then x<z }
  ;  NB: the order of middle elt first is crucial!
  ; 
;;; (defalias 'transitively-close-preposet 't-close-preorder-alt)


(deff t-close-digraph (P)
  "\
 Replaces an arbitrary relation by its transitive closure, the relation being 
 represented by the operator `past' on P.  The name refers to the fact
 that a relation can be thought of as a digraph with vertex-set P or 
 a ``pre-pseudo-order''with P as substrate. 
  "  
  (&localize P j reps)
  (varbind 
    j 0
    reps (o ceiling log_2 length P))
  (kwa 
   j from 1 to reps
   (image 
    out nil on P of
    (setf 
     (past $)
     (union-of (past $) (union-m (mapcar 'past (past $)))))))
  "Transitive closure complete")
  ;
  ;  alternate last line (see printout)
  ;  (past-of :set (pastof :elt $ :convention inclusive))
  ;
  ; Method is to successively, for each element x,
  ;  set the past of x to  U {past(y) | y <= x}
  ;
  ;; VERY SLOW.  Use only for cases where cycles might be present.
  ;
  ; Notice that transitivity requires reflexivity for members of cycles, 
  ; and so forming closure will of course put in such relations.
  ; In this sense, it might be easier to just define x<x always;
  ; on the other hand, the distinction between x<x and its opposite gives
  ; a simple signature for which elts belong to cycles!
  ;
  ; remark: could simplify code slightly if `prec' was defined reflexively.
  ; 
(defalias 'transitively-close-digraph          't-close-digraph)
(defalias 'transitively-close-pre-pseudo-order 't-close-digraph)



(deff count-chains-for (name)
  "\
 The argument should be a symbol naming a poset.
 We ASSUME that the poset is sorted upwards and that it has a unique
 mimimal element (its ``origin''), which is the first element in the list.
 We compute the number of chains joining this origin to each element of the
 poset and put the answer in `chains' for that element.  To deal with the 
 large integer problem in elisp we use floats. \
  " 
  (&localize name P j x max-reliable)
  (&bind-too j x)
  (varbind max-reliable 1e16)
  (varbind P (eval name))
  (if 
    (and
      (not (eq (get name 'sorted) 'upward)) 
      *carefully*
      (not (sortedp P)))
      (error "ERROR: Poset is not sorted upward for chain counting"))
  (put 
    (car P) 
    'chains
    (tcl-or-elisp 1 1.0))	   ; prepare the origin for the recursion
  (kwa j 
     from 1 upto (length P)	   ; start the recursion with second elt in P
     (setq x (nth j P))
     (put x 'chains (sum (image on (past x) of (get $ 'chains)))))
  (pushnew 'chains  (get name 'already-computed))
  (setq x  (sup (image of (get $ 'chains) on P)))
  (if (> x max-reliable)
      (error "too many chains counted to be reliable: %s" x))
  "Chains counted (but only ones emanating from the poset's first element)")
 ;;
 ;; Comments:  the 0th element will be the origin if P is sorted upwards    
 ;;
 ;;; hardly checked out. 
 ;;; the large integer problem is present again, so we float everything
 
(defalias 'count-chains-of-poset-inayoitwa 'count-chains-for)


(deff count-paths-for (name)
  "\
 The argument should be a symbol naming a poset.
 We ASSUME that the poset is sorted upwards and that it has a unique
 mimimal element (its ``origin''), which is the first element in the list.  
 We compute the number of paths joining the origin to each element of the
 poset and put the answer in `paths' for that element. \
  "
  (&localize name P j x)
  (&bind-too j x)
  (varbind max-reliable 1e16)
  (varbind P (eval name))
  (if 
    (and
     (not (eq (get name 'sorted) 'upward)) 
     *carefully*
     (not (sortedp P)))
      (error "ERROR: Poset is not sorted upward for path counting"))
  (put (car P) 
       'paths
       (tcl-or-elisp 1 1.0))		; prepare the origin for the recursion
  (kwa j from 1 upto (length P)		; then loop over subsequent elts
     (setq x (nth j P))
     (put x 'paths 
       (sum (image of (get $ 'paths) on (jana x)))))
  (pushnew 'paths  (get name 'already-computed))
  (setq x  (sup (image of (get $ 'paths) on P)))
  (if (> x max-reliable)
    (error "too many paths counted to be reliable: %s" x))
  "Paths counted (but only ones emanating from the poset's first element)")
   ;
   ; Comments 
   ;
   ; The 0th element will be the origin if P is sorted upwards    
   ;
   ; This fcn ASSUMES that its argument is the (name of the) global poset, not
   ; some subset thereof.  If not must rewrite it to use relative past. 
 
(defalias 'count-paths-of-poset-inayoitwa 'count-paths-for)


(deff count-spanning-trees-for (name)
  "\
 This fcn assumes that its argument is the (name of the) global poset, not
 some subset thereof.  If not, must rewrite it to use relative past.  It
 computes log_2 of the number of spanning-trees, and puts the answer in the
 plist of its argument.
 For convenience we assume that the first element of the poset is minimal."
  ;
  (&localize name P)
  (varbind P (eval name))
  (unless (null (past (car P))) 
    (error "ERROR: first element is not minimal (for tree-counting)"))
  (put
   name 'log_2-spanning-trees 
   (sum 
    (image 
     of (o log_2 length (past $))
     on (cdr P)))))
  ;
  ;;; hardly checked out.   
  ;; will now give error if poset is not originary
 
(defalias 'count-spanning-trees-of-poset-inayoitwa 'count-spanning-trees-for)


(deff count-spanning-hasse-trees-for (name) 
  "\
 This function assumes that `jana' will return the ``wazazi'' OF ITS
 ARGUMENT  (i.e. it does NOT intersect them with the actual poset, which might
 be a subposet). It computes log_2 of the number of ``spanning Hasse trees'',
 and puts the answer the plist of its argument.
 For convenience we assume that the first element of the poset is minimal." 
  ;
  (&localize name P)
  (varbind P (eval name))
  (unless (null (past (car P))) 
    (error "ERROR: first element is not minimal (for tree-counting)"))
  (put 
   name 'log_2-spanning-hasse-trees 
   (sum
    (image 
     of (o log_2 length (jana $))
     on (cdr P)))))
   ;
   ; Comments 
   ;
   ; This fcn ASSUMES that its argument is the (name of the) global poset, not
   ; some subset thereof.  If not must rewrite it to use relative past. 
   ;
 ;;; hardly checked out. 
 ;; will now give error if poset is not originary!
 
(defalias
  'count-spanning-hasse-trees-of-poset-inayoitwa
  'count-spanning-hasse-trees-for)


;_________________________________________________________________________

(defun level (x)   "gets level of x from its plist" (get x 'level))
 
(defsetf level (x) (v) `(setf (get ,x 'level) ,v))
 
;; these last two fairly useless, could inline them in next

(deff assign-levels-to-sorted-preorder (P)
  "\
 We compute the bottom-up levels and install them.  
 The substrate must be SORTED UPWARD.   \
  "
  (&localize  P)
  (image out nil on P
   of (put 
       $ 'level 
       (if (o not past $) 0
	 (o 1+ sup (mapcar #'level (past $))))))
  "Levels assigned, did you remember to sort first?")
  ;
  ; Notes
  ;
  ; `using (past $) seems to be faster than (jana $) for big posets, of course
  ; latter would be faster if jana's were already computed.
  ;
  ; In theory it could be faster to use in last line something like following
  ; rather than create a whole new list with the levels and then run thru it
  ; again to take sup:
  ;      (1+(loop for y in (past x) maximize (level y)))
  ; on the other hand, in elisp `sup' uses max which is a subr so probably
  ; faster. 
  ; A few experiments confirm no gain from the change (in TCL might though).
 
(deff compute-levels-for (name)
  "\
 The argument should be a symbol naming a poset.  We compute the bottom-up
 levels of its elements, and install them in their plists.  (We first sort 
 the poset if necessary then call `assign-levels-to-sorted-preorder'.)   
 We also install in the symbol the information that levels have been
 computed, for what it's worth.     \
  "
  (&localize name P)
  (varbind P (eval name))
  (if (not (eq (get name 'sorted) 'upward)) 
      (setq P (posort-upward P)))
  (assign-levels-to-sorted-preorder P)
  (pushnew 'levels  (get name 'already-computed)))           
  ;
  ;  This could be made faster if its links (majana) are already computed, and
  ;  if the level finder it calls could use the info.
  ;
  ;; Shouldn't we just use `posort-order-upward' since it's faster??

(defalias 'compute-levels-of-poset-inayoitwa 'compute-levels-for)



(deff disconnected-part (C A)
  "\
 (disconnected-part X Y) => that part of X which is not connected to the
 subset Y of X.  When Y is empty, this is all of X by definition. 
 When X is itself a suborder we interpret ``connected'' *within* X, not
 relative to some ambient poset that may contain X.
 Probably X and Y can be any preorders here, they needn't be transitive.  \
  "
  (&localize C A)
  (cond
   ((null C) nil)
   ((null A) C)
   (t 
    (disconnected-part
     (less C A)
     (less   
      (union-of (past-of :set A :rel C) (future-of :set A :rel C))  
      A)))))
    ;
    ; another strategy above would be to first copy the subposet, then not have
    ; to relativize to C
 
(deff connected-part (whole part)
  "\
 Let `whole' be an order and `part' a subset of it.  Then 
   (connected-part whole part) 
 is that part of `whole' which is connected to `part'.  
 See `disconnected-part'for more explanations. 
  "
  (&localize whole part)
  (less whole (disconnected-part whole part)))



(LUN     "create-relation" (H P)
(defmacro create-relation (&rest H)
  "\
 Input should be something like

      a < b < c   d < c  e < f < g < c

 Returned will be a list of the symbols given, and their pasts will have 
 been set as specified. (No transitive closure is done.)
 This version designed to be usable within compiled functions.
  "
  `(progn
     (let 
         ((H (quote ,H))
          (P nil)) 
      ;-----------------------------------------------------------------
      ;/Gather elements and erase their pasts (but not their properties)
      ;-----------------------------------------------------------------
       (loop 
          for $ in H 
          unless (eq '< $) 
          do (setq P (adjoin $ P)))
       (image of (setf (past $) nil) on P out nil)    
      ;---------------------------
      ;/Install pasts as specified
      ;---------------------------
       (loop
          for $ from 0 below (length H)
          if (eq '< (nth $ H))
          do (pushnew 
              (nth (1- $) H)
              (past (nth (1+ $) H))))
      ;--------------------
      ;/Return list of elts
      ;--------------------
       (reverse P)))))
    ;
    ; Note: the complications are to ensure that the desired side effects of
    ; setting the pasts occur at run time, not at compile time.  This is why
    ; the simpler version contained in "developing" can be used only
    ; interpretively
    ;
    ; We don't localize `$' since it is used within `image'
    ;
    ; Improvement: might want to reverse all the pasts at the end
    ;
(defalias 'make-relation 'create-relation)

(LUN     "create-order" (H P)
(defmacro create-order (&rest H)
  "\
 This is the same as `create-relation' except that transitive closure is taken
 to get an order.  Therefore the input should be a preorder (it should not
 contain cycles).  We don't check for mistakes in this. 
 Returned will be (the substrate of) the order.  \
  "
  `(progn
     (setq P (create-relation ,@H))
     (t-close-preorder P)
     P)))


(deff n-delete-from-poset (elt order)
  "\
 The arguments should be an element of an order and the order itself (meaning 
 its substrate).  We DESTRUCTIVELY delete the elt from the poset and 
 return the new substrate.  Comparison is done with `eq' (i.e. deletion is 
 done with `delq').
 REMEMBER: This function doesn't know about any symbol ``naming'' the poset
 (we'd need a `delete-from-poset-inayoitwa' for that).  Hence if `P' (say) is
 the symbolic name of the poset then you have to reset it by hand, ie we have 
 here a function like `delq', not a macro like `pop'.  \
  "  
  (&localize elt order y)
  (setq order (delq elt order))
  (loop
   for y in order do
   (setf
    (past y)
    (delq elt (past y))))
  order)

(deff copy-subposet (P)
  "\
 The argument should be a subset of the substrate of some poset.  
 We return the substrate of an isomorphic copy of this subposet whose
 elements are *uninterned* symbols with the *same names* as the originals.   
 So beware: Confusion will reign if you try to refer to the elements of the
 copy by name!   \
  " 
  (&localize P C x y)
 ;----------------------------------------------------------------------
 ;/make new substrate, with each old elt pointing to its new counterpart
 ;----------------------------------------------------------------------
  (varbind C (loop for x in P collect (o make-symbol symbol-name x)))
  (fbind counterpart (x) (get x 'temp-point))
  (loop 
   for x in P
   for y in C
   do (put x 'temp-point y))
 ;---------------
 ;/form new pasts
 ;---------------
  (loop 
   for x in P do
   (setf 
    (o past counterpart x)
    (mapcar #'counterpart (past-of :elt x :rel P))))
 ;---------------------
 ;/return new substrate
 ;---------------------
  C)
 ;
 ; Improvements and Notes
 ;
 ; we could add a keyword arg :anon which, if NIL, we use interned symbols with
 ; new names formed from the old names, like `e32*'
 ;
 ; we could assign collation keys to the new elements for future sorting
 ; purposes 
 ;
 ; This function could also be called `subposet' , however it is so only in
 ; the category-theoretic sense (since it is a copy).  
 
(defalias 'copy-poset 'copy-subposet)


(deff count-relations (arg &key ((:arg-is arg-is) 'substrate))
  "\
 The argument should be either the substrate of an order, or a symbol naming an
 order, as specified by the keyword arg :arg-is being respectively `substrate'
 (the default) or `name'.  We return the total number of relations and also
 install it in the plist of the argument in case :arg-is=name.    \
  "
  (&localize arg arg-is named nombre S R)
  (&bind-too nombre S named)
  (case arg-is
    ((substrate elements) (setq named nil S arg))
    ((name jina) (setq named t nombre arg S (eval nombre)))
    (otherwise 
     (error 
      (concat
       "argument type to `count-relations' is neither substrate nor name:"
       " %s")
      arg-is)))
  (varbind R (sum (image of (o length past $) on S)))
  (if named (put nombre 'relations R))
  R)

;; the following is kept only for historical reasons, should be retired
;;
(deff count-relations-for (name) 
  (count-relations name :arg-is 'name))

; (deff count-relations-for (name)
;   "\
;  The argument should be a symbol naming a poset.  We count the total number 
;  of relations and install it in the plist of the argument.  We return the
;  number of relations.         \
;   "
;   (&localize name)
;   (put 
;    name 'relations
;    (sum (image of (o length past $) on (eval name)))))


(deff count-links-for (name)
  "\
 The argument should be a symbol naming a poset.  We count the total number 
 of links and install it in the plist of the argument.
  "
  (&localize name)
  (put 
   name 'links
   (sum (image of (length (jana $)) on (eval name)))))



